// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/canonical_sampling/SimulatedTemperingMover.cc
/// @brief SimulatedTempering methods implemented
/// @author Oliver Lange ( oliver.lange@tum.de )


// Unit Headers
#include <protocols/canonical_sampling/SimulatedTempering.hh>
#include <protocols/canonical_sampling/SimulatedTemperingCreator.hh>


// protocols headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>


#include <protocols/jd2/util.hh>

// core headers
#include <basic/options/option_macros.hh>

#include <basic/Tracer.hh>


#include <core/types.hh>

// numeric headers
#include <numeric/random/random.hh>

// utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

// C++ Headers
#include <cmath>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <basic/options/keys/OptionKeys.hh> // AUTO IWYU For BooleanOptionKey, IntegerOptionKey, Opt...
#include <basic/options/option.hh> // AUTO IWYU For OptionCollection, option

// cmd-line options
OPT_2GRP_KEY( Integer, tempering, reweight, stride )
OPT_2GRP_KEY( Boolean, tempering, temp, jump )
OPT_2GRP_KEY( Real, tempering, temp, offset )

static basic::Tracer tr( "protocols.canonical_sampling.SimulatedTempering" );


bool protocols::canonical_sampling::SimulatedTempering::options_registered_( false );

//Mike: when you want to remove these Macros... leave them at least here as comment - since they provide documentation
void protocols::canonical_sampling::SimulatedTempering::register_options() {
	if ( !options_registered_ ) {
		options_registered_ = true;
		TemperingBase::register_options();
		NEW_OPT( tempering::temp::offset, "offset for score (effectively scales all weights)", 40 );
		NEW_OPT( tempering::reweight::stride, "every X trials update the weights distribution - 0 for no reweighting", 0 );
		NEW_OPT( tempering::temp::jump, "if true we can jump to any temperature instead of +/- 1 level", false );
	}
}

namespace protocols {
namespace canonical_sampling {
using namespace core;




SimulatedTempering::SimulatedTempering() {
	set_defaults();
}

SimulatedTempering::SimulatedTempering( SimulatedTempering const & ) = default;


/// @brief callback executed before any Monte Carlo trials
/// use to fill count_ with 0
void
SimulatedTempering::reset_raw_counter() {
	counts_.resize( 0 );
	counts_.assign( n_temp_levels(), 1 ); //initialize with 1 to avoid div by zero...
}

void
SimulatedTempering::initialize_simulation(
	pose::Pose & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
	core::Size cycle //default=0; non-zero if trajectory is restarted
) {
	Parent::initialize_simulation(pose, metropolis_hastings_mover,cycle);
	total_count_ = 0;
	reset_raw_counter();
	tr.Debug << std::setprecision(2);
	if ( weights_.size() != n_temp_levels() ) {
		weights_.clear();
		weighted_counts_.clear();
		for ( core::Size ct = 0; ct < n_temp_levels(); ++ct ) {
			weights_.push_back( 1.0 );
			weighted_counts_.push_back( 0 );
		}
	}
}

void
SimulatedTempering::reweight() {
	for ( core::Size i = 1; i <= counts_.size(); ++i ) {
		weighted_counts_[ i ] += 1.0 * counts_[ i ] / weights_[ i ];
	}
	//update weights...
	for ( core::Size i = 1; i <= weights_.size(); ++i ) {
		weights_[i] = 1.0 / weighted_counts_[i];
	}
	Real const norm_w ( weights_.back() ); //and normalize by last element...
	for ( core::Size i = 1; i <= weights_.size(); ++i ) {
		weights_[i] /= norm_w;
	}
}

void
SimulatedTempering::finalize_simulation(
	pose::Pose& pose,
	protocols::canonical_sampling::MetropolisHastingsMover const & mhm
) {
	finalize_simulation( mhm.output_name() );
	Parent::finalize_simulation( pose, mhm );
}

void
SimulatedTempering::finalize_simulation( std::string const& output_name ) {
	reweight();
	write_to_file( stats_file(), output_name, weighted_counts_ );
	reset_raw_counter();

}

core::Real
SimulatedTempering::temperature_move( core::Real score ) {
	check_temp_consistency();
	if ( !time_for_temp_move() ) return temperature();
	//temperature increase, decrease or wait?
	core::Size new_temp( current_temp() );
	core::Size const nlevels( n_temp_levels() );
	if ( temperature_jumps_ ) {
		new_temp=numeric::random::rg().random_range( 1, nlevels );
	} else {
		Real const r1( numeric::random::rg().uniform() );
		if ( r1 > 0.5+self_transition_*0.5 ) {
			++new_temp;
		} else if ( r1 < 0.5-self_transition_*0.5 ) {
			--new_temp;
		}
	}
	if ( new_temp > nlevels ) new_temp = nlevels;
	if ( new_temp < 1 ) new_temp = 1;

	//make the decision based on score and temperature ratio
	Real real_temp( temperature() );
	if ( new_temp!=current_temp() ) {
		Real const prefac( weights_[ new_temp ]/weights_[ current_temp() ] );
		Real const temp_ratio =
			(temperature() - temperature( new_temp ))/(temperature() * temperature( new_temp ));
		if ( numeric::random::rg().uniform() < std::min( 1.0, prefac*std::exp(-(score+score_offset_)*temp_ratio) ) ) {
			set_current_temp( new_temp );
			real_temp = temperature();
			tr.Debug << "set new temperature to level " << new_temp << " T=" << real_temp << std::endl;
		}
	}
	++counts_[ current_temp() ];
	++total_count_;
	if ( reweight_stride_ > 0 && total_count_ % reweight_stride_ == 0 ) {
		reweight();
		write_to_file( stats_file(), jd2::current_output_name(), weighted_counts_ );
		reset_raw_counter();
	}
	return real_temp;
}



protocols::moves::MoverOP
SimulatedTempering::clone() const
{
	return utility::pointer::make_shared< protocols::canonical_sampling::SimulatedTempering >(*this);
}

protocols::moves::MoverOP
SimulatedTempering::fresh_instance() const
{
	return utility::pointer::make_shared< SimulatedTempering >();
}

void
SimulatedTempering::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) {
	Parent::parse_my_tag( tag, data );
	//simple options
	score_offset_ = tag->getOption< Real >( "score_offset", 40.0 );
	temperature_jumps_ = tag->getOption< bool >( "temperature_jumps", false );
	reweight_stride_ = tag->getOption< core::Size >( "reweight_stride", core::Size(0) );
}

/// handling of options including command-line
void SimulatedTempering::set_defaults() {
	self_transition_ = 0.0;
	temperature_jumps_ = false;
	score_offset_ = 40;
}

/// @brief Assigns user specified values to primitive members using command line options
void SimulatedTempering::init_from_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;
	tr.Debug << "initialize from options..." << std::endl;
	Parent::init_from_options();
	temperature_jumps_ = option[ tempering::temp::jump ]();
	reweight_stride_ = option[ tempering::reweight::stride ]();
	score_offset_ = option[ tempering::temp::offset ]();
}

bool SimulatedTempering::initialize_from_file( std::string const& filename ) {
	utility::vector1< core::Real > temperatures;
	weights_.clear();

	utility::io::izstream in( filename );
	if ( !in.good() ) {
		tr.Error << "cannot open file " << filename << std::endl;
		return false;
	}
	std::string line;
	getline( in, line );
	std::istringstream line_stream( line );
	std::string tag;
	core::Size n_levels( 0 );

	line_stream >> tag >> n_levels >> score_offset_;

	bool line_format( false );
	if ( !line_stream.good() ) {
		tr.Error << "format not recognized in temperature file: " << filename << " at line " << line << std::endl;
		tr.Error << "excpected TEMPERING or TEMPERING_TABLE" << std::endl;
		return false;
	}

	if ( tag == "TEMPERING" ) {
		line_format=true;
	} else if ( tag == "TEMPERING_TABLE" ) {
		line_format=false;
	} else {
		tr.Error << "format not recognized in temperature file: " << filename << " at line " << line << std::endl;
		tr.Error << "excpected TEMPERING or TEMPERING_TABLE" << std::endl;
		return false;
	}

	if ( line_format ) {
		Real temp, weight, count, wcount;
		for ( core::Size ct=1; ct <= n_levels; ++ct ) {
			line_stream >> temp >> weight >> wcount >> count;
			temperatures.push_back( temp );
			weights_.push_back( weight );
			weighted_counts_.push_back( wcount );
			//ignore the counts
		}
		if ( !line_stream.good() ) {
			tr.Error << "format not recognized in temperature file: " << filename << " at line " << line << std::endl;
			tr.Error << "excpected TEMPERING N t1 w1 c1 t2 w2 c2 ... tN wN cN" << std::endl;
			return false;
		}
		return true;
	}

	// table format
	while ( getline( in, line ) ) {
		std::istringstream line_stream2( line );
		Real temp, weight, wcount;
		line_stream2 >> temp;
		if ( !line_stream2.good() ) tr.Error << "format error in temperature file: " << filename << " at line " << line << std::endl;
		line_stream2 >> weight;
		if ( !line_stream2.good() ) {
			tr.Warning << "no weights in temperature file: " << filename << " initialize with 1.0" << std::endl;
			weight = 1.0;
		}
		line_stream2 >> wcount;
		if ( !line_stream2.good() ) {
			wcount=1;
		}
		temperatures.push_back( temp );
		weights_.push_back( weight );
		weighted_counts_.push_back( wcount );
	}
	set_temperatures( temperatures );
	return true; //succesfully initialized
}

void SimulatedTempering::write_to_file( std::string const& file_in, std::string const& output_name, utility::vector1< Real > const& wcounts ) {

	//either write all these things to a single file, or write to <jobname>.file_in
	std::string file;
	if ( stats_silent_output() ) {
		file=file_in;
	} else {
		file=output_name+"."+file_in;
	}

	//open file
	utility::io::ozstream out( file, stats_line_output() ? std::ios::app : std::ios::out );

	//write
	if ( stats_line_output() ) { //line format
		out << "TEMPERING " << n_temp_levels() << " " << score_offset_ << " " << total_count_;
		for ( core::Size i=1; i <= n_temp_levels(); ++i ) {
			out << " " << temperature( i ) << " " << weights_[ i ] << " " << wcounts[ i ] << " " << counts_[ i ];
		}
		out << " " << output_name << std::endl;
	} else { //table format
		out << std::setw( 10 );
		out << "TEMPERING_TABLE " << n_temp_levels() << " " << score_offset_ << " " << " " << total_count_ << " " << output_name << std::endl;
		for ( core::Size i=1; i <= n_temp_levels(); ++i ) {
			out << temperature( i ) << " " << weights_[ i ] << " " << wcounts[ i ] << " " << counts_[ i ] << std::endl;
		}
	}
}

std::string SimulatedTempering::get_name() const {
	return mover_name();
}

std::string SimulatedTempering::mover_name() {
	return "SimulatedTempering";
}

void SimulatedTempering::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	TemperingBase::attributes_for_tempering_base( attlist, xsd );
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "score_offset", xsct_real, "Offset for score (scales all weights)", "40.0" )
		+ XMLSchemaAttribute::attribute_w_default( "temperature_jumps", xsct_rosetta_bool, "Jump to any temperature, not just by one level", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "reweight_stride", xsct_rosetta_bool, "How many trials between automatic reweighting", "false" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Single process switching stochastically between temperature levels", attlist );
}

std::string SimulatedTemperingCreator::keyname() const {
	return SimulatedTempering::mover_name();
}

protocols::moves::MoverOP
SimulatedTemperingCreator::create_mover() const {
	return utility::pointer::make_shared< SimulatedTempering >();
}

void SimulatedTemperingCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SimulatedTempering::provide_xml_schema( xsd );
}



} //moves
} //protocols

