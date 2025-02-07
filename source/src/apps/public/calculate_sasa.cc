// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief

// libRosetta headers
//#include <basic/options/option.hh>

#include <core/types.hh>


#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


#include <protocols/minimization_packing/MinMover.hh>

#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/options/option.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


// C++ headers
#include <iostream>
#include <string>

//Auto Headers

#include <basic/options/keys/in.OptionKeys.gen.hh>

//silly using/typedef

using namespace core;
using namespace core::id;
using namespace protocols::minimization_packing;
using namespace pose;
using namespace chemical;
using namespace utility;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void
calculate_sasa()
{
	using namespace pose;
	using namespace conformation;
	using namespace chemical;
	using namespace core::import_pose;
	using namespace core::scoring;


	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Pose pose;
	core::import_pose::pose_from_file( pose, option[ in::file::s ]().vector().front() );

	Pose start_pose = pose;

	std::string weights( "sasa_only" );
	ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( weights ) );

	Real sasa( (*score_fxn)( pose ) );
	std::cout << "Solvent Accessible Surface Area is: " << sasa << std::endl;

	return;

}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		//using namespace core;
		devel::init( argc, argv );

		calculate_sasa();
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
}
