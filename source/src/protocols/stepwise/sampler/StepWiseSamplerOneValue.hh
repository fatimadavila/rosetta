// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerOneValue.hh
/// @brief Generate rotamer for one value. Abstract base class.
/// @author Rhiju Das (rhiju@stanford.edu)


#ifndef INCLUDED_stepwise_sampler_StepWiseSamplerOneValue_HH
#define INCLUDED_stepwise_sampler_StepWiseSamplerOneValue_HH

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerOneValue.fwd.hh>

// Package headers
#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>

#include <utility/vector1.hh> // AUTO IWYU For vector1

typedef utility::vector1<core::Real> ValueList;

namespace protocols {
namespace stepwise {
namespace sampler {

class StepWiseSamplerOneValue : public StepWiseSamplerSized {
public:

	StepWiseSamplerOneValue();

	StepWiseSamplerOneValue(
		ValueList const & values
	);

	StepWiseSamplerOneValue(
		ValueList const & allowed_values,
		std::string const & name
	);

	~StepWiseSamplerOneValue() override;

	/// @brief Initialization
	void init() override {
		set_init( true );
		reset();
	}

	/// @brief Apply the current rotamer to pose
	void apply( core::pose::Pose & ) override{}

	/// @brief Apply the i-th rotamer to pose
	void apply( core::pose::Pose &, core::Size const ) override{}

	/// @brief Get the total number of rotamers in sampler
	core::Size size() const override {
		runtime_assert( is_init() );
		return values_.size();
	}

	/// @brief Get the value of current value
	virtual core::Real value() const {
		runtime_assert( is_init() );
		return values_[id_];
	}

	/// @brief Get the value of i-th value
	virtual core::Real value( core::Size const i ) const {
		runtime_assert( is_init() );
		return values_[i];
	}

	/// @brief Set the allowed values in sampler
	virtual void set_values( ValueList const & setting ) {
		set_and_reinit( values_, setting );
	}

	/// @brief Name of the class
	std::string get_name() const override;

	/// @brief Type of class (see enum in SamplerPlusPlusTypes.hh)
	toolbox::SamplerPlusPlusType type() const override { return toolbox::ONE_VALUE; }

protected: // will be available to derived classes.
	ValueList values_;

private:
	std::string const tag_;
};

} //sampler
} //stepwise
} //protocols

#endif
