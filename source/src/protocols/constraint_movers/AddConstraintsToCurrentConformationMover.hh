// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song

#ifndef INCLUDED_protocols_constraint_movers_AddConstraintsToCurrentConformationMover_hh
#define INCLUDED_protocols_constraint_movers_AddConstraintsToCurrentConformationMover_hh

#include <protocols/moves/Mover.hh>

#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

namespace protocols {
namespace constraint_movers {

class AddConstraintsToCurrentConformationMover : public moves::Mover {
public:
	enum class AtomSelector {
		ALL,
		CA_ONLY,
		BB_ONLY,
		SC_TIP_ONLY//only atoms affected by the final chi angle
	};

public:
	typedef core::pack::task::TaskFactoryOP TaskFactoryOP;

	AddConstraintsToCurrentConformationMover();
	~AddConstraintsToCurrentConformationMover() override;

	void
	apply( Pose & pose ) override;

	virtual core::scoring::constraints::ConstraintCOPs
	generate_constraints( Pose const & pose );


	void task_factory( TaskFactoryOP tf );
	void residue_selector( core::select::residue_selector::ResidueSelectorCOP selector );

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void
	parse_my_tag( TagCOP, basic::datacache::DataMap & ) override;

	AtomSelector atom_selector() const {
		return atom_sele_;
	}

	void set_atom_selector( AtomSelector const setting ) {
		atom_sele_ = setting;
	}

	bool CA_only() const { return atom_sele_ == AtomSelector::CA_ONLY; }
	bool bb_only() const { return atom_sele_ == AtomSelector::BB_ONLY; }
	bool sc_tip_only() const { return atom_sele_ == AtomSelector::SC_TIP_ONLY; }


	bool       & use_distance_cst() { return use_distance_cst_; }
	bool       & use_harmonic_func() { return use_harmonic_func_; }
	bool       & use_bounded_func() { return use_bounded_func_; }
	bool       & inter_chain() { return inter_chain_; }
	core::Real & cst_weight() { return cst_weight_; }
	core::Real & max_distance() { return max_distance_; }
	core::Real & coord_dev() { return coord_dev_; }
	core::Real & bound_width() { return bound_width_; }
	core::Size & min_seq_sep() { return min_seq_sep_; }
	bool       const & use_distance_cst() const { return use_distance_cst_; }
	bool       const & use_harmonic_func() const { return use_harmonic_func_; }
	bool       const & use_bounded_func() const { return use_bounded_func_; }
	bool       const & inter_chain() const { return inter_chain_; }
	core::Real const & cst_weight() const { return cst_weight_; }
	core::Real const & max_distance() const { return max_distance_; }
	core::Real const & coord_dev() const { return coord_dev_; }
	core::Real const & bound_width() const { return bound_width_; }
	core::Size const & min_seq_sep() const { return min_seq_sep_; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::Size
	find_best_anchor( core::pose::Pose const & pose ) const;

	/// @brief parse "task_operations" XML option (can be employed virtually by derived movers)
	virtual void parse_task_operations(
		TagCOP,
		basic::datacache::DataMap const &
	);

	core::scoring::constraints::ConstraintCOPs
	generate_coordinate_constraints(
		core::pose::Pose const & pose,
		core::select::residue_selector::ResidueSubset const & subset ) const;

	core::scoring::constraints::ConstraintCOPs
	generate_atom_pair_constraints(
		core::pose::Pose const & pose,
		core::select::residue_selector::ResidueSubset const & subset ) const;

private:
	bool use_distance_cst_;
	bool use_harmonic_func_;
	AtomSelector atom_sele_ = AtomSelector::CA_ONLY;
	bool inter_chain_;
	core::Real cst_weight_;
	core::Real max_distance_;
	core::Real coord_dev_;
	core::Real bound_width_;
	core::Size min_seq_sep_;

	core::select::residue_selector::ResidueSelectorCOP selector_;

	bool use_bounded_func_; //Based on the tag, are we using the Bound[ed]Func or the default for the constraint type (HarmonicFunc for CoordinateConstraint, SOGFunc for AtomPairConstraint)

};

} // moves
} // protocols

#endif
