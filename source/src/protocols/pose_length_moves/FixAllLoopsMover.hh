// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/pose_length_moves/FixAllLoopsMover.hh
/// @details connects chains using a very fast RMSD lookback. only works for chains <5 residues. Designed to make loops look within .4 RMSD to naturally occuring loops
/// @author TJ Brunette tjbrunette@gmail.com


#ifndef INCLUDED_protocols_pose_length_moves_FixAllLoopsMover_hh
#define INCLUDED_protocols_pose_length_moves_FixAllLoopsMover_hh


#include <protocols/moves/Mover.hh>

#include <protocols/pose_length_moves/FixAllLoopsMover.fwd.hh>

#include <core/pose/Pose.fwd.hh>

//#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <protocols/indexed_structure_store/SSHashedFragmentStore.fwd.hh>
// C++ Headers
#include <string>
// Utility Headers
#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>

namespace protocols {
namespace pose_length_moves {

class FixAllLoopsMover : public protocols::moves::Mover {
public:
	FixAllLoopsMover();
	void apply( Pose & pose ) override;
	moves::MoverOP clone() const override { return utility::pointer::make_shared< FixAllLoopsMover >( *this ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	protocols::loops::Loops get_loops(core::pose::Pose const & pose);
	int resAdjustmentStartLow_;
	int resAdjustmentStartHigh_;
	int resAdjustmentStopLow_;
	int resAdjustmentStopHigh_;
	int resAdjustmentStartLow_sheet_;
	int resAdjustmentStartHigh_sheet_;
	int resAdjustmentStopLow_sheet_;
	int resAdjustmentStopHigh_sheet_;
	core::Size loopLengthRangeLow_;
	core::Size loopLengthRangeHigh_;
	core::Real rmsThreshold_;
	core::Real max_vdw_change_;
	bool ideal_;
	bool reject_failed_loops_;
	core::Size firstResidue_;
	core::Size lastResidue_;
	std::string allowed_loop_abegos_;
	bool refix_loops_;
	std::string label_loop_;
	std::string fragment_store_path_;
	std::string fragment_store_format_;
	std::string fragment_store_compression_;
	core::Size numb_stubs_to_consider_;

	//time_t start_time_;

	protocols::indexed_structure_store::SSHashedFragmentStoreOP SSHashedFragmentStoreOP_;
};


} // pose_length_moves
} // protocols

#endif
