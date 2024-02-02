// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/bootcamp/fold_tree_from_ss.cc
/// @brief  Fold tree from secondary structure elements
/// @author fatimadavila (fatima.a.davila@gmail.com)

// Test headers

#include <protocols/match/upstream/ProteinSCSampler.hh>
#include <protocols/match/upstream/OriginalScaffoldBuildPoint.hh>

// Utility headers
#include <iostream>
#include <string>

/// Project headers
#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

// C++ headers

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>

namespace protocols {
namespace bootcamp  {

using namespace protocols::match;
using namespace protocols::match::upstream;

utility::vector1< std::pair< core::Size, core::Size > >
identify_secondary_structure_spans( 
    std::string const & ss_string
);

core::Size
find_middle_ss(
    core::Size beginning, 
    core::Size end
);

core::kinematics::FoldTree 
fold_tree_from_dssp_string(
    std::string input_dssp
);

utility::vector1< std::pair< core::Size, core::Size > > 
identify_gaps_spans(
    utility::vector1< std::pair< core::Size, core::Size > > ss_spans,core::Size last_residue
);

}
}