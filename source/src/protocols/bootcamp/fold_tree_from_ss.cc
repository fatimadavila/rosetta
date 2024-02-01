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
#include <string>
#include <iostream>

/// Project headers
#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <protocols/bootcamp/fold_tree_from_ss.hh>

// C++ headers

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>


namespace protocols {
namespace bootcamp  {

using namespace protocols::match;
using namespace protocols::match::upstream;

utility::vector1< std::pair< core::Size, core::Size > >
identify_secondary_structure_spans( std::string const & ss_string )
{
  utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries;
  core::Size strand_start = -1;
  for ( core::Size ii = 0; ii < ss_string.size(); ++ii ) {
    if ( ss_string[ ii ] == 'E' || ss_string[ ii ] == 'H'  ) {
      if ( int( strand_start ) == -1 ) {
        strand_start = ii;
      } else if ( ss_string[ii] != ss_string[strand_start] ) {
        ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
        strand_start = ii;
      }
    } else {
      if ( int( strand_start ) != -1 ) {
        ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
        strand_start = -1;
      }
    }
  }
  if ( int( strand_start ) != -1 ) {
    // last residue was part of a ss-eleemnt                                                                                                                                
    ss_boundaries.push_back( std::make_pair( strand_start+1, ss_string.size() ));
  }
  for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
    std::cout << "SS Element " << ii << " from residue "
      << ss_boundaries[ ii ].first << " to "
      << ss_boundaries[ ii ].second << std::endl;
  }
  return ss_boundaries;
}

core::Size find_middle_ss(core::Size beginning, core::Size end){
	core::Size middle = (end - beginning)/2 + beginning;
	return middle;
}

core::kinematics::FoldTree fold_tree_from_dssp_string(std::string input_dssp){
	core::kinematics::FoldTree ft;
	utility::vector1< std::pair< core::Size, core::Size > > input_vector = identify_secondary_structure_spans( input_dssp );
	core::Size in_beginning_1 = input_vector[1].first;
	core::Size in_end_1 = input_vector[1].second;
	core::Size middle_1 = find_middle_ss(in_beginning_1, in_end_1);
	ft.add_edge(middle_1, 1, core::kinematics::Edge::PEPTIDE);
	ft.add_edge(middle_1, in_end_1, core::kinematics::Edge::PEPTIDE);
	core::Size middle_2 = find_middle_ss(input_vector[1].second, input_vector[2].first);
	ft.add_edge(middle_1, middle_2, 1);
	ft.add_edge(middle_2, in_end_1 + 1, core::kinematics::Edge::PEPTIDE);
	core::Size end_loop_1 = input_vector[2].first -1;
	ft.add_edge(middle_2, end_loop_1, core::kinematics::Edge::PEPTIDE);
	core::Size jump_counter = 2;
	for(core::Size i = 2; i < input_vector.size(); ++i){
		core::Size in_beginning_2 = input_vector[i].first;
		core::Size in_end = input_vector[i].second;
		core::Size in_beginning_3 = input_vector[i + 1].first;
		core::Size middle = find_middle_ss(in_beginning_2, in_end);
		core::Size middle_3 = find_middle_ss(in_end, in_beginning_3);
		ft.add_edge(middle_1, middle, jump_counter);
		ft.add_edge(middle, in_beginning_2, core::kinematics::Edge::PEPTIDE);
		ft.add_edge(middle, in_end, core::kinematics::Edge::PEPTIDE);
		jump_counter++;
		ft.add_edge(middle_1, middle_3, jump_counter);
		ft.add_edge(middle_3, in_end + 1, core::kinematics::Edge::PEPTIDE);
		core::Size end_loop_2 = input_vector[i + 1].first - 1;
		ft.add_edge(middle_3, end_loop_2, core::kinematics::Edge::PEPTIDE);
		jump_counter++;
	}
	core::Size in_beginning_final = input_vector[input_vector.size()].first;
	core::Size in_end_final = input_vector[input_vector.size()].second;
	core::Size middle_final = find_middle_ss(in_beginning_final, in_end_final);
	core::Size very_end = input_dssp.size();
	ft.add_edge(middle_1, middle_final, jump_counter);
	ft.add_edge(middle_final, in_beginning_final, core::kinematics::Edge::PEPTIDE);
	ft.add_edge(middle_final, very_end, core::kinematics::Edge::PEPTIDE);
	return ft;
}

}
}
