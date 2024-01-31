// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/match/ProteinSCSampler.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

#include <protocols/match/upstream/ProteinSCSampler.hh>
#include <protocols/match/upstream/OriginalScaffoldBuildPoint.hh>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Utility headers
#include <string>

/// Project headers
#include <core/types.hh>

// C++ headers

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>


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

// --------------- Test Class --------------- //

class FoldTreeFromSSTests : public CxxTest::TestSuite {

public:

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// --------------- Test Cases --------------- //
	void test_hello_world() {
		std::cout << "This is a test that passed!" << std::endl;
		TS_ASSERT( true );
	}

	void test_case_1(){
		std::string test_string = "   EEEEE   HHHHHHHH  EEEEE   IGNOR EEEEEE   HHHHHHHHHHH  EEEEE  HHHH   ";
		utility::vector1< std::pair< core::Size, core::Size > > reference_vector = {{4,8},{12,19},{22,26},{36,41},{45,55},{58,62},{65,68}};
		utility::vector1< std::pair< core::Size, core::Size > > test_vector = identify_secondary_structure_spans( test_string );
		std::cout << "Checking that the two vectors in the first case are the same length\n";
		TS_ASSERT_EQUALS(reference_vector.size(), test_vector.size())
		for(core::Size i = 1; i <= reference_vector.size(); ++i){
			core::Size reference_lower = reference_vector[i].first;
			core::Size reference_upper = reference_vector[i].second;
			core::Size test_lower = test_vector[i].first;
			core::Size test_upper = test_vector[i].second;
			TS_ASSERT_EQUALS(reference_lower, test_lower);
			TS_ASSERT_EQUALS(reference_upper, test_upper);
		}
	}

	void test_case_2(){
		std::string test_string = "HHHHHHH   HHHHHHHHHHHH      HHHHHHHHHHHHEEEEEEEEEEHHHHHHH EEEEHHH ";
		utility::vector1< std::pair< core::Size, core::Size > > reference_vector = {{1,7},{11,22},{29,40},{41,50},{51,57},{59,62},{63,65}};
		utility::vector1< std::pair< core::Size, core::Size > > test_vector = identify_secondary_structure_spans( test_string );
		std::cout << "Checking that the two vectors in the first case are the same length\n";
		TS_ASSERT_EQUALS(reference_vector.size(), test_vector.size())
		for(core::Size i = 1; i <= reference_vector.size(); ++i){
			core::Size reference_lower = reference_vector[i].first;
			core::Size reference_upper = reference_vector[i].second;
			core::Size test_lower = test_vector[i].first;
			core::Size test_upper = test_vector[i].second;
			TS_ASSERT_EQUALS(reference_lower, test_lower);
			TS_ASSERT_EQUALS(reference_upper, test_upper);
		}
	}

	void test_case_3(){
		std::string test_string = "EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE";
		utility::vector1< std::pair< core::Size, core::Size > > reference_vector = {{1,9},{11,18},{20,28},{30,30},{32,36},{38,38},{40,40},{42,42},{44,51}};
		utility::vector1< std::pair< core::Size, core::Size > > test_vector = identify_secondary_structure_spans( test_string );
		std::cout << "Checking that the two vectors in the first case are the same length\n";
		TS_ASSERT_EQUALS(reference_vector.size(), test_vector.size())
		for(core::Size i = 1; i <= reference_vector.size(); ++i){
			core::Size reference_lower = reference_vector[i].first;
			core::Size reference_upper = reference_vector[i].second;
			core::Size test_lower = test_vector[i].first;
			core::Size test_upper = test_vector[i].second;
			TS_ASSERT_EQUALS(reference_lower, test_lower);
			TS_ASSERT_EQUALS(reference_upper, test_upper);
		}
	}

};
