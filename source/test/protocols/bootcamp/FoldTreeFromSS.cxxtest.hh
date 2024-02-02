// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author fatimadavila (fatima.a.davila@gmail.com)


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
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <protocols/bootcamp/fold_tree_from_ss.hh>


// C++ headers

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>


using namespace protocols::match;
using namespace protocols::match::upstream;

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

	void test_case_1(){
		std::string test_string = "   EEEEE   HHHHHHHH  EEEEE   IGNOR EEEEEE   HHHHHHHHHHH  EEEEE  HHHH   ";
		utility::vector1< std::pair< core::Size, core::Size > > reference_vector = {{4,8},{12,19},{22,26},{36,41},{45,55},{58,62},{65,68}};
		utility::vector1< std::pair< core::Size, core::Size > > test_vector = protocols::bootcamp::identify_secondary_structure_spans( test_string );
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
		utility::vector1< std::pair< core::Size, core::Size > > test_vector = protocols::bootcamp::identify_secondary_structure_spans( test_string );
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
		utility::vector1< std::pair< core::Size, core::Size > > test_vector = protocols::bootcamp::identify_secondary_structure_spans( test_string );
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

	void test_4(){
		std::cout << "\n\n\nDSSP TEST 1\n\n\n";
		std::string test_string = "   EEEEEEE    EEEEEEE         EEEEEEEEE    EEEEEEEEEE   HHHHHH         EEEEEEEEE         EEEEE     ";
		core::kinematics::FoldTree test_ft = protocols::bootcamp::fold_tree_from_dssp_string(test_string);
		std::cout << test_ft << std::endl;
		TS_ASSERT( test_ft.check_fold_tree() );
	}

	void test_5(){
		std::cout << "\n\n\nDSSP TEST 2\n\n\n";
		std::string test_string = "LEEEEEELLLLEEEEELLLLLEHHHHHHHHHHHHLLLHHHEEEEELLEELLLLLELHHHLLLLLLEEEEEELLLLL";
		core::kinematics::FoldTree test_ft = protocols::bootcamp::fold_tree_from_dssp_string(test_string);
		std::cout << test_ft << std::endl;
		TS_ASSERT( test_ft.check_fold_tree() );
	}
};
