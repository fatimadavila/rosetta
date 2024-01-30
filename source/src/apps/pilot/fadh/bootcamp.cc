// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <iostream>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <numeric/random/random.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

int 
main( int argc, char ** argv ) {
    
    std::cout << "Hello World!" << std::endl;

    devel::init( argc, argv );
    utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
    if ( filenames.size() > 0 ) {
        std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
        core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );
        core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();

        core::Real temp = 1.0;
        protocols::moves::MonteCarloOP monte_carlo(new protocols::moves::MonteCarlo(*mypose, *sfxn, temp));

        core::kinematics::MoveMap mm;
        mm.set_bb( true );
        mm.set_chi( true );
        core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
        core::optimization::AtomTreeMinimizer atm;

        // Speed up by copying pose outside of the loop
        core::pose::Pose copy_pose;

        for (int i = 0; i < 100; i++){
            double uniform_random_number = numeric::random::uniform();
            int N = mypose->size();
            core::Size randres = static_cast< core::Size > ( uniform_random_number * N + 1 );
            core::Real pert1 = numeric::random::gaussian();
            core::Real pert2 = numeric::random::gaussian();
            core::Real orig_phi = mypose->phi( randres );
            core::Real orig_psi = mypose->psi( randres );
            mypose->set_phi( randres, orig_phi + pert1 );
            mypose->set_psi( randres, orig_psi + pert2 );
            core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( *mypose );
            repack_task->restrict_to_repacking();
            core::pack::pack_rotamers( *mypose, *sfxn, repack_task );
            copy_pose = *mypose; // Copying pose to speed code up
            atm.run( copy_pose, mm, *sfxn, min_opts );
            *mypose = copy_pose;
            monte_carlo->boltzmann( *mypose );
            core::Real score = sfxn->score( *mypose );
            std::cout << "The score of your pose is: " << score << std::endl;
        }


    } 
    else {
        std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
        return 1;
    }


    return 0;
}
