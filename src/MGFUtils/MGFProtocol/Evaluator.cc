/// @file   MGFUtils/MGFProtocol/Evaluator.cc
///
/// @brief
/// @author Mario Garza-Fabre

#include <MGFUtils/MGFGlobal.hh>
#include <MGFUtils/MGFProtocol/Evaluator.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

// #include <core/scoring/ScoreFunction.hh>
// #include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <numeric/xyzVector.hh>

// #include <iostream>
// #include <fstream>

namespace MGFUtils {
namespace MGFProtocol {

/////////////////////////////////////////////////////////////////////////////
// Constructor
// Creates an evaluator according to the given identifier
Evaluator::Evaluator( std::string const & identifier ){

	// Set evaluator identifier and configure
	set_evaluator( identifier );

}
/////////////////////////////////////////////////////////////////////////////
// Sets and cofigure the evaluator according to the given identifier
void Evaluator::set_evaluator( std::string const & identifier ){

	// Set evaluator identifier and configure
	evaluator_identifier = identifier;
	create_score_function();

}
/////////////////////////////////////////////////////////////////////////////
// Constructs an evaluation function according to the provided evaluation identifier
void Evaluator::create_score_function( ){


	if( evaluator_identifier == "STAGE1" ){

		score_function = core::scoring::ScoreFunctionFactory::create_score_function( "score0" );

	}else if( evaluator_identifier == "STAGE2" ){

		score_function = core::scoring::ScoreFunctionFactory::create_score_function( "score1" );

	}else if( evaluator_identifier == "STAGE3_A" ){

		score_function = core::scoring::ScoreFunctionFactory::create_score_function( "score2" );

	}else if( evaluator_identifier == "STAGE3_B" ){

		score_function = core::scoring::ScoreFunctionFactory::create_score_function( "score5" );

	}else if( evaluator_identifier == "STAGE4" ){

		score_function = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );

	}else{

		std::cout << "\n\tERROR: Unrecognised evaluator identifier: " << evaluator_identifier << std::endl << std::endl;
		exit(1);

	} 

}
/////////////////////////////////////////////////////////////////////////////
// Evaluates the given solution, and stores the obtained scores by setting the corresponding solution's attributes
void Evaluator::evaluate( Solution & solution ){

	// Retrieve pose 
	core::pose::PoseOP p = solution.get_pose();

	// Evaluate pose encoded in solution
	(*score_function)( *p );

	// Set solution attributes
	solution.set_energy( p->energies().total_energy() );

}
/////////////////////////////////////////////////////////////////////////////
double Evaluator::compute_RMSD( Solution & s1, Solution & s2){
	return core::scoring::CA_rmsd( *(s1.get_pose()), *(s2.get_pose()) );
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

} // MGFProtocol
} // MGFUtils

