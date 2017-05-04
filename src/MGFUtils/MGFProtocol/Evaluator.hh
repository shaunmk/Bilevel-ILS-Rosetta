/// @file   MGFUtils/MGFProtocol/Evaluator.hh
///
/// @brief
/// @author Mario Garza-Fabre

#ifndef INCLUDED_MGFUtils_MGFProtocol_Evaluator_HH
#define INCLUDED_MGFUtils_MGFProtocol_Evaluator_HH

#include <MGFUtils/MGFProtocol/Evaluator.fwd.hh>
#include <MGFUtils/MGFProtocol/Solution.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>

// C++ headers
#include <string>

namespace MGFUtils {
namespace MGFProtocol {

	/// @brief
class Evaluator : public utility::pointer::ReferenceCount{

public: // Creation

	/// @brief Constructor 
  	Evaluator( std::string const & );

	/// @brief Destructor
	~Evaluator(){}

private: // Private Methods

	void create_score_function();	

public: // Public Methods
	
	void set_evaluator( std::string const & );
	std::string get_identifier(){ return evaluator_identifier; }
	void evaluate( Solution & );
	double compute_RMSD( Solution &, Solution & );

private: // Private attributes

	std::string evaluator_identifier;				// String identifying the current settings of the evaluator
	core::scoring::ScoreFunctionOP score_function;	// Evaluation function using by the evaluator

public: // Public attributes

}; // class Evaluator


} // MGFProtocol
} // MGFUtils


#endif // INCLUDED_MGFUtils_MGFProtocol_Evaluator_HH

