/// @file   MGFUtils/MGFProtocol/Solution.hh
///
/// @brief
/// @author Mario Garza-Fabre

#ifndef INCLUDED_MGFUtils_MGFProtocol_Solution_HH
#define INCLUDED_MGFUtils_MGFProtocol_Solution_HH

// Project forward headers
#include <MGFUtils/MGFProtocol/Solution.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh>
// #include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/types.hh>

// External library headers
// C++ headers
#include <string>
// #include <cstdlib> 
//#include <utility/vector1.hh>

namespace MGFUtils {
namespace MGFProtocol {

	/// @brief
class Solution : public utility::pointer::ReferenceCount{

public: // Creation

	/// @brief Constructor 1
  	Solution( );

	/// @brief Constructor 2
	Solution( core::pose::Pose const & );

	/// @brief Destructor
	~Solution(){}

private: // Private Methods

	void generate_extended_pose( core::pose::PoseOP, std::string const & );	

public: // Public Methods
	
	core::pose::PoseOP get_pose(){ return pose; }
	double get_energy(){ return total_energy; }
	void set_energy( double energy ){ total_energy = energy; }

	// SMK added these
	bool is_new() {return is_new_; }
	void set_newness(bool const newness) {is_new_ = newness;}

private: // Private attributes

	core::pose::PoseOP pose;		// Structure
	double total_energy;			// Energy	
	unsigned long evaluation;		// Evaluation at which the solution was found

	// SMK added this attribute. When outputting trajectories, we want to know if a solution
	// that has "made the cut" of retained solutions after a reduction operation is "new"
	// recommended only when Archive desired size == max size
	// Implementation of this system is incomplete 28/01/2016
	bool is_new_;

public: // Public attributes



}; // class Solution


} // MGFProtocol
} // MGFUtils


#endif // INCLUDED_MGFUtils_MGFProtocol_Solution_HH
