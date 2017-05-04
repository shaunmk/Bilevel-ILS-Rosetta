/// @file   MGFUtils/MGFProtocol/Solution.cc
///
/// @brief
/// @author Mario Garza-Fabre

#include <MGFUtils/MGFGlobal.hh>
#include <MGFUtils/MGFProtocol/Solution.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <core/scoring/Energies.hh>
#include <core/chemical/ChemicalManager.hh>

#include <numeric/xyzVector.hh>

// #include <iostream>
// #include <fstream>

namespace MGFUtils {
namespace MGFProtocol {

using core::Size;

/////////////////////////////////////////////////////////////////////////////
// Constructor 1 (no pose is given, but the protein sequence)
// Generate a new extended pose from the input protein sequence
Solution::Solution(){

	// Create new extended pose
	pose = core::pose::PoseOP( new core::pose::Pose() );
	generate_extended_pose( pose, MGFUtils::MGFGlobal::protein_sequence );

	// Initialize parameters
	evaluation = 0;
	total_energy = 0;

}
/////////////////////////////////////////////////////////////////////////////
// Constructor 2 (a pose is given)
// Create a copy of the given pose and use it to initialize the new instantiated solution
Solution::Solution( core::pose::Pose const & p ){

	// Create copy of the given pose
	// pose = p;
	pose = core::pose::PoseOP( new core::pose::Pose( p ) );

	// Initialize parameters
	evaluation = 0;
	total_energy = 0;
  
}
/////////////////////////////////////////////////////////////////////////////
// Constructor 3 (a Solution object is given)
// Copy constructor
// Solution::Solution( SolutionOP s ){

// 	std::cout << ">>>> Solution constructor 3" << std::endl;

	// // Set pose
	// pose = p;

	// // Initialize parameters
	// evaluation = 0;
	// total_energy = 0;
  
// }
/////////////////////////////////////////////////////////////////////////////
void Solution::generate_extended_pose( core::pose::PoseOP extended_pose, std::string const & sequence ){
	
	core::pose::make_pose_from_sequence(
		*extended_pose,
		sequence,
		*( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID ))
	);

	// make extended chain
	for ( Size pos = 1; pos <= extended_pose->total_residue(); pos++ ) {
		if ( ! extended_pose->residue(pos).is_protein() ) continue;
		extended_pose->set_phi( pos, -150 );
		extended_pose->set_psi( pos, 150);
		extended_pose->set_omega( pos, 180 );
	}

}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

} // MGFProtocol
} // MGFUtils

