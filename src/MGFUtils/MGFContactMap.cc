/// @file   MGFUtils/MGF_ContactMap.cc
///
/// @brief
/// @author Mario Garza-Fabre

#include <protocols/moves/MGFContactMap.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>

#include <iostream>
#include <fstream>

namespace MGFUtils {

using core::Size;
using std::string;
using utility::vector1;

/////////////////////////////////////////////////////////////////////////////
MGFContactMap::MGFContactMap( core::pose::Pose const & pose ){

  	computed = false; distance_cutoff = 8.0;
	// std::cout << ">>>> MGFContactMap constructor\n";	
	// std::cout << "\t total_residue() = " << pose.total_residue() << std::endl;

	compute_contact_map( pose );
}
/////////////////////////////////////////////////////////////////////////////
void MGFContactMap::set_pose( core::pose::Pose const & pose ){
	
	// std::cout << ">>>> MGFContactMap constructor\n";	
	// std::cout << "\t total_residue() = " << pose.total_residue() << std::endl;

	compute_contact_map( pose );
}
/////////////////////////////////////////////////////////////////////////////
void MGFContactMap::allocate_memory(){

	// distance_map = (double **)malloc( sizeof(double *) * nres );
	contact_map = (int **)malloc( sizeof(int *) * nres );
	for(core::Size i=0; i<nres; i++){
		// distance_map[i] = (double *)malloc( sizeof(double) * nres );
		contact_map[i] = (int *)malloc( sizeof(int) * nres );

		// Initialize (not required, values are always computed)
		// for(core::Size j=0; j<nres; j++) distance_map[i][j] = contact_map[i][j] = 0;
	}
}
/////////////////////////////////////////////////////////////////////////////
void MGFContactMap::free_memory(){

	if( computed ){

		// if(distance_map!=NULL){
		// 	for(core::Size i=0; i<nres; i++) free( distance_map[i] );
		// 	free(distance_map);
		// } 
		if(contact_map!=NULL){
			for(core::Size i=0; i<nres; i++) free( contact_map[i] );
			free(contact_map); 
		}

	}else{

		std::cerr << "\n\tWARNING: Cannot free memory because CM has not been computed!!!\n";

	}
}
/////////////////////////////////////////////////////////////////////////////
void MGFContactMap::update_distance_map( core::Size r1, core::Size r2, core::Real distance ){
	distance_map[r1-1][r2-1] = distance;
}
/////////////////////////////////////////////////////////////////////////////
void MGFContactMap::update_contact_map( core::Size r1, core::Size r2, core::Real distance ){
	contact_map[r1-1][r2-1] = ( distance <= distance_cutoff )? 1 : 0;
}
/////////////////////////////////////////////////////////////////////////////
double MGFContactMap::read_distance_map( core::Size r1, core::Size r2 ){
	return distance_map[r1-1][r2-1];
}
/////////////////////////////////////////////////////////////////////////////
int MGFContactMap::read_contact_map( core::Size r1, core::Size r2 ){
	return contact_map[r1-1][r2-1];
}
/////////////////////////////////////////////////////////////////////////////
void MGFContactMap::compute_contact_map( core::pose::Pose const & pose ){

	if( !computed ){

		// Set number of residues
		nres = pose.total_residue();

		// Allocate memory for data structures
		allocate_memory();

	}

	// Consider each pair of residues
	for(core::Size r1=1; r1<=nres; r1++){

		// Decide if computing distance from CA or CB atom
		//std::string target_atom1 =  pose.residue_type( r1 ).name1() == 'G'  ?   "CA" : "CB";
		std::string target_atom1 =  "CA";

		// Get coordinates of residue 1
		numeric::xyzVector<core::Real> coord_r1 = pose.residue( r1 ).atom( target_atom1 ).xyz();

		for(core::Size r2=r1+1; r2<=nres; r2++){

			// Decide if computing distance from CA or CB atom
			//std::string target_atom2 =  pose.residue_type( r2 ).name1() == 'G'  ?   "CA" : "CB";
			std::string target_atom2 =  "CA";

			// Get coordinates of residue 2 
			numeric::xyzVector<core::Real> coord_r2 = pose.residue( r2 ).atom( target_atom2 ).xyz();

			// Compute distance between residues
			core::Real distance = coord_r1.distance(coord_r2);

			// Update distance map
			// update_distance_map( r1, r2, distance );

			// Update contact map
			update_contact_map( r1, r2, distance );

		}
	}

	// Set flag
	computed = true;
}
/////////////////////////////////////////////////////////////////////////////
void MGFContactMap::print_distance_map(){

	if( computed ){

		for(core::Size r1=1; r1<=nres; r1++){
			for(core::Size r2=1; r2<=nres; r2++){
				std::cout << read_distance_map( r1, r2 ) << "\t";
			}
			std::cout << std::endl;
		}

	}else{

		std::cerr << "\n\tWARNING: Cannot print DM because it has not been computed!!!\n";

	}
}
/////////////////////////////////////////////////////////////////////////////
void MGFContactMap::print_contact_map(){

	if( computed ){

		for(core::Size r1=1; r1<=nres; r1++){
			for(core::Size r2=1; r2<=nres; r2++){
				std::cout << read_contact_map( r1, r2 ) << "\t";
			}
			std::cout << std::endl;
		}

	}else{

		std::cerr << "\n\tWARNING: Cannot print CM because it has not been computed!!!\n";

	}
}
/////////////////////////////////////////////////////////////////////////////
void MGFContactMap::print_distance_map_tofile( std::string filename ){

	if( computed ){

		std::ofstream file;
		file.open( filename.c_str() );
		for(core::Size r1=1; r1<=nres; r1++){
			for(core::Size r2=1; r2<=nres; r2++){
				file << read_distance_map( r1, r2 ) << "\t";
			}
			file << std::endl;
		}
		file.close();

	}else{

		std::cerr << "\n\tWARNING: Cannot print DM because it has not been computed!!!\n";

	}

}
/////////////////////////////////////////////////////////////////////////////
void MGFContactMap::print_contact_map_tofile( std::string filename ){

	if( computed ){

		std::ofstream file;
		file.open( filename.c_str() );
		for(core::Size r1=1; r1<=nres; r1++){
			for(core::Size r2=1; r2<=nres; r2++){
				file << read_contact_map( r1, r2 ) << "\t";
			}
			file << std::endl;
		}
		file.close();

	}else{

		std::cerr << "\n\tWARNING: Cannot print CM because it has not been computed!!!\n";

	}

}
/////////////////////////////////////////////////////////////////////////////
void MGFContactMap::print_distance_map_tofile_compact( std::string filename ){

	if( computed ){

		std::ofstream file;
		file.open( filename.c_str() );
		for(core::Size r1=1; r1<=nres; r1++){
			for(core::Size r2=r1+1; r2<=nres; r2++){
				file << read_distance_map( r1, r2 ) << std::endl;
			}
		}
		file.close();

	}else{

		std::cerr << "\n\tWARNING: Cannot print DM because it has not been computed!!!\n";

	}

}
/////////////////////////////////////////////////////////////////////////////
void MGFContactMap::print_contact_map_tofile_compact( std::string filename ){

	if( computed ){

		std::ofstream file;
		file.open( filename.c_str() );
		for(core::Size r1=1; r1<=nres; r1++){
			for(core::Size r2=r1+1; r2<=nres; r2++){
				file << read_contact_map( r1, r2 ) << std::endl;
			}
		}
		file.close();

	}else{

		std::cerr << "\n\tWARNING: Cannot print CM because it has not been computed!!!\n";

	}

}
/////////////////////////////////////////////////////////////////////////////
void MGFContactMap::print_contact_map_tofile_compact2( std::string filename ){

	if( computed ){

		std::ofstream file;
		file.open( filename.c_str() );
		for(core::Size r1=1; r1<=nres; r1++){
			for(core::Size r2=r1+1; r2<=nres; r2++){
				if( read_contact_map( r1, r2 ) == 1 ) file << r1 << "\t" << r2 << std::endl;
			}
		}
		file.close();

	}else{

		std::cerr << "\n\tWARNING: Cannot print CM because it has not been computed!!!\n";

	}

}
/////////////////////////////////////////////////////////////////////////////
int MGFContactMap::hamming_contact_map( MGFContactMap &other ){

	// Initialize
	int hamming_distance = 0;

	if( computed && other.is_computed() ){

		// Compute Hamming distance with respect to the given CM
		for(core::Size r1=1; r1<=nres; r1++){
			for(core::Size r2=r1+1; r2<=nres; r2++){

				// Increase HD if contact maps differ w.r.t. these residues
				if( read_contact_map( r1, r2 ) != other.read_contact_map( r1, r2 )) hamming_distance++;

			}
		}

	}else{

		std::cerr << "\n\tWARNING: Cannot compute Hamming Distance because either one of the maps or both have not been computed!!!\n";

	}

	return hamming_distance;
}
/////////////////////////////////////////////////////////////////////////////
double MGFContactMap::euclidean_distance_map( MGFContactMap &other){

	// Initialize
	int euclidean_distance = 0;

	if( computed && other.is_computed() ){

		// Compute Euclidean distance with respect to the given DM
		for(core::Size r1=1; r1<=nres; r1++){
			for(core::Size r2=r1+1; r2<=nres; r2++){

				// Increase ED 
				euclidean_distance += pow( read_distance_map( r1, r2 ) - other.read_distance_map( r1, r2 ), 2);

			}
		}

	}else{

		std::cerr << "\n\tWARNING: Cannot compute Euclidean Distance because either one of the maps or both have not been computed!!!\n";

	}

	return sqrt(euclidean_distance);
}


} // MGFUtils

