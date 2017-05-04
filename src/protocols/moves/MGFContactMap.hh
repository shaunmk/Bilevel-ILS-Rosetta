/// @file   MGFUtils/MGFContactMap.hh
///
/// @brief
/// @author Mario Garza-Fabre

#ifndef INCLUDED_MGFUtils_MGFContactMap_HH
#define INCLUDED_MGFUtils_MGFContactMap_HH

// Project forward headers
#include <protocols/moves/MGFContactMap.fwd.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

// External library headers
// C++ headers
#include <string>
#include <cstdlib> 
//#include <utility/vector1.hh>

namespace MGFUtils {

	/// @brief
class MGFContactMap{

// public: // Types
// private: // Types
// public: // Constants
// private: // Constants
public: // Creation

	/// @brief Constructor 1
  MGFContactMap( ){ computed =false;distance_cutoff = 8.0;}

	/// @brief Constructor 2
	MGFContactMap( core::pose::Pose const & pose );

	/// @brief Destructor
	~MGFContactMap(){}

	/// @brief Copy constructor
	//MGFContactMap( MGFContactMap const & src ){ *this = src; }

private: // Private Methods

	void compute_contact_map( core::pose::Pose const & pose );
	void allocate_memory();
	void update_distance_map( core::Size r1, core::Size r2, core::Real distance );
	void update_contact_map( core::Size r1, core::Size r2, core::Real distance );
	double read_distance_map( core::Size r1, core::Size r2 );
	int read_contact_map( core::Size r1, core::Size r2 );
	bool is_computed() const { return computed; }

public: // Public Methods
	
	// 
	void set_pose( core::pose::Pose const & pose );

	// Memory manipulation
	void free_memory();

	// Printing results
	void print_distance_map();
	void print_contact_map();
	void print_distance_map_tofile( std::string filename );
	void print_contact_map_tofile( std::string filename );
	void print_distance_map_tofile_compact( std::string filename );
	void print_contact_map_tofile_compact( std::string filename );
	void print_contact_map_tofile_compact2( std::string filename );

	// Calculating distances
	int hamming_contact_map( MGFContactMap & );
	double euclidean_distance_map( MGFContactMap & );


private: // Fields
	
	bool computed;

	// Number of residues (map size)
	core::Size nres;

	// Distance map
	double **distance_map;

	// Contact map
  	int **contact_map;

	// Distance cutoff
	core::Real distance_cutoff;

}; // StandaloneClass


} // MGFUtils


#endif // INCLUDED_MGFUtils_MGFContactMap_HH
