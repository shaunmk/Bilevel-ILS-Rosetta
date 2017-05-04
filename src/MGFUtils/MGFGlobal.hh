/// @file   MGFUtils/MGFGlobal.hh
///
/// @brief
/// @author Mario Garza-Fabre

#ifndef INCLUDED_MGFUtils_MGFGlobal_HH
#define INCLUDED_MGFUtils_MGFGlobal_HH

// Project forward headers
// #include <MGFUtils/MGFGlobal.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>

// #include <core/scoring/rms_util.hh>
#include <MGFUtils/MGFContactMap.fwd.hh>

// External library headers
// C++ headers
#include <string>
#include <cstdlib> 
#include <utility/vector1.hh>
#include <utility/vector0.hh>

namespace MGFUtils {
typedef core::Real Real;
class MGFGlobal{


public: // Creation

	/// @brief Constructor
	// MGFGlobal( ){}

	/// @brief Destructor
	// ~MGFGlobal(){}

public: // Public Methods
	
	static std::string get_filename_base( std::string const stage, int abiter);
	static void write_pose_output( core::pose::Pose &, std::string const, int abiter );
	static void write_pose_pdb( core::pose::Pose &, std::string const, int abiter );
	static bool isH( char res );
	static void configure_loop_information( core::pose::Pose & p );
	static bool is_loop( core::Size pos );

	static void set_stage( int s ){ 

		stage = s; 

		if( s == 1 || s == 2 || s == 3 ){

			window_size = 9;
			total_frag = 25;

		}else if( s == 4 ){

			window_size = 3;
			total_frag = 200;

		}else{

			window_size = 0;
			total_frag = 0;
		}

	}


	static void initialize( core::pose::Pose &, core::pose::PoseOP const &, std::string const );
	static void load_configure_strategy( core::pose::Pose &, core::pose::PoseOP const &,std::string const );
//	static int random_integer( int min, int max );
//	static void shuffle(utility::vector0< int > & vector, int size);

	static void extend_pose_singleloopresidue( core::pose::Pose & );
	static void extend_pose_singleloop( core::pose::Pose & );
	static void extend_pose_loops( core::pose::Pose & );
	static void extend_pose( core::pose::Pose & );
	static void configure_window_loop_correspondence();
	
	static void configure_ss_chunks( core::pose::Pose & );
	static int get_ss_chunk( int pos );

	static bool accept_window( core::Size window );
	static bool accept_residue_change( core::Size pos );

	static void clear_candidate_changes( );
	static void add_candidate_change( core::Size pos );
	static void reset_applied_changes( );
	static void apply_candidate_changes( );
	static bool stop_stage1( );

	// Other previously used
	// static int get_stage( ){ return stage; }
	// static void report_frame_fragment( core::Size const frame, core::Size const frag );
	// static void abego( core::pose::Pose & p );
	static void add_to_features_file( core::pose::Pose &, std::string const );
	// static void add_to_features_file_fullatom( core::pose::Pose &, std::string const );

	static void quicksortWID( utility::vector0< double > & value, long int inf, long int sup, int order );
	static void quicksort( utility::vector0< double > & value, utility::vector0< int > & ids, long int inf, long int sup, int order);

	static double calculate_rg_score( core::pose::Pose const & pose );
	void print_current_archive();

	// SMK
	template <class T>
	Real const median ( std::vector < T > & values );

public: 

	// Problem instance
	static core::pose::PoseOP native_pose;	// Native structure
	static std::string protein_sequence;	// Input protein sequence

	// General information
	static int stage;												// Stage number
	static int window_size;											// Size of window/fragment
	static int total_frag;											// Total number of fragments in the library

	static int total_residues;										// Total number of residues in the sequence
	static int total_windows_9;										// Total number of 9-length windows
	static int total_windows_3;										// Total number of 3-length windows

	static utility::vector1< char > ss;								// Secondary structure of each residue
	static utility::vector1< bool > isloop;							// Binary array indicating whether each residue is in a loop region or not
	static int total_loop_residues;									// Total number of residues in the sequence
	static int total_loop_regions;									// Total number of loop regions
	static utility::vector1< int > loop_start;						// Array containing the ids of residues representing each chunk of secondary structure
	static utility::vector1< int > loop_end;						// Array containing the ids of residues representing each chunk of secondary structure
	
	static utility::vector1< core::Real > window_loop_9;			// Percentage of residues in each window that are believed to lie in a loop region (9-length windows)
	static utility::vector1< core::Real > window_loop_9_normalized;	// Percentage of residues in each window that are believed to lie in a loop region (9-length windows)
	static utility::vector1< core::Real > window_loop_3;			// Percentage of residues in each window that are believed to lie in a loop region (3-length windows)
	static utility::vector1< core::Real > window_loop_3_normalized;	// Percentage of residues in each window that are believed to lie in a loop region (3-length windows)

	static utility::vector1< int > ss_chunk_centre;							// Array containing the ids of residues representing each chunk of secondary structure
	static utility::vector1< int > ss_chunk_start;							// Array containing the ids of residues representing each chunk of secondary structure
	static utility::vector1< int > ss_chunk_end;							// Array containing the ids of residues representing each chunk of secondary structure
	static utility::vector1< int > ss_chunk_residue;							// Array containing the ids of residues representing each chunk of secondary structure
	static int total_ss_chunks;									// Total number of sec. struct. chunks in the chain


	// Search strategy information/parameters
	static std::string strategy_input;								// Given input
	static std::string strategy;									// Name of strategy
	static int ABINITIO_ITERATIONS;									// Number of iterations of abinitio procedure
	static int ABITER;												// Current iteration of abinitio procedure
	static utility::vector1< core::Size > candidate_changes; 		// List of residue ids that were affected by the candidate move
	static utility::vector1< bool > applied_changes; 				// List of residue ids that were affected by the candidate move
	static int total_loops_changed;									// Current iteration of abinitio procedure
	static int total_residues_changed;								// Current iteration of abinitio procedure
	static core::Real prob_change_nonloopres;					// Probability of changing non-loop residues when applying a fragment insertion
	static bool to_change_nonloopres;							// Flag to indicate if the selected insertion will be allowed to affect non-loop residues, based on prob_change_nonloopres

	static bool terminate_bilevel_by_trial_count;	//SMK this is used to terminate the Bilevel optimiser once a certain number of "valid" L attempts have been made.

}; // MGFGlobal




} // MGFUtils


#endif // INCLUDED_MGFUtils_MGFGlobal_HH
