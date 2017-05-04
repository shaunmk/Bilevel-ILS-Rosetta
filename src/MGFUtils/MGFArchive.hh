/// @file   MGFUtils/MGFArchive.hh
///
/// @brief
/// @author Mario Garza-Fabre

#ifndef INCLUDED_MGFUtils_MGFArchive_HH
#define INCLUDED_MGFUtils_MGFArchive_HH

// Project forward headers
#include <MGFUtils/MGFArchive.fwd.hh>

#include <MGFUtils/MGFGlobal.hh>
#include <MGFUtils/MGFProtocol/Solution.hh>
// #include <MGFUtils/MGFProtocol/Solution.fwd.hh>
#include <MGFUtils/MGFProtocol/Evaluator.hh>
// #include <MGFUtils/MGFProtocol/Evaluator.fwd.hh>
#include <numeric/random/random.hh>


// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// External library headers
// C++ headers
#include <iostream>
// #include <string>
// #include <cstdlib> 
#include <utility/vector0.hh>

namespace MGFUtils {

using std::string;
using core::pose::PoseOP;
using core::pose::Pose;
using MGFUtils::MGFProtocol::SolutionOP;
using MGFUtils::MGFProtocol::Solution;
using MGFUtils::MGFProtocol::EvaluatorOP;
using MGFUtils::MGFProtocol::Evaluator;

/// @brief
class MGFArchive : public utility::pointer::ReferenceCount{

typedef utility::vector0< SolutionOP >::iterator SolutionIterator;

public: // Creation

	/// @brief Constructor 1
  	MGFArchive( int, int, string, string tag="" );

	/// @brief Constructor 2
	MGFArchive( int, string, string tag="" );

	/// @brief Destructor
	~MGFArchive(){}

private: // Private Methods

	void initialize( int, int, string );
	void configure( );
	void add( SolutionOP s );

	// --------------------------------
	// Archive reductors (private methods)
	// --------------------------------

	void reduce_archive_size_ENERGY( int );	

	void reduce_archive_size_ENERGY_DA( int );	
	bool is_duplicate( SolutionOP const &, utility::vector0< SolutionOP > & );
	
	void reduce_archive_size_DEGRADATION( int );	
	void degrade_unselected_neighbours( SolutionOP &, utility::vector0< SolutionOP > &, utility::vector0< bool > &, utility::vector0< int > & );

	void reduce_archive_size_DIVERSITY( int );

	void reduce_archive_size_COMMA( int );

	void reduce_archive_size_SR( int );
	void reduce_archive_size_SR_CM( int );
	// SMK are these archivers better than one that chooses one structure based on energy, and all others at random?
	void reduce_archive_size_ELITIST_RANDOM( int );

	utility::vector0< int > get_CM( Pose const &);
	int hamming_distance( utility::vector0< int > &, utility::vector0< int > & );

	void reduce_archive_size_CLUSTERING( int );

	// SMK output some metrics after each reduction step
	void assess_solutions_and_output();

	// --------------------------------

public: // Public Methods

	int 	get_current_size(){ /*std::cout << "get_current_size() was called! size is " << current_size << std::endl; */return current_size; }
	int 	get_desired_size(){ return desired_size; }
	void 	set_desired_size( int size ){ desired_size = size; }
	int 	get_max_size(){ return max_size; }
	void 	set_max_size( int size ){ 		
		// max_size needs to be greater than or equal to desired_size
		max_size = ( size >= desired_size ) ? size : desired_size; 
	}

	void set_evaluator( string const &, bool const );

	void add_pose( Pose const & );
	void add_extended_pose( );
	void add_existing_solution( SolutionOP s );

	//SMK to add a pose to a given position in the archive. Must do bounds checking
	void add_pose_in_position(Pose const & pose, int pos);

	// void add_solution( Solution const & );
	PoseOP get_member_pose( int );
	SolutionOP get_member( int );

	void reduce_archive();
	void reduce_archive( int );

	void dump_scored_archive_pdbs(core::scoring::ScoreFunction const & scorefxn, std::string info="");
	void replace_evaluator(string stage_tag);
	void print_current_archive();

	// --------------------------------
	// Archive reductors (public methods)
	// --------------------------------

	void set_DUPLICATE_RMSD_CUTOFF( double val ){ DUPLICATE_RMSD_CUTOFF = val; }
	double get_DUPLICATE_RMSD_CUTOFF( ){ return DUPLICATE_RMSD_CUTOFF; }

	void set_NICHE_RADIUS( double val ){ NICHE_RADIUS = val; }
	double get_NICHE_RADIUS( ){ return NICHE_RADIUS; }

	void set_SR_PROB( double val ){ SR_PROB = val; }
	double get_SR_PROB( ){ return SR_PROB; }

	const string& get_tag() const {
		return tag_;
	}

	core::Size getReductorCalls() const {
		return n_reductor_calls;
	}

	void setReductorCalls(core::Size reductorCalls) {
		n_reductor_calls = reductorCalls;
	}

	// SMK
	// we have the size of the retained set following the last reduction step.
	// if this equals 10 (the number of solutions finally output),
	// 	the last step corresponded either to a reduction triggered by a change in scorefxn or by the end of stage3.
	// We want to know what fraction of the solutions retained using the previous scorefxn are retained in the current step.
	// I also want to know this for all other reductions.
	void assess_frac_retained_from_previous_reduction(int size_after_current_reduction, const utility::vector0<int>& indices);

	// --------------------------------

private: // Private attributes

	core::Size n_reductor_calls;

	int desired_size; 	// Size of bounded archive
	int max_size;		// To allow the archive to grow beyond the desired size (this allows to reduce/process the archive less frequenlty);
	int current_size;	// Current size of the archive (this is a variable size archive)

	utility::vector0< SolutionOP > archive;		// Main storage
	utility::vector0< SolutionOP > auxiliar;	// Auxiliar, temporal storage to be used during reduction step

	string tag_;				// A tag which appears on output files. Set in ctors.
	string archive_type;		// Type of archive to be used
	EvaluatorOP evaluator;		// Evaluator to be use to assess the quality/features of solutions in the archive
	void (MGFUtils::MGFArchive::*reductor)( int );	// Pointer to archive reduction function

	int size_after_last_reduction; // SMK the size of the archive following the reduction step just prior to the current one.

	// --------------------------------
	// Parameters used for particular archive reductors
	// --------------------------------
	double DUPLICATE_RMSD_CUTOFF;		// Maximum RMSD value for considering two structures to be duplicates (ENERGY_DA reductor)
	double NICHE_RADIUS;				// Maximum RMSD value for considering two structures to be neighbours, members of the same niche (DEGRADATION reductor)
	double SR_PROB;						// Probability for stochastic ranking-based procedure
//	double RG_CUTOFF;					// Cutoff of acceptable RG values
	// --------------------------------
public: // Public attributes


}; // MGFArchive


} // MGFUtils


#endif // INCLUDED_MGFUtils_MGFArchive_HH
