/// @file   MGFUtils/MGFGlobal.cc
///
/// @brief
/// @author Mario Garza-Fabre

#include <MGFUtils/MGFGlobal.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/types.hh>
#include <numeric/random/random.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> 
#include <basic/options/keys/run.OptionKeys.gen.hh> 
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

// #include <protocols/abinitio/ResolutionSwitcher.hh>
// #include <protocols/moves/Mover.hh>

#include <numeric/xyzVector.hh>
#include <core/scoring/rms_util.hh>
// #include <core/scoring/rms_util.tmpl.hh>
// #include <core/scoring/TMscore.hh>

#include <core/util/ABEGOManager.hh>
#include <core/util/ABEGOManager.fwd.hh>

#include <MGFUtils/MGFContactMap.hh>

#include <iostream>
#include <fstream>

namespace MGFUtils {

using core::Size;
using std::string;
using utility::vector1;
using core::Vector;

/////////////////////////////////////////////////////////////////////////////
// Initialization/memory allocation for static variables

core::pose::PoseOP MGFGlobal::native_pose=NULL;
std::string MGFGlobal::protein_sequence="";
int MGFGlobal::stage = 0;
int MGFGlobal::window_size = 0;
int MGFGlobal::total_frag = 0;
// bool MGFGlobal::ss_computed = false;
utility::vector1< char > MGFGlobal::ss(0);
utility::vector1< bool > MGFGlobal::isloop(0);
std::string MGFGlobal::strategy_input = "";	
std::string MGFGlobal::strategy = "";	
int MGFGlobal::ABINITIO_ITERATIONS = 0;
int MGFGlobal::ABITER = 0;
int MGFGlobal::total_loop_residues = 0;
int MGFGlobal::total_loop_regions = 0;
utility::vector1< int > MGFGlobal::loop_start(0);	
utility::vector1< int > MGFGlobal::loop_end(0);	
int MGFGlobal::total_residues = 0;
int MGFGlobal::total_windows_9 = 0;
int MGFGlobal::total_windows_3 = 0;
utility::vector1< core::Real > MGFGlobal::window_loop_9(0);
utility::vector1< core::Real > MGFGlobal::window_loop_9_normalized(0);
utility::vector1< core::Real > MGFGlobal::window_loop_3(0);
utility::vector1< core::Real > MGFGlobal::window_loop_3_normalized(0);
utility::vector1< core::Size > MGFGlobal::candidate_changes(0);
utility::vector1< bool > MGFGlobal::applied_changes(0);
int MGFGlobal::total_loops_changed = 0;
int MGFGlobal::total_residues_changed = 0;
core::Real MGFGlobal::prob_change_nonloopres = 1.0;
bool MGFGlobal::to_change_nonloopres = false;
utility::vector1< int > MGFGlobal::ss_chunk_centre(0);	
utility::vector1< int > MGFGlobal::ss_chunk_start(0);	
utility::vector1< int > MGFGlobal::ss_chunk_end(0);	
utility::vector1< int > MGFGlobal::ss_chunk_residue(0);	
int MGFGlobal::total_ss_chunks = 0;
bool MGFGlobal::terminate_bilevel_by_trial_count = false; // SMK

/////////////////////////////////////////////////////////////////////////////
void MGFGlobal::initialize( core::pose::Pose & p, core::pose::PoseOP const & native, std::string const seq ){

	native_pose = native;
	protein_sequence = seq;

	total_residues = p.total_residue();
	total_windows_9 = total_residues - 8;
	total_windows_3 = total_residues - 2;

	configure_loop_information( p );
	configure_window_loop_correspondence();
	configure_ss_chunks( p );
}
/////////////////////////////////////////////////////////////////////////////
// Reads specified strategy in input parameters and configures settings
/*
void MGFGlobal::load_configure_strategy( core::pose::Pose & p, core::pose::PoseOP const & native, std::string const seq ){

	// Initialize
	initialize( p, native, seq );

	// ------------------------------------------------
	// Save input string
	// ------------------------------------------------
	strategy_input = basic::options::option[ basic::options::OptionKeys::abinitio::MGF_strategy ];

	// ------------------------------------------------
	// Parse strategy name and parameters
	// ------------------------------------------------
	std::stringstream sstream( strategy_input );	

	// Read strategy identifier (first given string/word)	
	std::getline(sstream, strategy, '*');
	// std::cout << "Configuring strategy: " << strategy << std::endl;	

	// Read parameters and configure
	if( strategy == "ROSETTA" ){

		// No more parameters to be read, just load default settings
		ABINITIO_ITERATIONS = 1;

	}else if( strategy.substr(0,12) == "ROSETTA_ILSA" || strategy.substr(0,12) == "ROSETTA_ILSB" || strategy.substr(0,12) == "ROSETTA_ILSC" ){

		// std::cout << sstream.str().empty() << std::endl;	
		// std::cout << sstream.eof() << std::endl;	
		// std::cout << sstream.rdbuf()->in_avail() << std::endl;	

		if( !(sstream >> ABINITIO_ITERATIONS) ){
			std::cout << "\n\tERROR: Iterations number required for strategy " << strategy << std::endl << std::endl;
			exit(1);
		}	
		std::string wg; std::getline(sstream, wg, '*'); // just to deal with the inserted *
		
		if( strategy.substr(0,12) == "ROSETTA_ILSB" || strategy.substr(0,12) == "ROSETTA_ILSC" ){
			
			if( !(sstream >> prob_change_nonloopres) ){
				std::cout << "\n\tERROR: Unable to read prob_change_nonloopres from input parameters " << strategy << std::endl << std::endl;
				exit(1);
			}	
			// std::string wg; std::getline(sstream, wg, '*'); // just to deal with the inserted *
		}

	}else{

		std::cout << "\n\tERROR: Unrecognised search strategy: " << strategy << std::endl << std::endl;
		exit(1);

	}

	//~ std::cout << "Configuration succesful" << std::endl;		
	//~ std::cout << "\tABINITIO_ITERS: " << ABINITIO_ITERATIONS << std::endl;		
	//~ std::cout << "\tprob_change_nonloopres: " << prob_change_nonloopres << std::endl;	
	//~ exit(1);
}
*/
/////////////////////////////////////////////////////////////////////////////
void MGFGlobal::write_pose_pdb( core::pose::Pose & p, std::string const stage, int abiter ){

	std::string fn = get_filename_base( stage, abiter ) + ".pdb" ;
	p.dump_pdb( fn.c_str() );

}
/////////////////////////////////////////////////////////////////////////////
/*
void MGFGlobal::write_pose_output( core::pose::Pose & p, std::string const stage, int abiter ){

	std::ofstream file;	

	// Generate filename base
	std::string filename_base = get_filename_base( stage, abiter );

	// Write energies file
	// std::cout << "\nEnergies:\n" << std::endl;	
	std::string fn = filename_base + ".energies";
	file.open( fn.c_str() );
	core::scoring::ScoreFunctionOP sf; 
	sf = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
	(*sf)( p ); 	
	file << "TOTAL\t" << p.energies().total_energy() << std::endl;
	//~ file << "ENV\t" << std::scientific<< (*sf).score_by_scoretype(p, core::scoring::env, true) << std::endl;
	file << "ENV\t" << (*sf).score_by_scoretype(p, core::scoring::env, true) << std::endl;
	file << "PAIR\t" << (*sf).score_by_scoretype(p, core::scoring::pair, true) << std::endl;
	file << "CBETA\t" << (*sf).score_by_scoretype(p, core::scoring::cbeta, true) << std::endl;
	file << "VDW\t" << (*sf).score_by_scoretype(p, core::scoring::vdw, true) << "\t" << std::endl;
	file << "RG\t" << (*sf).score_by_scoretype(p, core::scoring::rg, true) << "\t" << std::endl;
	file << "CENPACK\t" << (*sf).score_by_scoretype(p, core::scoring::cenpack, true) << std::endl;
	file << "HSPAIR\t" << (*sf).score_by_scoretype(p, core::scoring::hs_pair, true) << std::endl;
	file << "SSPAIR\t" << (*sf).score_by_scoretype(p, core::scoring::ss_pair, true) << std::endl;
	file << "RSIGMA\t" << (*sf).score_by_scoretype(p, core::scoring::rsigma, true) << std::endl;
	file << "SHEET\t" << (*sf).score_by_scoretype(p, core::scoring::sheet, true) << std::endl;
	// Applied Weights
	// file << "\n*** Weights" << std::endl;
	// file << "ENV\t" << (*sf).get_weight(core::scoring::env) << std::endl;
	// file << "PAIR\t" << (*sf).get_weight(core::scoring::pair) << std::endl;
	// file << "CBETA\t" << (*sf).get_weight(core::scoring::cbeta) << std::endl;
	// file << "VDW\t" << (*sf).get_weight(core::scoring::vdw) << std::endl;
	// file << "RG\t" << (*sf).get_weight(core::scoring::rg) << std::endl;
	// file << "CENPACK\t" << (*sf).get_weight(core::scoring::cenpack) << std::endl;
	// file << "HSPAIR\t" << (*sf).get_weight(core::scoring::hs_pair) << std::endl;
	// file << "SSPAIR\t" << (*sf).get_weight(core::scoring::ss_pair) << std::endl;
	// file << "RSIGMA\t" << (*sf).get_weight(core::scoring::rsigma) << std::endl;
	// file << "SHEET\t" << (*sf).get_weight(core::scoring::sheet) << std::endl;
	file.close();

	// Write PDB file
	// std::cout << "\nPDB:\n" << std::endl;	
	fn = filename_base + ".pdb" ;
	p.dump_pdb( fn.c_str() );

	// Write Torsion angles file
	fn = filename_base + ".torsion";
	file.open( fn.c_str() );
	for ( Size pos = 1; pos <= p.total_residue(); pos++ ) {
		if ( ! p.residue(pos).is_protein() ) continue;
		file << p.phi( pos ) << "\t" <<  p.psi( pos ) << "\t" << p.omega( pos ) << "\t" << p.pdbid( pos ) << "\t" << p.res( pos ) << std::endl;
	}
	file.close();


	// Write Contact map and Distance map files
	MGFUtils::MGFContactMap C( p );
	// std::cout << "\nDistance map:\n" << std::endl;
	// C.print_distance_map();
	//C.print_distance_map_tofile_compact( filename_base + ".dm" );
	// std::cout << "\nContact map:\n" << std::endl;
	// C.print_contact_map();
	// C.print_contact_map_tofile_compact( filename_base + ".cm" );
	C.print_contact_map_tofile_compact2( filename_base + ".cm" );
	C.free_memory();
	
	// Write ss chunk distance file
	fn = filename_base + ".ss_chunk_dist";
	file.open( fn.c_str() );
	int ssc1, ssc2;
	numeric::xyzVector<core::Real> coord_r1, coord_r2;
	for( ssc1=1; ssc1 < total_ss_chunks; ssc1++ ){
		for( ssc2=ssc1+1; ssc2 <= total_ss_chunks; ssc2++ ){
			
			coord_r1 = p.residue( ss_chunk_centre[ ssc1 ] ).atom( "CA" ).xyz();
			coord_r2 = p.residue( ss_chunk_centre[ ssc2 ] ).atom( "CA" ).xyz();
			file << coord_r1.distance(coord_r2) << std::endl;
			
		}
	}	
	file.close();

}
*/
/////////////////////////////////////////////////////////////////////////////
std::string MGFGlobal::get_filename_base( std::string const stage, int abiter ){

	// Read parameters and define path and filename base
	std::stringstream convert;
	int seed = basic::options::option[basic::options::OptionKeys::run::jran].value();
	convert << ((seed - 1) * ABINITIO_ITERATIONS + abiter);
	
	std::string seed_used = convert.str();
	std::string out_path = basic::options::option[basic::options::OptionKeys::out::path::all].value();

	// std::cout << "Out_path: " << out_path << std::endl;
	// std::cout << "rfind: " << out_path.rfind("/run_") << std::endl;
	// std::cout << "substr: " << out_path.substr(0, out_path.rfind("/run_")) << std::endl;

	std::string sub_path = out_path.substr(0, out_path.rfind("/run_"));
	if( sub_path != out_path ){
		out_path = sub_path + "/run_" + seed_used + "/";
	}

	// std::cout << "Out_path: " << out_path << std::endl;

	std::string filename_base = out_path + stage + "_" + seed_used;

	return filename_base;
}
/////////////////////////////////////////////////////////////////////////////
// Decide if alterations are allowed at a given position of the chain,
// based on the current search strategy
bool MGFGlobal::accept_residue_change( core::Size pos ){

	if( strategy.substr(0,12) == "ROSETTA_ILSA" && ABITER > 1 && stage == 1 ){

		// In stage 1, for ABITER>1, ROSETTA_ILSA only alter loop-residues
		// if( !is_loop( pos ) ) return false; 
		return is_loop( pos ); 
		
	}else if( strategy.substr(0,12) == "ROSETTA_ILSB" && ABITER > 1 && stage < 4 ){
		
		// If stage 1, ROSETTA_ILSB only alter loop-residues
		// If stage 2 or 3, ROSETTA_ILSB accepts changes in loop residues, and non-loop residues based on decision previously made ( to_change_nonloopres )
		return ( stage == 1) ? is_loop( pos ) : is_loop( pos ) || to_change_nonloopres;
		
	}else if( strategy.substr(0,12) == "ROSETTA_ILSC" && ABITER > 1 && stage < 4 ){
		
		// If stage 1, ROSETTA_ILSC only alter loop-residues
		// If stage 2 or 3, ROSETTA_ILSC accepts changes in loop residues, and non-loop residues based on decision previously made ( to_change_nonloopres )
		return ( stage == 1) ? is_loop( pos ) : is_loop( pos ) || to_change_nonloopres;
		
	}

	return true;
}
/////////////////////////////////////////////////////////////////////////////
// Decide if insertions in a given window are acceptable based on
// the search strategy currently being applied
/*
bool MGFGlobal::accept_window( core::Size window ){

	if( strategy.substr(0,12) == "ROSETTA_ILSA" && ABITER > 1 && stage == 1 ){

		// In stage 1, for ABITER>1, ROSETTA_ILSA only applies insertions in windows where at least 1 residue is believed to be loop
		if( window_loop_9[ window ] < 0.01 ) return false; 
		// Using < 0.01 to avoid using == 0.0 (to prevent problems with precision). 
		// Value for 1 loop-residue in the window is 0.1111 (all values below this are equivalent to 0.0 for these purposes)

		// If all loop residues in this window were already changed, then reject the window
		bool dec = false; // Assume window is rejected
		for(int i=0; i<9 && !dec; i++){ 
			if( is_loop( window+i ) && !applied_changes[ window+i ] ) dec = true;
		}
		// std::cout << "All already changed!!! ( " << total_loop_residues << " vs " << total_loops_changed << ")\n";
		if( !dec ) return false; // No unchanged loops were found
		
		// AT THIS POINT THE WINDOW IS ACCEPTED
		
	}else if( strategy.substr(0,12) == "ROSETTA_ILSB" && ABITER > 1 && stage < 4 ){
		
		// ROSETTA_ILSB only applies insertions in windows where at least 1 residue is believed to be loop (Stages 1,2,3)
		if( window_loop_9[ window ] < 0.01 ) return false;
		
		// In stage 1, reject the window if all loop residues were already changed
		if( stage == 1){
			bool dec = false; // Assume window is rejected
			for(int i=0; i<9 && !dec; i++){ 
				if( is_loop( window+i ) && !applied_changes[ window+i ] ) dec = true;
			}
			if( !dec ) return false; // No unchanged loops were found
		}		
		
		// AT THIS POINT THE WINDOW IS ACCEPTED
		
		// Decide, based on the given probability, if changes at non-loop positions are to be accepted
		to_change_nonloopres = (  prob_change_nonloopres < 1.0 && numeric::random::rg().uniform() >= prob_change_nonloopres ) ? false : true;
		
	}else if( strategy.substr(0,12) == "ROSETTA_ILSC" && ABITER > 1 && stage < 4 ){
		
		// ROSETTA_ILSC applies insertions at a given window on the basis of a probability that corresponds to the proportion of loop residues in the window
		if( window_loop_9_normalized[ window ] < 1.0 && numeric::random::rg().uniform() >= window_loop_9_normalized[ window ] ) return false;
		
		// In stage 1, reject the window if all loop residues were already changed
		if( stage == 1){
			bool dec = false; // Assume window is rejected
			for(int i=0; i<9 && !dec; i++){ 
				if( is_loop( window+i ) && !applied_changes[ window+i ] ) dec = true;
			}
			if( !dec ) return false; // No unchanged loops were found
		}	
		
		// AT THIS POINT THE WINDOW IS ACCEPTED
		
		// Decide, based on the given probability, if changes at non-loop positions are to be accepted
		to_change_nonloopres = (  prob_change_nonloopres < 1.0 && numeric::random::rg().uniform() >= prob_change_nonloopres ) ? false : true;
	}

	return true;
}
*/
/////////////////////////////////////////////////////////////////////////////
int MGFGlobal::get_ss_chunk( int pos ){

	return ss_chunk_residue[ pos ];
}
/////////////////////////////////////////////////////////////////////////////
void MGFGlobal::configure_ss_chunks( core::pose::Pose & p ){
	
	// Determine total number of ss chunks
	total_ss_chunks = 0;
	ss_chunk_residue.assign( total_residues, 0 );
	core::Size start, end;
	for(int i=1; i <= total_residues; i++){

		if( ss[i] != 'L' ){

			// A new ss chunk was identified
			total_ss_chunks++;
		
			// Start of the ss chunk identified
			// start = i++;
			start = i;
			
			// Identify size of chunk
			for( ; i <= total_residues && ss[i] == ss[start]; i++ ) ss_chunk_residue[i] = total_ss_chunks;
				
			// End of the ss chunk identified
			end = --i;
			
			// Allocate memory for storing chunk info			
			ss_chunk_centre.resize( total_ss_chunks );
			ss_chunk_start.resize( total_ss_chunks );
			ss_chunk_end.resize( total_ss_chunks );

			// Compute and store representative residue info 
			ss_chunk_centre[ total_ss_chunks ] = start + (int)(  (end-start)/2.0 );
			ss_chunk_start[ total_ss_chunks ] = start;
			ss_chunk_end[ total_ss_chunks ] = end;
			
		}
	}
	
	/////////////////////////////////////////////////////////////////////////////
	// Showing information of SS chunks
	/////////////////////////////////////////////////////////////////////////////

	/*
	std::cout << "Total number of chunks: " << total_ss_chunks << "\n";
	std::cout << "Chunk starts: " << ss_chunk_start << std::endl;
	std::cout << "Chunk ends:: " << ss_chunk_end << std::endl;
	std::cout << "Chunk centres: " << ss_chunk_centre << std::endl << std::endl;

	// Distances between elements from different ss chunks
	int ssc1, ssc2;
	numeric::xyzVector<core::Real> coord_r1, coord_r2;
	for( ssc1=1; ssc1 < total_ss_chunks; ssc1++ ){
		for( ssc2=ssc1+1; ssc2 <= total_ss_chunks; ssc2++ ){
			
			std::cout << "Distances between SS chunks: " << ssc1 << "[ " << ss_chunk_start[ ssc1 ] << " - " << ss_chunk_end[ ssc1 ] << " ] , " << ssc2 << "[ " << ss_chunk_start[ ssc2 ] << " - " << ss_chunk_end[ ssc2 ] << " ] " << std::endl;

			core::Real min=100000, max=0;
			int min1,min2, max1, max2;

			for(int i = ss_chunk_start[ ssc1 ]; i <= ss_chunk_end[ ssc1 ]; i++){
				for(int j = ss_chunk_start[ ssc2 ]; j <= ss_chunk_end[ ssc2 ]; j++){
					coord_r1 = p.residue( i ).atom( "CA" ).xyz();
					coord_r2 = p.residue( j ).atom( "CA" ).xyz();

					core::Real distance = coord_r1.distance(coord_r2);
					if( distance > max ) { max = distance; max1=i; max2=j; }
					if( distance < min ) { min = distance; min1=i; min2=j; }
				}
			}			

			std::cout << "\tMin: d(" << min1 << " , " << min2 << ") = " << min << "\n";
			std::cout << "\tMax: d(" << max1 << " , " << max2 << ") = " << max << "\n";

		}
	}
	*/
}	
/////////////////////////////////////////////////////////////////////////////
void MGFGlobal::configure_window_loop_correspondence(){

	core::Real max=0.0;

	// 9-length windows
	window_loop_9.assign( total_windows_9, 0.0 );
	window_loop_9_normalized.assign( total_windows_9, 0.0 );
	for(int wpos=1; wpos<=total_windows_9; wpos++){
		for(int i=0; i<9; i++){ 
			if( is_loop( wpos+i ) )	window_loop_9[ wpos ] += 1.0;
		}
		window_loop_9[ wpos ] /= 9.0;

		if( window_loop_9[ wpos ] > max ) max = window_loop_9[ wpos ];
	}

	// std::cout << "window_loop_9:" << std::endl;
	// std::cout << window_loop_9 << std::endl;

	// Normalize to [0,1] range
	for(int wpos=1; wpos<=total_windows_9; wpos++) window_loop_9_normalized[ wpos ] = window_loop_9[ wpos ]/max;

	// std::cout << "window_loop_9 (normalized):" << std::endl;
	// std::cout << window_loop_9_normalized << std::endl;

	// 3-length windows 
	window_loop_3.assign( total_windows_3, 0.0 );
	window_loop_3_normalized.assign( total_windows_3, 0.0 );
	max=0.0;
	for(int wpos=1; wpos<=total_windows_3; wpos++){
		for(int i=0; i<3; i++){ 
			if( is_loop( wpos+i ) )	window_loop_3[ wpos ] += 1.0;
		}
		window_loop_3[ wpos ] /= 3.0;

		if( window_loop_3[ wpos ] > max ) max = window_loop_3[ wpos ];
	}
	// std::cout << "window_loop_3:" << std::endl;
	// std::cout << window_loop_3 << std::endl;

		// Normalize to [0,1] range
	for(int wpos=1; wpos<=total_windows_3; wpos++) window_loop_3_normalized[ wpos ] = window_loop_3[ wpos ]/max;

	// std::cout << "window_loop_3 (normalized):" << std::endl;
	// std::cout << window_loop_3_normalized << std::endl;

}
/////////////////////////////////////////////////////////////////////////////
// Retrieve information of secondary structure from psipred file
void MGFGlobal::configure_loop_information( core::pose::Pose & p ){

	total_loop_residues = 0;
	total_loop_regions = 0;
	// int cont=0;

	// Read secondary structure info
	ss = core::pose::read_psipred_ss2_file( p );	

	// Determine whether a residue is in a loop region (discard terminal loops)
	// VERSION 2: Save information about number of loop regions, as well as start and end point of each of them
	isloop.assign(p.total_residue(), false); // Assume not loop	
	int start, end;
	for(int i=1; i <= total_residues; i++){

		if( ss[i] == 'L' ){		
		
			// Start of the loop identified
			start = i++;
			
			// Identify size of loop
			for( ; i <= total_residues && ss[i] == 'L'; i++ );
				
			// End of the loop identified
			end = --i;
			
			// Discard terminal loops (left and right tails)
			if( start == 1 || end == total_residues ) continue;

			// A new loop was identified
			total_loop_regions++;

			// Save start and end points
			loop_start.resize( total_loop_regions );
			loop_end.resize( total_loop_regions );
			loop_start[ total_loop_regions ] = start;
			loop_end[ total_loop_regions ] = end;

			// Mark all residues in this region as loops
			for( int j=start; j<=end; j++){

				total_loop_residues++;
				isloop[ j ] = true;

			}

			// std::cout << total_loop_regions << ": [ " << start << " , " << end << " ]\n";
		}
	}

	/*// Determine whether a residue is in a loop region (discard terminal loops)	
	// VERSION 1: Working, but more information was required and a different process is implemented in VERSION 2
	isloop.assign(p.total_residue(), false); // Assume not loop
	for(core::Size i=1; i<=p.total_residue(); i++){

		if( ss[i] == 'L' ){

			// Discard terminal loops
			bool is_left_tail = true;
			for(core::Size j = i-1; j >= 1 && is_left_tail; j--){
				if( ss[j] != 'L' ) is_left_tail = false;
			}
			bool is_right_tail = true;
			for(core::Size j = i+1; j <= p.total_residue() && is_right_tail; j++){
				if( ss[j] != 'L' ) is_right_tail = false;
			}

			// If not left tail, and not right tail, then LOOP!!!
			isloop[ i ] = !is_left_tail && !is_right_tail;

			// Update total loops
			if( isloop[ i ] ) total_loop_residues++;

			// if( isloop[ i ] ){

			// 	std::cout << i << "\t" << total_loop_residues << "\t" << ++cont << " (phi) \t" << ++cont << " (psi) \t" << ++cont << " (omg)" << std::endl;			

			// }
		}
	}*/
}
/////////////////////////////////////////////////////////////////////////////
bool MGFGlobal::is_loop( core::Size pos ){	

	// Determine if the given residue is in a loop region
	return isloop[ pos ];
}
/////////////////////////////////////////////////////////////////////////////
/*int MGFGlobal::random_integer( int min, int max ){

	return ( static_cast< int >( numeric::random::rg().uniform() * (max-min+1) ) + 1 );

}*/
/////////////////////////////////////////////////////////////////////////////
// Changes torsion angles to an extended-pose configuration at LOOP REGIONS
/*void MGFGlobal::extend_pose_singleloopresidue( core::pose::Pose & p ){
	
	// Decide which loop and residue to extend at random
	int loop = random_integer( 1, total_loop_regions );
	int pos = random_integer( loop_start[loop], loop_end[loop] );

	p.set_phi( pos, -150 );
	p.set_psi( pos, 150);
	p.set_omega( pos, 180 );
}*/
/////////////////////////////////////////////////////////////////////////////
// Changes torsion angles to an extended-pose configuration at LOOP REGIONS
/*void MGFGlobal::extend_pose_singleloop( core::pose::Pose & p ){

	// Decide which loop to extend at random
	int loop = random_integer( 1, total_loop_regions );

	for ( int pos = loop_start[loop]; pos <= loop_end[loop]; pos++ ) {
			p.set_phi( pos, -150 );
			p.set_psi( pos, 150);
			p.set_omega( pos, 180 );
	}	
}*/
/////////////////////////////////////////////////////////////////////////////
// Changes torsion angles to an extended-pose configuration at LOOP REGIONS
void MGFGlobal::extend_pose_loops( core::pose::Pose & p ){
	for ( Size pos = 1; pos <= p.total_residue(); pos++ ) {
		if ( ! p.residue(pos).is_protein() ) continue;

		if( is_loop( pos ) ){			
			p.set_phi( pos, -150 );
			p.set_psi( pos, 150);
			p.set_omega( pos, 180 );
		}
	}	
}
/////////////////////////////////////////////////////////////////////////////
// Changes all torsion angles to generate a fully extended pose
void MGFGlobal::extend_pose( core::pose::Pose & p ){
	for ( Size pos = 1; pos <= p.total_residue(); pos++ ) {
		// if ( ! p.residue(pos).is_protein() ) continue;
		p.set_phi( pos, -150 );
		p.set_psi( pos, 150);
		p.set_omega( pos, 180 );
	}
}
/////////////////////////////////////////////////////////////////////////////
void MGFGlobal::clear_candidate_changes( ){

	// if( strategy.substr(0,12) == "ROSETTA_ILSA" && ABITER > 1 && stage == 1 ){
		candidate_changes.clear();
	// }
}
/////////////////////////////////////////////////////////////////////////////
void MGFGlobal::add_candidate_change( core::Size pos ){

	//~ if( strategy.substr(0,12) == "ROSETTA_ILSA" && ABITER > 1 && stage == 1 ){
	if( strategy.substr(0,11) == "ROSETTA_ILS" && ABITER > 1 && stage == 1 ){
		candidate_changes.push_back( pos );
	}
}
/////////////////////////////////////////////////////////////////////////////
void MGFGlobal::reset_applied_changes( ){

	candidate_changes.clear();
	total_loops_changed = 0;
	total_residues_changed = 0;
	applied_changes.assign(total_residues, false); // nothing has changed

}
/////////////////////////////////////////////////////////////////////////////
void MGFGlobal::apply_candidate_changes( ){

	//~ if( strategy.substr(0,12) == "ROSETTA_ILSA" && ABITER > 1 && stage == 1 ){
	if( strategy.substr(0,11) == "ROSETTA_ILS" && ABITER > 1 && stage == 1 ){

		for ( Size i = 1; i <= candidate_changes.size(); i++ ){
			Size pos = candidate_changes[i];

			if( !applied_changes[ pos ] ){

				total_residues_changed++;
				applied_changes[ pos ] = true;

				if( is_loop( pos ) ){
					total_loops_changed++;
				}
			}
		}
	}
}
/////////////////////////////////////////////////////////////////////////////
bool MGFGlobal::stop_stage1(){

	//~ if( strategy.substr(0,12) == "ROSETTA_ILSA" && ABITER > 1){
	if( strategy.substr(0,11) == "ROSETTA_ILS" && ABITER > 1){

		if( total_loop_residues <= total_loops_changed ){
			// std::cout << "Stopping stage1!!!!\n";
			return true;
		}
	}

	return false;
}
/////////////////////////////////////////////////////////////////////////////
bool MGFGlobal::isH( char res ){

	if(
		res == 'A' ||
		res == 'C' ||
		res == 'G' ||
		res == 'I' ||
		res == 'L' ||
		res == 'M' ||
		res == 'F' ||
		res == 'P' ||
		res == 'W' ||
		res == 'V' 
	)
		return true;

	return false;
}
// /////////////////////////////////////////////////////////////////////////////
// void MGFGlobal::add_to_features_file_fullatom( core::pose::Pose & p, std::string const stage ){

// 	std::ofstream file;		
// 	std::string filename_base = get_filename_base( stage, 1 );
// 	std::string fn;
	
// 	/////////////////////////////////////////////////////////////////////////////
// 	// ENERGIES (full-atom)
// 	/////////////////////////////////////////////////////////////////////////////

// 	core::scoring::ScoreFunctionOP sf; 
// 	sf = core::scoring::get_score_function();
// 	(*sf)( p );

// 	std::cout  << (*sf).score_by_scoretype(p, core::scoring::hbond_sr_bb, true) << "\n";
// 	std::cout  << (*sf).score_by_scoretype(p, core::scoring::hbond_lr_bb, true) << "\n";
// 	std::cout  << (*sf).score_by_scoretype(p, core::scoring::hbond_bb_sc, true) << "\n";
// 	std::cout  << (*sf).score_by_scoretype(p, core::scoring::hbond_sc, true) << "\n";
// 	std::cout  << (*sf).score_by_scoretype(p, core::scoring::pro_close, true) << "\n";
// 	std::cout  << (*sf).score_by_scoretype(p, core::scoring::dslf_fa13, true) << "\n";
// 	std::cout  << (*sf).score_by_scoretype(p, core::scoring::ref, true) << "\n";
// 	std::cout  << (*sf).score_by_scoretype(p, core::scoring::fa_atr, true) << "\n";
// 	std::cout  << (*sf).score_by_scoretype(p, core::scoring::fa_sol, true) << "\n";
// 	std::cout  << (*sf).score_by_scoretype(p, core::scoring::fa_elec, true) << "\n";
// 	std::cout  << (*sf).score_by_scoretype(p, core::scoring::fa_dun, true) << "\n";
// 	std::cout  << (*sf).score_by_scoretype(p, core::scoring::omega, true) << "\n";
// 	std::cout  << (*sf).score_by_scoretype(p, core::scoring::fa_rep, true) << "\n";
// 	std::cout  << (*sf).score_by_scoretype(p, core::scoring::p_aa_pp, true) << "\n";
// 	std::cout  << (*sf).score_by_scoretype(p, core::scoring::rama, true) << "\n";
// 	std::cout  << (*sf).score_by_scoretype(p, core::scoring::fa_intra_rep, true) << "\n";


// 	std::cout << "Total energy: " << p.energies().total_energy() << "\n\n";

// 	// std::cout << "ENERGIES\n\n";
// 	// p.energies().total_energies().print();

// 	// std::cout << "WEIGHTS\n\n";
// 	// (*sf).weights().print(); // Full list of weights


// 	// core::pose::Pose p_copy( p );
// 	// core::pose::Pose const centroid_pose ( p_copy );
// 	// protocols::abinitio::ResolutionSwitcher res_switch( centroid_pose, false, true, true );
// 	// res_switch.apply( p_copy );
// }
// /////////////////////////////////////////////////////////////////////////////
/*
void MGFGlobal::add_to_features_file( core::pose::Pose & p, std::string const stage ){

	// Generate filename base and create/open file
	std::ofstream file;		
	std::string filename_base = get_filename_base( stage, 1 );
	std::string fn;

	/////////////////////////////////////////////////////////////////////////////
	// ENERGIES (low-res)
	/////////////////////////////////////////////////////////////////////////////

	fn = filename_base + ".features.energies";
	// file.open( fn.c_str(), std::ofstream::out | std::ofstream::app );
	file.open( fn.c_str() );

	core::scoring::ScoreFunctionOP sf; 
	sf = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
	(*sf)( p ); 

	file << p.energies().total_energy() << "\t";
//	file << (*sf).score_by_scoretype(p, core::scoring::env, true) << "\t";
//	file << (*sf).score_by_scoretype(p, core::scoring::pair, true) << "\t";
//	file << (*sf).score_by_scoretype(p, core::scoring::cbeta, true) << "\t";
//	file << (*sf).score_by_scoretype(p, core::scoring::vdw, true) << "\t" << "\t";
//	file << (*sf).score_by_scoretype(p, core::scoring::rg, true) << "\t" << "\t";
//	file << (*sf).score_by_scoretype(p, core::scoring::cenpack, true) << "\t";
//	file << (*sf).score_by_scoretype(p, core::scoring::hs_pair, true) << "\t";
//	file << (*sf).score_by_scoretype(p, core::scoring::ss_pair, true) << "\t";
//	file << (*sf).score_by_scoretype(p, core::scoring::rsigma, true) << "\t";
//	file << (*sf).score_by_scoretype(p, core::scoring::sheet, true) << "\t";

	file << std::endl;
	file.close();

	/////////////////////////////////////////////////////////////////////////////
	// CONTACTS: Summary | All contacts | SS chunk contacts
	/////////////////////////////////////////////////////////////////////////////

	fn = filename_base + ".features.contacts_all";
	// file.open( fn.c_str(), std::ofstream::out | std::ofstream::app );
	file.open( fn.c_str() );

	fn = filename_base + ".features.contacts_sschunk";
	// file.open( fn.c_str(), std::ofstream::out | std::ofstream::app );
	std::ofstream file_ss;
	file_ss.open( fn.c_str() );

	core::Real distance_cutoff = 8.0;
	core::Size nres = p.total_residue();

	int TOTAL_CONTACTS = 0;
	int LOCAL_CONTACTS = 0;
	int NONLOCAL_CONTACTS = 0;
	int HH_CONTACTS = 0;
	int HP_CONTACTS = 0;
	int PP_CONTACTS = 0;

	// Consider each pair of residues
	core::Size r1;
	core::Size r2;
	for(r1=1; r1<nres; r1++){

		// Get coordinates of residue 1
		numeric::xyzVector<core::Real> coord_r1 = p.residue( r1 ).atom( "CA" ).xyz();

		for(r2=r1+1; r2<=nres; r2++){

			// Get coordinates of residue 2 
			numeric::xyzVector<core::Real> coord_r2 = p.residue( r2 ).atom( "CA" ).xyz();

			// Compute distance between residues
			core::Real distance = coord_r1.distance(coord_r2);

			int is_contact = 0;

			// r1 and r2 are in contact
			if( distance <= distance_cutoff ){

				is_contact = 1;

				// Total number of contacts
				TOTAL_CONTACTS++;

				// Local contacts and Non-local contacts
				if (r2 - r1 <= 5){
					LOCAL_CONTACTS++;
				}else{
					NONLOCAL_CONTACTS++;
				}

				// Hydrophobicity-based contacts
				bool r1_H = isH( p.residue_type( r1 ).name1() );
				bool r2_H = isH( p.residue_type( r2 ).name1() );

				if( r1_H && r2_H ){

					HH_CONTACTS++;

				}else if( !r1_H && !r2_H ){

					PP_CONTACTS++;
					
				}else{

					HP_CONTACTS++;
				}

			}

			// Store info for "all contacts" file
			file << is_contact << "\t";

			// Store info for "ss chunk contacts" file
			if( get_ss_chunk(r1)>0 && get_ss_chunk(r2)>0 && get_ss_chunk(r1) != get_ss_chunk(r2) ){
				file_ss << is_contact << "\t";
			}

		}
	}

	file << std::endl;
	file.close();
	file_ss << std::endl;
	file_ss.close();


	// Summary file
	fn = filename_base + ".features.contacts_summary";
	// file.open( fn.c_str(), std::ofstream::out | std::ofstream::app );
	file.open( fn.c_str() );

	file << TOTAL_CONTACTS << "\t";
	file << LOCAL_CONTACTS << "\t";
	file << NONLOCAL_CONTACTS << "\t";

	file << std::endl;
	file.close();

	/////////////////////////////////////////////////////////////////////////////
	// HP MODEL-BASED
	/////////////////////////////////////////////////////////////////////////////

	fn = filename_base + ".features.hpmodel";
	// file.open( fn.c_str(), std::ofstream::out | std::ofstream::app );
	file.open( fn.c_str() );

	file << HH_CONTACTS << "\t";
	file << PP_CONTACTS << "\t";
	file << HP_CONTACTS << "\t";

	// Radius of gyration of H and P residues
	double RgH=0.0, RgP=0.0;
	double meanHcoord[3], meanPcoord[3];
	int contH = 0;
	int contP = 0;

	// Compute mean coordinates of H and P residues	
	meanHcoord[0]=meanPcoord[0]=0.0;	
	meanHcoord[1]=meanPcoord[1]=0.0;	
	meanHcoord[2]=meanPcoord[2]=0.0;	

	for(r1=1; r1<=nres; r1++){

		bool H = isH( p.residue_type( r1 ).name1() );
		numeric::xyzVector<core::Real> coord = p.residue( r1 ).atom( "CA" ).xyz();

		if( H ){

			meanHcoord[0] += coord.x();
			meanHcoord[1] += coord.y();
			meanHcoord[2] += coord.z();
			contH++;

		}else{

			meanPcoord[0] += coord.x();
			meanPcoord[1] += coord.y();
			meanPcoord[2] += coord.z();
			contP++;
		}
	}

	meanHcoord[0] /= contH;
	meanHcoord[1] /= contH;
	meanHcoord[2] /= contH;

	meanPcoord[0] /= contP;
	meanPcoord[1] /= contP;
	meanPcoord[2] /= contP;

	// Compute RgH and RgP
	for(r1=1; r1<=nres; r1++){

		bool H = isH( p.residue_type( r1 ).name1() );
		numeric::xyzVector<core::Real> coord = p.residue( r1 ).atom( "CA" ).xyz();

		if( H ){

			RgH += pow((coord.x() - meanHcoord[0]), 2);	
			RgH += pow((coord.y() - meanHcoord[1]), 2);	
			RgH += pow((coord.z() - meanHcoord[2]), 2);	

		}else{

			RgP += pow((coord.x() - meanPcoord[0]), 2);	
			RgP += pow((coord.y() - meanPcoord[1]), 2);	
			RgP += pow((coord.z() - meanPcoord[2]), 2);	

		}
	}
	
	RgH = sqrt(RgH/(double)contH);
	RgP = sqrt(RgP/(double)contP);

	// Report difference if H residues are more exposed, 0 otherwise
	file << ((RgH > RgP)? RgH-RgP : 0.0 )<< "\t";

	file << std::endl;
	file.close();

	/////////////////////////////////////////////////////////////////////////////
	// DISTANCES
	/////////////////////////////////////////////////////////////////////////////
	
	// Distance between L and R terminus
	fn = filename_base + ".features.distances_terminus";
	// file.open( fn.c_str(), std::ofstream::out | std::ofstream::app );
	file.open( fn.c_str() );

	numeric::xyzVector<core::Real> coord_r1 = p.residue( 1 ).atom( "CA" ).xyz();
	numeric::xyzVector<core::Real> coord_r2 = p.residue( nres ).atom( "CA" ).xyz();
	file << coord_r1.distance(coord_r2) << "\t";

	file << std::endl;
	file.close();

	// Pairwise distances between ss chunk representatives (centre)
	fn = filename_base + ".features.distances_sschunk_centre";
	// file.open( fn.c_str(), std::ofstream::out | std::ofstream::app );
	file.open( fn.c_str() );

	int ssc1, ssc2;
	for( ssc1=1; ssc1 < total_ss_chunks; ssc1++ ){
		for( ssc2=ssc1+1; ssc2 <= total_ss_chunks; ssc2++ ){
			
			coord_r1 = p.residue( ss_chunk_centre[ ssc1 ] ).atom( "CA" ).xyz();
			coord_r2 = p.residue( ss_chunk_centre[ ssc2 ] ).atom( "CA" ).xyz();
			file << coord_r1.distance(coord_r2) << "\t";
			
		}
	}
	
	file << std::endl;
	file.close();


	// Pairwise distances between ss chunk representatives (extreme points)
	fn = filename_base + ".features.distances_sschunk_extreme";
	// file.open( fn.c_str(), std::ofstream::out | std::ofstream::app );
	file.open( fn.c_str() );

	// int ssc1, ssc2;
	for( ssc1=1; ssc1 < total_ss_chunks; ssc1++ ){
		for( ssc2=ssc1+1; ssc2 <= total_ss_chunks; ssc2++ ){

			// Start of SSC1 - Start SSC2
			coord_r1 = p.residue( ss_chunk_start[ ssc1 ] ).atom( "CA" ).xyz();
			coord_r2 = p.residue( ss_chunk_start[ ssc2 ] ).atom( "CA" ).xyz();
			file << coord_r1.distance(coord_r2) << "\t";

			// Start of SSC1 - End SSC2
			coord_r1 = p.residue( ss_chunk_start[ ssc1 ] ).atom( "CA" ).xyz();
			coord_r2 = p.residue( ss_chunk_end[ ssc2 ] ).atom( "CA" ).xyz();
			file << coord_r1.distance(coord_r2) << "\t";

			// End of SSC1 - Start SSC2
			coord_r1 = p.residue( ss_chunk_end[ ssc1 ] ).atom( "CA" ).xyz();
			coord_r2 = p.residue( ss_chunk_start[ ssc2 ] ).atom( "CA" ).xyz();
			file << coord_r1.distance(coord_r2) << "\t";

			// End of SSC1 - End SSC2
			coord_r1 = p.residue( ss_chunk_end[ ssc1 ] ).atom( "CA" ).xyz();
			coord_r2 = p.residue( ss_chunk_end[ ssc2 ] ).atom( "CA" ).xyz();
			file << coord_r1.distance(coord_r2) << "\t";			
		}
	}
	
	file << std::endl;
	file.close();

	/////////////////////////////////////////////////////////////////////////////
	// TORSION ANGLES
	/////////////////////////////////////////////////////////////////////////////

	fn = filename_base + ".features.torsion";
	// file.open( fn.c_str(), std::ofstream::out | std::ofstream::app );
	file.open( fn.c_str() );

	for ( Size pos = 1; pos <= p.total_residue(); pos++ ) {
		if ( !p.residue(pos).is_protein() ) continue;
		file << p.phi( pos ) << "\t" <<  p.psi( pos ) << "\t" << p.omega( pos ) << "\t";
	}

	file << std::endl;
	file.close();

	/////////////////////////////////////////////////////////////////////////////
	// TORSION ANGLES (LOOP REGIONS)
	/////////////////////////////////////////////////////////////////////////////

	fn = filename_base + ".features.torsionloops";
	// file.open( fn.c_str(), std::ofstream::out | std::ofstream::app );
	file.open( fn.c_str() );

	for ( Size pos = 1; pos <= p.total_residue(); pos++ ) {
		if ( !p.residue(pos).is_protein() ) continue;

		if( is_loop( pos ) )
			file << p.phi( pos ) << "\t" <<  p.psi( pos ) << "\t" << p.omega( pos ) << "\t";
	}

	file << std::endl;
	file.close();

	/////////////////////////////////////////////////////////////////////////////
	// ABEGO 
	/////////////////////////////////////////////////////////////////////////////

	fn = filename_base + ".features.abego";
	// file.open( fn.c_str(), std::ofstream::out | std::ofstream::app );
	file.open( fn.c_str() );

	core::util::ABEGOManagerOP ABEGO( new core::util::ABEGOManager());

	for ( Size pos = 1; pos <= p.total_residue(); pos++ ) {
		if ( !p.residue(pos).is_protein() ) continue;
		file << ABEGO->torsion2index( p.phi( pos ), p.psi( pos ), p.omega( pos ) ) << "\t";
	}

	file << std::endl;
	file.close();

	/////////////////////////////////////////////////////////////////////////////
	// ABEGO (LOOP REGIONS)
	/////////////////////////////////////////////////////////////////////////////

	fn = filename_base + ".features.abegoloops";
	// file.open( fn.c_str(), std::ofstream::out | std::ofstream::app );
	file.open( fn.c_str() );

	for ( Size pos = 1; pos <= p.total_residue(); pos++ ) {
		if ( !p.residue(pos).is_protein() ) continue;

		if( is_loop( pos ) )
			file << ABEGO->torsion2index( p.phi( pos ), p.psi( pos ), p.omega( pos ) ) << "\t";
	}

	file << std::endl;
	file.close();

	/////////////////////////////////////////////////////////////////////////////
	// FRAGMENT IDS 
	/////////////////////////////////////////////////////////////////////////////

	fn = filename_base + ".features.fragid";
	// file.open( fn.c_str(), std::ofstream::out | std::ofstream::app );
	file.open( fn.c_str() );

	for ( Size pos = 1; pos <= p.total_residue(); pos++ ) {
		if ( !p.residue(pos).is_protein() ) continue;
		file << p.get_mgf_fragid( pos ) << "\t";
	}

	file << std::endl;
	file.close();

	/////////////////////////////////////////////////////////////////////////////
	// FRAGMENT IDS (LOOP REGIONS)
	/////////////////////////////////////////////////////////////////////////////

	fn = filename_base + ".features.fragidloops";
	// file.open( fn.c_str(), std::ofstream::out | std::ofstream::app );
	file.open( fn.c_str() );

	for ( Size pos = 1; pos <= p.total_residue(); pos++ ) {
		if ( !p.residue(pos).is_protein() ) continue;

		if( is_loop( pos ) )
			file << p.get_mgf_fragid( pos ) << "\t";
	}

	file << std::endl;
	file.close();

	/////////////////////////////////////////////////////////////////////////////
	// RMSD to native
	/////////////////////////////////////////////////////////////////////////////

	fn = filename_base + ".features.rmsd";
	// file.open( fn.c_str(), std::ofstream::out | std::ofstream::app );
	file.open( fn.c_str() );

	file << core::scoring::CA_rmsd( p, *native_pose ) << "\t";

	file << std::endl;
	file.close();
}
*/
/////////////////////////////////////////////////////////////////////////////
/*void MGFGlobal::abego( core::pose::Pose & p ){

	std::cout << "Analising ABEGO!!!!" << std::endl;

	// MGFUtils::MGFStrategyConfigurationOP ptr( new MGFUtils::MGFStrategyConfiguration( basic::options::option[ basic::options::OptionKeys::abinitio::MGF_strategy ] ) );
	core::util::ABEGOManagerOP ABEGO( new core::util::ABEGOManager());
	std::cout << ABEGO->get_abego_string( ABEGO->get_symbols(p) ) << std::endl;

}*/
/////////////////////////////////////////////////////////////////////////////
/*void MGFGlobal::report_frame_fragment( Size frame, Size frag ){

	std::string filename;

	// std::cout << "STAGE " << stage << ": using fragment " << frag << " of frame " << frame << std::endl;
	if(stage == 1){
		filename = get_filename_base( "stage1" ) + ".frame_frag";
	}else if(stage == 2){
		filename = get_filename_base( "stage2" ) + ".frame_frag";
	}else if(stage == 3){
		filename = get_filename_base( "stage3" ) + ".frame_frag";
	}else if(stage == 4){
		filename = get_filename_base( "stage4" ) + ".frame_frag";
	}

	std::ofstream file;	
	file.open( filename.c_str(), std::ofstream::out | std::ofstream::app );
	file << frame << "\t" << frag << std::endl;
	file.close();

	if( stage>=1 && stage<=3){

		filename = get_filename_base( "all" ) + ".frame_frag";
		file.open( filename.c_str(), std::ofstream::out | std::ofstream::app );
		file << frame << "\t" << frag << std::endl;
		file.close();

	}	
}*/
/////////////////////////////////////////////////////////////////////////////
// Receives the vector of values and the elements' ids
// order=0: Ascending    order=1: Descending
void MGFGlobal::quicksort( utility::vector0< double > & value, utility::vector0< int > & ids, long int inf, long int sup, int order ){

	if (inf >= sup)   return;

	double elem_div = value[sup], temp;
	int tmp, i = inf - 1, j = sup, cont = 1;

	while (cont){
		if(order==0){
			while (i< j && value[++i] < elem_div);
			while (j>i && value[--j] > elem_div);
		}else{		
			while (i< j && value[++i] > elem_div);
			while (j>i && value[--j] < elem_div);
		}
		
		if (i < j){
			temp = value[i];
			value[i] = value[j];
			value[j] = temp;
			tmp=ids[i];
			ids[i]=ids[j];
			ids[j]=tmp;
		}else
			cont = 0;
	}	

	temp = value[i];
	value[i] = value[sup];
	value[sup] = temp;	
	tmp=ids[i];
	ids[i]=ids[sup];
	ids[sup]=tmp;
	quicksort (value, ids, inf, i - 1, order);
	quicksort (value, ids, i + 1, sup, order);
}
/////////////////////////////////////////////////////////////////////////////
// Does not receive the ids vector
// order=0: Ascending    order=1: Descending
void MGFGlobal::quicksortWID( utility::vector0< double > & value, long int inf, long int sup, int order ){ 
	
	if (inf >= sup)   return;

	double elem_div = value[sup], temp;
	int i = inf - 1, j = sup, cont = 1;

	while (cont){
		if(order==0){
			while (i< j && value[++i] < elem_div);
			while (j>i && value[--j] > elem_div);
		}else{		
			while (i< j && value[++i] > elem_div);
			while (j>i && value[--j] < elem_div);
		}
		
		if (i < j){
			temp = value[i];
			value[i] = value[j];
			value[j] = temp;
		}else
			cont = 0;
	}	

	temp = value[i];
	value[i] = value[sup];
	value[sup] = temp;	
	quicksortWID (value, inf, i - 1, order);
	quicksortWID (value, i + 1, sup, order);
}
/////////////////////////////////////////////////////////////////////////////
// Shuffles (generates a permutation of) the elements of a given integer vector
//void MGFGlobal::shuffle(utility::vector0< int > & vector, int size){
//	int i, temp, r;
//	for(i=0; i<size; i++){
//		r = random_integer(i, size-1);
//		temp = vector[i];
//		vector[i] = vector[r];
//		vector[r] = temp;
//	}
//}
/////////////////////////////////////////////////////////////////////////////
double MGFGlobal::calculate_rg_score( core::pose::Pose const & pose ){

	Size const nres( pose.total_residue() );
	Size nres_counted=0;

	///////////////////////////////////////
	//
	// RG SCORE

	// calculate center of mass
	Vector center_of_mass( 0, 0, 0 );
	for ( Size i = 1; i <= nres; ++i ) {
		// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
		if ( pose.residue(i).has_variant_type( core::chemical::REPLONLY ) ) continue;
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;

		Vector const v( pose.residue(i).nbr_atom_xyz() );
		center_of_mass += v;
		nres_counted++;
	}
	center_of_mass /= nres_counted;

	// calculate RG based on distance from center of mass
	double rg_score = 0;
	for ( Size i = 1; i <= nres; ++i ) {
		// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
		if ( pose.residue(i).has_variant_type( core::chemical::REPLONLY ) ) continue;
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;

		Vector const v( pose.residue(i).nbr_atom_xyz() );
		rg_score += v.distance_squared( center_of_mass );
	}

	// This definition of rg differs with the conventional definition which
	// divides by nres and not by nres-1.  For the sake of matching r++, it's
	// being left at nres-1 for now, but is a candidate for change in the near
	// future.
	rg_score /= (nres_counted - 1);

	return sqrt( rg_score );

}
/////////////////////////////////////////////////////////////////////////////
// SMK get median of values
template < class T >
Real const MGFGlobal::median ( std::vector < T > & values ){

	Size const length = values.size();

	if (length == 0) {
		std::cerr  << "WARNING: In MGFGlobal::median : length of supplied vector is zero!" << std::endl;
		return 0.0;
	} else if ( length == 1){
		return values.at(0);
	}

	// possible to avoid sorting the complete vector. See std::nth_element on cppreference.com
	std::sort (values.begin(), values.end());

	Size const mid_point = (length / 2);

	// if size of vector is odd, get middle element
	// if it is even, get the mean of the two middle elements
	if (length % 2 == 0){
		// zero-indexed!
		return Real( values.at(mid_point-1) + values.at(mid_point) ) / 2.0;
	} else {
		return values.at(mid_point); // e.g. if length=5, mid_point = Size(5/2) = 2. values[2] is the third element.
	}

}


} // MGFUtils
