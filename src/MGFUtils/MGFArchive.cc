/// @file   MGFUtils/MGFArchive.cc
///
/// @brief
/// @author Mario Garza-Fabre

#include <MGFUtils/MGFArchive.hh>
#include <numeric/random/random.hh>
#include <core/pose/util.hh> //SMK
#include <utility/string_util.hh> //SMK
#include <fstream> //SMK
#include <iostream> //SMK

namespace MGFUtils {

using core::Size;
using std::string;
using core::pose::PoseOP;
using core::pose::Pose;
using MGFUtils::MGFProtocol::SolutionOP;
using MGFUtils::MGFProtocol::Solution;
using MGFUtils::MGFProtocol::EvaluatorOP;
using MGFUtils::MGFProtocol::Evaluator;

// using utility::vector1;
static numeric::random::RandomGenerator archive_RNG(35674); // SMK

/////////////////////////////////////////////////////////////////////////////
// Constructor 1 - both desired_size and max_size are provided
MGFArchive::MGFArchive( int desired, int max, string type, string tag ){
	tag_=tag;
	initialize( desired, max, type );
}
/////////////////////////////////////////////////////////////////////////////
// Constructor 2 - only desired_size is provided, max_size is set to 2*desired_size
MGFArchive::MGFArchive( int desired, string type, string tag){
	tag_=tag;
	initialize( desired, 2*desired, type ); // Default value for max_size
}
/////////////////////////////////////////////////////////////////////////////
// Sets main general parameters
void MGFArchive::initialize( int desired, int max, string type ){

	// Set archive parameters and initialize
	desired_size = desired;
	max_size = ( max >= desired_size ) ? max : desired_size; // max_size needs to be greater than or equal to desired_size
	current_size = 0;
	archive.clear();
	auxiliar.clear();
	archive_type = type;

	DUPLICATE_RMSD_CUTOFF = 1.0;
	NICHE_RADIUS = 1.0;
	SR_PROB = 0.5;
	//	RG_CUTOFF = 0.0;
	n_reductor_calls=0;

	size_after_last_reduction = 0; //SMK

	// Configure archiver
	configure();

}
/////////////////////////////////////////////////////////////////////////////
// Sets the evaluator, reductor, and specific parameters of the implemented archivers
void MGFArchive::configure( ){

	if( archive_type == "ENERGY" ){
		// Only store lowest-energy structures.
		// Should give identical results as SMKarchive with LowEnergyCriterion.
		evaluator = EvaluatorOP( new Evaluator("STAGE4") ); // Default evaluator
		reductor = &MGFUtils::MGFArchive::reduce_archive_size_ENERGY;

	}else if( archive_type == "ENERGY_DA" ){
		// Low-energy criterion, with Duplicate Avoidance.
		// Duplication is defined by RMSD: structures within DUPLICATE_RMSD_CUTOFF are duplicates.
		// This parameter should be defined just after creating an object of this archive type.
		evaluator = EvaluatorOP( new Evaluator("STAGE4") ); // Default evaluator
		reductor = &MGFUtils::MGFArchive::reduce_archive_size_ENERGY_DA;

	}else if( archive_type == "DEGRADATION" ){
		// Solutions within NICHE_RADIUS are degraded (i.e. rated worse). Idea is to promote diversity
		// The NICHE_RADIUS parameter should be defined just after creating an object of this archive type.
		evaluator = EvaluatorOP( new Evaluator("STAGE4") ); // Default evaluator
		reductor = &MGFUtils::MGFArchive::reduce_archive_size_DEGRADATION;

	}else if( archive_type == "DIVERSITY" ){
		// Only structural diversity (RMSD) as a criterion
		evaluator = EvaluatorOP( new Evaluator("STAGE4") ); // Default evaluator
		reductor = &MGFUtils::MGFArchive::reduce_archive_size_DIVERSITY;

	}else if( archive_type == "COMMA" ){

		evaluator = EvaluatorOP( new Evaluator("STAGE4") ); // Default evaluator
		reductor = &MGFUtils::MGFArchive::reduce_archive_size_COMMA;

	}else if( archive_type == "SR" ){
		// Stochastic Ranking between RMSD (optionally including low RG), and energy
		evaluator = EvaluatorOP( new Evaluator("STAGE4") ); // Default evaluator
		reductor = &MGFUtils::MGFArchive::reduce_archive_size_SR;

	}else if( archive_type == "SRCM" ){
		// Stochastic Ranking between Hamming distance (contact maps) and energy
		evaluator = EvaluatorOP( new Evaluator("STAGE4") ); // Default evaluator
		reductor = &MGFUtils::MGFArchive::reduce_archive_size_SR_CM;

	}else if( archive_type == "ELITIST_RANDOM" ){
		// Elitist step (energy), followed by random decisions
		evaluator = EvaluatorOP( new Evaluator("STAGE4") ); // Default evaluator
		reductor = &MGFUtils::MGFArchive::reduce_archive_size_ELITIST_RANDOM;

	}else{

		std::cerr << "\n\tERROR: Unrecognised archive type: " << archive_type << std::endl << std::endl;
		exit(1);

	}

}
/////////////////////////////////////////////////////////////////////////////
// Sets a new evaluator, and re-evaluates all archive members if the new evaluator is different to the previous one
void MGFArchive::set_evaluator( string const & new_evaluator, bool const reevaluate ){

	// if( evaluator->get_evaluator() != new_evaluator ){

	// Set new evaluator
	evaluator->set_evaluator( new_evaluator );

	if( reevaluate ){

		// Evaluate all archive members
		for(int i=0; i<current_size; i++){

			evaluator->evaluate( *archive[ i ] );

		}
	}

	// }
}
/////////////////////////////////////////////////////////////////////////////
// Returns the archive member at the given index
// TODO use archive.at(). does bounds checking internally
SolutionOP MGFArchive::get_member( int index ){

	if( index < current_size ){

		return archive[ index ];
	}else{
		utility_exit_with_message("Accesing to archive member with index out of range!!! " );
	}

	return NULL;
}
/////////////////////////////////////////////////////////////////////////////
// Returns the pose of archive member at the given index
// TODO use archive.at(). does bounds checking internally
PoseOP MGFArchive::get_member_pose( int index ){

	if( index < current_size ){

		return archive[ index ]->get_pose();

	}else{
		utility_exit_with_message("Accesing to archive member with index out of range!!! " );
	}

	return NULL;
}
/////////////////////////////////////////////////////////////////////////////
void MGFArchive::add_pose( Pose const & p ){

	// Create a copy of the given pose and include it in the archive
	add( SolutionOP( new Solution( p ) ) );

}
/////////////////////////////////////////////////////////////////////////////
void MGFArchive::add_extended_pose( ){

	// Create new solution with an extended pose, and include it in the archive
	add( SolutionOP( new Solution( ) ) );

}
/////////////////////////////////////////////////////////////////////////////
void MGFArchive::add_existing_solution( SolutionOP solution ){

	// Simply add the given solution pointer
	add( solution );

}
/////////////////////////////////////////////////////////////////////////////
void MGFArchive::add_pose_in_position(Pose const & pose, int pos){
	//	SolutionOP sol = SolutionOP(new Solution(pose));

	archive.at(pos) = SolutionOP(new Solution(pose));
}

/////////////////////////////////////////////////////////////////////////////
void MGFArchive::add( SolutionOP solution ){

	solution->set_newness(true); //SMK every solution is "new" when it has just been added to the archive

	// Add solution to the archive
	archive.push_back( solution );

	// Evaluate the inserted solution. This calls scorefxn as well as sets the internal energy variable.
	evaluator->evaluate( *archive[ current_size ] );

	// Increase archive size
	current_size++;

	// Is it time to reduce/prune the archive?
	if( current_size > max_size ) reduce_archive();

}
/////////////////////////////////////////////////////////////////////////////
void MGFArchive::reduce_archive(){
	( this->*reductor )( desired_size );

}
/////////////////////////////////////////////////////////////////////////////
void MGFArchive::reduce_archive( int size ){
	( this->*reductor )( size );
	//	this->assess_solutions_and_output();
}
/////////////////////////////////////////////////////////////////////////////
// Reduces the archive to the given size based on ENERGY
void MGFArchive::reduce_archive_size_ENERGY( int required_size ){

	// Only process the archive if reduction is required
	if( required_size < current_size ){
		n_reductor_calls++;
		// std::cout << "\nENERGY-based reductor. Current size: " << current_size << "   Required size: " << required_size << std::endl;

		// Sort solutions based on energy value
		utility::vector0< int > index( current_size ); //int index[ current_size ];
		utility::vector0< double > value( current_size ); //double value[ current_size ];			
		for(int i=0; i<current_size; i++){
			index[i] = i;
			value[i] = archive[i]->get_energy();
		}
		//		std::cout << "Printing state of MGF energy archive before sort:" << std::endl;
		//		this->print_current_archive();
		MGFUtils::MGFGlobal::quicksort(value, index, 0, current_size-1, 0);

		//		std::cout << "Printing state of MGF energy archive after sort:" << std::endl;
		//		this->print_current_archive();

//		assess_frac_retained_from_previous_reduction(required_size, index);// SMK

		// Fill the new archive with the best half of the solutions
		for(int i=0; i < required_size; i++){
			auxiliar.push_back( archive[ index[i] ] );
		}

		// Empty old archive, and replace it with the new one
		archive.clear();
		archive.swap( auxiliar );

		//		std::cout << "ENERGY archive pruned to size " << required_size
		//						<< ", current_size = " << current_size<< ", archive.size() = "<< archive.size()
		//						<< ". set current_size=required_size." << std::endl;

		// Set current size of the archive
		current_size = required_size;
		//		std::cout << "Printing state of MGF energy archive after reduction:" << std::endl;
		//		this->print_current_archive();

		// SMK check newness of each solution in the new archive
		// set newness to false, and output as a traj point

		// SMK update this var
		size_after_last_reduction = required_size;

//		this->assess_solutions_and_output(); //SMK
	}

}
/////////////////////////////////////////////////////////////////////////////
// Reduces the archive to the given size based on ENERGY, DUPLICATE AVOIDANCE
void MGFArchive::reduce_archive_size_ENERGY_DA( int required_size ){

	// Only process the archive if reduction is required
	if( required_size < current_size ){
		n_reductor_calls++;
		// std::cout << "\nENERGY_DA-based reductor. Current size: " << current_size << "   Required size: " << required_size << std::endl;

		// Sort solutions based on energy value
		utility::vector0< int > index( current_size ); //int index[ current_size ];
		utility::vector0< double > value( current_size ); //double value[ current_size ];			
		for(int i=0; i<current_size; i++){
			index[i] = i;
			value[i] = archive[i]->get_energy();
		}		
		MGFUtils::MGFGlobal::quicksort(value, index, 0, current_size-1, 0);

		// Fill the new archive with the best half of nonduplicate solutions
		utility::vector0< int > duplicate(0);
		for(int i=0; i < current_size && auxiliar.size() < (unsigned long)required_size; i++){

			// To avoid duplicate solutions, 
			// each new candidate is to be compared to all solutions previously included in the archive.
			// If the candidate is identified to be a duplicate (is within a given RMSD to another solution),
			// then we keep the solution with the lowest energy value (the one already included, therefore candidate is discarded)
			if( !is_duplicate( archive[ index[i] ], auxiliar ) ){
				auxiliar.push_back( archive[ index[i] ] );
			}else{
				duplicate.push_back( index[i] );
			}

		}

//		assess_frac_retained_from_previous_reduction(required_size, index);// SMK

		// If the number of 'unique' (nonduplicated) solutions is lower than the number of required solutions
		// then include the best duplicates to the archive (this allows to output the required number of solutions)
		for(int i=0; auxiliar.size() < (unsigned long)required_size; i++){

			auxiliar.push_back( archive[ duplicate[i] ] );

		}

		// Empty old archive, and replace it with the new one
		archive.clear();
		archive.swap( auxiliar );

		// Set current size of the archive
		// Because of duplicates removal, the new archive can have less than 'required_size' solutions (this is no longer a problen, already addressed)
		//		std::cout << "ENERGY_DA archive pruned to size " << required_size
		//				<< ", current_size (" << current_size<< ") will be set to "<< archive.size() << std::endl;
		current_size = archive.size(); 

		// SMK update this var
		size_after_last_reduction = required_size;

//		this->assess_solutions_and_output(); //SMK

	}

}

bool MGFArchive::is_duplicate( SolutionOP const & candidate, utility::vector0< SolutionOP > & solution_set ){

	for(int i=0, size=solution_set.size() ; i<size; i++){

		if( evaluator->compute_RMSD( *candidate, *solution_set[i] ) <= DUPLICATE_RMSD_CUTOFF ){ 
			return true;
		}
	}

	return false;
}
/////////////////////////////////////////////////////////////////////////////
// Reduces the archive to the given size based on NEIGHBOUR DEGRADATION
void MGFArchive::reduce_archive_size_DEGRADATION( int required_size ){

	// Only process the archive if reduction is required
	if( required_size < current_size ){
		n_reductor_calls++;
		// std::cout << "DEGRADE-based reductor. Current size: " << current_size << "   Required size: " << required_size << std::endl;

		// Initialize all archive members as nonselected, and nondegraded
		utility::vector0< bool > selected( current_size, false ); 
		utility::vector0< int > degraded( current_size, 0 ); 

		// Main cycle
		int degrading_level = 0;	
		do{
			// Identify available (nondegraded, unselected) individuals from the current archive
			int n_available = 0;
			utility::vector0< int > index(0);
			utility::vector0< double > value(0);
			for(int i=0; i<current_size; i++){

				if( !selected[i] && degraded[i] <= degrading_level ){
					n_available++;
					index.push_back( i );
					value.push_back( archive[i]->get_energy() );
				}
			}
			// Discrimitate among available individuals
			if( n_available > 0 ){

				// Sort available solutions on the basis of energy 
				MGFUtils::MGFGlobal::quicksort(value, index, 0, n_available-1, 0);

				// Select solutions iteratively starting from the best
				// Degrade solutions within a certain distance from each selected solution
				for(int i=0; i<n_available && auxiliar.size() < (unsigned long)required_size; i++){

					// Select if it is within the acceptable degradation level
					if( degraded[ index[i] ] <= degrading_level ){

						// Include solution in the archive, and mark as selected
						auxiliar.push_back( archive[ index[i] ] );
						selected[ index[i] ] = true;

						// Degrade all neighbouring unselected solutions 
						degrade_unselected_neighbours( archive[ index[i] ], archive, selected, degraded );
					}
				}
			}
			// Increase acceptable degradation level
			degrading_level++;

		}while( auxiliar.size() < (unsigned long)required_size ); 

		// SMK TODO How to use index for assess_frac_retained_from_previous_reduction()?

		// Empty old archive, and replace it with the new one
		archive.clear();
		archive.swap( auxiliar );

		//		std::cout << "DEGRADATION archive pruned to size " << required_size
		//					<< ", current_size = " << current_size<< ", archive.size() = "<< archive.size()
		//					<< ". set current_size=required_size." << std::endl;
		// Set current size of the archive
		current_size = required_size;

		// SMK update this var
		size_after_last_reduction = required_size;

//		this->assess_solutions_and_output(); //SMK
	}

}	

void MGFArchive::degrade_unselected_neighbours( SolutionOP & solution, utility::vector0< SolutionOP > & solution_set, utility::vector0< bool > & selected, utility::vector0< int > & degraded ){

	for(int i=0, size=solution_set.size() ; i<size; i++){

		// Only degrade nonselected individuals
		if( !selected[ i ] ){

			// If the solution is within NICHE_RADIUS distance from the reference solution, degrade it!!!
			if( evaluator->compute_RMSD( *solution, *solution_set[ i ] ) <= NICHE_RADIUS ){ 
				degraded[ i ]++;
			}
		}
	}
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// Reduces the archive to the given size based on ENERGY
void MGFArchive::reduce_archive_size_DIVERSITY( int required_size ){

	// Only process the archive if reduction is required
	if( required_size < current_size ){
		n_reductor_calls++;
		// std::cout << "\nDIVERSITY-based reductor. Current size: " << current_size << "   Required size: " << required_size << std::endl;

		// For each solution in the archive, compute the minimum RMSD value from this solution to all other solutions
		// Identify also the solution with the lowest energy value
		int min_energy_idx = 0;
		utility::vector0< double > min_rmsd( current_size, 100000.0 );
		utility::vector0< int > index( current_size ); //int index[ current_size ];
		for(int i=0; i<current_size; i++){

			// Locate lowest-energy solution
			if( archive[ i ]->get_energy() < archive[ min_energy_idx ]->get_energy() ) min_energy_idx = i;

			// Compute min RMSD for each solution
			for(int j=i+1; j<current_size; j++){

				double rmsd = evaluator->compute_RMSD( *archive[ i ], *archive[ j ] );
				if( rmsd < min_rmsd[ i ] ) min_rmsd[ i ] = rmsd;
				if( rmsd < min_rmsd[ j ] ) min_rmsd[ j ] = rmsd;
			}
			// Save index for sorting
			index[ i ] = i;
		}

		// Sort solutions (descending order) based on the obtained minimum rmsd values
		// This tends to prefer selection of solutions in less crowded areas
		MGFUtils::MGFGlobal::quicksort( min_rmsd, index, 0, current_size-1, 1 );

		// ELITIST STEP: add the lowest-energy solution to the new archive
		auxiliar.push_back( archive[ min_energy_idx ] );

		// SMK TODO the index vector is handled differently in this reductor;
		// need to think about what needs to be done in order to call assess_frac_retained_from_previous_reduction().

		// DIVERSITY STEP: fill new archive with the most diversified solutions
		for(int i=0; i < current_size && auxiliar.size() < (unsigned long)required_size; i++){

			// Verify not to include the lowest-energy solution again
			if( index[i] != min_energy_idx ){
				auxiliar.push_back( archive[ index[i] ] );
			}
		}

		// Empty old archive, and replace it with the new one
		archive.clear();
		archive.swap( auxiliar );

		//		std::cout << "DIVERSITY archive pruned to size " << required_size
		//								<< ", current_size = " << current_size<< ", archive.size() = "<< archive.size()
		//								<< ". set current_size=required_size." << std::endl;

		// Set current size of the archive
		current_size = required_size;

		// SMK update this var
		size_after_last_reduction = required_size;

//		this->assess_solutions_and_output(); //SMK
	}
}
/////////////////////////////////////////////////////////////////////////////
// Reduces the archive based on the (m,l)-selection: keep second half of individuals (children)
void MGFArchive::reduce_archive_size_COMMA( int required_size ){

	// Only process the archive if reduction is required
	if( current_size == 2*desired_size ){
		n_reductor_calls++;
		// std::cout << "\nCOMMA-selection reductor. Current size: " << current_size << "   Required size: " << required_size << std::endl;

		// Preserve only second half of the individuals
		for(int i=desired_size; i<2*desired_size; i++){
			auxiliar.push_back( archive[ i ] );
		}

		//SMK TODO No index system here; how to use assess_frac_retained_from_previous_reduction()?
		// Empty old archive, and replace it with the new one
		archive.clear();
		archive.swap( auxiliar );

		// Set current size of the archive
		current_size = desired_size;

		// SMK update this var
		size_after_last_reduction = required_size;

		if( (unsigned long)current_size != archive.size() ){

			std::cerr << "\n\tERROR: Wrong assumption in comma selection on the size of archive after reduction!!!" << std::endl << std::endl;
			exit(1);

		}

//		this->assess_solutions_and_output(); //SMK

	}else{

		std::cerr << "\n\tERROR: COMMA reductor applied, but current_size is " << current_size << " and desired_size is " << desired_size << std::endl << std::endl;
		exit(1);

	}

}
/////////////////////////////////////////////////////////////////////////////
// Reduces the archive based on the (m,l)-selection: keep second half of individuals (children)
void MGFArchive::reduce_archive_size_SR( int required_size ){

	//	std::cout << "SR archive reductor called on " << this->get_tag() << "; SR_PROB is " << this->get_SR_PROB() <<std::endl;
	// Only process the archive if reduction is required
	if( required_size < current_size ){
		n_reductor_calls++;
		// std::cout << "\nSR reductor (p=" << SR_PROB <<  "). Current size: " << current_size << "   Required size: " << required_size << std::endl;

		// For each solution in the archive, compute the minimum RMSD value from this solution to all other solutions
		// Identify also the solution with the lowest energy value
		// Finally, compute the RG for each solution
		int min_energy_idx = 0;
		utility::vector0< double > min_rmsd( current_size, 100000.0 );
		//		utility::vector0< double > rg( current_size, 0.0 );
		utility::vector0< int > index( current_size ); //int index[ current_size ];
		//		RG_CUTOFF = 0.0;
		for(int i=0; i<current_size && this->SR_PROB<1.0; i++){ //IF SR==1.0 this block can be skipped

			// Locate lowest-energy solution
			if( archive[ i ]->get_energy() < archive[ min_energy_idx ]->get_energy() ) min_energy_idx = i;

			// Compute min RMSD for each solution
			for(int j=i+1; j<current_size; j++){

				double rmsd = evaluator->compute_RMSD( *archive[ i ], *archive[ j ] );
				if( rmsd < min_rmsd[ i ] ) min_rmsd[ i ] = rmsd;
				if( rmsd < min_rmsd[ j ] ) min_rmsd[ j ] = rmsd;

			}

			// Save index for sorting
			index[ i ] = i;

			// Compute RG
			//			rg[ i ] = MGFUtils::MGFGlobal::calculate_rg_score( *(archive[ i ]->get_pose()) );
			//			RG_CUTOFF += rg[ i ];

		}
		//		RG_CUTOFF /= current_size;


		//		std::cout << "\nSR reductor (p=" << SR_PROB <<  "). Current size: " << current_size << "   Required size: " << required_size << ". RG_CUTOFF = " << RG_CUTOFF << std::endl;


		// ELITIST STEP: make the best energy solution be at the beginning of the index list
		// Later sorting step will preserve this first position without change
		index[ min_energy_idx ] = 0;
		index[ 0 ] = min_energy_idx;

		// Stochastic ranking-based sorting process
		for(int i=0,j=0; i<current_size; i++){

			bool swap_occurred;
			for(j=1, swap_occurred = false; j<current_size-1; j++){ // Skip j=0 (elitism)

				bool swap = false;

				if( archive[index[ j ]]->get_energy() == archive[index[ j+1 ]]->get_energy() ){

					if( min_rmsd[index[ j ]] == min_rmsd[index[ j+1 ]] ){

						// Decide randomly if they are equal based on both criteria
						//						swap = ( numeric::random::rg().uniform() <= 0.5 ) ? true : false;
						swap = ( archive_RNG.uniform() <= 0.5 ) ? true : false;
					}else if( min_rmsd[index[ j ]] < min_rmsd[index[ j+1 ]] ){
						swap = true;
					}

				}else if( min_rmsd[index[ j ]] == min_rmsd[index[ j+1 ]] ){

					if( archive[index[ j ]]->get_energy() > archive[index[ j+1 ]]->get_energy() ) swap = true;

				}else{

					//					if( numeric::random::rg().uniform() <= SR_PROB ){
					if( archive_RNG.uniform() <= SR_PROB ){
						// Compare based only on objective value
						if( archive[index[ j ]]->get_energy() > archive[index[ j+1 ]]->get_energy() ) swap = true;

					}else{

						// Compare based only on diversity
						if( min_rmsd[index[ j ]] < min_rmsd[index[ j+1 ]] ) swap = true;

						// Compare based on diversity, but including constraints to prevent diverse but poor quality solutions (high RG value) from being considered
						//						if( min_rmsd[index[ j ]] < min_rmsd[index[ j+1 ]] && ( rg[index[ j+1 ]] <= RG_CUTOFF || rg[index[ j ]] > rg[index[ j+1 ]] ) )  swap = true;

					}

				}

				// Swap individuals
				if( swap ){

					// Swap index (ranking position) of solutions
					int temp = index[ j ];
					index[ j ] = index[ j+1 ];
					index[ j+1 ] = temp;
					swap_occurred = true;

				}

			}

			// Terminate process if no swap ocurred
			if( !swap_occurred ) break;
		}

//		assess_frac_retained_from_previous_reduction(required_size, index);// SMK

		// Fill new archive with solutions at the top of the list
		for(int i=0; i < current_size && auxiliar.size() < (unsigned long)required_size; i++){

			auxiliar.push_back( archive[ index[i] ] );

		}

		// Empty old archive, and replace it with the new one
		archive.clear();
		archive.swap( auxiliar );

		// Set current size of the archive
		current_size = required_size;

		// SMK update this var
		size_after_last_reduction = required_size;

//		this->assess_solutions_and_output(); //SMK

	}

}

/////////////////////////////////////////////////////////////////////////////
// SR-based reductor, using HDCM as diversity measure
void MGFArchive::reduce_archive_size_SR_CM( int required_size ){

	//	std::cout << "SRCM archive reductor called on " << this->get_tag() << "; SR_PROB is " << this->get_SR_PROB() <<std::endl;

	// Only process the archive if reduction is required
	if( required_size < current_size ){
		n_reductor_calls++;
		// std::cout << "\nSR_CM reductor (p=" << SR_PROB <<  "). Current size: " << current_size << "   Required size: " << required_size << std::endl;

		// For each solution in the archive, compute the minimum HDCM value from this solution to all other solutions
		// Identify also the solution with the lowest energy value
		// Finally, compute the RG for each solution
		int min_energy_idx = 0;
		utility::vector0< int > min_hdcm( current_size, 100000 );
		//		utility::vector0< double > rg( current_size, 0.0 );
		utility::vector0< int > index( current_size ); //int index[ current_size ];
		//		RG_CUTOFF = 0.0;
		for(int i=0; i<current_size && this->SR_PROB<1.0; i++){ //if SR_PROB==1.0 this block can be skipped.

			// Locate lowest-energy solution
			if( archive[ i ]->get_energy() < archive[ min_energy_idx ]->get_energy() ) min_energy_idx = i;

			// Compute CM of i-th solution
			utility::vector0< int > CM_i = get_CM( *(archive[ i ]->get_pose()) );

			// Compute min HDCM from each solution
			for(int j=i+1; j<current_size; j++){

				// Compute CM of j-th solution
				utility::vector0< int > CM_j = get_CM( *(archive[ j ]->get_pose()) );

				// Compute HDCM
				int hdcm = hamming_distance( CM_i, CM_j );

				// Update min values
				if( hdcm < min_hdcm[ i ] ) min_hdcm[ i ] = hdcm;
				if( hdcm < min_hdcm[ j ] ) min_hdcm[ j ] = hdcm;

			}

			// Save index for sorting
			index[ i ] = i;

			// Compute RG
			//			rg[ i ] = MGFUtils::MGFGlobal::calculate_rg_score( *(archive[ i ]->get_pose()) );
			//			RG_CUTOFF += rg[ i ];

		}
		//		RG_CUTOFF /= current_size;

		// std::cout << "\nSR_CM reductor (p=" << SR_PROB <<  "). Current size: " << current_size << "   Required size: " << required_size << ". RG_CUTOFF = " << RG_CUTOFF << std::endl;

		// ELITIST STEP: make the best energy solution be at the beginning of the index list
		// Later sorting step will preserve this first position without change
		index[ min_energy_idx ] = 0;
		index[ 0 ] = min_energy_idx;

		// Stochastic ranking-based sorting process
		for(int i=0,j=0; i<current_size; i++){

			bool swap_occurred;
			for(j=1, swap_occurred = false; j<current_size-1; j++){ // Skip j=0 (elitism)

				bool swap = false;

				if( archive[index[ j ]]->get_energy() == archive[index[ j+1 ]]->get_energy() ){

					if( min_hdcm[index[ j ]] == min_hdcm[index[ j+1 ]] ){

						// Decide randomly if they are equal based on both criteria
						swap = ( archive_RNG.uniform() <= 0.5 ) ? true : false;

					}else if( min_hdcm[index[ j ]] < min_hdcm[index[ j+1 ]] ){
						swap = true;
					}

				}else if( min_hdcm[index[ j ]] == min_hdcm[index[ j+1 ]] ){

					if( archive[index[ j ]]->get_energy() > archive[index[ j+1 ]]->get_energy() ) swap = true;

				}else{

					if( archive_RNG.uniform() <= SR_PROB ){

						// Compare based only on objective value
						if( archive[index[ j ]]->get_energy() > archive[index[ j+1 ]]->get_energy() ) swap = true;

					}else{

						// Compare based only on diversity
						if( min_hdcm[index[ j ]] < min_hdcm[index[ j+1 ]] ) swap = true;

						// Compare based on diversity, but including constraints to prevent diverse but poor quality solutions (high RG value) from being considered
						//						if( min_hdcm[index[ j ]] < min_hdcm[index[ j+1 ]] && ( rg[index[ j+1 ]] <= RG_CUTOFF || rg[index[ j ]] > rg[index[ j+1 ]] ) )  swap = true;

					}

				}

				// Swap individuals
				if( swap ){

					// Swap index (ranking position) of solutions
					int temp = index[ j ];
					index[ j ] = index[ j+1 ];
					index[ j+1 ] = temp;
					swap_occurred = true;

				}

			}

			// Terminate process if no swap ocurred
			if( !swap_occurred ) break;
		}

//		assess_frac_retained_from_previous_reduction(required_size, index);// SMK

		// Fill new archive with solutions at the top of the list
		for(int i=0; i < current_size && auxiliar.size() < (unsigned long)required_size; i++){

			auxiliar.push_back( archive[ index[i] ] );

		}

		// Empty old archive, and replace it with the new one
		archive.clear();
		archive.swap( auxiliar );

		// Set current size of the archive
		current_size = required_size;

		// SMK update this var
		size_after_last_reduction = required_size;

//		this->assess_solutions_and_output(); //SMK

	}

}
//////////////////////
// SMK
void
MGFArchive::reduce_archive_size_ELITIST_RANDOM( int required_size ){
	// Only process the archive if reduction is required
	if( required_size < current_size ){

		n_reductor_calls++;

		utility::vector0< int > index( current_size );

		// same elitist step as in SR
		int min_energy_idx = 0;

		for(int i=0; i<current_size; i++){
			// Locate lowest-energy solution
			if( archive[ i ]->get_energy() < archive[ min_energy_idx ]->get_energy() ) min_energy_idx = i;
			index[ i ] = i;
		}
		index[ min_energy_idx ] = 0;
		index[ 0 ] = min_energy_idx;

		// now fill the remainder of the selected set with a random selection of the remaining solutions
		int last_index = index.size() - 1;
		int selection, temp;
//		std::cout << "Erandom: choosing between indices 1 and " << last_index << std::endl;
		for (int i=1; i < required_size; i++){
			// choose an index in the rest of the list of indices
			selection = archive_RNG.random_range(i,last_index);
//			std::cout << "selected " << selection << " with value " << index[selection] << std::endl;
			//move that to the position 'i', and move what was in position 'i' to position 'selection'
			temp = index[i];
			index[i] = index[selection];
			index[selection] = temp;
		}

//		assess_frac_retained_from_previous_reduction(required_size, index);// SMK

		// Fill new archive with solutions at the top of the list
		for(int i=0; i < required_size; i++){
			auxiliar.push_back( archive[ index[i] ] );
		}

		// Empty old archive, and replace it with the new one
		archive.clear();
		archive.swap( auxiliar );

		// Set current size of the archive
		current_size = required_size;

		// SMK update this var
		size_after_last_reduction = required_size;

//		this->assess_solutions_and_output(); //SMK
	}
}
//////////////////////
utility::vector0< int > MGFArchive::get_CM( Pose const & p ){

	int nres = p.total_residue();
	utility::vector0< int > CM(0);

	// Consider each pair of residues
	for(int r1=1; r1<=nres; r1++){

		// Get coordinates of residue 1
		numeric::xyzVector<core::Real> coord_r1 = p.residue( r1 ).atom( "CA" ).xyz();

		for(int r2=r1+1; r2<=nres; r2++){

			// Get coordinates of residue 2
			numeric::xyzVector<core::Real> coord_r2 = p.residue( r2 ).atom( "CA" ).xyz();

			// Update CM
			CM.push_back( ( coord_r1.distance(coord_r2) <= 8.0 )? 1 : 0 );

		}
	}

	return CM;
}
/////////////////////
int MGFArchive::hamming_distance( utility::vector0< int > & a, utility::vector0< int > & b ){

	// Validate size
	if( a.size() != b.size() ){

		std::cout << "\n\tERROR: Computing Hamming distance between vectors of different size!!!" << std::endl << std::endl;
		exit(1);

	}

	// Compute HD
	int hd = 0;
	for(int i=0, size=a.size(); i<size; i++){

		if( a[i] != b[i] ) hd++;

	}

	return hd;

}

/////////////////////////////////////////////////////////////////////////////
void MGFArchive::print_current_archive(){

	std::cout << "Archive size: " << current_size << std::endl;

	for(int i=0; i<current_size; i++){
		std::cout << "\t" << i << ": " << archive[i]->get_energy() << std::endl;
	}

}
/////////////////////////////////////////////////////////////////////////////
void
MGFArchive::dump_scored_archive_pdbs(core::scoring::ScoreFunction const & scorefxn, std::string info){
	// tag_ stores the type of archiving strategy used. Use info to add details like stage number

	//	for(SolutionIterator it = this->archive.begin(); it != this->archive.end(); ++it){
	for (int i = 1; i <= archive.size(); ++i){
		// core::pose::PoseOP p = archive.get_member_pose( i );
		// std::cout << "\t\t" << i << ": " << (*score_stage4_)(*p) << std::endl;
		// std::cout << "\t\t" << i << ": " << this->get_member( i )->get_energy() << std::endl;
		std::stringstream filename;
		filename << "ArchiveStruct" << info << tag_ << "_struct" << i <<".PDB";
		this->get_member_pose(i-1)->dump_scored_pdb(filename.str(),scorefxn);

		//		(*it)->pose->dump_scored_pdb( "Archive_struct_"+info+tag_+".PDB" , scorefxn);
	}
}

void
MGFArchive::replace_evaluator(string stage_tag) {
	evaluator = EvaluatorOP( new Evaluator(stage_tag) );
	for (SolutionIterator it = archive.begin(); it != archive.end(); ++it){
		evaluator->evaluate(*(*it));
	}
}


// SMK
void
MGFArchive::assess_solutions_and_output(){

	if (MGFGlobal::stage == 2 || MGFGlobal::stage == 3){

		//use tag_
		std::ofstream stats_rmsds, stats_HDs;
		std::string prefilename = "archive_metrics" + tag_;

		std::string rmsdfilename = prefilename + "_stg" + utility::to_string(MGFGlobal::stage) + "_RMSD.txt";
		std::string HDfilename = prefilename + "_stg" + utility::to_string(MGFGlobal::stage) + "_HD.txt";

		stats_rmsds.open (rmsdfilename.c_str(), std::ios::app);
		stats_HDs.open (HDfilename.c_str(), std::ios::app);

		//	std::vector < Real > rmsds;
		//	std::vector < int > HDs; // Hamming distances;

		//		stats_rmsds << MGFGlobal::stage << '\t';
		//		stats_HDs << MGFGlobal::stage << '\t';

		// pairwise comparisons
		for (int i = 0; i < current_size; ++i){

			//get pose and contact map
			Pose const p_i = *get_member_pose(i);
			utility::vector0< int > cm_i = get_CM(p_i);

			for (int j = i+1; j < current_size; ++j){

				//get pose and contact map
				Pose const p_j = *get_member_pose(j);
				utility::vector0< int > cm_j = get_CM(p_j);

				stats_rmsds << core::scoring::CA_rmsd( p_i, p_j) << '\t';
				stats_HDs << hamming_distance( cm_i, cm_j ) << '\t';
			}
		}

		//one line per reduction step
		stats_rmsds << std::endl;
		stats_HDs << std::endl;

		stats_rmsds.close();
		stats_HDs.close();
	}
}

//SMK
void MGFArchive::assess_frac_retained_from_previous_reduction(int size_after_current_reduction, const utility::vector0<int>& indices) {

	// This assumes that the solutions retained in the prior reduction step will have indices in the range [1, size_after_last_reduction].
	// this is only true because the retained solutions from the last step are at the top of the list when the current reduction step begins.
	// If e.g. a randomisation step is added before the actual reduction is done, this function won't work as intended.

	int n_retained_from_last_reduction = 0;
	double fraction = 0.0;
	if (size_after_last_reduction > 0) {
		for (int i = 0; i < size_after_current_reduction; i++) {
			// have to check all required_size solutions
			if (indices[i] < size_after_last_reduction)
				// index starts from 0
				n_retained_from_last_reduction++;
		}
		// calculate fraction
		fraction = double(n_retained_from_last_reduction) / double(size_after_last_reduction);
	}
	//write to file
	std::string fname = "frac_archived_solutions_retained" + get_tag() + ".txt";
	std::ofstream frac_retained_file(fname.c_str(), std::ios::app);
	frac_retained_file << n_retained_from_last_reduction << '\t'
			<< size_after_last_reduction << '\t' << fraction << '\t'
			<< MGFUtils::MGFGlobal::stage << std::endl;
	frac_retained_file.close();
}


} // MGFUtils
