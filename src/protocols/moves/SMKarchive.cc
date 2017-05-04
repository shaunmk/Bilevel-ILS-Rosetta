// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/subproject/SMKarchive.cc
///
/// @brief
/// @author	Shaun Kandathil

// Project forward headers
#include <protocols/moves/SMKarchive.hh>

// type headers
#include <core/types.hh>

// Project headers
#include <core/pose/Pose.hh>

// External library headers

// Utility headers
#include <basic/Tracer.hh>

// C++ headers
// Operating system headers

namespace protocols {
namespace moves {

static basic::Tracer trace("protocols.moves.SMKarchive");

SMKarchive::SMKarchive() :
	archive_size(10)
{
	trace << "created archive!" << std::endl;
}

SMKarchive&
SMKarchive::operator = (SMKarchive const& ar){
	//		this->resize(ar.archive_size);

	archive_size = ar.archive_size;

//	for (unsigned it=0; it != archive_size; ++it){
//		this->structures_.push_back(ar.structures_.at(it));
//	}

	this->structures_ = ar.structures_;

	//		for (ConditionIterator cit = ar.conditions_.begin(); cit != ar.conditions_.end(); ++cit){
	//			this->push_back_condition(*cit);
	//		}

	this->conditions_ = ar.conditions_;

	return *this;
}

SMKarchive::~SMKarchive(){};

int
SMKarchive::add_to_archive(Pose& pose){
	if (this->structures_.size() < archive_size){
		this->structures_.push_back(pose);
		trace << "Adding Lmin to initial archive (" << this->structures_.size() << "/" << archive_size << ")" << std::endl;
		return 1;
	}

	// if the archive is full, then the list of SMKarchiveConditions has to be evaluated.
	// as a first trial, we add something iff all conditions return true.
	// Eventually we will want to use some multiobjective way to do this

//	SMKarchive::ArchiveIterator str_to_replace; 	// will (and should!) only be dereferenced if it has been assigned to something
	// You can't check for NULL

	unsigned str_to_replace; //this is now simply an index number

	for (ConditionIterator cond_iterator=conditions_.begin(); cond_iterator!=conditions_.end(); cond_iterator++){

		str_to_replace = 9999; //default for debugging

		// (*cond_iterator) is an OP
		if( (*cond_iterator)->evaluate(pose, *this, str_to_replace) ) {

			// debug output is in the evaluate() method for each condition

			structures_.at(str_to_replace) = pose;
			print_archive_info(trace);
		}
	}

	return 0;
}

void
SMKarchive::print_archive_info(std::ostream& os){
	os << "CURRENT ARCHIVE STATE" << std::endl;
	os << "structure\tscore" << std::endl; // maybe add rmsd to native?
	for (unsigned i=1; i <= this->structures_.size(); i++){
		os << i << "\t\t\t" << this->structures_.at(i-1).energies().total_energy() << std::endl;
	}
}

void
SMKarchive::rescore_archive(ScoreFunction const & scorefxn) {
	if (this->structures_.empty())
		return;

	trace << "Rescore triggered. Printing state of archive pre-rescore:" << std::endl;
	print_archive_info(trace);

	for(ArchiveIterator it=this->structures_.begin(); it!=this->structures_.end(); it++){
		scorefxn(*it);
	}

	trace<< "Archive rescored; new state:" << std::endl;
	print_archive_info(trace);
}

void
SMKarchive::dump_archive_pdbs(std::string tag){
//	char filenumber[5];

	for(unsigned i=0; i<this->structures_.size(); i++){
		std::stringstream filename;
//		std::sprintf(filenumber, "%04u", i+1); //u is for unsigned
		filename << "ArchiveStruct" << tag << "_" << i+1 << ".pdb";
		this->structures_.at(i).dump_pdb(filename.str());
		trace << "Outputting archive structure " << filename.str() << std::endl;
	}
}

void
SMKarchive::dump_scored_archive_pdbs(ScoreFunction const & scorefxn, std::string tag){
//	char filenumber[5];

	for(unsigned i=0; i<this->structures_.size(); i++){
		std::stringstream filename;
//		std::sprintf(filenumber, "%04u", i+1); //u is for unsigned
		filename << "ArchiveStruct" << tag << "_" << i+1 << ".pdb";
		this->structures_.at(i).dump_scored_pdb(filename.str(), scorefxn, "test");
		trace << "Outputting archive structure " << filename.str() << std::endl;
	}
}


/////////////////////////////////////////////////////////////////////////
// LOW ENERGY CONDITION

bool
LowEnergyCondition::evaluate (Pose& pose, SMKarchive& archive, unsigned& index){

	bool result=false;
	trace << "Evaluating LowEnergyCondition" << std::endl;
	trace << "Printing archive state before sorting:" << std::endl;
	archive.print_archive_info(trace);

	std::sort(archive.structures_.begin(), archive.structures_.end(), EC);
	trace << "Printing archive state after sorting:" << std::endl;
	archive.print_archive_info(trace);

	SMKarchive::ArchiveIterator max_energy_pose = std::max_element(archive.structures_.begin(), archive.structures_.end(), EC);

	if (pose.energies().total_energy() < max_energy_pose->energies().total_energy()){
		result=true;
		index = std::distance(archive.structures_.begin(), max_energy_pose);
	}

	//debug messages
	if (result){
		trace << "Lmin should be added to archive at index " << index <<"; energy of current pose = " << pose.energies().total_energy() << std::endl;
		archive.print_archive_info(trace);
		trace << "Energy of pose being replaced = " << archive.structures_.at(index).energies().total_energy() << std::endl;
	} else {
		trace << "Lmin should NOT be added to archive; energy of current pose = " << pose.energies().total_energy() << std::endl;
		trace << "Index is " << index << std::endl;
		archive.print_archive_info(trace);
	}

	return result;
}

//////////////////////////////////////////////////////////////////////////
// DIVERSITY CONDITION

void
DiversityCondition::calc_pairwise_Hdistances(SMKarchive& archive){
	unsigned const dim = archive.archive_size;
	//
	//	check_mtx_allocate_status(dim);
	//	//		for (unsigned i = 0; i < dim; i++){
	//	//			for (unsigned j=0; j < dim; j++){
	//	//				Hdist_mtx[i][j] = get_bcmapdistance(archive.at(i), archive.at(j));
	//	//			}
	//	//		}
	//
	//	//compute all contact maps and store
	//	for (SMKarchive::ArchiveIterator it = archive.begin(); it != archive.end(); ++it){
	//		contact_maps.push_back(MGFUtils::MGFContactMap(*it));
	//	}
	//
	//	//get Hdist values
	//	for (unsigned i=0; i<contact_maps.size(); i++){
	//		for (unsigned j=0; j<contact_maps.size(); j++){
	//			Hdist_mtx[i][j] = contact_maps.at(i).hamming_contact_map(contact_maps.at(j));
	//		}
	//	}

	if (!Hdistances.empty())
		Hdistances.clear();

	for(unsigned i=0; i<dim; ++i) {
		for (unsigned j=i+1; j<dim; ++j) //skip self-distances
			this->Hdistances.push_back(contact_maps.at(i).hamming_contact_map(contact_maps.at(j)));
	}

	//check
	unsigned const n_elems = (dim*(dim+1))/2 - dim;
	runtime_assert(Hdistances.size()==n_elems);

}

bool
DiversityCondition::evaluate(Pose& pose, SMKarchive& archive, unsigned& index) {
	trace << "Evaluating Diversity condition (dummy)" << std::endl;
	return false;
}

void
DiversityCondition::check_mtx_allocate_status(unsigned dim){
	if (!Hdist_mtx_allocated){
		Hdist_mtx = new unsigned* [dim];
		for (unsigned i =0; i < dim; i++){
			Hdist_mtx[i] = new unsigned[dim];
			for (unsigned j=0; j<dim; j++){
				Hdist_mtx[i][j] = 0;
			}
		}
		Hdist_mtx_allocated = true;
	}
}

} // namespace moves
} // namespace protocols
