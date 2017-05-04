// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/subproject/SMKarchive.hh
///
/// @brief	An archive of best solutions
/// @author Shaun Kandathil


#ifndef INCLUDED_protocols_moves_SMKarchive_HH
#define INCLUDED_protocols_moves_SMKarchive_HH

// type headers
#include <core/types.hh>

// Project forward headers
#include <protocols/moves/SMKarchive.fwd.hh>

// Project headers
#include <utility/pointer/ReferenceCount.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <protocols/moves/MGFContactMap.hh>
// External library headers


// C++ headers
#include <algorithm>

// Operating system headers


// Forward declarations




namespace protocols {
namespace moves {

/*
	A vector of Pose objects, copied by value, representing the set of "good" solutions seen during the optimisation
	What exactly is meant by "good" is controlled by SMKarchiveCondition objects.
	Deriving from STL containers is generally a no-no since they don't have virtual destructors.
	see http://stackoverflow.com/questions/14089088/are-stl-containers-designed-to-allow-inheritance?lq=1
*/
class SMKarchive {
public:
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef std::vector < core::pose::Pose >::iterator ArchiveIterator;
	typedef std::vector < SMKarchiveConditionOP >::iterator ConditionIterator;

	SMKarchive();

	SMKarchive(SMKarchive const & ar){
		*this = ar;
	}

	SMKarchive& operator = (SMKarchive const& ar);

	~SMKarchive();

	//add a condition. uses vector::push_back so order is important
	void push_back_condition(SMKarchiveConditionOP condition){
		conditions_.push_back(condition);
	}

	//add a structure to the archive if conditions are met
	int add_to_archive(Pose& pose);

	//output the structures in the archive
	void dump_archive_pdbs(std::string tag = "");

	//output the structures in the archive
	void dump_scored_archive_pdbs(ScoreFunction const & scorefxn, std::string tag = "");

	//for debugging: show info about the current state of the archive
	void print_archive_info(std::ostream& os);

	// Use when scorefxn changes
	void rescore_archive(ScoreFunction const & scorefxn);

	std::vector < core::pose::Pose > structures_;
	std::vector < SMKarchiveConditionOP > conditions_;

	unsigned archive_size;
};


/*
	Archive conditions dictate whether a given structure should be added to the archive.
	They are added to SMKarchive::conditions_ using the push_back_condition() function there.
	The order in which the conditions are added is important; the condition added first is evaluated first.
*/

class SMKarchiveCondition : public utility::pointer::ReferenceCount {
public:
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;

	SMKarchiveCondition(){};
	virtual ~SMKarchiveCondition(){};

	// overloads of this method should return two things:
	// the returned bool indicates whether a replacement should actually occur. true=replace, false=don't replace
	// index: the (0-based) index of the element that has been selected, according to the criterion.
	virtual bool evaluate(Pose& pose, SMKarchive& archive, unsigned& index)=0;

};

class LowEnergyCondition : public SMKarchiveCondition {
public:
	typedef core::Real Real;
	LowEnergyCondition(){};
	~LowEnergyCondition(){};

	// this will be used correctly whether you are looking for max_element or min_element. Neat!
	// should not be <= operator.
	struct EnergyComparator {
		bool operator () (Pose i, Pose j){
			return ( i.energies().total_energy() < j.energies().total_energy() );
		}
	} EC;

	// if the energy of the structure being considered is lower than the highest energy in the archive, replace it.
	// if the archive has fewer than archive_size elements, just push_back the pose to the archive. This is done by the archive automatically
	bool evaluate (Pose& pose, SMKarchive& archive, unsigned& index);

};

//How to implement?
class DiversityCondition : public SMKarchiveCondition {

public:
//	using namespace MGFUtils;
	DiversityCondition():
		Hdist_mtx_allocated(false)
	{};

	~DiversityCondition(){};

	bool evaluate (Pose& pose, SMKarchive& archive, unsigned& index);

	void calc_pairwise_Hdistances(SMKarchive& archive);
	unsigned do_each_replace_and_test(SMKarchive& archive);

	void check_mtx_allocate_status(unsigned dim);

	bool Hdist_mtx_allocated;
	unsigned **Hdist_mtx;

	std::vector < MGFUtils::MGFContactMap > contact_maps;
	//since we are only concerned with the minimum Hdist (or some other statistic)
	// just puch_back values in the upper triangle of the pairwise matrix.
	// then we can use any function that operates on a range (std;:max, std::min etc).
	std::vector < int > Hdistances;
};

} // namespace moves
} // namespace protocols


#endif // INCLUDED_protocols_moves_SMKarchive_HH
