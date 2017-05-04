// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/subproject/SMKarchive.fwd.hh
///
/// @brief
/// @author Shaun Kandathil


#ifndef INCLUDED_protocols_moves_SMKarchive_fwd_HH
#define INCLUDED_protocols_moves_SMKarchive_fwd_HH

// Project headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace moves {

class SMKarchive;
class SMKarchiveCondition;

//typedef utility::pointer::owning_ptr < SMKarchive > SMKarchiveOP;
typedef utility::pointer::owning_ptr < SMKarchiveCondition > SMKarchiveConditionOP;


} // namespace moves
} // namespace protocols


#endif // INCLUDED_protocols_moves_SMKarchive_fwd_HH
