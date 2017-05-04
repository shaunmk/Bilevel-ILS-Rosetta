/// @file   MGFUtils/MGFArchive.fwd.hh
///
/// @brief
/// @author Mario Garza-Fabre

#ifndef INCLUDED_MGFUtils_MGFArchive_FWD_HH
#define INCLUDED_MGFUtils_MGFArchive_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace MGFUtils {

class MGFArchive;
typedef utility::pointer::owning_ptr< MGFArchive > MGFArchiveOP;
typedef utility::pointer::owning_ptr< MGFArchive const > MGFArchiveCOP;

} // MGFUtils

#endif // INCLUDED_MGFUtils_MGFArchive_FWD_HH
