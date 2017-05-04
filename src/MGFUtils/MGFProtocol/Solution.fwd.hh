/// @file   MGFUtils/MGFProtocol/Solution.fwd.hh
///
/// @brief
/// @author Mario Garza-Fabre

#ifndef INCLUDED_MGFUtils_MGFProtocol_Solution_FWD_HH
#define INCLUDED_MGFUtils_MGFProtocol_Solution_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace MGFUtils {
namespace MGFProtocol {

class Solution;
typedef utility::pointer::owning_ptr< Solution > SolutionOP;
typedef utility::pointer::owning_ptr< Solution const > SolutionCOP;

} // MGFProtocol
} // MGFUtils

#endif // INCLUDED_MGFUtils_MGFProtocol_Solution_FWD_HH
