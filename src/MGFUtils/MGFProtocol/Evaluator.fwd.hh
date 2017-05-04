/// @file   MGFUtils/MGFProtocol/Evaluator.fwd.hh
///
/// @brief
/// @author Mario Garza-Fabre

#ifndef INCLUDED_MGFUtils_MGFProtocol_Evaluator_FWD_HH
#define INCLUDED_MGFUtils_MGFProtocol_Evaluator_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace MGFUtils {
namespace MGFProtocol {

class Evaluator;
typedef utility::pointer::owning_ptr< Evaluator > EvaluatorOP;
typedef utility::pointer::owning_ptr< Evaluator const > EvaluatorCOP;

} // MGFProtocol
} // MGFUtils

#endif // INCLUDED_MGFUtils_MGFProtocol_Evaluator_FWD_HH
