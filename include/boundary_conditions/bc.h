#ifndef BC_H
#define BC_H

#include <concepts>
#include "common.h"

template <typename BC_type, typename Model>
concept BoundaryCondition = requires(BC_type& bc, Model& model) {
	{ bc.init(model) } -> std::same_as<void>;

	{ bc.apply_bc(model) } -> std::same_as<void>;
};

// 	template <typename SETUP_T, typename MODEL_T>
// concept HasInternalBoundary = Setup<SETUP_T,MODEL_T> &&
// 	requires(SETUP_T& setup, const MODEL_T& model) {	//, StepInfo info
// 		{
// 			setup.diagnose(model)
// 		} -> std::same_as<void>;	// setup.diagnose(model, info)
// 	};


#endif