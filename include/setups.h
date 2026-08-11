#ifndef SETUPS_H
#define SETUPS_H

// #include <array>
// #include <memory>
// #include <vector>

// #include <cassert>
// #include <cmath>

#include <concepts>
#include <type_traits>

#include "common.h"

template <typename SETUP_T, typename MODEL_T>
concept Setup = requires(SETUP_T setup, MODEL_T model) {
	{ setup.initialize(model) } -> std::same_as<void>;
};

template <typename SETUP_T, typename MODEL_T>
concept HasDiagnostics = Setup<SETUP_T,MODEL_T> &&
	requires(SETUP_T& setup, const MODEL_T& model) {	//, StepInfo info
		{
			setup.diagnose(model)
		} -> std::same_as<void>;	// setup.diagnose(model, info)
	};

#endif