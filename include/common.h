#ifndef COMMON_H
#define COMMON_H

#include <omp.h>
#include <cassert>

#include <algorithm>
#include <memory>
#include <vector>
#include <array>


// filesystem and io
#include <fstream>
// #include <print>
#include <sys/stat.h>
#include <filesystem>
#include <iostream>

//
#include <concepts>
#include <type_traits>

// commmon math functions
#include <cstdlib>
#include <numbers>
#include <cmath>


// some mathematical constants
constexpr double PI = std::numbers::pi_v<double>;
constexpr double INV_PI = std::numbers::inv_pi_v<double>;
constexpr double INV_SQRT_PI = std::numbers::inv_sqrtpi_v<double>;

constexpr double SQRT_2 = std::numbers::sqrt2_v<double>;

constexpr double SQRT_3 = std::numbers::sqrt3_v<double>;
constexpr double INV_SQRT_3 = std::numbers::inv_sqrt3_v<double>;

constexpr double LOG10_E = std::numbers::log10e_v<double>;
constexpr double LOG2_E = std::numbers::log2e_v<double>;

constexpr double LN_2 = std::numbers::ln2_v<double>;
constexpr double LN_10 = std::numbers::ln10_v<double>;

// basic timer interface
namespace timer {
using namespace std::chrono_literals;
using Time = std::chrono::steady_clock;
// using Time = std::chrono::high_resolution_clock;
// can use ::system_clock which is not guaranteed to be monotonic.
// hence, prefer ::steady_clock.

using DoubleSeconds = std::chrono::duration<double>;
using FloatTimePoint = std::chrono::time_point<Time, DoubleSeconds>;

inline FloatTimePoint GetCurrentTime() {
	return Time::now();
}
}	 // namespace timer


#endif