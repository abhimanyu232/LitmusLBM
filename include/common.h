#ifndef COMMON_H
#define COMMON_H

#include <cassert>

#include <algorithm>
#include <array>
#include <memory>
#include <vector>

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
#include <cmath>
#include <cstdlib>
#include <numbers>

// std time
#include <chrono>

//!!! todo: integrate to easily switch precision and save memory
// type aliases for floating point and index types
using float_type = double;
using index_type = uint32_t;

enum class Axis : uint32_t {
	X = 0,
	Y = 1,
	Z = 2,
};

// user-defined literal for specifying float_type
constexpr float_type operator""_fp(long double v) noexcept {
	return static_cast<float_type>(v);
}

// some mathematical constants
constexpr float_type PI = std::numbers::pi_v<float_type>;
constexpr float_type INV_PI = std::numbers::inv_pi_v<float_type>;
constexpr float_type INV_SQRT_PI = std::numbers::inv_sqrtpi_v<float_type>;

constexpr float_type SQRT_2 = std::numbers::sqrt2_v<float_type>;

constexpr float_type SQRT_3 = std::numbers::sqrt3_v<float_type>;
constexpr float_type INV_SQRT_3 = std::numbers::inv_sqrt3_v<float_type>;

constexpr float_type LOG10_E = std::numbers::log10e_v<float_type>;
constexpr float_type LOG2_E = std::numbers::log2e_v<float_type>;

constexpr float_type LN_2 = std::numbers::ln2_v<float_type>;
constexpr float_type LN_10 = std::numbers::ln10_v<float_type>;

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