#ifndef LATTICES_H
#define LATTICES_H

#include <array>
#include <memory>
#include <vector>

#include <cassert>
#include <cmath>

#include <concepts>
#include <type_traits>

template <typename L>
concept LatticeType = requires {
	// Compile-time constants exist and have sensible types
	requires std::is_integral_v<decltype(L::DIM)>;
	requires std::is_integral_v<decltype(L::Q)>;
	requires std::is_floating_point_v<decltype(L::Cs)>;
	// Static arrays exist with correct sizes (catches DIM*Q mismatch)
	requires std::tuple_size_v<std::remove_cvref_t<decltype(L::velocities)>> ==
						 static_cast<std::size_t>(L::DIM* L::Q);
	requires std::tuple_size_v<std::remove_cvref_t<decltype(L::weights)>> ==
						 static_cast<std::size_t>(L::Q);
};

// D2Q9 lattice implementation
class D2Q9 {
 public:
	static constexpr int DIM = 2;
	static constexpr int Q = 9;
	// static constexpr T0 = 1. / 3.0;
	static constexpr double Cs = 1.0 / std::sqrt(3.0);

	static constexpr std::array<int, DIM * Q> velocities = {
		// x components
		0, 1, 0, -1, 0, 1, -1, -1, 1,
		// y components
		0, 0, 1, 0, -1, 1, 1, -1, -1};

	static constexpr std::array<int, Q> reflected_index = {0, 3, 4, 1, 2,
																												 7, 8, 6, 5};

	static constexpr std::array<double, Q> weights = {
		4.0 / 9.0,	1.0 / 9.0,	1.0 / 9.0,	1.0 / 9.0, 1.0 / 9.0,
		1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};

 public:
	// static getters

	const int getQ() const { return Q; }

	// double getT0() const { return T0; }

	const double getLatticeSpeedofSound() const { return Cs; }

	const std::array<int, DIM * Q>& getVelocities() const { return velocities; }

	const std::array<double, Q>& getWeights() const { return weights; }

	static const std::string getName() { return "D2Q9"; }
};

// D2Q5 lattice implementation
class D2Q5 {
 public:
	static constexpr int DIM = 2;
	static constexpr int Q = 5;
	static constexpr std::array<int, DIM * Q> velocities = {	// x components
		0, 1, 0, -1, 0,
		// y components
		0, 0, 1, 0, -1};
	static constexpr std::array<int, Q> reflected_index = {0, 3, 4, 1, 2};
	static constexpr std::array<double, Q> weights = {
		(1.0 / 3.0), (1.0 / 6.0), (1.0 / 6.0), (1.0 / 6.0), (1.0 / 6.0)};

	// static constexpr T0 = 1. / 3.0;
	static constexpr double Cs = 1. / std::sqrt(3.0);

 public:
	// static getters
	const int getQ() const { return Q; }

	// static double getT0() const { return T0; }

	const double getLatticeSpeedofSound() const { return Cs; }

	const std::array<int, DIM * Q>& getVelocities() const { return velocities; }

	const std::array<double, Q>& getWeights() const { return weights; }

	static const std::string getName() { return "D2Q5"; }
};

// D3Q19 lattice implementation
class D3Q19 {
 public:
	static constexpr int DIM = 3;
	static constexpr int Q = 19;
	static constexpr double Cs = 1. / std::sqrt(3.0);

	static constexpr std::array<int, DIM * Q> velocities = {
		// x components
		0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1, 0, 0, 0, 0,
		// y components
		0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1,
		// z components
		0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1};

	static constexpr std::array<double, Q> weights = {
		// Rest particle
		1.0 / 3.0,

		// Face neighbors (6)
		1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0,

		// Edge neighbors (12)
		1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
		1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};

 public:
	// static getters
	const int getQ() const { return Q; }

	const double getLatticeSpeedofSound() const { return Cs; }

	const std::array<int, DIM * Q>& getVelocities() const { return velocities; }

	const std::array<double, Q>& getWeights() const { return weights; }

	static const std::string getName() { return "D3Q19"; }
};

// D3Q27 lattice implementation
class D3Q27 {
 public:
	static constexpr int DIM = 3;
	static constexpr int Q = 27;
	static constexpr double Cs = 1. / std::sqrt(3.0);
	// Initialize D3Q27 velocities: {cx, cy, cz}
	// Order: rest particle, face neighbors, edge neighbors, corner neighbors
	static constexpr std::array<int, DIM * Q> velocities = {
		// x components
		0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1, 0, 0, 0, 0, 1, 1, 1, 1,
		-1, -1, -1, -1,
		// y components
		0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1,
		1, 1, -1, -1,
		// z components
		0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1,
		1, -1, 1, -1};

	static constexpr std::array<double, Q> weights = {
		// Rest particle
		8.0 / 27.0,

		// Face neighbors (6)
		2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0,

		// Edge neighbors (12)
		1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,
		1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,

		// Corner neighbors (8)
		1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0,
		1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0};

 public:
	// static getters

	const int getQ() const { return Q; }

	const double getLatticeSpeedofSound() const { return Cs; }

	const std::array<int, DIM * Q>& getVelocities() const { return velocities; }

	const std::array<double, Q>& getWeights() const { return weights; }

	static const std::string getName() { return "D3Q27"; }
};


#endif