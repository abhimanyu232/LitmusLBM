#include <omp.h>
#include <cassert>
#include <cmath>

#include <algorithm>
#include <array>
#include <memory>
#include <vector>

// filesystem and io
#include <sys/stat.h>
#include <filesystem>
#include <fstream>

#include <concepts>
#include <type_traits>

#include <iostream>
// #include <print>

#include <chrono>

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

// CRTP base class for simulation models
template <typename Derived, LatticeType LATTICE>
class Models {
 public:
	// virtual ~Models() = default;

	void init() { static_cast<Derived*>(this)->init(); }

	void step() { static_cast<Derived*>(this)->step(); }

	void computeVorticity() { static_cast<Derived*>(this)->computeVorticity(); }

	// Output functions
	void saveVelocityField(const std::string& filename) {
		static_cast<Derived*>(this)->saveVelocityField(filename);
	}

	void saveVorticityField(const std::string& filename) {
		static_cast<Derived*>(this)->saveVorticityField(filename);
	}

	void saveVelocitySlice2D(int index, int axis, const std::string& filename) {
		static_cast<Derived*>(this)->saveVelocitySlice2D(index, axis, filename);
	};

	// Getters for simulation state
	std::string getModelName() const {
		return static_cast<const Derived*>(this)->getModelName();
	}

	int getCurrentStep() const {
		return static_cast<const Derived*>(this)->getCurrentStep();
	}

	double getCurrentTime() const {
		return static_cast<const Derived*>(this)->getCurrentTime();
	}

	// Configuration functions
	void setViscosity(double nu) {
		static_cast<Derived*>(this)->setViscosity(nu);
	}

	void setRelaxationTime(double tau) {
		static_cast<Derived*>(this)->setRelaxationTime(tau);
	}

	double getViscosity() const {
		return static_cast<const Derived*>(this)->getViscosity();
	}

	double getRelaxationTime() const {
		return static_cast<const Derived*>(this)->getRelaxationTime();
	}

	// Initialization helpers for external setups
	int getTotalNodes() const {
		return static_cast<const Derived*>(this)->getTotalNodes();
	}

	int getQ() const { return static_cast<const Derived*>(this)->getQ(); }

	std::array<int, LATTICE::DIM> getLatticeSize() const {
		return static_cast<const Derived*>(this)->getLatticeSize();
	}

	std::array<int, LATTICE::DIM> getPositionFromIndex(int idx) const {
		return static_cast<const Derived*>(this)->getPositionFromIndex(idx);
	}

	void setRhoAtIndex(int idx, double rho) {
		static_cast<Derived*>(this)->setRhoAtIndex(idx, rho);
	}

	void setVelocityAtIndex(int idx, const std::array<double, LATTICE::DIM>& u) {
		static_cast<Derived*>(this)->setVelocityAtIndex(idx, u);
	}

	void setFAt(int k, int idx, double value) {
		static_cast<Derived*>(this)->setFAt(k, idx, value);
	}

	double computeEquilibriumForInit(
		int k, double rho_val, const std::array<double, LATTICE::DIM>& u) const {
		return static_cast<const Derived*>(this)->computeEquilibriumForInit(
			k, rho_val, u);
	}

	const std::array<double, LATTICE::DIM>& getVelocityAtIndex(int idx) const {
		return static_cast<const Derived*>(this)->getVelocityAtIndex(idx);
	}

	const double getRhoAtIndex(int idx) const {
		return static_cast<const Derived*>(this)->getRhoAtIndex(idx);
	}
};

// Standard LBM implementation inheriting from Models using CRTP
template <LatticeType LATTICE>
class LBM : public Models<LBM<LATTICE>, LATTICE> {
	// template <int DIM, int Q>
 private:
	static constexpr int DIM = LATTICE::DIM;
	static constexpr int Q = LATTICE::Q;

	std::array<int, DIM> domain_size;
	int total_nodes;
	double tau, nu;
	int current_step;
	double current_time;
	double dt;

	// Fields
	// Flattened storage arrays
	//!!! why are these dynamic vectors, the depend only on above static variables.
	std::vector<double> f, f_new, rho;
	std::vector<std::array<double, DIM>> velocity;
	std::vector<double> vorticity;

	// Note: In 3D, vorticity would be a vector field

	// todo: precompute and store as an unordered map perhaps. or an ordered one whichever is better.
	// then O(1) lookup at runtime.
	// inline Helper function to convert multi-dimensional position to flat index
	int getNodeIndex(const std::array<int, DIM>& pos) const {
		int index = 0;
		for (int stride = 1, d = 0; d < DIM; ++d) {
			index += pos[d] * stride;
			stride *= domain_size[d];
		}
		return index;
	}

	// TODO: AoS or SoA : which is better?
	// inline Helper function to get distribution function index
	int getDistIndex(const int q, const std::array<int, DIM>& pos) const {
		// [AoS] // !!! stride with q =>  1
		// p1f1 p1f2 p1f3 ... p1fQ, p2f1 p2f2 p2f3 ... p2fQ, ... => P1(f1,f2,f3,...,fQ), P2(f1,f2,f3,...,fQ), ...
		return getNodeIndex(pos) * Q + q;

		// for collision (equilibrium calc): AoS is better, since we need all f_i sequentially for each node

		// [SoA] 	// !!! stride with q => total_nodes
		// p1f1 p2f1 p3f1 ... p(total_nodes)f1, p1f2 p2f2 p3f2 ... p(total_nodes)f2, ... =>
		// return q * total_nodes + getNodeIndex(pos);
		// for streaming: SoA is better, since we need neighbouring nodes for each f_i
	}

	// todo: precompute at initialisation and store in an array.
	//  inline Convert flat index back to multi-dimensional position
	std::array<int, DIM> getPosition(int flat_index) const {
		std::array<int, DIM> pos;
		for (int d = 0; d < DIM; ++d) {
			pos[d] = flat_index % domain_size[d];
			flat_index /= domain_size[d];
		}
		return pos;	 // todo: NRVO using std::move
	}

 private:
	double computeEquilibrium(int k, double rho_val,
														const std::array<double, DIM>& u_vel) const {
		const auto& v = LATTICE::velocities;
		const auto& w = LATTICE::weights;

		double cu = 0.0;
		double usqr = 0.0;
		for (int d = 0; d < DIM; ++d) {
			size_t cIdx = k + d * Q;	// !!! cache : strided access
			cu += v[cIdx] * u_vel[d];
			// cu += v[k][d] * u[d];
			usqr += u_vel[d] * u_vel[d];
		}
		return w[k] * rho_val * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * usqr);
	}

	void computeMacroscopic() {
		const auto& v = LATTICE::velocities;
#pragma omp parallel for
		for (int i = 0; i < total_nodes; ++i) {
			// reset
			double rho_acc = 0.0;
			std::array<double, DIM> u_acc{};

			// sum over directions
			for (int k = 0; k < Q; ++k) {
				// [SoA]
				// int idx = k * total_nodes + i; //!!! change access pattern as above in getDistIndex()

				// [AoS]
				int idx = k + i * Q;	// stride over k = 1
				// here prefer to have sequential access in k (over populations)

				double fk = this->f[idx];
				rho_acc += fk;
				for (int d = 0; d < DIM; ++d) {
					size_t cIdx = k + d * Q;
					u_acc[d] += v[cIdx] * fk;
					// u_acc[d] += v[k][d] * fk;
				}
			}
			// normalize
			rho[i] = rho_acc;
			if (rho_acc > 0.0) {
				for (int d = 0; d < DIM; ++d)
					velocity[i][d] = u_acc[d] / rho_acc;
			} else [[unlikely]] {
				// /!!! todo: throw error here, rho <= 0 is not physical behaviour
				std::cerr << " negative density \n ";
				// return;
			}
		}
	}

	void collide() {
#pragma omp parallel for
		for (int i = 0; i < total_nodes; ++i) {
			for (int k = 0; k < Q; ++k) {
				// [AoS] // stride over k = 1
				int idx = k + i * Q;
				//!!! access pattern: [p0] 0, 1, 2, ... Q-1, [p1] 0, 1, 2, ... Q-1, [p2] 0, 1, 2, ... Q-1, ...

				// [SoA]
				// int idx = k * total_nodes + i;  [SoA]
				//!!! ^^^ results in non-sequential access: 0,total_nodes,2*total_nodes,...,1,1+total_nodes,1+2*total_nodes...,2,...

				double feq = computeEquilibrium(k, rho[i], velocity[i]);
				f_new[idx] = f[idx] - (f[idx] - feq) / tau;
			}
		}
	}

	void stream() {

		const auto& v = LATTICE::velocities;
// stream each population along its direction
#pragma omp parallel for
		for (int i = 0; i < total_nodes; ++i) {
			auto pos = getPosition(i);
			for (int k = 0; k < Q; ++k) {
				// !!! todo:  precompute neighbours or pointers to neighbours possible? instead of computing at runtime
				std::array<int, DIM> npos;
				for (int d = 0; d < DIM; ++d) {
					// applies the periodic boundary conditions, wraps at the domain boundary
					size_t cIdx = k + d * Q;
					npos[d] = (pos[d] + v[cIdx] + domain_size[d]) % domain_size[d];
					// npos[d] = (pos[d] + v[k][d] + lattice_size[d]) % lattice_size[d];
				}

				// [AoS]
				int dst = k + getNodeIndex(npos) * Q;
				int src = k + i * Q;

				// [SoA]
				// int dst = k * total_nodes + getNodeIndex(npos);
				// int src = k * total_nodes + i;

				f[dst] = f_new[src];
				// todo: replace with inplace streaming,
				//!!! check race conditions, none exist currently
				//!!! possible with inplace streaming
			}
		}
	}

 public:
	LBM(const std::array<int, DIM>& size, double viscosity)
			: domain_size(size),
				nu(viscosity),
				current_step(0),
				current_time(0.0),
				dt(1.0) {
		//todo: error checking here, check for tau > 0.5 (nu < 0) for stability
		tau = 3.0 * nu + 0.5;
		if (tau <= 0.5) {
			std::cerr << "Error: tau must be > 0.5 for stability (nu > 0). Got tau="
								<< tau << std::endl;
			std::exit(EXIT_FAILURE);
		}

		// Calculate total number of nodes
		total_nodes = 1;
		for (int d = 0; d < DIM; ++d) {
			total_nodes *= domain_size[d];
		}

		// Initialize storage
		f.resize(Q * total_nodes);
		f_new.resize(Q * total_nodes);
		rho.resize(total_nodes);
		velocity.resize(total_nodes);
		vorticity.resize(total_nodes);

		std::cout << "Initialized " << LATTICE::getName() << " LBM with "
							<< total_nodes << " nodes" << std::endl;
	}

	// todo: should also ideally call the setup->init() function here
	void init() {
		current_step = 0;
		current_time = 0.0;
	}

	// todo:!!! inline this function for better everything, lots of reuse here.
	void step() {

		/**
		 *  computeMacroscopic();
		*/
		const auto v = LATTICE::velocities;
		const auto w = LATTICE::weights;

#pragma omp parallel for
		for (int i = 0; i < total_nodes; ++i) {
			// reset
			double rho_acc = 0.0;
			std::array<double, DIM> u_acc{};

			// sum over directions
			for (int k = 0; k < Q; ++k) {
				// [SoA]
				// int idx = k * total_nodes + i; //!!! change access pattern as above in getDistIndex()

				// [AoS]
				int idx = k + i * Q;	// stride over k = 1
				// here prefer to have sequential access in k (over populations)

				double fk = this->f[idx];
				rho_acc += fk;
				for (int d = 0; d < DIM; ++d) {
					size_t cIdx = k + d * Q;
					u_acc[d] += v[cIdx] * fk;
					// u_acc[d] += v[k][d] * fk;
				}
			}
			// normalize
			rho[i] = rho_acc;
			if (rho_acc > 0.0) {
				for (int d = 0; d < DIM; ++d)
					velocity[i][d] = u_acc[d] / rho_acc;
			} else [[unlikely]] {
				// /!!! todo: throw error here, rho <= 0 is not physical behaviour
				std::cerr << " negative density \n ";
				// return;
			}
		}
		/**
		 *  computeMacroscopic();
		*/

		/**
		 * collide();
		*/
#pragma omp parallel for
		for (int i = 0; i < total_nodes; ++i) {

			// std::array<int,DIM> usqr{0.};
			double usqr{0.};
			for (size_t d = 0; d < DIM; ++d) {	// dot(u,u)
				usqr += velocity[i][d] * velocity[i][d];
			}

			for (int k = 0; k < Q; ++k) {
				// [AoS] // stride over k = 1
				int idx = k + i * Q;
				//!!! access pattern: [p0] 0, 1, 2, ... Q-1, [p1] 0, 1, 2, ... Q-1, [p2] 0, 1, 2, ... Q-1, ...

				// [SoA]
				// int idx = k * total_nodes + i;  [SoA]
				//!!! ^^^ results in non-sequential access: 0,total_nodes,2*total_nodes,...,1,1+total_nodes,1+2*total_nodes...,2,...

				// double feq = computeEquilibrium(k, rho[i], velocity[i]);
				/**
				 * computeEquilibrium
				 */
				// const auto& v = LATTICE::velocities;
				// const auto& w = LATTICE::weights;

				double cu = 0.0;
				// double usqr = 0.0;
				for (int d = 0; d < DIM; ++d) {
					size_t cIdx = k + d * Q;	// !!! cache : strided access
					cu += v[cIdx] * velocity[i][d];
					// cu += v[k][d] * u[d];
					// usqr += velocity[i][d] * velocity[i][d];
				}
				// return w[k] * rho_val * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * usqr);
				/**
				 * computeEquilibrium
				*/
				double feq =
					w[k] * rho[i] * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * usqr);

				f_new[idx] = f[idx] - (f[idx] - feq) / tau;
			}
		}
		/**
		 *  collide();
		*/

		/**
		 *   stream();
		*/
		// const auto& v = LATTICE::velocities;
// stream each population along its direction
#pragma omp parallel for
		for (int i = 0; i < total_nodes; ++i) {
			auto pos = getPosition(i);	// precompute
			for (int k = 0; k < Q; ++k) {
				// !!! todo:  precompute neighbours or pointers to neighbours possible? instead of computing at runtime
				std::array<int, DIM> npos;
				for (int d = 0; d < DIM;
						 ++d) {	 // precompute this loop, avoid the math inside.
					// applies the periodic boundary conditions, wraps at the domain boundary
					size_t cIdx = k + d * Q;
					npos[d] = (pos[d] + v[cIdx] + domain_size[d]) % domain_size[d];
					// npos[d] = (pos[d] + v[k][d] + lattice_size[d]) % lattice_size[d];
				}

				// [AoS]
				int dst = k + getNodeIndex(npos) *
												Q;	//todo: preocompute the flat index to vector mapping
				int src = k + i * Q;

				// [SoA]
				// int dst = k * total_nodes + getNodeIndex(npos);
				// int src = k * total_nodes + i;

				f[dst] = f_new[src];
				// todo: replace with inplace streaming,
				//!!! check race conditions, none exist currently
				//!!! possible with inplace streaming
			}
		}
		/**
		 *   stream();
		*/

		current_step++;
		current_time += dt;
	}

	void computeVorticity() {
		static_assert(DIM == 2,
									"Vorticity computation currently only implemented for 2D");

#pragma omp parallel for
		for (int i = 0; i < total_nodes; ++i) {
			auto pos = getPosition(i);

			// Neighbor positions with periodic boundaries
			std::array<int, 2> pos_xp = {(pos[0] + 1) % domain_size[0], pos[1]};
			std::array<int, 2> pos_xm = {
				(pos[0] - 1 + domain_size[0]) % domain_size[0], pos[1]};
			std::array<int, 2> pos_yp = {pos[0], (pos[1] + 1) % domain_size[1]};
			std::array<int, 2> pos_ym = {
				pos[0], (pos[1] - 1 + domain_size[1]) % domain_size[1]};

			int idx_xp = getNodeIndex(pos_xp);
			int idx_xm = getNodeIndex(pos_xm);
			int idx_yp = getNodeIndex(pos_yp);
			int idx_ym = getNodeIndex(pos_ym);

			// Compute derivatives using central differences
			double duydx = (velocity[idx_xp][1] - velocity[idx_xm][1]) * 0.5;
			double duxdy = (velocity[idx_yp][0] - velocity[idx_ym][0]) * 0.5;

			vorticity[i] = duydx - duxdy;
		}
	}

	// todo: handle 3D  as well
	void saveVelocityField(const std::string& filename) {
		static_assert(DIM == 2, "File output currently only implemented for 2D");

		std::ofstream outFile(filename);
		for (int y = 0; y < domain_size[1]; ++y) {
			for (int x = 0; x < domain_size[0]; ++x) {
				std::array<int, 2> pos = {x, y};
				int idx = getNodeIndex(pos);

				double vel_mag = 0.0;
				for (int d = 0; d < DIM; ++d) {
					vel_mag += velocity[idx][d] * velocity[idx][d];
				}
				vel_mag = std::sqrt(vel_mag);

				outFile << vel_mag << " ";
			}
			outFile << "\n";
		}
		outFile.close();
	}

	void saveVelocitySlice2D(const int slice_location, const int slice_axis,
													 const std::string& filename) {

		if (slice_axis < 0 || slice_axis > 2) [[unlikely]] {
			std::cerr << "Error: Slice Axis out of bounds, please choose 0,1,2. ";
			return;
		}
		std::ofstream outFile(filename);
		std::array<int, 2> pos = {};

		if (slice_axis == 0) {
			pos[1] = slice_location;
		} else
			pos[0] = slice_location;

		for (int xi = 0; xi < domain_size[slice_axis]; ++xi) {
			pos[slice_axis] = xi;
			int idx = getNodeIndex(pos);

			// save both components of velocity
			double normal_pos =
				static_cast<double>(xi) / static_cast<double>(domain_size[0] - 1);
			outFile << normal_pos << '\t' << velocity[idx][0] << '\t'
							<< velocity[idx][1] << std::endl;
		}
		outFile << "\n";
		outFile.close();
	}

	// todo: handle 3D  as well
	void saveVorticityField(const std::string& filename) {
		static_assert(DIM == 2, "File output currently only implemented for 2D");

		std::ofstream outFile(filename);
		for (int y = 0; y < domain_size[1]; ++y) {
			for (int x = 0; x < domain_size[0]; ++x) {
				std::array<int, 2> pos = {x, y};
				int idx = getNodeIndex(pos);
				outFile << vorticity[idx] << " ";
			}
			outFile << "\n";
		}
		outFile.close();
	}

	std::string getModelName() const {
		return "Standard LBM (" + LATTICE::getName() + ")";
	}

	int getCurrentStep() const { return current_step; }

	double getCurrentTime() const { return current_time; }

	void setViscosity(double viscosity) {
		nu = viscosity;
		tau = 3.0 * nu + 0.5;
	}

	void setRelaxationTime(double relaxation_time) {
		tau = relaxation_time;
		nu = (tau - 0.5) / 3.0;
	}

	double getViscosity() const { return nu; }

	double getRelaxationTime() const { return tau; }

	// Initialization helpers for external setups
	int getTotalNodes() const { return total_nodes; }

	int getQ() const { return Q; }

	double getLatticeSpeedofSound() const { return LATTICE::Cs; }

	std::array<int, DIM> getLatticeSize() const { return domain_size; }

	std::array<int, DIM> getPositionFromIndex(int flat_index) const {
		return getPosition(flat_index);
	}

	// class member functions are default inline
	void setRhoAtIndex(int idx, double rho_val) { rho[idx] = rho_val; }

	const double getRhoAtIndex(int idx) const { return rho[idx]; }

	const std::array<double, DIM>& getVelocityAtIndex(int idx) const {
		return velocity[idx];
	}

	void setVelocityAtIndex(int idx, const std::array<double, DIM>& u) {
		velocity[idx] = u;
	}

	void setFAt(int k, int idx, double value) {
		// [SoA]
		// f[k * total_nodes + idx] = value;

		// [AoS]
		f[k + (Q * idx)] = value;
	}

	double computeEquilibriumForInit(int k, double rho_val,
																	 const std::array<double, DIM>& u) const {
		return computeEquilibrium(k, rho_val, u);
	}
};

// todo: need class members, initial values of macro quantities etc, useful for calculations later.
//!!! todo!!! maybe should derive from the model and not have it as a member.
//!!! that simplifies the initialisation.
template <typename Derived, typename ModelType>
class Setups {
 public:
	ModelType& model;

	Setups(ModelType& m) : model(m) {}

	void initialize() { static_cast<Derived*>(this)->initialize(); }

	// not sure if it is a good idea doing this here
	// std::array<double, DIM> u0{0.};
	// double rho0 = 1.0;
	// double T0{model.getLatticeSpeedofSound()*model.getLatticeSpeedofSound()};
	// double p0 = rho0*T0;
};

// todo: need class members, initial values of macro quantities etc, useful for calculations later.
template <typename ModelType, LatticeType LATTICE>
class TaylorGreenVortex
		: public Setups<TaylorGreenVortex<ModelType, LATTICE>, ModelType> {
 public:
	using Setups<TaylorGreenVortex<ModelType, LATTICE>, ModelType>::model;

	explicit TaylorGreenVortex(ModelType& m)
			: Setups<TaylorGreenVortex<ModelType, LATTICE>, ModelType>(m) {

		initialize();
	}

	void initialize() {
		const auto size =
			model.getLatticeSize();	 // !!! assert that size[0] == size[1] == size[2]
		const int total = model.getTotalNodes();
		const int Q = model.getQ();
		// set-up from Kallikounis Thesis : Re= 50 = [u0*L/nu]
		double rho0 = 1.;	 // valid for standard LB
		double u0 = 0.0287;
		float reynolds = 50;
		model.setViscosity(u0 * size[0] / reynolds);	// kinematic viscosity
		std::cout << " Viscosity Set : nu = " << model.getViscosity() << std::endl;
		for (int i = 0; i < total; ++i) {
			auto pos = model.getPositionFromIndex(i);
			double x_pos = static_cast<double>(pos[0]) / (size[0] - 1);
			double y_pos = static_cast<double>(pos[1]) / (size[1] - 1);
			std::array<double, 2> u;
			u[0] = -u0 * cos(2 * M_PI * x_pos) * sin(2 * M_PI * y_pos);
			u[1] = u0 * sin(2 * M_PI * x_pos) * cos(2 * M_PI * y_pos);
			model.setVelocityAtIndex(i, u);

			double density = rho0 - 0.25 * (u0 * u0) *
																(cos(2 * M_PI * x_pos) + cos(2 * M_PI * y_pos));
			model.setRhoAtIndex(i, density);

			for (int k = 0; k < Q; ++k) {
				double feq = model.computeEquilibriumForInit(k, density, u);
				model.setFAt(k, i, feq);
			}
		}
	}	 // initialize

	void compute_error(double time) {

		// exact solution
		const int total_nodes = model.getTotalNodes();
		const auto size = model.getLatticeSize();

		// also compute a 2D slice along x-axis at the middle-ish plane.
		int slice_y = static_cast<int>(size[0] / 2);
		const std::string filename_slice =
			"slice2D" + std::to_string(slice_y) + "_" +
			std::to_string(static_cast<int>(time)) + ".dat";
		std::ofstream ofile_slice(filename_slice);
		ofile_slice << "Xpos" << '\t' << "Ux_{Exact}" << ' ' << "Ux_{LBM}" << '\t'
								<< "Uy_{Exact}" << ' ' << "Uy_{LBM}" << std::endl;

		double lattice_ref_temp =
			model.getLatticeSpeedofSound() * model.getLatticeSpeedofSound();
		double p0 = 1.0 / lattice_ref_temp;	 // valid for standard LB
		double u0 = 0.0287;
		// float reynolds = 50;
		double nu = model.getViscosity();	 // kinematic viscosity

		const double k_squared = 4 * (M_PI * M_PI) / (size[0] * size[0]);
		const double exp_K_time = std::exp(-2 * k_squared * nu * time);

		double ux_sq_err_runner{0.}, uy_sq_err_runner{0.};
		double ux_exact_runner{0.}, uy_exact_runner{0.};
		double L2_err_ux_total{0.}, L2_err_uy_total{0.};

		double slice_ux_sq_err_runner{0.}, slice_uy_sq_err_runner{0.};
		double slice_ux_exact_runner{0.}, slice_uy_exact_runner{0.};
		double L2_err_ux_slice{0.}, L2_err_uy_slice{0.};

		for (int i = 0; i < total_nodes; ++i) {
			auto u_lbm = model.getVelocityAtIndex(i);
			double ux_exact{0.}, uy_exact{0.};
			double p_exact{0.};

			auto pos = model.getPositionFromIndex(i);
			double x_pos = static_cast<double>(pos[0]) / (size[0] - 1);
			double y_pos = static_cast<double>(pos[1]) / (size[1] - 1);

			ux_exact =
				-u0 * cos(2 * M_PI * x_pos) * sin(2 * M_PI * y_pos) * exp_K_time;
			uy_exact =
				u0 * sin(2 * M_PI * x_pos) * cos(2 * M_PI * y_pos) * exp_K_time;

			p_exact = p0 - 0.25 * (u0 * u0) *
											 (cos(2 * M_PI * x_pos) + cos(2 * M_PI * y_pos)) *
											 exp_K_time;

			double ux_exact_sq = std::pow(ux_exact, 2);
			double uy_exact_sq = std::pow(uy_exact, 2);

			double local_ux_sq_err = std::pow(ux_exact - u_lbm[0], 2);
			double local_uy_sq_err = std::pow(uy_exact - u_lbm[1], 2);

			// calculate L2 vel error over slice
			if (pos[1] == slice_y) {
				slice_ux_sq_err_runner +=
					local_ux_sq_err;	// std::pow(ux_exact - u_lbm[0], 2);
				slice_uy_sq_err_runner +=
					local_uy_sq_err;	// std::pow(uy_exact - u_lbm[1], 2);
				slice_ux_exact_runner += ux_exact_sq;
				slice_uy_exact_runner += uy_exact_sq;

				// write slice to file
				double normal_pos =
					static_cast<double>(pos[0]) / static_cast<double>(size[0] - 1);
				ofile_slice << normal_pos << '\t' << ux_exact << ' ' << u_lbm[0] << '\t'
										<< uy_exact << ' ' << u_lbm[1] << std::endl;
				// ofile_slice << uy_exact << " " << normal_pos << std::endl;
			}

			ux_exact_runner += ux_exact_sq;
			uy_exact_runner += uy_exact_sq;

			ux_sq_err_runner += local_ux_sq_err;	// std::pow(ux_exact - u_lbm[0], 2);
			uy_sq_err_runner += local_uy_sq_err;	// std::pow(uy_exact - u_lbm[1], 2);
		}

		L2_err_ux_total = std::sqrt(ux_sq_err_runner / ux_exact_runner);
		L2_err_uy_total = std::sqrt(uy_sq_err_runner / uy_exact_runner);

		L2_err_ux_slice = std::sqrt(slice_ux_sq_err_runner / slice_ux_exact_runner);
		L2_err_uy_slice = std::sqrt(slice_uy_sq_err_runner / slice_uy_exact_runner);

		// // plot x-slice at y=slice_y
		// model.saveVelocitySlice2D(
		// 	(slice_y), 0, "vel_slice_y=" + std::to_string(slice_y) + ".txt");

		std::cout
			<< "Velocity RMS Error (Relative) over domain, (Ux_Err,Uy_Err) = ("
			<< L2_err_ux_total << ", " << L2_err_uy_total << ")" << std::endl;

		std::cout << "Velocity RMS Error (Relative) over slice, (Ux_Err,Uy_Err) = ("
							<< L2_err_ux_slice << ", " << L2_err_uy_slice << ")" << std::endl;
	}
};

// todo: need class members, initial values of macro quantities etc, useful for calculations later.
template <typename ModelType, LatticeType LATTICE>
class DoublePeriodicShearLayer
		: public Setups<DoublePeriodicShearLayer<ModelType, LATTICE>, ModelType> {
 public:
	using Setups<DoublePeriodicShearLayer<ModelType, LATTICE>, ModelType>::model;

	explicit DoublePeriodicShearLayer(ModelType& m)
			: Setups<DoublePeriodicShearLayer<ModelType, LATTICE>, ModelType>(m) {
		initialize();
	}

	void initialize() {
		const auto size = model.getLatticeSize();
		const int total = model.getTotalNodes();
		const int Q = model.getQ();
		for (int i = 0; i < total; ++i) {
			auto pos = model.getPositionFromIndex(i);
			double x_pos = static_cast<double>(pos[0]) / (size[0] - 1);
			double y_pos = static_cast<double>(pos[1]) / (size[1] - 1);
			std::array<double, 2> u;

			// setup: https://www.researchgate.net/publication/335029866_Pseudoentropic_derivation_of_the_regularized_lattice_Boltzmann_method
			// https://www.researchgate.net/publication/266088157_Gibbs'_principle_for_the_lattice-kinetic_theory_of_fluid_dynamics
			double u0 = 0.04;
			double perturbation = 0.05 * u0;
			double lambda = 80;	 // 20;
			float reynolds = 30000;
			// only stable ish for size 200x200.
			// Blows up at 160x160
			// !!! need entropic etc for higher reynolds numbers.
			model.setViscosity(u0 * size[0] / reynolds);
			if (y_pos <= 0.5) {
				u[0] = u0 * tanh((y_pos - 0.25) * lambda);
			} else {
				u[0] = u0 * tanh((0.75 - y_pos) * lambda);
			}
			u[1] = perturbation * sin(2 * M_PI * (x_pos + 0.25));

			model.setVelocityAtIndex(i, u);
			model.setRhoAtIndex(i, 1.0);
			for (int k = 0; k < Q; ++k) {
				double feq = model.computeEquilibriumForInit(k, 1.0, u);
				model.setFAt(k, i, feq);
			}
		}
	}
};

int main() {
	// Set up OpenMP
	std::cout << "Using " << omp_get_max_threads() << " OpenMP threads"
						<< std::endl;

	// Create output directory
	std::string result_directory_name = "lbm_ai";
	std::filesystem::path dir(result_directory_name);
	if (!std::filesystem::exists(
				dir)) {	 // check if directory exists, if yes, do nothing
		if (std::filesystem::create_directory(dir)) {	 // default permissions
			std::cout << " Results directory created : " << dir << std::endl;
		} else {
			std::cerr << "Error: Failed to create directory " << dir << std::endl;
		}
	} else {
		std::cout << " Results directory already exists : " << dir << std::endl;
	}

	// Simulation parameters
	constexpr size_t DIM = 2;
	constexpr double scale = 25.;
	constexpr size_t LX = 500;	//16 * scale;
	constexpr size_t LY = 500;	//16 * scale;
	// timestep to do validation plots and error checks.
	bool SAVE_FILE = true;
	constexpr size_t SAVE_INTERVAL = 1250;

	// constexpr size_t NSTEPS = 6000;
	constexpr int checkpoint = 50 * scale;
	constexpr size_t NSTEPS = checkpoint + 1;

	std::array<int, DIM> domain_size = {LX, LY};

	//!!! todo: make it so that the viscosity is set by the setup, currently confusing
	constexpr double visc = 0.01;	 // can be overridden by the setup

	// Use abstract Models interface with lattice
	// std::unique_ptr<Models<DIM>> model =
	// 	std::make_unique<LBM<D2Q9>>(lattice_size, visc);

	// typedef typename LBM<D2Q9>> model_type;

	// Instantiate specific model
	auto model = std::make_unique<LBM<D2Q9>>(domain_size, visc);

	// Initialize the model :
	model->init();
	//!!! ideally the setup class should also do the domain init here.
	// !!! very strange way of doing things?
	// !!! we locally initialize the model by instantiating a setup object USING the model
	// shouldnt we probably use the setup class to initialize the model instead
	// Choose setup and initialize
	int test_case = 1;	// 1: Taylor-Green, 2: Double Periodic Shear Layer
											// if (test_case == 1) {
	std::cout << "Initializing Taylor-Green Vortex..." << std::endl;

	TaylorGreenVortex<LBM<D2Q9>, D2Q9> setup(*model);
	// } else {
	// 	std::cout << "Initializing Double Periodic Shear Layer..." << std::endl;
	// DoublePeriodicShearLayer<LBM<D2Q9>, D2Q9> setup(*model);

	// Mass Conservation Check // cover with if-def DEBUG block
	double initial_mass = 0.;
	for (int i = 0; i < model->getTotalNodes(); ++i) {
		initial_mass += model->getRhoAtIndex(i);
	}
	// Mass Conservation Check // cover with if-def DEBUG block

	std::cout << "Running " << model->getModelName() << " simulation..."
						<< std::endl;

	auto time_start = timer::GetCurrentTime();
	for (size_t step = 0; step < NSTEPS; ++step) {
		// setup.compute_error(model->getCurrentTime());

		model->step();

		// Save fields every SAVE_INTERVAL steps
		if (SAVE_FILE && (step % SAVE_INTERVAL == 0)) {
			model->saveVelocityField(result_directory_name + "/velocity_" +
															 std::to_string(step) + ".txt");

			model->computeVorticity();
			model->saveVorticityField(result_directory_name + "/vorticity_" +
																std::to_string(step) + ".txt");

			std::cout << "Step " << step << " (t=" << model->getCurrentTime()
								<< ") completed" << std::endl;
		}

		// run Taylor-Green Vortex validation, probably instead do a guard block here
		double time_non_dim = model->getCurrentTime() * 0.0287 / LX;
		//!!! magic number, should come from setup->u_init;
		if (test_case == 1 && (model->getCurrentStep() % checkpoint == 0)) {
			std::cout << "Step " << step << " (t=" << model->getCurrentTime()
								<< ") completed" << "; non Dim Time = " << time_non_dim
								<< std::endl;
			setup.compute_error(model->getCurrentTime());
		}
	}
	auto time_end = timer::GetCurrentTime();

	// Mass Conservation Check // cover with if-def DEBUG block
	double final_mass = 0.;
	for (int i = 0; i < model->getTotalNodes(); ++i) {
		final_mass += model->getRhoAtIndex(i);
	}

	std::cout << "Mass Conservation Check. Initial Mass = " << initial_mass
						<< std::endl;
	std::cout << "Mass Conservation Check. Total Mass = " << final_mass
						<< std::endl;
	// Mass Conservation Check // cover with if-def DEBUG block

	std::cout << "Simulation completed. Final step: " << model->getCurrentStep()
						<< ", Final time: " << model->getCurrentTime() << std::endl;

	auto elapsed_time = time_end - time_start;
	std::cout << "total time elapsed : " << elapsed_time << std::endl;

	return 0;
}