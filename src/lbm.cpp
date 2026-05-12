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

// LBM includes
#include "../include/common.h"
#include "../include/lattices.h"
#include "../include/mesh.h"

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

	void saveVelocitySlice2D(size_t index, size_t axis,
													 const std::string& filename) {
		static_cast<Derived*>(this)->saveVelocitySlice2D(index, axis, filename);
	};

	// Getters for simulation state
	std::string getModelName() const {
		return static_cast<const Derived*>(this)->getModelName();
	}

	size_t getCurrentStep() const {
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
	size_t getTotalNodes() const {
		return static_cast<const Derived*>(this)->getTotalNodes();
	}

	size_t getQ() const { return static_cast<const Derived*>(this)->getQ(); }

	std::array<size_t, LATTICE::DIM> getLatticeSize() const {
		return static_cast<const Derived*>(this)->getLatticeSize();
	}

	std::array<size_t, LATTICE::DIM> getPositionFromIndex(size_t idx) const {
		return static_cast<const Derived*>(this)->getPositionFromIndex(idx);
	}

	void setRhoAtIndex(size_t idx, double rho) {
		static_cast<Derived*>(this)->setRhoAtIndex(idx, rho);
	}

	void setVelocityAtIndex(size_t idx,
													const std::array<double, LATTICE::DIM>& u) {
		static_cast<Derived*>(this)->setVelocityAtIndex(idx, u);
	}

	void setFAt(size_t k, size_t idx, double value) {
		static_cast<Derived*>(this)->setFAt(k, idx, value);
	}

	double computeEquilibriumForInit(
		size_t k, double rho_val, const std::array<double, LATTICE::DIM>& u) const {
		return static_cast<const Derived*>(this)->computeEquilibriumForInit(
			k, rho_val, u);
	}

	const std::array<double, LATTICE::DIM>& getVelocityAtIndex(size_t idx) const {
		return static_cast<const Derived*>(this)->getVelocityAtIndex(idx);
	}

	const double getRhoAtIndex(size_t idx) const {
		return static_cast<const Derived*>(this)->getRhoAtIndex(idx);
	}
};

// Standard LBM implementation inheriting from Models using CRTP
template <LatticeType LATTICE>
class LBM : public Models<LBM<LATTICE>, LATTICE> {
 private:
	// todo: make these short or uint_8t
	static constexpr size_t DIM = LATTICE::DIM;
	static constexpr size_t Q = LATTICE::Q;

	std::array<size_t, DIM> domain_size;

	std::unique_ptr<Mesh<LATTICE::DIM, LATTICE>> mesh;

	size_t total_nodes;
	double tau, nu;
	size_t current_step;
	double current_time;
	double dt;

	// Fields
	// Flattened storage arrays
	//!!! why are these dynamic vectors, the depend only on above static variables.
	std::vector<double> f, f_new, rho;
	std::vector<std::array<double, DIM>> velocity;
	std::vector<double> vorticity;
	// Note: In 3D, vorticity would be a vector field

 private:
	double computeEquilibrium(size_t k, double rho_val,
														const std::array<double, DIM>& u_vel) const {
		const auto& v = LATTICE::velocities;
		const auto& w = LATTICE::weights;

		double cu = 0.0;
		double usqr = 0.0;
		for (size_t d = 0; d < DIM; ++d) {
			size_t cIdx = k + d * Q;	// !!! cache : strided access
			cu += v[cIdx] * u_vel[d];
			// cu += v[k][d] * u[d];
			usqr += u_vel[d] * u_vel[d];
		}
		return w[k] * rho_val * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * usqr);
	}

// 	void computeMacroscopic() {
// 		const auto& v = LATTICE::velocities;
// #pragma omp parallel for
// 		for (size_t i = 0; i < mesh->total_nodes; ++i) {
// 			// reset
// 			double rho_acc = 0.0;
// 			std::array<double, DIM> u_acc{};

// 			// sum over directions
// 			for (size_t k = 0; k < Q; ++k) {
// 				// [SoA]
// 				// size_t idx = k * total_nodes + i; //!!! change access pattern as above in getDistIndex()

// 				// [AoS]
// 				size_t idx = k + i * Q;	 // stride over k = 1
// 				// here prefer to have sequential access in k (over populations)

// 				double fk = this->f[idx];
// 				rho_acc += fk;
// 				for (size_t d = 0; d < DIM; ++d) {
// 					size_t cIdx = k + d * Q;
// 					u_acc[d] += v[cIdx] * fk;
// 					// u_acc[d] += v[k][d] * fk;
// 				}
// 			}
// 			// normalize
// 			rho[i] = rho_acc;
// 			if (rho_acc > 0.0) {
// 				for (size_t d = 0; d < DIM; ++d)
// 					velocity[i][d] = u_acc[d] / rho_acc;
// 			} else [[unlikely]] {
// 				// /!!! todo: throw error here, rho <= 0 is not physical behaviour
// 				std::cerr << " negative density \n ";
// 				// return;
// 			}
// 		}
// 	}

// 	void collide() {
// #pragma omp parallel for
// 		for (size_t i = 0; i < mesh->total_nodes; ++i) {
// 			for (size_t k = 0; k < Q; ++k) {
// 				// [AoS] // stride over k = 1
// 				size_t idx = k + i * Q;
// 				//!!! access pattern: [p0] 0, 1, 2, ... Q-1, [p1] 0, 1, 2, ... Q-1, [p2] 0, 1, 2, ... Q-1, ...

// 				// [SoA]
// 				// size_t idx = k * total_nodes + i;  [SoA]
// 				//!!! ^^^ results in non-sequential access: 0,total_nodes,2*total_nodes,...,1,1+total_nodes,1+2*total_nodes...,2,...

// 				double feq = computeEquilibrium(k, rho[i], velocity[i]);
// 				f_new[idx] = f[idx] - (f[idx] - feq) / tau;
// 			}
// 		}
// 	}

// 	void stream() {
// 		const auto& v = LATTICE::velocities;
// // stream each population along its direction
// #pragma omp parallel for
// 		for (size_t i = 0; i < mesh->total_nodes; ++i) {
// 			// auto pos = mesh->nodes[i];
// 			for (size_t k = 0; k < Q; ++k) {
// 				// !!! todo:  precompute neighbours or pointers to neighbours possible? instead of computing at runtime
// 				// std::array<int, DIM> npos; // todo: check for negatives at the boundary.
// 				// for (size_t d = 0; d < DIM; ++d) {
// 				// 	// applies the periodic boundary conditions, wraps at the domain boundary
// 				// 	size_t cIdx = k + d * Q;
// 				// 	npos[d] = (pos[d] + v[cIdx] + domain_size[d]) % domain_size[d];
// 				// 	// npos[d] = (pos[d] + v[k][d] + lattice_size[d]) % lattice_size[d];
// 				// }

// 				// [AoS]
// 				size_t dst =
// 					mesh->neighbour_index[i][k];	// k + mesh->getNodeIndex(npos) * Q;
// 				size_t src = k + i * Q;

// 				// [SoA]
// 				// size_t dst = k * total_nodes + mesh->getNodeIndex(npos);
// 				// size_t src = k * total_nodes + i;

// 				f[dst] = f_new[src];
// 				// todo: replace with inplace streaming,
// 				//!!! check race conditions, none exist currently
// 				//!!! possible with inplace streaming
// 			}
// 		}
// 	}

 public:
	LBM(const std::array<size_t, DIM>& size, double viscosity)
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
		for (size_t d = 0; d < DIM; ++d) {
			total_nodes *= domain_size[d];
		}

		mesh = std::make_unique<Mesh<LATTICE::DIM, LATTICE>>(domain_size);
		// mesh->writeToFile();

		// Initialize storage
		f.resize(Q * mesh->total_nodes);
		f_new.resize(Q * mesh->total_nodes);
		rho.resize(mesh->total_nodes);
		velocity.resize(mesh->total_nodes);
		vorticity.resize(mesh->total_nodes);

		std::cout << "Initialized " << LATTICE::getName() << " LBM with "
							<< mesh->total_nodes << " nodes" << std::endl;
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
		const auto& v = LATTICE::velocities;
		const auto& w = LATTICE::weights;

#pragma omp parallel for
		for (size_t i = 0; i < mesh->total_nodes; ++i) {
			// reset
			double rho_acc = 0.0;
			std::array<double, DIM> u_acc{};

			// sum over directions
			for (size_t k = 0; k < Q; ++k) {
				// [SoA]
				// size_t idx = k * total_nodes + i; //!!! change access pattern as above in getDistIndex()

				// [AoS]
				size_t idx = k + i * Q;	 // stride over k = 1
				// here prefer to have sequential access in k (over populations)

				double fk = this->f[idx];
				rho_acc += fk;
				for (size_t d = 0; d < DIM; ++d) {
					size_t cIdx = k + d * Q;
					u_acc[d] += v[cIdx] * fk;
					// u_acc[d] += v[k][d] * fk;
				}
			}
			// normalize
			rho[i] = rho_acc;
			if (rho_acc > 0.0) {
				for (size_t d = 0; d < DIM; ++d)
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
		for (size_t i = 0; i < mesh->total_nodes; ++i) {

			double usqr{0.};
			for (size_t d = 0; d < DIM; ++d) {	// dot(u,u)
				usqr += velocity[i][d] * velocity[i][d];
			}

			for (size_t k = 0; k < Q; ++k) {
				// [AoS] // stride over k = 1
				size_t idx = k + i * Q;
				//!!! access pattern: [p0] 0, 1, 2, ... Q-1, [p1] 0, 1, 2, ... Q-1, [p2] 0, 1, 2, ... Q-1, ...

				// [SoA]
				// size_t idx = k * total_nodes + i;  [SoA]
				//!!! ^^^ results in non-sequential access: 0,total_nodes,2*total_nodes,...,1,1+total_nodes,1+2*total_nodes...,2,...

				// double feq = computeEquilibrium(k, rho[i], velocity[i]);
				/**
				 * computeEquilibrium
				 */
				// const auto& v = LATTICE::velocities;
				// const auto& w = LATTICE::weights;

				double cu = 0.0;
				// double usqr = 0.0;
				for (size_t d = 0; d < DIM; ++d) {
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
		for (size_t i = 0; i < mesh->total_nodes; ++i) {
			// auto pos = mesh->nodes[i];	// precompute
			for (size_t k = 0; k < Q; ++k) {
				// !!! todo:  precompute neighbours or pointers to neighbours possible? instead of computing at runtime
				// std::array<size_t, DIM> npos; // todo: check for negative at boundary if not periodic bc
				// for (size_t d = 0; d < DIM;
				// 		 ++d) {	 // precompute this loop, avoid the math inside.
				// 	// applies the periodic boundary conditions, wraps at the domain boundary
				// 	size_t cIdx = k + d * Q;
				// 	npos[d] = (pos[d] + v[cIdx] + domain_size[d]) % domain_size[d];
				// 	// npos[d] = (pos[d] + v[k][d] + lattice_size[d]) % lattice_size[d];
				// }

				// [AoS]
				size_t dst =
					mesh->neighbour_index[i][k];	// k + mesh->getNodeIndex(npos) * Q;
				size_t src = k + i * Q;

				// [SoA]
				// size_t dst = k * total_nodes + mesh->getNodeIndex(npos);
				// size_t src = k * total_nodes + i;

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
		for (size_t i = 0; i < mesh->total_nodes; ++i) {
			auto pos = mesh->nodes[i];

			// todo: check for negatives at boundary
			// Neighbor positions with periodic boundaries
			// replace with mesh->nodes[mesh->neighbour_index[i][1,2,3,4]]
			std::array<size_t, 2> pos_xp = {(pos[0] + 1) % domain_size[0], pos[1]};
			std::array<size_t, 2> pos_xm = {
				(pos[0] - 1 + domain_size[0]) % domain_size[0], pos[1]};
			std::array<size_t, 2> pos_yp = {pos[0], (pos[1] + 1) % domain_size[1]};
			std::array<size_t, 2> pos_ym = {
				pos[0], (pos[1] - 1 + domain_size[1]) % domain_size[1]};

			size_t idx_xp = mesh->getNodeIndex(pos_xp);
			size_t idx_xm = mesh->getNodeIndex(pos_xm);
			size_t idx_yp = mesh->getNodeIndex(pos_yp);
			size_t idx_ym = mesh->getNodeIndex(pos_ym);

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
		for (size_t y = 0; y < domain_size[1]; ++y) {
			for (size_t x = 0; x < domain_size[0]; ++x) {
				std::array<size_t, 2> pos = {x, y};
				size_t idx = mesh->getNodeIndex(pos);

				double vel_mag = 0.0;
				for (size_t d = 0; d < DIM; ++d) {
					vel_mag += velocity[idx][d] * velocity[idx][d];
				}
				vel_mag = std::sqrt(vel_mag);

				outFile << vel_mag << " ";
			}
			outFile << "\n";
		}
		outFile.close();
	}

	void saveVelocitySlice2D(const size_t slice_location, const size_t slice_axis,
													 const std::string& filename) {

		if (slice_axis < 0 || slice_axis > 2) [[unlikely]] {
			std::cerr << "Error: Slice Axis out of bounds, please choose 0,1,2. ";
			return;
		}
		std::ofstream outFile(filename);
		std::array<size_t, 2> pos = {};

		if (slice_axis == 0) {
			pos[1] = slice_location;
		} else
			pos[0] = slice_location;

		for (size_t xi = 0; xi < domain_size[slice_axis]; ++xi) {
			pos[slice_axis] = xi;
			size_t idx = mesh->getNodeIndex(pos);

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
		for (size_t y = 0; y < domain_size[1]; ++y) {
			for (size_t x = 0; x < domain_size[0]; ++x) {
				std::array<size_t, 2> pos = {x, y};
				size_t idx = mesh->getNodeIndex(pos);
				outFile << vorticity[idx] << " ";
			}
			outFile << "\n";
		}
		outFile.close();
	}

	std::string getModelName() const {
		return "Standard LBM (" + LATTICE::getName() + ")";
	}

	size_t getCurrentStep() const { return current_step; }

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
	size_t getTotalNodes() const { return mesh->total_nodes; }

	size_t getQ() const { return Q; }

	double getLatticeSpeedofSound() const { return LATTICE::Cs; }

	std::array<size_t, DIM> getLatticeSize() const { return domain_size; }

	std::array<size_t, DIM> getPositionFromIndex(size_t flat_index) const {
		return mesh->nodes[flat_index];
	}

	// class member functions are default inline
	void setRhoAtIndex(size_t idx, double rho_val) { rho[idx] = rho_val; }

	const double getRhoAtIndex(size_t idx) const { return rho[idx]; }

	const std::array<double, DIM>& getVelocityAtIndex(size_t idx) const {
		return velocity[idx];
	}

	void setVelocityAtIndex(size_t idx, const std::array<double, DIM>& u) {
		velocity[idx] = u;
	}

	void setFAt(size_t k, size_t idx, double value) {
		// [SoA]
		// f[k * total_nodes + idx] = value;

		// [AoS]
		f[k + (Q * idx)] = value;
	}

	// todo: precompute and store as an unordered map perhaps. or an ordered one whichever is better.
	// then O(1) lookup at runtime.
	// inline Helper function to convert multi-dimensional position to flat index
	// !!! moved to mesh.h
	// size_t getNodeIndex(std::array<size_t, DIM> pos) const {
	// 	size_t index = 0;
	// 	for (size_t stride = 1, d = 0; d < DIM; ++d) {
	// 		index += pos[d] * stride;
	// 		stride *= domain_size[d];
	// 	}
	// 	return index;
	// }

	// TODO: AoS or SoA : which is better?
	// inline Helper function to get distribution function index
	size_t getDistIndex(size_t q, std::array<size_t, DIM> pos) const {
		// [AoS] // !!! stride with q =>  1
		// p1f1 p1f2 p1f3 ... p1fQ, p2f1 p2f2 p2f3 ... p2fQ, ... => P1(f1,f2,f3,...,fQ), P2(f1,f2,f3,...,fQ), ...
		return mesh->getNodeIndex(pos) * Q + q;

		// for collision (equilibrium calc): AoS is better, since we need all f_i sequentially for each node

		// [SoA] 	// !!! stride with q => total_nodes
		// p1f1 p2f1 p3f1 ... p(total_nodes)f1, p1f2 p2f2 p3f2 ... p(total_nodes)f2, ... =>
		// return q * total_nodes + mesh->getNodeIndex(pos);
		// for streaming: SoA is better, since we need neighbouring nodes for each f_i
	}

	double computeEquilibriumForInit(size_t k, double rho_val,
																	 const std::array<double, DIM>& u) const {
		return computeEquilibrium(k, rho_val, u);
	}
};

// todo: need class members, initial values of macro quantities etc, useful for calculations later.
//!!! maybe should derive from the model and not have it as a member.
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
		const auto size = model.getLatticeSize();
		const size_t total = model.getTotalNodes();
		const size_t Q = model.getQ();
		// set-up from Kallikounis Thesis : Re= 50 = [u0*L/nu]
		double rho0 = 1.;	 // valid for standard LB
		double u0 = 0.0287;
		float reynolds = 50;
		model.setViscosity(u0 * size[0] / reynolds);	// kinematic viscosity
		std::cout << " Viscosity Set : nu = " << model.getViscosity() << std::endl;
		for (size_t i = 0; i < total; ++i) {
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

			for (size_t k = 0; k < Q; ++k) {
				double feq = model.computeEquilibriumForInit(k, density, u);
				model.setFAt(k, i, feq);
			}
		}
	}	 // initialize

	void compute_error(double time) {

		// exact solution
		const size_t total_nodes = model.getTotalNodes();
		const auto size = model.getLatticeSize();

		// also compute a 2D slice along x-axis at the middle-ish plane.
		size_t slice_y = static_cast<int>(size[0] / 2);
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

		for (size_t i = 0; i < total_nodes; ++i) {
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
				slice_ux_sq_err_runner += local_ux_sq_err;
				slice_uy_sq_err_runner += local_uy_sq_err;
				slice_ux_exact_runner += ux_exact_sq;
				slice_uy_exact_runner += uy_exact_sq;

				// write slice to file
				double normal_pos =
					static_cast<double>(pos[0]) / static_cast<double>(size[0] - 1);
				ofile_slice << normal_pos << '\t' << ux_exact << ' ' << u_lbm[0] << '\t'
										<< uy_exact << ' ' << u_lbm[1] << std::endl;
			}

			ux_exact_runner += ux_exact_sq;
			uy_exact_runner += uy_exact_sq;

			ux_sq_err_runner += local_ux_sq_err;
			uy_sq_err_runner += local_uy_sq_err;
		}

		L2_err_ux_total = std::sqrt(ux_sq_err_runner / ux_exact_runner);
		L2_err_uy_total = std::sqrt(uy_sq_err_runner / uy_exact_runner);

		L2_err_ux_slice = std::sqrt(slice_ux_sq_err_runner / slice_ux_exact_runner);
		L2_err_uy_slice = std::sqrt(slice_uy_sq_err_runner / slice_uy_exact_runner);

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
		const size_t total = model.getTotalNodes();
		const size_t Q = model.getQ();
		for (size_t i = 0; i < total; ++i) {
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
			// NOTE: only stable ish for size 200x200.
			// NOTE: Blows up at 160x160
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
			for (size_t k = 0; k < Q; ++k) {
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
	std::string result_directory_name = "result_fields";
	std::filesystem::path dir(result_directory_name);
	if (!std::filesystem::exists(dir)) {
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
	constexpr size_t checkpoint = 50 * scale;
	constexpr size_t NSTEPS = checkpoint + 1;

	std::array<size_t, DIM> domain_size = {LX, LY};

	//!!! todo: make it so that the viscosity is set by the setup, currently confusing
	constexpr double visc = 0.01;	 // can be overridden by the setup

	// Instantiate specific model

	// typedef typename LBM<D2Q9>> model_type;
	auto model = std::make_unique<LBM<D2Q9>>(domain_size, visc);

	// Initialize the model :
	model->init();
	//!!! ideally the setup class should also do the domain init here.
	// !!! very strange way of doing things?
	// !!! we locally initialize the model by instantiating a setup object USING the model
	// shouldnt we probably use the setup class to initialize the model instead

	// Choose setup and initialize
	int test_case = 1;	// Taylor-Green Vortex
	// if (test_case == 1) {
	std::cout << "Initializing Taylor-Green Vortex..." << std::endl;
	TaylorGreenVortex<LBM<D2Q9>, D2Q9> setup(*model);
	// } 
	// else 
	// {
	// int test_case = 2;  Double Periodic Shear Layer
	// 	std::cout << "Initializing Double Periodic Shear Layer..." << std::endl;
	// DoublePeriodicShearLayer<LBM<D2Q9>, D2Q9> setup(*model);
	// }
	
	// Mass Conservation Check // cover with if-def DEBUG block
	double initial_mass = 0.;
	for (size_t i = 0; i < model->getTotalNodes(); ++i) {
		initial_mass += model->getRhoAtIndex(i);
	}
	// Mass Conservation Check // cover with if-def DEBUG block

	std::cout << "Running " << model->getModelName() << " simulation..."
						<< std::endl;

	auto time_start = timer::GetCurrentTime();
	for (size_t step = 0; step < NSTEPS; ++step) {
		// setup.compute_error(model->getCurrentTime());

		model->step();

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
	for (size_t i = 0; i < model->getTotalNodes(); ++i) {
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