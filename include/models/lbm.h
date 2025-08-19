#ifndef LBM_H
#define LBM_H

#include "models.h"

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

	// Flattened storage arrays
	std::vector<double> f, f_new, rho;
	std::vector<std::array<double, DIM>> velocity;
	std::vector<double> vorticity;

	// Note: In 3D, vorticity would be a vector field

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
														const std::array<double, DIM>& u) const {
		const auto& v = LATTICE::velocities;
		const auto& w = LATTICE::weights;

		double cu = 0.0;
		double usqr = 0.0;
		for (int d = 0; d < DIM; ++d) {
			size_t cIdx = k + d * Q;	// !!! cache : strided access
			cu += v[cIdx] * u[d];
			// cu += v[k][d] * u[d];
			usqr += u[d] * u[d];
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

	void step() {
		computeMacroscopic();
		collide();
		stream();
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


#endif