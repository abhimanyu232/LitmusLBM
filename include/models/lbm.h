#ifndef LBM_H
#define LBM_H

#include "models.h"


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

				// Layout agnostic index based on LayoutPolicy
				// !!! idx  = LAYOUT::getIndex(k,i,Q,total_nodes);

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
				std::exit(EXIT_FAILURE);
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

				// Layout agnostic index based on LayoutPolicy
				// !!! idx  = LAYOUT::getIndex(k,i,Q,total_nodes);

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
// !!! convert to pull streaming
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

				// [AoS] : PUSH
				// size_t dst =
				// 	mesh->next_neighbour_index[i][k];	// k + mesh->getNodeIndex(npos) * Q;
				// size_t src = k + i * Q;
				// f[dst] = f_new[src];

				// Layout agnostic index based on LayoutPolicy : PUSH
				// dst unchanged, depends on push or pull
				// !!! src  = LAYOUT::getIndex(k,i,Q,total_nodes);

				// [AoS] : PULL
				size_t src =
					mesh
						->prev_neighbour_index[i][k];	 // k + mesh->getNodeIndex(npos) * Q;
				size_t dst = k + i * Q;
				f[dst] = f_new[src];

				// Layout agnostic index based on LayoutPolicy : PULL
				// source unchanged, depends on push or pull
				// !!! dst  = LAYOUT::getIndex(k,i,Q,total_nodes);

				// [SoA] : PUSH
				// size_t dst = k * total_nodes + mesh->getNodeIndex(npos);
				// size_t src = k * total_nodes + i;

				// f[dst] = f_new[src];
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

	void setFAt(size_t k, size_t idx, double value) {
		// [SoA]
		// f[k * total_nodes + idx] = value;

		// [AoS]
		f[k + (Q * idx)] = value;

		// Layout agnostic index
		// !!! f[LAYOUT::getIndex(k,idx,LATTICE::Q,total_nodes)];
	}

	void computeVorticity() {
		static_assert(DIM == 2,
									"Vorticity computation currently only implemented for 2D");

#pragma omp parallel for
		for (size_t i = 0; i < mesh->total_nodes; ++i) {
			auto pos = mesh->nodes[i];

			// todo: check for negatives at boundary
			// Neighbor positions with periodic boundaries
			// replace with mesh->nodes[mesh->next_neighbour_index[i][1,2,3,4]]
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

	// TODO: AoS or SoA : which is better?
	// inline Helper function to get distribution function index
	size_t getDistIndex(size_t q, std::array<size_t, DIM> pos) const {

		// Layout agnostic index
		// for collision (equilibrium calc): AoS is better, since we need all f_i sequentially for each node
		// for streaming: SoA is better, since we need neighbouring nodes for each f_i
		// !!! next_idx  = LAYOUT::getIndex(q,mesh->getNodeIndex(pos),Q,total_nodes);

		// [AoS] // !!! stride with q =>  1
		// p1f1 p1f2 p1f3 ... p1fQ, p2f1 p2f2 p2f3 ... p2fQ, ... => P1(f1,f2,f3,...,fQ), P2(f1,f2,f3,...,fQ), ...
		return mesh->getNodeIndex(pos) * Q + q;

		// [SoA] 	// !!! stride with q => total_nodes
		// p1f1 p2f1 p3f1 ... p(total_nodes)f1, p1f2 p2f2 p3f2 ... p(total_nodes)f2, ... =>
		// return q * total_nodes + mesh->getNodeIndex(pos);
	}

	double computeEquilibriumForInit(size_t k, double rho_val,
																	 const std::array<double, DIM>& u) const {
		return computeEquilibrium(k, rho_val, u);
	}
};

#endif