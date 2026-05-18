#ifndef LBM_H
#define LBM_H

#include "../policy/layout_policy.h"
#include "models.h"

// Standard LBM implementation inheriting from Models using CRTP
template <LatticeType LATTICE, LayoutPolicy LAYOUT>
class LBM : public Models<LBM<LATTICE, LAYOUT>, LATTICE> {
 private:
	// todo: make these short or uint_8t
	static constexpr index_type DIM = LATTICE::DIM;
	static constexpr index_type Q = LATTICE::Q;

	std::array<index_type, DIM> domain_size;

	//!!! todo: (perf) : causes pointer dereferencing in the hot loop
	// consider storing it as a simple Mesh member. // !!! (profile)
	std::unique_ptr<Mesh<LATTICE::DIM, LATTICE, LAYOUT>> mesh;

	index_type total_nodes;
	float_type tau, nu;
	size_t current_step;
	double current_time;
	double dt;

	// Fields
	// Flattened storage arrays
	//!!! why are these dynamic vectors, the depend only on above static variables.
	std::vector<float_type> f, f_new, rho;
	std::vector<std::array<float_type, DIM>> velocity;
	std::vector<float_type> vorticity;
	// Note: In 3D, vorticity would be a vector field

 private:
	float_type computeEquilibrium(index_type k, float_type rho_val,
														const std::array<float_type, DIM>& u_vel) const {
		const auto& v = LATTICE::velocities;
		const auto& w = LATTICE::weights;

		float_type cu = 0.0;
		float_type usqr = 0.0;
		for (index_type d = 0; d < DIM; ++d) {
			index_type cIdx = k + d * Q;	// !!! cache : strided access
			cu += v[cIdx] * u_vel[d];
			// cu += v[k][d] * u[d];
			usqr += u_vel[d] * u_vel[d];
		}
		return w[k] * rho_val * (1.0_fp + 3.0_fp * cu + 4.5_fp * cu * cu - 1.5_fp * usqr);
	}

 public:
	LBM(const std::array<index_type, DIM>& size, float_type viscosity)
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
		for (index_type d = 0; d < DIM; ++d) {
			total_nodes *= domain_size[d];
		}

		mesh = std::make_unique<Mesh<LATTICE::DIM, LATTICE, LAYOUT>>(domain_size);
		// mesh->writeToFile();

		// Initialize storage
		// !!! todo: (mem) Q * total_nodes should fit inside uint32_t for domain < 500
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
		for (index_type i = 0; i < mesh->total_nodes; ++i) {
			// reset
			float_type rho_acc = 0.0;
			std::array<float_type, DIM> u_acc{};

			for (index_type k = 0; k < Q; ++k) {
				// note: [AoS] v [SoA]
				// here prefer to have sequential access in k (over populations)

				// Layout agnostic index based on LayoutPolicy
				index_type idx = LAYOUT::getIndex(k, i, Q, total_nodes);

				float_type fk = this->f[idx];
				rho_acc += fk;
				for (index_type d = 0; d < DIM; ++d) {
					index_type cIdx = k + d * Q;
					u_acc[d] += v[cIdx] * fk;	 // v[k][d]
				}
			}
			// normalize
			rho[i] = rho_acc;
			if (rho_acc > 0.0) {
				for (index_type d = 0; d < DIM; ++d)
					velocity[i][d] = u_acc[d] / rho_acc;
			} else [[unlikely]] {
				std::cerr << " negative density \n ";
				// std::exit(EXIT_FAILURE);
				// todo: save state, and exit after the loop ends.
				// todo: example: bool sim_error = true; later check sim_error == true and then std::exit
				// std::exit inside parallel worker region is UB for openmp
			}
		}
		/**
		 *  end computeMacroscopic();
		*/

		/**
		 * collide();
		*/
#pragma omp parallel for
		for (index_type i = 0; i < mesh->total_nodes; ++i) {

			float_type usqr{0.};
			for (index_type d = 0; d < DIM; ++d) {	// dot(u,u)
				usqr += velocity[i][d] * velocity[i][d];
			}

			for (index_type k = 0; k < Q; ++k) {

				// Layout agnostic index based on LayoutPolicy
				index_type idx = LAYOUT::getIndex(k, i, Q, total_nodes);

				/**
				 * computeEquilibrium
				*/
				float_type cu{0.0};
				for (index_type d = 0; d < DIM; ++d) {
					index_type cIdx = k + d * Q;
					cu += v[cIdx] * velocity[i][d];
				}
				float_type feq =
					w[k] * rho[i] * (1.0_fp + 3.0_fp * cu + 4.5_fp * cu * cu - 1.5_fp * usqr);
				/**
				 * end computeEquilibrium
				*/

				this->f_new[idx] = this->f[idx] - (this->f[idx] - feq) / tau;
			}
		}
		/**
		 *  end collide();
		*/

		/**
		 *   stream();
		*/
#pragma omp parallel for
		for (index_type i = 0; i < mesh->total_nodes; ++i) {
			for (index_type k = 0; k < Q; ++k) {

				// Layout agnostic index based on LayoutPolicy
				// PULL : source pre-computed, destination using policy
				index_type src = mesh->prev_neighbour_index[i][k];
				index_type dst = LAYOUT::getIndex(k, i, Q, total_nodes);
				f[dst] = f_new[src];

				// PUSH : destination pre-computed, source using policy
				// index_type src =  LAYOUT::getIndex(k, i, Q, total_nodes); // k * total_nodes + i;
				// index_type dst =  mesh->next_neighbour_index[i][k];;
				// f[dst] = f_new[src];

				// todo: (mem) (2) inplace streaming, AA pattern
				//!!! check race conditions, none exist currently
			}
		}
		/**
		 *  end stream();
		*/
		current_step++;
		current_time += dt;
	}

	void setFAt(index_type k, index_type idx, float_type value) {

		// Layout agnostic index
		f[LAYOUT::getIndex(k, idx, LATTICE::Q, total_nodes)] = value;
	}

	void computeVorticity() {
		static_assert(DIM == 2,
									"Vorticity computation currently only implemented for 2D");

#pragma omp parallel for
		for (index_type i = 0; i < mesh->total_nodes; ++i) {
			auto pos = mesh->nodes[i];

			// todo: check for negatives at boundary
			// Neighbor positions with periodic boundaries
			// replace with mesh->nodes[mesh->next_neighbour_index[i][1,2,3,4]]
			std::array<index_type, 2> pos_xp = {(pos[0] + 1) % domain_size[0],
																					pos[1]};
			std::array<index_type, 2> pos_xm = {
				(pos[0] - 1 + domain_size[0]) % domain_size[0], pos[1]};
			std::array<index_type, 2> pos_yp = {pos[0],
																					(pos[1] + 1) % domain_size[1]};
			std::array<index_type, 2> pos_ym = {
				pos[0], (pos[1] - 1 + domain_size[1]) % domain_size[1]};

			index_type idx_xp = mesh->getNodeIndex(pos_xp);
			index_type idx_xm = mesh->getNodeIndex(pos_xm);
			index_type idx_yp = mesh->getNodeIndex(pos_yp);
			index_type idx_ym = mesh->getNodeIndex(pos_ym);

			// Compute derivatives using central differences
			float_type duydx = (velocity[idx_xp][1] - velocity[idx_xm][1]) * 0.5;
			float_type duxdy = (velocity[idx_yp][0] - velocity[idx_ym][0]) * 0.5;

			vorticity[i] = duydx - duxdy;
		}
	}

	// todo: handle 3D  as well
	void saveVelocityField(const std::string& filename) {
		static_assert(DIM == 2, "File output currently only implemented for 2D");

		std::ofstream outFile(filename);
		for (index_type y = 0; y < domain_size[1]; ++y) {
			for (index_type x = 0; x < domain_size[0]; ++x) {
				std::array<index_type, 2> pos = {x, y};
				index_type idx = mesh->getNodeIndex(pos);

				float_type vel_mag = 0.0;
				for (index_type d = 0; d < DIM; ++d) {
					vel_mag += velocity[idx][d] * velocity[idx][d];
				}
				vel_mag = std::sqrt(vel_mag);

				outFile << vel_mag << " ";
			}
			outFile << "\n";
		}
		outFile.close();
	}

	void saveVelocitySlice2D(const index_type slice_location,
													 const index_type slice_axis,
													 const std::string& filename) {

		if (slice_axis > 2) [[unlikely]] {
			std::cerr << "Error: Slice Axis out of bounds, please choose 0,1,2. ";
			return;
		}
		std::ofstream outFile(filename);
		std::array<index_type, 2> pos = {};

		if (slice_axis == 0) {
			pos[1] = slice_location;
		} else
			pos[0] = slice_location;

		for (index_type xi = 0; xi < domain_size[slice_axis]; ++xi) {
			pos[slice_axis] = xi;
			index_type idx = mesh->getNodeIndex(pos);

			// save both components of velocity
			float_type normal_pos =
				static_cast<float_type>(xi) / static_cast<float_type>(domain_size[0] - 1);
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
		for (index_type y = 0; y < domain_size[1]; ++y) {
			for (index_type x = 0; x < domain_size[0]; ++x) {
				std::array<index_type, 2> pos = {x, y};
				index_type idx = mesh->getNodeIndex(pos);
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

	void setViscosity(float_type viscosity) {
		nu = viscosity;
		tau = 3.0_fp * nu + 0.5_fp;
	}

	void setRelaxationTime(float_type relaxation_time) {
		tau = relaxation_time;
		nu = (tau - 0.5_fp) / 3.0_fp;
	}

	float_type getViscosity() const { return nu; }

	float_type getRelaxationTime() const { return tau; }

	// Initialization helpers for external setups
	index_type getTotalNodes() const { return mesh->total_nodes; }

	index_type getQ() const { return Q; }

	float_type getLatticeSpeedofSound() const { return LATTICE::Cs; }

	std::array<index_type, DIM> getLatticeSize() const { return domain_size; }

	std::array<index_type, DIM> getPositionFromIndex(
		index_type flat_index) const {
		return mesh->nodes[flat_index];
	}

	// class member functions are default inline
	void setRhoAtIndex(index_type idx, float_type rho_val) { rho[idx] = rho_val; }

	float_type getRhoAtIndex(index_type idx) const { return rho[idx]; }

	const std::array<float_type, DIM>& getVelocityAtIndex(index_type idx) const {
		return velocity[idx];
	}

	void setVelocityAtIndex(index_type idx, const std::array<float_type, DIM>& u) {
		velocity[idx] = u;
	}

	float_type computeEquilibriumForInit(index_type k, float_type rho_val,
																	 const std::array<float_type, DIM>& u) const {
		return computeEquilibrium(k, rho_val, u);
	}
};

#endif