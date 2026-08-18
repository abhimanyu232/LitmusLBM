#ifndef LBM_H
#define LBM_H

#include <format>
#include "../policy/layout_policy.h"
#include "models.h"

//todo: possibly avoided or improved using concepts
// Standard LBM implementation inheriting from Models using CRTP
template <
	LatticeType Lattice,
	LayoutPolicy
		Layout>	 // possibly here: BC_type BoundaryCondition = NoBC (default)
class LBM {
 public:
	// todo: make these short or uint_8t
	static constexpr index_type Dim = Lattice::Dim;
	static constexpr index_type Q = Lattice::Q;
	static constexpr float_type Cs = Lattice::Cs;
	std::array<index_type, Dim> domain_size;

	// Fields
	// Flattened storage arrays
	//!!! why are these dynamic vectors, the depend only on above static variables.
	std::vector<float_type> f, f_new, rho;
	// velocity is default 3D irrespective of spatial dim : for VTK compatibility
	std::vector<std::array<float_type, 3>> velocity = {};

 private:
	index_type total_nodes;
	float_type nu = 0.1_fp, tau = 0.8_fp, beta = 0.625_fp ;
	// float_type gamma = 2_fp, PrandtlNo = 1, Cv = 1 ;
	index_type current_step = 0;
	double current_time = 0;
	double dt = 1.0;

	//!!! todo: (perf) : causes pointer dereferencing in the hot loop
	// consider storing it as a simple Mesh member. // !!! (profile)
	// std::unique_ptr<Mesh<Lattice::Dim, Lattice, Layout>> mesh;
	Mesh<Lattice::Dim, Lattice, Layout> mesh;

 public:
	// Note: In 3D, vorticity would be a vector field

	LBM(const std::array<index_type, Dim>& _size,
			float_type _dt)	 //, float_type viscosity)
			: domain_size(_size), dt(_dt), mesh(_size) {
		// Calculate total number of nodes
		total_nodes = 1;
		for (index_type d = 0; d < Dim; ++d) {
			total_nodes *= domain_size[d];
		}

		// mesh = std::make_unique<Mesh<Lattice::Dim, Lattice, Layout>>(domain_size);
		// mesh.writeToFile();

		// Initialize fields storage
		// !!! todo: (mem) Q * total_nodes should fit inside uint32_t for domain < 500
		f.resize(Q * total_nodes);
		f_new.resize(Q * total_nodes);
		rho.resize(total_nodes);
		// energy.resize(total_nodes);
		velocity.resize(total_nodes);

		// todo: now initisalize the fields using the SETUP::init() class call

		std::cout << "Initialized " << Lattice::getName() << " LBM with "
							<< total_nodes << " nodes" << '\n';
	}

	// todo: should also ideally call the setup->init() function here
	void init() {
		current_step = 0;
		current_time = 0.0;

		//todo: error checking here, check for tau > 0.5 (nu < 0) for stability
		tau = 3.0 * nu + 0.5;
		beta = 1 / (2 * tau);
		if (tau <= 0.5) {
			std::cerr << "Error: tau must be > 0.5 for stability (nu > 0). Got tau="
								<< tau << '\n';
			std::exit(EXIT_FAILURE);
		}

		// should compute the relaxation etc from the given viscosity etc
		// must be called after the setup(eg tg-vortex) has initialized fields
		//!!! call SETUP::init(); // this should be a call to taylor_green.init() etc
		//!!! should be called by constructor of the lbm class
	}

	[[nodiscard]] float_type computeTotalMass() {
		float_type mass_counter = 0.;
#pragma omp parallel for reduction(+ : mass_counter)
		for (index_type idx = 0; idx < total_nodes; ++idx) {
			mass_counter += rho[idx];
		}
		return mass_counter;
	}

	float_type computeFEquilibrium(index_type k, float_type rho_val,
																 const std::array<float_type, 3>& u_vel) const {

		const auto& v = Lattice::velocities;
		const auto& w = Lattice::weights;

		float_type cu = 0.0;
		float_type usqr = 0.0;
		for (index_type d = 0; d < Dim; ++d) {
			index_type cIdx = k + d * Q;	// !!! cache : strided access
			cu += v[cIdx] * u_vel[d];
			// cu += v[k][d] * u[d];
			usqr += u_vel[d] * u_vel[d];
		}
		return w[k] * rho_val *
					 (1.0_fp + 3.0_fp * cu + 4.5_fp * cu * cu - 1.5_fp * usqr);
	}

	void step() {
		const auto& v = Lattice::velocities;
		const auto& w = Lattice::weights;

#pragma omp parallel
		{
			/**
		 *   stream();
		*/
#pragma omp for schedule(static)
			for (index_type i = 0; i < total_nodes; ++i) {
				// reset
				float_type rho_acc = 0.0;
				std::array<float_type, Dim> u_acc{};
				// float_type energy_acc = 0.0;

				for (index_type k = 0; k < Q; ++k) {
					// Layout agnostic index based on LayoutPolicy
					// PULL : source pre-computed, destination using policy
					// !!! if src is missing (boundary), -> boundary condition enters here:
					// !!! bb=reflection, immersed=stream as usual, grads=equilibrium for missing pops
					index_type src = mesh.prev_neighbour_index[i][k];
					index_type dst = Layout::getIndex(k, i, Q, total_nodes);
					f[dst] =
						f_new[src];	 // !!! for grads: reconstruct missing = f_grads_eq[src]

					float_type fk = f[dst];
					// macro computation
					// float_type temp_v2;
					rho_acc += fk;
					for (index_type d = 0; d < Dim; ++d) {
						index_type cIdx = k + d * Q;
						u_acc[d] += v[cIdx] * fk;	 // v[k][d]
																			 // temp_v2 += v[cIdx] * v[cIdx];
					}
					// todo: (mem) (2) inplace streaming, AA pattern
				}
				// normalize
				rho[i] = rho_acc;
				if (rho_acc > 0.0) {
					for (index_type d = 0; d < Dim; ++d) {
						velocity[i][d] = u_acc[d] / rho_acc;
					}
					// energy[i] = energy_acc;
					// temperature = (2 * energy_acc - (dot(velocity[i], velocity[i]))) / Dim;
				} else [[unlikely]] {
					std::cerr << " negative density \n ";
					// std::exit(EXIT_FAILURE);
					// todo: save state, and exit after the loop ends.
					// todo: example: bool sim_error = true; later check sim_error == true and then std::exit
					// std::exit inside parallel worker region is UB for openmp
				}
			}
			/**
		 *  end stream();
		*/
			// //!!! todo: apply boundary condition
			// if constexpr (BoundaryCondition::type != NoBC::type) {
			// 	bc.apply(model);
			// }

			/**
		 * collide();
		*/
#pragma omp for schedule(static)
			for (index_type i = 0; i < total_nodes; ++i) {

				float_type usqr{0.};
				for (index_type d = 0; d < Dim; ++d) {	// dot(u,u)
					usqr += velocity[i][d] * velocity[i][d];
				}

				for (index_type k = 0; k < Q; ++k) {

					// Layout agnostic index based on LayoutPolicy
					index_type idx = Layout::getIndex(k, i, Q, total_nodes);

					/**
				 * computeEquilibrium
         */
					float_type cu{0.0};
					for (index_type d = 0; d < Dim; ++d) {
						index_type cIdx = k + d * Q;
						cu += v[cIdx] * velocity[i][d];
					}

					/**
           * - Standard second-order Hermite expansion of the Maxwell-Boltzmann distribution
           * - Only first 2 Hermite moments are matched (ρ and ρu)
           * - stress tensor P_{ij} = Σ c_{k,i} c_{k,j} f_k^eq  correct only up to O(Ma²) due to the truncation
           * - purely isothermal
           */
					float_type feq =
						w[k] * rho[i] *
						(1.0_fp + 3.0_fp * cu + 4.5_fp * cu * cu - 1.5_fp * usqr);

					/**
          * end computeEquilibrium
  				*/

					/** 
           * - collision ≡ f[k] − (f[k]−f_eq)/τ 
           * - single relaxation rate
           * - higher-order moments relax at same rate, i.e no control over numerical hyperviscosity
           * - for 2D, Shear viscosity ν = cs²·(τ − 0.5) == bulk viscosity ζ
					 * - for 3D, Shear viscosity ν =  2/3 * ζ, /bulk viscosity 
           * - Pr = 1 ; since viscosity = thermal diffusivity
					 * - BGK cannot independently control bulk viscosity
           */
					this->f_new[idx] = this->f[idx] - 2 * beta * (this->f[idx] - feq);
				}
			}
			/**
		 *  end collide();
		 */
		}

		current_step++;
		current_time += dt;
	}

	inline void setFAt(index_type k, index_type idx, float_type value) {
		// Layout agnostic index
		f_new[Layout::getIndex(k, idx, Lattice::Q, total_nodes)] = value;
	}

	template <FieldWriter Writer_t>
	void write_fields(Writer_t& writer) const {
		writer.write_scalar_double_dataset("Density", rho);
		writer.write_vector_double_dataset("Velocity", velocity);
	}

// 	void computeVorticity() {
// 		static_assert(Dim == 2,
// 									"Vorticity computation currently only implemented for 2D");

// #pragma omp parallel for
// 		for (index_type i = 0; i < total_nodes; ++i) {
// 			auto pos = mesh.nodes[i];

// 			// todo: check for negatives at boundary
// 			// Neighbor positions with periodic boundaries
// 			// replace with mesh.nodes[mesh.next_neighbour_index[i][1,2,3,4]]
// 			std::array<index_type, 2> pos_xp = {(pos[0] + 1) % domain_size[0],
// 																					pos[1]};
// 			std::array<index_type, 2> pos_xm = {
// 				(pos[0] - 1 + domain_size[0]) % domain_size[0], pos[1]};
// 			std::array<index_type, 2> pos_yp = {pos[0],
// 																					(pos[1] + 1) % domain_size[1]};
// 			std::array<index_type, 2> pos_ym = {
// 				pos[0], (pos[1] - 1 + domain_size[1]) % domain_size[1]};

// 			index_type idx_xp = mesh.getNodeIndex(pos_xp);
// 			index_type idx_xm = mesh.getNodeIndex(pos_xm);
// 			index_type idx_yp = mesh.getNodeIndex(pos_yp);
// 			index_type idx_ym = mesh.getNodeIndex(pos_ym);

// 			// Compute derivatives using central differences
// 			float_type duydx = (velocity[idx_xp][1] - velocity[idx_xm][1]) * 0.5;
// 			float_type duxdy = (velocity[idx_yp][0] - velocity[idx_ym][0]) * 0.5;

// 			vorticity[i] = duydx - duxdy;
// 		}
// 	}

	// todo: handle 3D  as well
	// void saveVelocityField(const std::string& filename) {
	// 	static_assert(Dim == 2, "File output currently only implemented for 2D");

	// 	std::ofstream outFile(filename);
	// 	for (index_type y = 0; y < domain_size[1]; ++y) {
	// 		for (index_type x = 0; x < domain_size[0]; ++x) {
	// 			std::array<index_type, 2> pos = {x, y};
	// 			index_type idx = mesh.getNodeIndex(pos);

	// 			float_type vel_mag = 0.0;
	// 			for (index_type d = 0; d < Dim; ++d) {
	// 				vel_mag += velocity[idx][d] * velocity[idx][d];
	// 			}
	// 			vel_mag = std::sqrt(vel_mag);

	// 			outFile << vel_mag << " ";
	// 		}
	// 		outFile << "\n";
	// 	}
	// 	outFile.close();
	// }

	// void saveVelocitySlice2D(const index_type slice_location,
	// 												 const index_type slice_axis,
	// 												 const std::string& filename,
	// 												 const index_type z_pos = 0) {

	// 	if (slice_axis >= Dim) [[unlikely]] {
	// 		std::cerr << "Error: Slice Axis out of bounds, please choose 0,1,2. ";
	// 		return;
	// 	}
	// 	std::ofstream outFile(filename);
	// 	std::array<index_type, Dim> pos = {};

	// 	if (slice_axis == 0) {
	// 		pos[1] = slice_location;
	// 	} else
	// 		pos[0] = slice_location;

	// 	if (Dim == 2)
	// 		pos[2] = z_pos;

	// 	for (index_type xi = 0; xi < domain_size[slice_axis]; ++xi) {

	// 		pos[slice_axis] = xi;
	// 		index_type idx = mesh.getNodeIndex(pos);

	// 		// save both components of velocity
	// 		float_type normal_pos = static_cast<float_type>(xi) /
	// 														static_cast<float_type>(domain_size[0] - 1);
	// 		outFile << normal_pos << '\t' << velocity[idx][0] << '\t'
	// 						<< velocity[idx][1] << '\n';
	// 	}
	// 	outFile << "\n";
	// 	outFile.close();
	// }

	// // todo: handle 3D  as well
	// void saveVorticityField(const std::string& filename) {
	// 	static_assert(Dim == 2, "File output currently only implemented for 2D");

	// 	std::ofstream outFile(filename);
	// 	for (index_type y = 0; y < domain_size[1]; ++y) {
	// 		for (index_type x = 0; x < domain_size[0]; ++x) {
	// 			std::array<index_type, 2> pos = {x, y};
	// 			index_type idx = mesh.getNodeIndex(pos);
	// 			outFile << vorticity[idx] << " ";
	// 		}
	// 		outFile << "\n";
	// 	}
	// 	outFile.close();
	// }

	inline std::string getModelName() const {
		return std::format("Standard LBM ({})", Lattice::getName());
	}

	inline index_type getCurrentStep() const { return current_step; }

	inline double getCurrentTime() const { return current_time; }

	void setViscosity(float_type viscosity) {
		nu = viscosity;
		tau = 3.0_fp * nu + 0.5_fp;
		beta = 1 / (2 * tau);
	}

	void setRelaxationTime(float_type relaxation_time) {
		tau = relaxation_time;
		beta = 1 / (2 * tau);
		nu = (tau - 0.5_fp) / 3.0_fp;
	}

	inline float_type getViscosity() const { return nu; }

	inline float_type getRelaxationTime() const { return tau; }

	// Initialization helpers for external setups
	inline index_type getTotalNodes() const { return total_nodes; }

	constexpr index_type getF_Q() const { return Q; }

	constexpr float_type getFLatticeSpeedofSound() const { return Lattice::Cs; }

	inline std::array<index_type, Dim> getLatticeSize() const {
		return domain_size;
	}

	inline std::array<index_type, Dim> getPositionFromIndex(
		index_type flat_index) const {
		return mesh.nodes[flat_index];
	}

	// class member functions are default inline
	inline void setRhoAtIndex(index_type idx, float_type rho_val) {
		rho[idx] = rho_val;
	}

	inline float_type getRhoAtIndex(index_type idx) const { return rho[idx]; }

	inline const std::array<float_type, 3>& getVelocityAtIndex(
		index_type idx) const {
		return velocity[idx];
	}

	inline void setVelocityAtIndex(index_type idx,
																 const std::array<float_type, 3>& u) {
		velocity[idx] = u;
	}
};

#endif
