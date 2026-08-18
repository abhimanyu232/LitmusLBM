#pragma once

#include <omp.h>
// #include <print>

// LBM includes
#include "common.h"

// LBM includes
#include "mesh.h"
#include "models/lbm.h"

// todo: need class members, initial values of macro quantities etc, useful for calculations later.
// possibly here: BC_type BoundaryCondition = NoBC
template <typename ModelType, typename RunConfig, typename WriterType>
requires IsothermalModel<ModelType,WriterType>
class TaylorGreenVortex_2D {
 public:
	// set-up from Kallikounis Thesis : Re= 50 = [u0*L/nu]
	float_type reynolds = 50;	 // = u0*L/nu
	float_type u0 = 0.0287;
	float_type rho0 = 1.;

	TaylorGreenVortex_2D() = default;

	TaylorGreenVortex_2D(const RunConfig& config)
			: reynolds(config.reynolds), u0(config.u0), rho0(config.rho0) {}

	TaylorGreenVortex_2D(float_type _reynolds, float_type _u0, float_type _rho0)
			: reynolds(_reynolds), u0(_u0), rho0(_rho0) {}

	void initialize(ModelType& model) {
		const auto size = model.getLatticeSize();
		const index_type total_nodes = model.getTotalNodes();
		//!!! why am I calling a runtime function for compile time info???
		const index_type Q = model.getF_Q();

		// set model kinematic viscosity
		model.setViscosity(u0 * size[0] / reynolds);
		std::cout << " Viscosity Set : nu = " << model.getViscosity() << '\n';

		float_type cs2_inv =
			1 / (model.getFLatticeSpeedofSound() * model.getFLatticeSpeedofSound());

		for (index_type i = 0; i < total_nodes; ++i) {
			auto pos = model.getPositionFromIndex(i);
			float_type x_pos = static_cast<float_type>(pos[0]) / (size[0]);	 // - 1)
			float_type y_pos = static_cast<float_type>(pos[1]) / (size[1]);	 // - 1)
			std::array<float_type, 3> u{};
			u[0] = -u0 * cos(2 * PI * x_pos) * sin(2 * PI * y_pos);
			u[1] = u0 * sin(2 * PI * x_pos) * cos(2 * PI * y_pos);
			model.setVelocityAtIndex(i, u);

			// p = cs2 * rho  //  p_amp = rho0 * u0 * u0 / 4;
			float_type density = rho0 - 0.25 * cs2_inv * (u0 * u0) *
																		(cos(4 * PI * x_pos) + cos(4 * PI * y_pos));
			model.setRhoAtIndex(i, density);

			for (index_type k = 0; k < Q; ++k) {
				float_type feq = model.computeFEquilibrium(k, density, u);
				model.setFAt(k, i, feq);
			}
		}
	}

	void diagnose(const ModelType& model) { compute_error(model); };

	// void apply_bc(const ModelType& model) {
	// 	if constexpr (BoundaryCondition::type != NoBC::type) {
	// 		bc.apply(model);
	// 	}
	// }

 private:
	void compute_error(const ModelType& model, Axis slice_axis = Axis::X,
										 float_type slice_norm_position = 0.125) {

		assert(slice_axis != Axis::Z &&
					 "Slicing along the Z-axis is not supported in 2D mode.");

		float_type time = model.getCurrentTime();
		const index_type total_nodes = model.getTotalNodes();
		const auto size = model.getLatticeSize();

		// also compute a 2D slice along x-axis at the 1/8 length of y-axis.
		index_type slice_position = static_cast<int>(
			size[static_cast<std::size_t>(Axis::Y)] * slice_norm_position);

		const std::string filename_slice =
			"slice2D_error" + std::to_string(slice_position) + "_" +
			std::to_string(static_cast<int>(time)) + ".dat";
		std::ofstream ofile_slice(filename_slice);
		ofile_slice << "Xpos" << '\t' << "Ux_{Exact}" << ' ' << "Ux_{LBM}" << '\t'
								<< "Uy_{Exact}" << ' ' << "Uy_{LBM}" << '\n';

		// float_type u0 = 0.0287;
		// float_type rho0 = 1.;	 // valid for standard LB
		float_type cs2 =
			model.getFLatticeSpeedofSound() * model.getFLatticeSpeedofSound();
		float_type p0 = cs2 * rho0;
		float_type p_amp = rho0 * u0 * u0 / 4;
		// float_type reynolds = 50;
		float_type nu = model.getViscosity();	 // kinematic viscosity

		const float_type k_squared = 4 * (PI * PI) / (size[0] * size[0]);
		const float_type exp_K_time = std::exp(-2 * k_squared * nu * time);

		float_type ux_sq_err_runner{0.}, uy_sq_err_runner{0.};
		float_type ux_exact_runner{0.}, uy_exact_runner{0.};
		float_type L2_err_ux_total{0.}, L2_err_uy_total{0.};

		float_type slice_ux_sq_err_runner{0.}, slice_uy_sq_err_runner{0.};
		float_type slice_ux_exact_runner{0.}, slice_uy_exact_runner{0.};
		float_type L2_err_ux_slice{0.}, L2_err_uy_slice{0.};

		// #pragma omp parallel for schedule(static)
		for (index_type i = 0; i < total_nodes; ++i) {
			auto u_lbm = model.getVelocityAtIndex(i);
			float_type ux_exact, uy_exact;
			float_type p_exact;

			auto pos = model.getPositionFromIndex(i);
			float_type x_pos = static_cast<float_type>(pos[0]) / (size[0]);	 // - 1)
			float_type y_pos = static_cast<float_type>(pos[1]) / (size[1]);	 // - 1)

			ux_exact = -u0 * cos(2 * PI * x_pos) * sin(2 * PI * y_pos) * exp_K_time;
			uy_exact = u0 * sin(2 * PI * x_pos) * cos(2 * PI * y_pos) * exp_K_time;

			p_exact = p0 - p_amp * (cos(4 * PI * x_pos) + cos(4 * PI * y_pos)) *
											 exp_K_time * exp_K_time;

			float_type ux_exact_sq = std::pow(ux_exact, 2);
			float_type uy_exact_sq = std::pow(uy_exact, 2);

			float_type local_ux_sq_err = std::pow(ux_exact - u_lbm[0], 2);
			float_type local_uy_sq_err = std::pow(uy_exact - u_lbm[1], 2);

			// calculate L2 vel error over slice
			if (pos[static_cast<std::size_t>(Axis::Y)] ==
					slice_position) {	 // pos[1] = y position
				slice_ux_sq_err_runner += local_ux_sq_err;
				slice_uy_sq_err_runner += local_uy_sq_err;
				slice_ux_exact_runner += ux_exact_sq;
				slice_uy_exact_runner += uy_exact_sq;

				// write slice to file
				float_type normal_pos =
					static_cast<float_type>(pos[0]) / static_cast<float_type>(size[0]);
				ofile_slice << normal_pos << '\t' << ux_exact << ' ' << u_lbm[0] << '\t'
										<< uy_exact << ' ' << u_lbm[1] << '\n';
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

		std::cout
			<< "Velocity RMS Error (Relative) over domain, (Ux_Err,Uy_Err) = ("
			<< L2_err_ux_total << ", " << L2_err_uy_total << ")" << '\n';

		std::cout << "Velocity RMS Error (Relative) over slice, (Ux_Err,Uy_Err) = ("
							<< L2_err_ux_slice << ", " << L2_err_uy_slice << ")" << '\n';
	}
};

// runner main