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
class TaylorGreenVortex_3D {
 public:
	// set-up from Kallikounis Thesis : Re= 50 = [u0*L/nu]
	float_type reynolds = 50;	 // = u0*L/nu
	float_type u0 = 0.0287;
	float_type rho0 = 1.;

	TaylorGreenVortex_3D(const RunConfig& config)
			: reynolds(config.reynolds), u0(config.u0), rho0(config.rho0) {}

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
			// float_type z_pos = static_cast<float_type>(pos[2]) / (size[2]);	 // - 1)

			std::array<float_type, 3> u{};
			u[0] = -u0 * cos(2 * PI * x_pos) * sin(2 * PI * y_pos);
			u[1] = u0 * sin(2 * PI * x_pos) * cos(2 * PI * y_pos);
			u[2] = 0;

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

	// void diagnose(const ModelType& model) { compute_error(model); };

	// void apply_bc(const ModelType& model) {
	// 	if constexpr (BoundaryCondition::type != NoBC::type) {
	// 		bc.apply(model);
	// 	}
	// }

 private:
	void compute_error(const ModelType& model, Axis slice_axis = Axis::X,
										 float_type slice_norm_position = 0.125) {}
};