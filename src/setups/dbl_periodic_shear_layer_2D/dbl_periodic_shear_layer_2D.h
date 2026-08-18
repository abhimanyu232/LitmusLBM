#pragma once

#include <omp.h>
// #include <print>

// LBM includes
#include "common.h"

// LBM includes
#include "mesh.h"
#include "models/lbm.h"

template <typename ModelType, typename RunConfig, typename WriterType>
requires IsothermalModel<ModelType,WriterType>
class DoublePeriodicShearLayer {
 public:
	// setup: https://www.researchgate.net/publication/335029866_Pseudoentropic_derivation_of_the_regularized_lattice_Boltzmann_method
	// https://www.researchgate.net/publication/266088157_Gibbs'_principle_for_the_lattice-kinetic_theory_of_fluid_dynamics
	float_type reynolds = 30000.;
	float_type u0 = 0.04;
	float_type rho0 = 1;
	float_type lambda = 80;	 // 20;
	float_type perturbation = 0.05 * u0;

	DoublePeriodicShearLayer() = default;

	DoublePeriodicShearLayer(const RunConfig& config)
			: reynolds(config.reynolds),
				u0(config.u0),
				rho0(config.rho0),
				lambda(config.lambda) {
		perturbation = 0.05 * u0;
	}

	// diagnostics(config.output_directory){}

	DoublePeriodicShearLayer(float_type _reynolds, float_type _u0,
													 float_type _rho0, float_type _lambda)
			: reynolds(_reynolds), u0(_u0), rho0(_rho0), lambda(_lambda) {
		perturbation = 0.05 * u0;
	}

	void initialize(ModelType& model) {
		const auto size = model.getLatticeSize();
		const index_type total_nodes = model.getTotalNodes();
		const index_type Q = model.getF_Q();
		//!!! why am I calling a runtime function for compile time info???
		model.setViscosity(u0 * size[0] / reynolds);
		std::cout << " Viscosity Set : nu = " << model.getViscosity() << '\n';

		for (index_type i = 0; i < total_nodes; ++i) {
			auto pos = model.getPositionFromIndex(i);
			float_type x_pos = static_cast<float_type>(pos[0]) / (size[0]);
			float_type y_pos = static_cast<float_type>(pos[1]) / (size[1]);
			std::array<float_type, 3> u{};

			// NOTE: only stable ish for size 200x200.
			// NOTE: Blows up at 160x160
			// !!! need entropic etc for higher reynolds numbers.
			if (y_pos <= 0.5) {
				u[0] = u0 * tanh((y_pos - 0.25) * lambda);
			} else {
				u[0] = u0 * tanh((0.75 - y_pos) * lambda);
			}
			u[1] = perturbation * sin(2 * PI * (x_pos + 0.25));

			model.setVelocityAtIndex(i, u);
			model.setRhoAtIndex(i, 1.0);
			for (index_type k = 0; k < Q; ++k) {
				float_type feq = model.computeFEquilibrium(k, 1.0, u);
				model.setFAt(k, i, feq);
			}
		}
	}

	// void diagnose(const ModelType& model) {};

 private:
	void compute_entropy(const ModelType& model) {}
};