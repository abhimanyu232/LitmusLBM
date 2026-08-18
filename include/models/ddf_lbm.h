#ifndef DDF_LBM_H
#define DDF_LBM_H

#include <format>
#include "../policy/layout_policy.h"
#include "models.h"

template <LatticeType F_Lattice, LatticeType G_Lattice, LayoutPolicy Layout>
class DDF {

	private:
	//!!! todo: (perf) : causes pointer dereferencing in the hot loop
	// consider storing it as a simple Mesh member. // !!! (profile)
	Mesh<Lattice::Dim, Lattice, Layout> mesh;

	float_type nu, alpha;	 // viscosity, thermal diffusivity
	float_type omega_f, omega_g, beta_f, beta_g;
	// omega_f = 2.0 / (2.0*(nu/dt)/T + 1.0)    # for f-population
	// omega_g = 2.0 / (2.0*(alpha/dt)/T + 1.0) # for g-population (thermal)
	float_type gamma, PrandtlNo, Cv;
	// omega_g/omega_f = α / ν = 1/Pr
	index_type total_nodes;
	size_t current_step;
	double current_time;
	double dt;

	// Fields
	// Flattened storage arrays
	//!!! why are these dynamic vectors, the depend only on above static variables.

 public:
	// todo: make these short or uint_8t
	static constexpr index_type Dim = Lattice::Dim;
	static constexpr index_type Q = Lattice::Q;
	std::array<index_type, Dim> domain_size;


	std::vector<float_type> f, f_new, rho;
	std::vector<float_type> g, g_new, energy;
	std::vector<std::array<float_type, 3>> velocity;

	// std::vector<float_type> moments_f, moments_g;
};

#endif