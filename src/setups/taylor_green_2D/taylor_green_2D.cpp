#include "taylor_green_2D.h"
#include "runner.h"

// generic helper for sim runner
template <index_type Dim>
struct RunConfig {
	// Model Parameters
	std::array<index_type, Dim> domain_size = {500, 500};
	index_type max_steps = 1250;

	float_type char_length = 1;
	float_type dt = 1;

	// Simulation parameters
	float_type reynolds = 50;	 // = u0*L/nu
	float_type u0 = 0.0287;
	float_type rho0 = 1.;

	// Output and Diagnostics Parameters
	bool save_file = true;
	index_type save_interval = 1250;

	std::string case_name = "taylor_green_vortex_2d";
	std::string output_directory = "tg_result_fields";

	bool print_out = true;
	index_type print_interval = 100;

	bool enable_diagnostic = true;
	index_type diagnostic_interval = max_steps - 1;
};

// for each setup-header
class SimTraits {
 public:
	using layout_policy_t = AoSLayout;
	using lattice_t = D2Q9;
	using config_t = RunConfig<lattice_t::Dim>;
	using model_t = LBM<lattice_t, layout_policy_t>;
	using setup_t = TaylorGreenVortex_2D<model_t,config_t>;
};

int main() {
	std::cout << "Using " << omp_get_max_threads() << " OpenMP threads"
						<< '\n';

	SimTraits::config_t tg_config;							 // read config
	SimRunner<SimTraits> simulation(tg_config);	 // construct runner from config

	simulation.init();
	simulation.run();
	simulation.clean_up();

	return 0;
}