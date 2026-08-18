#include "taylor_green_2D.h"
#include "writer/hdf5_writer.h"
#include "runner.h"

// generic helper for sim runner
template <index_type Dim>
struct RunConfig {
	// Model Parameters
	static constexpr index_type dim = Dim;
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
	index_type save_interval = 125;

	std::string case_name = "tg_2D";
	std::string output_directory = "tg2D_results";

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
	// using bc_t = IBB_BC;
	using model_t = LBM<lattice_t, layout_policy_t>;
	using config_t = RunConfig<lattice_t::Dim>;
	using writer_t = Hdf5VTKWriter<config_t::dim, config_t>;
	using setup_t = TaylorGreenVortex_2D<model_t, config_t,writer_t>;	// bc_t
};

int main() {
	std::cout << "Using " << omp_get_max_threads() << " OpenMP threads" << '\n';

	try {
		SimTraits::config_t tg_config;							 // read config
		SimRunner<SimTraits> simulation(tg_config);	 // construct runner from config

		simulation.init();
		simulation.run();
		simulation.clean_up();
	} catch (const std::exception& e) {
		std::cerr << "Error encountered: " << e.what() << "\n";
		return EXIT_FAILURE;
	}

	return 0;
}