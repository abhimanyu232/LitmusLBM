#ifndef RUNNER_H
#define RUNNER_H

#include <format>
#include "common.h"
#include "setups.h"

// generic simulation runner
template <typename Traits>
struct SimRunner {
	using model_t = typename Traits::model_t;
	using config_t = typename Traits::config_t;
	using writer_t = typename Traits::writer_t;
	using setup_t = typename Traits::setup_t;

	config_t config;
	model_t model;
	setup_t setup;	// default values defined in taylor_green2d.h

	SimRunner(config_t _config)
			: config(_config),
				model(_config.domain_size, _config.dt),
				setup(_config) {
		// refactor model ctor to take a sim_params object to initalilze the domain size write freq etc.
		// refactor model ctor to take a setups object ( say of type TG),
		// it then calls setup.init() from model.init()
		static_assert(KineticModel<model_t, writer_t>);
	}

	void init() {
		setup.initialize(model);	// sets the fields and viscosity
		model.init();							// sets the model parameters : relxation etc

		// create folder to save results
		if (config.save_file) {
			std::filesystem::path dir(config.output_directory);
			if (!std::filesystem::exists(dir)) {
				if (std::filesystem::create_directory(dir)) {	 // default permissions
					std::cout << " Results directory created : " << dir << '\n';
				} else {
					std::cerr << "Error: Failed to create directory " << dir << '\n';
				}
			} else {
				std::cout << " Results directory already exists : " << dir << '\n';
			}
		}
	}

	// advance simulation
	void run() {
		std::cout << "Running " << model.getModelName() << " simulation..." << '\n';
		// mass conservation check
		float_type initial_mass = model.computeTotalMass();

		auto time_start = timer::GetCurrentTime();
		for (size_t step = 0; step < config.max_steps; ++step) {

			// stream-collide
			model.step();

			// save and write fields
			if (config.save_file && (step % config.save_interval == 0)) {

				std::filesystem::path save_file = config.output_directory + "/" +
																					config.case_name +
																					std::format("_{:08d}", step);

				writer_t writer(save_file, config, static_cast<std::uint64_t>(step),
												model.getCurrentTime());
				model.write_fields(writer);

				std::cout << "Step " << step << " (t=" << model.getCurrentTime()
									<< ") completed" << '\n';
			}

			if constexpr (HasDiagnostics<setup_t, model_t>) {
				if (config.enable_diagnostic && config.diagnostic_interval != 0 &&
						step % config.diagnostic_interval == 0) {
					setup.diagnose(model);	// diagnose(model, info)
				}
			}
		}
		auto time_end = timer::GetCurrentTime();
		// mass conservation check
		float_type final_mass = model.computeTotalMass();

		std::cout << "Mass Conservation Check. Mass Error = "
							<< std::fabs(final_mass - initial_mass) / initial_mass << '\n';
		// todo: call additional diagnostics.post_run()

		std::cout << "Simulation completed. Final step: " << model.getCurrentStep()
							<< ", Final time: " << model.getCurrentTime() << '\n';

		auto elapsed_time = time_end - time_start;
		std::cout << "total time elapsed : " << elapsed_time << '\n';
	}

	void clean_up() {
		// construct optional diagnostics/output observers
		// invoke observers at configured times
	}
};

#endif