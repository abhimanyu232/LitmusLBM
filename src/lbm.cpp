#include <omp.h>
// #include <print>

// LBM includes
#include "../include/common.h"
#include "../include/mesh.h"
#include "../include/models/lbm.h"
#include "../include/policy/layout_policy.h"

// todo: need class members, initial values of macro quantities etc, useful for calculations later.
//!!! maybe should derive from the model and not have it as a member.
//!!! that simplifies the initialisation.
// concept ModelType requires { requires model.step(); }
template <typename Derived, typename ModelType>
class Setups {
 public:
	ModelType& model;

	Setups(ModelType& m) : model(m) {}

	void initialize() { static_cast<Derived*>(this)->initialize(); }

	// not sure if it is a good idea doing this here
	// std::array<double, DIM> u0{0.};
	// double rho0 = 1.0;
	// double T0{model.getLatticeSpeedofSound()*model.getLatticeSpeedofSound()};
	// double p0 = rho0*T0;
};

// todo: need class members, initial values of macro quantities etc, useful for calculations later.
template <IsothermalModel ModelType>
class TaylorGreenVortex
		: public Setups<TaylorGreenVortex<ModelType>, ModelType> {

 public:
	using Setups<TaylorGreenVortex<ModelType>, ModelType>::model;

	explicit TaylorGreenVortex(ModelType& m)
			: Setups<TaylorGreenVortex<ModelType>, ModelType>(m) {
		initialize();
	}

	void initialize() {
		const auto size = model.getLatticeSize();
		const index_type total_nodes = model.getTotalNodes();
		const index_type Q = model.getF_Q(); //!!! why am I calling a runtime function for compile time info???
		// set-up from Kallikounis Thesis : Re= 50 = [u0*L/nu]
		float_type rho0 = 1.;	 // valid for standard LB
		float_type u0 = 0.0287;
		float_type reynolds = 50.;
		model.setViscosity(u0 * size[0] / reynolds);	// kinematic viscosity
		std::cout << " Viscosity Set : nu = " << model.getViscosity() << std::endl;
		for (index_type i = 0; i < total_nodes; ++i) {
			auto pos = model.getPositionFromIndex(i);
			float_type x_pos = static_cast<float_type>(pos[0]) / (size[0] - 1);
			float_type y_pos = static_cast<float_type>(pos[1]) / (size[1] - 1);
			std::array<float_type, 2> u;
			u[0] = -u0 * cos(2 * PI * x_pos) * sin(2 * PI * y_pos);
			u[1] = u0 * sin(2 * PI * x_pos) * cos(2 * PI * y_pos);
			model.setVelocityAtIndex(i, u);

			float_type density =
				rho0 - 0.25 * (u0 * u0) * (cos(2 * PI * x_pos) + cos(2 * PI * y_pos));
			model.setRhoAtIndex(i, density);

			for (index_type k = 0; k < Q; ++k) {
				float_type feq = model.computeFEquilibrium(k, density, u);
				model.setFAt(k, i, feq);
			}
		}
	}

	void compute_error(double time) {

		const index_type total_nodes = model.getTotalNodes();
		const auto size = model.getLatticeSize();

		// also compute a 2D slice along x-axis at the middle-ish plane.
		index_type slice_y = static_cast<int>(size[0] / 2);
		const std::string filename_slice =
			"slice2D" + std::to_string(slice_y) + "_" +
			std::to_string(static_cast<int>(time)) + ".dat";
		std::ofstream ofile_slice(filename_slice);
		ofile_slice << "Xpos" << '\t' << "Ux_{Exact}" << ' ' << "Ux_{LBM}" << '\t'
								<< "Uy_{Exact}" << ' ' << "Uy_{LBM}" << std::endl;

		float_type lattice_ref_temp =
			model.getFLatticeSpeedofSound() * model.getFLatticeSpeedofSound();
		float_type p0 = 1.0 / lattice_ref_temp;	 // valid for standard LB
		float_type u0 = 0.0287;
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

#pragma omp parallel for schedule(static)
		for (index_type i = 0; i < total_nodes; ++i) {
			auto u_lbm = model.getVelocityAtIndex(i);
			float_type ux_exact, uy_exact;
			float_type p_exact;

			auto pos = model.getPositionFromIndex(i);
			float_type x_pos = static_cast<float_type>(pos[0]) / (size[0] - 1);
			float_type y_pos = static_cast<float_type>(pos[1]) / (size[1] - 1);

			ux_exact = -u0 * cos(2 * PI * x_pos) * sin(2 * PI * y_pos) * exp_K_time;
			uy_exact = u0 * sin(2 * PI * x_pos) * cos(2 * PI * y_pos) * exp_K_time;

			p_exact = p0 - 0.25 * (u0 * u0) *
											 (cos(2 * PI * x_pos) + cos(2 * PI * y_pos)) * exp_K_time;

			float_type ux_exact_sq = std::pow(ux_exact, 2);
			float_type uy_exact_sq = std::pow(uy_exact, 2);

			float_type local_ux_sq_err = std::pow(ux_exact - u_lbm[0], 2);
			float_type local_uy_sq_err = std::pow(uy_exact - u_lbm[1], 2);

			// calculate L2 vel error over slice
			if (pos[1] == slice_y) {
				slice_ux_sq_err_runner += local_ux_sq_err;
				slice_uy_sq_err_runner += local_uy_sq_err;
				slice_ux_exact_runner += ux_exact_sq;
				slice_uy_exact_runner += uy_exact_sq;

				// write slice to file
				float_type normal_pos = static_cast<float_type>(pos[0]) /
																static_cast<float_type>(size[0] - 1);
				ofile_slice << normal_pos << '\t' << ux_exact << ' ' << u_lbm[0] << '\t'
										<< uy_exact << ' ' << u_lbm[1] << std::endl;
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

		// model.saveVelocitySlice2D(
		// 	(slice_y), 0, "vel_slice_y=" + std::to_string(slice_y) + ".txt");

		std::cout
			<< "Velocity RMS Error (Relative) over domain, (Ux_Err,Uy_Err) = ("
			<< L2_err_ux_total << ", " << L2_err_uy_total << ")" << std::endl;

		std::cout << "Velocity RMS Error (Relative) over slice, (Ux_Err,Uy_Err) = ("
							<< L2_err_ux_slice << ", " << L2_err_uy_slice << ")" << std::endl;
	}
};

// todo: need class members, initial values of macro quantities etc, useful for calculations later.
template <IsothermalModel ModelType>
class DoublePeriodicShearLayer
		: public Setups<DoublePeriodicShearLayer<ModelType>, ModelType> {

 public:
	using Setups<DoublePeriodicShearLayer<ModelType>, ModelType>::model;

	explicit DoublePeriodicShearLayer(ModelType& m)
			: Setups<DoublePeriodicShearLayer<ModelType>, ModelType>(m) {
		initialize();
	}

	void initialize() {
		const auto size = model.getLatticeSize();
		const index_type total_nodes = model.getTotalNodes();
		const index_type Q = model.getF_Q(); //!!! why am I calling a runtime function for compile time info???
		for (index_type i = 0; i < total_nodes; ++i) {
			auto pos = model.getPositionFromIndex(i);
			float_type x_pos = static_cast<float_type>(pos[0]) / (size[0] - 1);
			float_type y_pos = static_cast<float_type>(pos[1]) / (size[1] - 1);
			std::array<float_type, 2> u;

			// setup: https://www.researchgate.net/publication/335029866_Pseudoentropic_derivation_of_the_regularized_lattice_Boltzmann_method
			// https://www.researchgate.net/publication/266088157_Gibbs'_principle_for_the_lattice-kinetic_theory_of_fluid_dynamics
			float_type u0 = 0.04;
			float_type perturbation = 0.05 * u0;
			float_type lambda = 80;	 // 20;
			float_type reynolds = 30000.;
			// NOTE: only stable ish for size 200x200.
			// NOTE: Blows up at 160x160
			// !!! need entropic etc for higher reynolds numbers.
			model.setViscosity(u0 * size[0] / reynolds);
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
};

int main() {
	std::cout << "Using " << omp_get_max_threads() << " OpenMP threads"
						<< std::endl;

	// Create output directory
	std::string result_directory_name = "result_fields";
	std::filesystem::path dir(result_directory_name);
	if (!std::filesystem::exists(dir)) {
		if (std::filesystem::create_directory(dir)) {	 // default permissions
			std::cout << " Results directory created : " << dir << std::endl;
		} else {
			std::cerr << "Error: Failed to create directory " << dir << std::endl;
		}
	} else {
		std::cout << " Results directory already exists : " << dir << std::endl;
	}

	// Simulation parameters
	constexpr index_type DIM = 2;
	constexpr float_type scale = 25.;
	constexpr index_type LX = 500;	//16 * scale;
	constexpr index_type LY = 500;	//16 * scale;
	// timestep to do validation plots and error checks.
	bool SAVE_FILE = true;
	constexpr index_type SAVE_INTERVAL = 1250;

	// constexpr index_type NSTEPS = 6000;
	constexpr index_type checkpoint = 50 * scale;
	constexpr index_type NSTEPS = checkpoint + 1;

	std::array<index_type, DIM> domain_size = {LX, LY};

	//!!! todo: make it so that the viscosity is set by the setup, currently confusing
	constexpr float_type visc = 0.01;	 // can be overridden by the setup

	// Choose between structure of array or array of structs for population storage
	//AoSLayout;	// SoALayout;
	using layout_type = AoSLayout;
	using model_type = LBM<D2Q9, layout_type>;
	auto model = std::make_unique<model_type>(domain_size, visc);

	// Initialize the model :
	model->init();
	//!!! ideally the setup class should also do the domain init here.
	// !!! very strange way of doing things?
	// !!! we locally initialize the model by instantiating a setup object USING the model
	// shouldnt we probably use the setup class to initialize the model instead

	// Choose setup and initialize
	int test_case = 1;	// Taylor-Green Vortex
	// if (test_case == 1) {
	std::cout << "Initializing Taylor-Green Vortex..." << std::endl;
	TaylorGreenVortex<model_type> setup(*model);
	// }
	// else
	// {
	// int test_case = 2;  Double Periodic Shear Layer
	// 	std::cout << "Initializing Double Periodic Shear Layer..." << std::endl;
	// DoublePeriodicShearLayer<model_type> setup(*model);
	// }

	// Mass Conservation Check // cover with if-def DEBUG block
	float_type initial_mass = 0.;
	for (index_type i = 0; i < model->getTotalNodes(); ++i) {
		initial_mass += model->getRhoAtIndex(i);
	}
	// Mass Conservation Check // cover with if-def DEBUG block

	std::cout << "Running " << model->getModelName() << " simulation..."
						<< std::endl;

	auto time_start = timer::GetCurrentTime();
	for (size_t step = 0; step < NSTEPS; ++step) {
		// setup.compute_error(model->getCurrentTime());

		model->step();

		if (SAVE_FILE && (step % SAVE_INTERVAL == 0)) {
			model->saveVelocityField(result_directory_name + "/velocity_" +
															 std::to_string(step) + ".txt");

			model->computeVorticity();
			model->saveVorticityField(result_directory_name + "/vorticity_" +
																std::to_string(step) + ".txt");

			std::cout << "Step " << step << " (t=" << model->getCurrentTime()
								<< ") completed" << std::endl;
		}

		// run Taylor-Green Vortex validation, probably instead do a guard block here
		float_type time_non_dim = model->getCurrentTime() * 0.0287 / LX;
		//!!! magic number, should come from setup->u_init;
		if (test_case == 1 && (model->getCurrentStep() % checkpoint == 0)) {
			std::cout << "Step " << step << " (t=" << model->getCurrentTime()
								<< ") completed" << "; non Dim Time = " << time_non_dim
								<< std::endl;
			setup.compute_error(model->getCurrentTime());
		}
	}
	auto time_end = timer::GetCurrentTime();

	// Mass Conservation Check // cover with if-def DEBUG block
	float_type final_mass = 0.;
	for (index_type i = 0; i < model->getTotalNodes(); ++i) {
		final_mass += model->getRhoAtIndex(i);
	}
	std::cout << "Mass Conservation Check. Initial Mass = " << initial_mass
						<< std::endl;
	std::cout << "Mass Conservation Check. Total Mass = " << final_mass
						<< std::endl;
	// Mass Conservation Check // cover with if-def DEBUG block

	std::cout << "Simulation completed. Final step: " << model->getCurrentStep()
						<< ", Final time: " << model->getCurrentTime() << std::endl;

	auto elapsed_time = time_end - time_start;
	std::cout << "total time elapsed : " << elapsed_time << std::endl;

	return 0;
}