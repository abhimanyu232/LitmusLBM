
#ifndef MODELS_H
#define MODELS_H

#include "../lattices.h"

template <typename ModelType>
concept KineticModel =
	requires(ModelType& m, index_type idx, float_type val,
					 const std::array<float_type, ModelType::Dim>& u) {
		// --- Lifecycle ---
		{ m.init() } -> std::same_as<void>;
		{ m.step() } -> std::same_as<void>;
		{ m.computeVorticity() } -> std::same_as<void>;
		{ m.computeTotalMass() } -> std::same_as<float_type>;

		// --- State queries ---
		{ m.getModelName() } -> std::convertible_to<std::string>;
		{ m.getCurrentStep() } -> std::convertible_to<index_type>;
		{ m.getCurrentTime() } -> std::convertible_to<double>;

		// --- Output ---
		{
			m.saveVelocityField(std::declval<const std::string&>())
		} -> std::same_as<void>;
		{
			m.saveVorticityField(std::declval<const std::string&>())
		} -> std::same_as<void>;
		{
			m.saveVelocitySlice2D(idx, idx, std::declval<const std::string&>())
		} -> std::same_as<void>;

		// --- Spatial grid ---
		{ m.getTotalNodes() } -> std::same_as<index_type>;
		{
			m.getLatticeSize()
		} -> std::same_as<std::array<index_type, ModelType::Dim>>;
		{
			m.getPositionFromIndex(idx)
		} -> std::same_as<std::array<index_type, ModelType::Dim>>;

		// --- Hydrodynamic fields ---
		{ m.setRhoAtIndex(idx, val) } -> std::same_as<void>;
		{ m.getRhoAtIndex(idx) } -> std::same_as<float_type>;
		{ m.setVelocityAtIndex(idx, u) } -> std::same_as<void>;
		{
			m.getVelocityAtIndex(idx)
		} -> std::same_as<const std::array<float_type, ModelType::Dim>&>;
	};

template <typename ModelType>
concept IsothermalModel =
	KineticModel<ModelType> &&
	requires(ModelType& m, index_type idx, float_type val,
					 const std::array<float_type, ModelType::Dim>& u) {
		// --- Single Population interface ---
		{ m.getF_Q() } -> std::same_as<index_type>;
		{ m.setFAt(idx, idx, val) } -> std::same_as<void>;
		{ m.computeFEquilibrium(idx, val, u) } -> std::same_as<float_type>;

		// --- Lattice property ---
		{ m.getFLatticeSpeedofSound() } -> std::same_as<float_type>;

		// --- Transport coefficients ---
		{ m.setViscosity(val) } -> std::same_as<void>;
		{ m.getViscosity() } -> std::same_as<float_type>;
		{ m.setRelaxationTime(val) } -> std::same_as<void>;
		{ m.getRelaxationTime() } -> std::same_as<float_type>;
	};

template <typename ModelType>
concept ThermalModel =
	IsothermalModel<ModelType> &&
	requires(ModelType& m, index_type idx, float_type val,
					 const std::array<float_type, ModelType::Dim>& u) {
		// --- Thermal population interface ---
		{ m.getG_Q() } -> std::same_as<index_type>;
		{ m.setGAt(idx, idx, val) } -> std::same_as<void>;
		{ m.computeGEquilibrium(idx, val, u) } -> std::same_as<float_type>;

		// --- Thermal lattice property ---
		{ m.getGLatticeSpeedofSound() } -> std::same_as<float_type>;

		// --- Thermal transport coefficients ---
		{ m.setThermalDiffusivity(val) } -> std::same_as<void>;
		{ m.getThermalDiffusivity() } -> std::same_as<float_type>;
		{ m.setPrandtlNumber(val) } -> std::same_as<void>;
		{ m.getPrandtlNumber() } -> std::same_as<float_type>;

		// --- Thermal fields ---
		{ m.setTemperatureAtIndex(idx, val) } -> std::same_as<void>;
		{ m.getTemperatureAtIndex(idx) } -> std::same_as<float_type>;
		{ m.setEnergyAtIndex(idx, val) } -> std::same_as<void>;
		{ m.getEnergyAtIndex(idx) } -> std::same_as<float_type>;
	};

#endif
