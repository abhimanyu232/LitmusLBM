
#ifndef MODELS_H
#define MODELS_H

#include "../lattices.h"

// CRTP base class for simulation models
template <typename Derived, LatticeType LATTICE>
class Models {
 public:
	// virtual ~Models() = default;

	void init() { static_cast<Derived*>(this)->init(); }

	void step() { static_cast<Derived*>(this)->step(); }

	void computeVorticity() { static_cast<Derived*>(this)->computeVorticity(); }

	// Output functions
	void saveVelocityField(const std::string& filename) {
		static_cast<Derived*>(this)->saveVelocityField(filename);
	}

	void saveVorticityField(const std::string& filename) {
		static_cast<Derived*>(this)->saveVorticityField(filename);
	}

	void saveVelocitySlice2D(index_type index, index_type axis,
													 const std::string& filename) {
		static_cast<Derived*>(this)->saveVelocitySlice2D(index, axis, filename);
	};

	// Getters for simulation state
	std::string getModelName() const {
		return static_cast<const Derived*>(this)->getModelName();
	}

	size_t getCurrentStep() const {
		return static_cast<const Derived*>(this)->getCurrentStep();
	}

	double getCurrentTime() const {
		return static_cast<const Derived*>(this)->getCurrentTime();
	}

	// Configuration functions
	void setViscosity(float_type nu) {
		static_cast<Derived*>(this)->setViscosity(nu);
	}

	void setRelaxationTime(float_type tau) {
		static_cast<Derived*>(this)->setRelaxationTime(tau);
	}

	float_type getViscosity() const {
		return static_cast<const Derived*>(this)->getViscosity();
	}

	float_type getRelaxationTime() const {
		return static_cast<const Derived*>(this)->getRelaxationTime();
	}

	// Initialization helpers for external setups
	index_type getTotalNodes() const {
		return static_cast<const Derived*>(this)->getTotalNodes();
	}

	index_type getQ() const { return static_cast<const Derived*>(this)->getQ(); }

	std::array<index_type, LATTICE::DIM> getLatticeSize() const {
		return static_cast<const Derived*>(this)->getLatticeSize();
	}

	std::array<index_type, LATTICE::DIM> getPositionFromIndex(
		index_type idx) const {
		return static_cast<const Derived*>(this)->getPositionFromIndex(idx);
	}

	void setRhoAtIndex(index_type idx, float_type rho) {
		static_cast<Derived*>(this)->setRhoAtIndex(idx, rho);
	}

	void setVelocityAtIndex(index_type idx,
													const std::array<float_type, LATTICE::DIM>& u) {
		static_cast<Derived*>(this)->setVelocityAtIndex(idx, u);
	}

	void setFAt(index_type k, index_type idx, float_type value) {
		static_cast<Derived*>(this)->setFAt(k, idx, value);
	}

	float_type computeEquilibriumForInit(
		index_type k, float_type rho_val,
		const std::array<float_type, LATTICE::DIM>& u) const {
		return static_cast<const Derived*>(this)->computeEquilibriumForInit(
			k, rho_val, u);
	}

	const std::array<float_type, LATTICE::DIM>& getVelocityAtIndex(
		index_type idx) const {
		return static_cast<const Derived*>(this)->getVelocityAtIndex(idx);
	}

	float_type getRhoAtIndex(index_type idx) const {
		return static_cast<const Derived*>(this)->getRhoAtIndex(idx);
	}
};

#endif