
#ifndef MODELS_H
#define MODELS_H

#include "../lattices.h"

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

	void saveVelocitySlice2D(int index, int axis, const std::string& filename) {
		static_cast<Derived*>(this)->saveVelocitySlice2D(index, axis, filename);
	};

	// Getters for simulation state
	std::string getModelName() const {
		return static_cast<const Derived*>(this)->getModelName();
	}

	int getCurrentStep() const {
		return static_cast<const Derived*>(this)->getCurrentStep();
	}

	double getCurrentTime() const {
		return static_cast<const Derived*>(this)->getCurrentTime();
	}

	// Configuration functions
	void setViscosity(double nu) {
		static_cast<Derived*>(this)->setViscosity(nu);
	}

	void setRelaxationTime(double tau) {
		static_cast<Derived*>(this)->setRelaxationTime(tau);
	}

	double getViscosity() const {
		return static_cast<const Derived*>(this)->getViscosity();
	}

	double getRelaxationTime() const {
		return static_cast<const Derived*>(this)->getRelaxationTime();
	}

	// Initialization helpers for external setups
	int getTotalNodes() const {
		return static_cast<const Derived*>(this)->getTotalNodes();
	}

	int getQ() const { return static_cast<const Derived*>(this)->getQ(); }

	std::array<int, LATTICE::DIM> getLatticeSize() const {
		return static_cast<const Derived*>(this)->getLatticeSize();
	}

	std::array<int, LATTICE::DIM> getPositionFromIndex(int idx) const {
		return static_cast<const Derived*>(this)->getPositionFromIndex(idx);
	}

	void setRhoAtIndex(int idx, double rho) {
		static_cast<Derived*>(this)->setRhoAtIndex(idx, rho);
	}

	void setVelocityAtIndex(int idx, const std::array<double, LATTICE::DIM>& u) {
		static_cast<Derived*>(this)->setVelocityAtIndex(idx, u);
	}

	void setFAt(int k, int idx, double value) {
		static_cast<Derived*>(this)->setFAt(k, idx, value);
	}

	double computeEquilibriumForInit(
		int k, double rho_val, const std::array<double, LATTICE::DIM>& u) const {
		return static_cast<const Derived*>(this)->computeEquilibriumForInit(
			k, rho_val, u);
	}

	const std::array<double, LATTICE::DIM>& getVelocityAtIndex(int idx) const {
		return static_cast<const Derived*>(this)->getVelocityAtIndex(idx);
	}

	const double getRhoAtIndex(int idx) const {
		return static_cast<const Derived*>(this)->getRhoAtIndex(idx);
	}
};


#endif