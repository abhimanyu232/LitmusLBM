#ifndef FIELDS_H
#define FIELDS_H

#include "common.h"
#include "lattices.h"

template <typename T>
class ScalarField {
 private:
	size_t size_;
	std::vector<T> data_;

 public:
	explicit ScalarField(size_t size) : size_(size), data_(size) {}

	T& operator[](size_t index) { return data_[index]; }

	const T& operator[](size_t index) const { return data_[index]; }
};

template <typename T, size_t DIM>
class VectorField {
 private:
	size_t size_;
	std::vector<std::array<T, DIM>> data_;

 public:
	explicit VectorField(size_t size) : size_(size), data_(size) {}

	std::array<T, DIM>& operator[](size_t index) { return data_[index]; }

	const std::array<T, DIM>& operator[](size_t index) const { return data_[index]; }
};

template <typename T, LatticeType LATTICE>
class Populations {
 private:
	static constexpr size_t Q = LATTICE::Q;
	static constexpr size_t DIM = LATTICE::DIM;

	// f0_i .... f0_Q-1, f_new_i , ... f_new_Q-1
	size_t domain_size_;
	std::vector<std::array<T, LATTICE::Q>> data_;

 public:
	Populations(size_t size) : domain_size_(size), data_(size * 2) {}

	// compute moments
	// here avoid copy using std::move if possible.
	// std::array<T,SIZE>& computeRhoField(){  };
	// computeRhoat(size_t index) { }

	static constexpr enum POP_ID { OLD = 0, NEW = 1 };

	// access: Populations<double,D2Q9> f; f[OLD,0,..Q] and f[NEW,0,...,Q]
	// 	access population vector at given node index
	std::array<T, LATTICE::Q>& operator()(enum POP_ID, size_t index) {
		return data_[(POP_ID * domain_size_) + index];
	}

	const std::array<T, LATTICE::Q>& operator()(enum POP_ID, size_t index) const {
		return data_[(POP_ID * domain_size_) + index];
	}

	// access individual population at k at given node index
	T& operator()(enum POP_ID, size_t index, size_t k) {
		return data_[(POP_ID * domain_size_) + index][k];
	}

	const T& operator()(enum POP_ID, size_t index, size_t k) const {
		return data_[(POP_ID * domain_size_) + index][k];
	}
};

template <size_t SIZE, LatticeType LATTICE>
class FieldSet {};

#endif