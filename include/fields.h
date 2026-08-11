#ifndef FIELDS_H
#define FIELDS_H

#include "common.h"
#include "lattices.h"

template <typename T>
class ScalarField {
 private:
	index_type size_;
	std::vector<T> data_;

 public:
	explicit ScalarField(index_type size) : size_(size), data_(size) {}

	T& operator[](index_type index) { return data_[index]; }

	const T& operator[](index_type index) const { return data_[index]; }
};

template <typename T, index_type Dim>
class VectorField {
 private:
	index_type size_;
	std::vector<std::array<T, Dim>> data_;

 public:
	explicit VectorField(index_type size) : size_(size), data_(size) {}

	std::array<T, Dim>& operator[](index_type index) { return data_[index]; }

	const std::array<T, Dim>& operator[](index_type index) const {
		return data_[index];
	}
};

template <typename T, LatticeType Lattice>
class Populations {
 private:
	static constexpr index_type Q = Lattice::Q;
	static constexpr index_type Dim = Lattice::Dim;

	// f0_i .... f0_Q-1, f_new_i , ... f_new_Q-1
	index_type domain_size_;
	std::vector<std::array<T, Lattice::Q>> data_;

 public:
	Populations(index_type size) : domain_size_(size), data_(size * 2) {}

	// compute moments
	// here avoid copy using std::move if possible.
	// std::array<T,SIZE>& computeRhoField(){  };
	// computeRhoat(index_type index) { }

	static constexpr enum POP_ID { OLD = 0, NEW = 1 };

	// access: Populations<float_type,D2Q9> f; f[OLD,0,..Q] and f[NEW,0,...,Q]
	// 	access population vector at given node index
	std::array<T, Lattice::Q>& operator()(enum POP_ID, index_type index) {
		return data_[(POP_ID * domain_size_) + index];
	}

	const std::array<T, Lattice::Q>& operator()(enum POP_ID, index_type index) const {
		return data_[(POP_ID * domain_size_) + index];
	}

	// access individual population at k at given node index
	T& operator()(enum POP_ID, index_type index, index_type k) {
		return data_[(POP_ID * domain_size_) + index][k];
	}

	const T& operator()(enum POP_ID, index_type index, index_type k) const {
		return data_[(POP_ID * domain_size_) + index][k];
	}
};

template <index_type SIZE, LatticeType Lattice>
class FieldSet {};

#endif