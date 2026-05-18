#ifndef LAYOUT_POLICY_H
#define LAYOUT_POLICY_H

#include <concepts>
#include <cstddef>
#include <string>
#include <string_view>
#include <type_traits>

// static policy to set the storage format for the populations.
// choice between policies at compile time
template <typename POLICY>
concept LayoutPolicy = requires(size_t pop_idx, size_t node_idx,
																size_t POP_SIZE, size_t DOMAIN_SIZE) {
	{
		POLICY::getIndex(pop_idx, node_idx, POP_SIZE, DOMAIN_SIZE)
	} -> std::same_as<size_t>;
	{ POLICY::getPopStride(POP_SIZE, DOMAIN_SIZE) } -> std::same_as<size_t>;
	{ POLICY::getNodeStride(POP_SIZE, DOMAIN_SIZE) } -> std::same_as<size_t>;
	{ POLICY::getPolicyName() } -> std::convertible_to<std::string_view>;
};

// Population Storage Policy: Structure of Arrays
// Storage:	f0_p0, f0_p1, ..., f0_pN, f1_p0, ...
struct SoALayout {
	static constexpr std::string getPolicyName() {
		return std::string{"Structure of Arrays"};
	}

	static constexpr size_t getIndex(size_t pop_idx, size_t node_idx,
																	 size_t POP_SIZE, size_t DOMAIN_SIZE) {
		return pop_idx * DOMAIN_SIZE + node_idx;
	}

	static constexpr size_t getPopStride(size_t POP_SIZE, size_t DOMAIN_SIZE) {
		return DOMAIN_SIZE;
	}

	static constexpr size_t getNodeStride(size_t POP_SIZE, size_t DOMAIN_SIZE) {
		return 1;
	}
};

// Population Storage Policy: Array of Structures
// Storage:	f0_p0, f1_p0, ..., f8_p0, f0_p1, ...
struct AoSLayout {
	static constexpr std::string getPolicyName() {
		return std::string{"Array of Structures"};
	}

	static constexpr size_t getIndex(size_t pop_idx, size_t node_idx,
																	 size_t POP_SIZE, size_t DOMAIN_SIZE) {
		return pop_idx + node_idx * POP_SIZE;
	}

	static constexpr size_t getPopStride(size_t POP_SIZE, size_t DOMAIN_SIZE) {
		return 1;
	}

	static constexpr size_t getNodeStride(size_t POP_SIZE, size_t DOMAIN_SIZE) {
		return POP_SIZE;
	}
};

// !!! todo: Array of Struct of Array.
// Effectively a blocked SoA layout, middle ground between GPU and CPU,
// block sizes can be tuned to simd registerr or cuda warp width
// example block_size=8 (AVX2) ; each block = block_size*Q = 8 nodes × 9 pops = 72 doubles
// Block 0: f0_p0..f0_p7, f1_p0..f1_p7, ..., f8_p0..f8_p7;
// Block 1: f0_p8..f0_p15, f1_p8..f1_p15, ...
// !!! todo: Array of Struct of Array.
// template <size_t BLOCK_SIZE = 8>
// struct AoSoA {
// 	static constexpr std::string getPolicyName() {
// 		return std::string{"Array of Structures of Arrays (Block Size: "} +
// 					 std::to_string(BLOCK_SIZE) + ")";
// 	}

// 	static constexpr size_t getIndex(size_t pop_idx, size_t node_idx,
// 																	 size_t POP_SIZE, size_t DOMAIN_SIZE) {
// 		const size_t block_id = node_idx / BLOCK_SIZE;
// 		const size_t local_idx = node_idx % BLOCK_SIZE;
// 		return block_id * (POP_SIZE * BLOCK_SIZE) + (pop_idx * BLOCK_SIZE) +
// 					 local_idx;
// 	}

// 	static constexpr size_t getPopStride(size_t POP_SIZE, size_t DOMAIN_SIZE) {
// 		return BLOCK_SIZE;
// 	}

// 	static constexpr size_t getNodeStride(size_t POP_SIZE, size_t DOMAIN_SIZE) {
// 		// Note: Node stride is non-uniform across block boundaries in AoSoA
// 		return 1;
// 	}
// };

#endif