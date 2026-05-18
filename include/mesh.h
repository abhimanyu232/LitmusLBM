#ifndef MESH_H
#define MESH_H

#include "common.h"
#include "lattices.h"
#include "policy/layout_policy.h"

template <size_t DIM, LatticeType LATTICE, LayoutPolicy LAYOUT>
struct Mesh {};

// 2D Implementation of Mesh
template <LatticeType LATTICE, LayoutPolicy LAYOUT>
struct Mesh<2, LATTICE, LAYOUT> {
	// members
	static constexpr size_t DIM = 2;
	std::array<size_t, 2> domain_size;
	size_t total_nodes{1};

	// ideally this should only have references to nodes to save on memory
	// but if we store it together, might be better for cache.
	//!!! todo: (perf) (memory) can skip on storing this in memory and only computing positions when needed.
	// in 2D, 500x500 domain, size_t == 4MB of memory
	// in 3D, 500x500x500 domain, size_t == 3 GB of memory 
	// for D3Q27, upto (541×541×541) fits inside uint32_t, only need size_t(unint64_t) for larger domains.
	// for D3Q19, upto (609x609x609) fits inside uint32_t, only need size_t(unint64_t) for larger domains.
	std::vector<std::array<size_t, 2>> nodes;

	// in 2D, D2Q9, 500x500, size_t = 18MB per neighbour index array
	// !!! in 3D, D3Q27, 500x500x500, size_t = 27 GB of memory per neighbour index array 
	std::vector<std::array<size_t, LATTICE::Q>> next_neighbour_index;
	std::vector<std::array<size_t, LATTICE::Q>> prev_neighbour_index;

	std::vector<std::array<size_t, 2>> inlet_nodes;
	std::vector<std::array<size_t, 2>> outlet_nodes;
	std::vector<std::array<size_t, 2>> top_boundary_nodes;
	std::vector<std::array<size_t, 2>> bottom_boundary_nodes;

	std::vector<size_t> inlet_nodes_index;
	std::vector<size_t> outlet_nodes_index;
	std::vector<size_t> top_nodes_index;
	std::vector<size_t> bottom_nodes_index;

	// !!! todo: why is this a double
	std::vector<std::array<double, 2>> wall_nodes;
	std::vector<size_t> wall_nodes_index;

 public:
	// Ctors
	// standard constructor without a body/geometry
	Mesh(std::array<size_t, 2> dimensions) : domain_size(dimensions) {

		static_assert(LATTICE::DIM == 2,
									"ERROR: 2D mesh initialized with 3D lattice \n");

		total_nodes = 1;
		for (size_t d = 0; d < 2; ++d)
			total_nodes *= dimensions[d];

		std::cout << "total nodes: " << total_nodes << std::endl;

		nodes.resize(total_nodes);
		next_neighbour_index.resize(total_nodes);
		prev_neighbour_index.resize(total_nodes);

		// boundary nodes
		inlet_nodes.reserve(dimensions[1]);
		outlet_nodes.reserve(dimensions[1]);
		top_boundary_nodes.reserve(dimensions[0]);
		bottom_boundary_nodes.reserve(dimensions[0]);

		for (size_t i = 0; i < total_nodes; ++i) {

			std::array<size_t, 2> pos;

			size_t tmp_flat_index = i;
			for (size_t d = 0; d < 2; ++d) {
				pos[d] = tmp_flat_index % domain_size[d];
				tmp_flat_index /= domain_size[d];
			}
			nodes[i] = pos;

			// find and save boundary nodes
			// corners belong to the top/bottom, not the inlet/outlet
			// todo: check if c++ algorithms class has somethign that will do it functionally
			if (pos[1] == 0)	// bottom  : y = 0
			{
				bottom_boundary_nodes.push_back(pos);
				bottom_nodes_index.push_back(i);
			} else if (pos[1] == domain_size[1] - 1)	// top : y = Ly -1
			{
				top_boundary_nodes.push_back(pos);
				top_nodes_index.push_back(i);
			} else if (pos[0] == 0)	 // inlet : x = 0
			{
				inlet_nodes.push_back(pos);
				inlet_nodes_index.push_back(i);
			} else if (pos[0] == domain_size[0] - 1)	// outlet : x = Lx -1
			{
				outlet_nodes.push_back(pos);
				outlet_nodes_index.push_back(i);
			}
		}

		compute_neighbour_indices();
	}

	// constructor for when theres a body inside.
	template <typename GEOMETRY>
	Mesh(std::array<size_t, 2> dimensions, GEOMETRY* body) : Mesh(dimensions) {
		std::cerr << " Mesh with bodies inside not implemented yet" << std::endl;
		std::exit(EXIT_FAILURE);

		// then do the geometry detection logic
		// detect_wall_nodes(body);
	}

	// compute the flat index from the position vector
	size_t getNodeIndex(const std::array<size_t, 2>& pos) const {
		size_t index = 0;
		for (size_t stride = 1, d = 0; d < 2; ++d) {
			index += pos[d] * stride;
			stride *= domain_size[d];
		}
		return index;
	}

	void compute_neighbour_indices();

	void writeToFile() {
		std::ofstream mesh_file("mesh.dat");
		std::ofstream inlet_file("inlet.dat");
		std::ofstream outlet_file("outlet.dat");
		std::ofstream top_file("top.dat");
		std::ofstream bottom_file("bottom.dat");

		for (size_t index = 0; index < total_nodes; ++index) {
			mesh_file << nodes[index][0] << "\t" << nodes[index][1] << "\n";
		}

		for (size_t index = 0; index < inlet_nodes.size(); ++index) {
			inlet_file << inlet_nodes[index][0] << "\t" << inlet_nodes[index][1]
								 << "\n";
		}
		for (size_t index = 0; index < outlet_nodes.size(); ++index) {
			outlet_file << outlet_nodes[index][0] << "\t" << outlet_nodes[index][1]
									<< "\n";
		}
		for (size_t index = 0; index < top_boundary_nodes.size(); ++index) {
			top_file << top_boundary_nodes[index][0] << "\t"
							 << top_boundary_nodes[index][1] << "\n";
		}
		for (size_t index = 0; index < bottom_boundary_nodes.size(); ++index) {
			bottom_file << bottom_boundary_nodes[index][0] << "\t"
									<< bottom_boundary_nodes[index][1] << "\n";
		}

		for (size_t index = 0; index < total_nodes; ++index) {
			mesh_file << nodes[index][0] << "\t" << nodes[index][1] << "\n";
		}
	}
};

template <LatticeType LATTICE, LayoutPolicy LAYOUT>
void Mesh<2, LATTICE, LAYOUT>::compute_neighbour_indices() {
	const auto& v = LATTICE::velocities;
	for (size_t idx = 0; idx < total_nodes; ++idx) {
		std::array<size_t, 2> pos = nodes[idx];

		// compute each neighbour : total Q neighbours
		for (size_t k = 0; k < LATTICE::Q; ++k) {
			std::array<size_t, 2> next_nb_pos, prev_nb_pos;
			for (size_t d = 0; d < 2; ++d) {
				size_t cIdx = k + d * LATTICE::Q;

				// !!! relies on unsigned integer underflow
				// note: here int v[cIdx] is implicitly converted to size_t resulting in underflow. 
				// note: math works as intended after wrapping, uint over/underflow is defined behaviour. 
				// !!! valid as long as (pos[d] + v[cIdx] + domain_size[d]) < (2^64  for uint64_t) and (2^32 for uint32_t)  
				next_nb_pos[d] = (pos[d] + v[cIdx] + domain_size[d]) % domain_size[d];
				prev_nb_pos[d] = ((pos[d]) - v[cIdx] + domain_size[d]) % domain_size[d];
				
				// note: safer pattern: avoids unsigned over/underflow, 
				// note: all math is signed, cast to unsigned in the end
				// auto tmp = static_cast<int>(pos[d]) - v[cIdx];
				// prev_nb_pos[d] = static_cast<index_t>( (tmp % static_cast<int>(domain_size[d])) + static_cast<int>(domain_size[d])) %
				// 								 static_cast<int>(domain_size[d]);
			}

			// Layout agnostic index
			size_t next_idx =
			LAYOUT::getIndex(k, getNodeIndex(next_nb_pos), LATTICE::Q, total_nodes);
			size_t prev_idx =
				LAYOUT::getIndex(k, getNodeIndex(prev_nb_pos), LATTICE::Q, total_nodes);


			next_neighbour_index[idx][k] = next_idx;
			prev_neighbour_index[idx][k] = prev_idx;
		}
	}
}

// template <>
// Mesh<3> {

// only for 3D
// std::vector<std::array<double, 3>> left_boundary_nodes;
// std::vector<std::array<double, 3>> right_boundary_nodes;
// std::vector<size_t> right_nodes_index;
// std::vector<size_t> left_nodes_index;

// if (DIM == 3) {
// 	left_boundary_nodes.reserve(dimensions[0]);
// 	right_boundary_nodes.reserve(dimensions[1]);
// }

// if (DIM == 3) {			// direction: looking inwards from the inlet
// 	if (pos[2] == 0)	// right
// 	{
// 		right_nodes.push_back(pos);
// 		right_nodes_index.push_back(i);
// 	}
// 	if (pos[2] == domain_size[2] - 1)	 // left : z = Lz -1
// 	{
// 		left_nodes.push_back(pos);
// 		left_nodes_index.push_back(i);
// 	}
// }
// };

#endif