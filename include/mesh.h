#ifndef MESH_H
#define MESH_H

#include "common.h"
#include "lattices.h"

template <size_t DIM, LatticeType LATTICE>
struct Mesh {};

// 2D Implementation of Mesh
template <LatticeType LATTICE>
struct Mesh<2, LATTICE> {

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
		neighbour_index.resize(total_nodes);

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
			if (pos[1] == 0)	// bottom
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

	// compute the flat index from the position vector
	size_t getNodeIndex(std::array<size_t, 2>& pos) const {
		size_t index = 0;
		for (size_t stride = 1, d = 0; d < 2; ++d) {
			index += pos[d] * stride;
			stride *= domain_size[d];
		}
		return index;
	}

	// members
	static constexpr size_t DIM = 2;
	std::array<size_t, 2> domain_size;
	size_t total_nodes{1};

	// ideally this should only have references to nodes to save on memory
	// but if we store it together, might be better for cache.
	std::vector<std::array<size_t, 2>> nodes;
	std::vector<std::array<size_t, LATTICE::Q>> neighbour_index;

	std::vector<std::array<size_t, 2>> inlet_nodes;
	std::vector<std::array<size_t, 2>> outlet_nodes;
	std::vector<std::array<size_t, 2>> top_boundary_nodes;
	std::vector<std::array<size_t, 2>> bottom_boundary_nodes;

	std::vector<size_t> inlet_nodes_index;
	std::vector<size_t> outlet_nodes_index;
	std::vector<size_t> top_nodes_index;
	std::vector<size_t> bottom_nodes_index;

	std::vector<std::array<double, 2>> wall_nodes;
	std::vector<size_t> wall_nodes_index;
};

template <LatticeType LATTICE>
void Mesh<2, LATTICE>::compute_neighbour_indices() {
	const auto& v = LATTICE::velocities;
	// const auto w = LATTICE::weights;
	for (size_t idx = 0; idx < total_nodes; ++idx) {
		std::array<size_t, 2> pos = nodes[idx];

		// compute each neighbour : total Q neighbours
		for (size_t k = 0; k < LATTICE::Q; ++k) {
			// todo: check for negatives at the boundary.
			std::array<size_t, 2> nb_pos;
			for (size_t d = 0; d < 2; ++d) {
				size_t cIdx = k + d * LATTICE::Q;
				nb_pos[d] = (pos[d] + v[cIdx] + domain_size[d]) % domain_size[d];
				// nb_pos[d] = (pos[d] + v[k][d] + lattice_size[d]) % lattice_size[d];
			}
			// [AoS]
			size_t dst_idx = k + getNodeIndex(nb_pos) * LATTICE::Q;
			// size_t src = k + i * Q;

			// [SoA]
			// size_t dst = k * total_nodes + getNodeIndex(npos);
			// size_t src = k * total_nodes + i;

			neighbour_index[idx][k] = dst_idx;
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