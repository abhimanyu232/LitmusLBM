#ifndef HDF5_WRITER_H
#define HDF5_WRITER_H

#include <hdf5.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>

#include <common.h>

#include "writer.h"

template <size_t Dim, typename Config_t>
class Hdf5VTKWriter {
 public:
	// Hdf5Writer() = delete;

	Hdf5VTKWriter(const std::filesystem::path& filename, const Config_t& config,
								std::uint64_t step, double time) {

		static_assert(Dim == 2 || Dim == 3,
									"Hdf5VTKWriter only supports 2D and 3D");

		// HDF5 order is [z,y,x], while VTK extent order is [x,y,z].
		// VTK extent uses x-min,x-max,y-min,y-max,z-min,z-max.
		std::array<int, 6> whole_extent = {};
		if constexpr (Dim == 2) {
			nx = config.domain_size[0];
			ny = config.domain_size[1];

			scalar_shape = {ny, nx};									// (1,ny,nx)
			vector_shape = {ny, nx, vtk_vector_dim};	// (1,ny,nx,3)
			whole_extent = {
				0, static_cast<int>(nx - 1), 0, static_cast<int>(ny - 1), 0, 0};

		} else {
			nx = config.domain_size[0];
			ny = config.domain_size[1];
			nz = config.domain_size[2];

			scalar_shape = {nz, ny, nx};
			vector_shape = {nz, ny, nx, vtk_vector_dim};
			whole_extent = {0, static_cast<int>(nx - 1), 0, static_cast<int>(ny - 1),
											0, static_cast<int>(nz - 1)};
		}

		//!!! todo fix format here
		// using std::filesystem
		std::filesystem::path temp_fname = filename;
		temp_fname.replace_extension(".vtkhdf");
		// char* temp_fname = _fname + "." + fformat;

		file = H5Fcreate(temp_fname.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
										 H5P_DEFAULT);
		if (file < 0) {
			hdf5_error("H5Fcreate");
			throw std::runtime_error("Failed to create HDF5 file");
			// std::exit(EXIT_FAILURE);
			// MPI_Abort(comm, EXIT_FAILURE);
		}

		vtkhdf = H5Gcreate2(file, "VTKHDF", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		if (vtkhdf < 0) {
			hdf5_error("create /VTKHDF");
			H5Fclose(file);
			throw std::runtime_error("Failed to create /VTKHDF group");
		}

		point_data =
			H5Gcreate2(vtkhdf, "PointData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (point_data < 0) {
			hdf5_error("create /VTKHDF/PointData");
			H5Gclose(vtkhdf);
			H5Fclose(file);
			throw std::runtime_error("Failed to create /VTKHDF/PointData group");
		}

		// hid_t cell_data =
		// 	H5Gcreate2(vtkhdf, "CellData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		// hid_t field_data =
		// 	H5Gcreate2(vtkhdf, "FieldData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		bool success = true;

		// VTKHDF dataset type: ImageData
		success &= write_string_attribute(vtkhdf, "Type", "ImageData");

		// VTKHDF version. Match this to the ParaView/VTK reader.
		const std::int64_t version[2] = {2, 0};
		success &= write_vector_attribute(vtkhdf, "Version", H5T_STD_I64LE,
																			H5T_NATIVE_INT64, 2, version);

		success &= write_vector_attribute(vtkhdf, "WholeExtent", H5T_STD_I32LE,
																			H5T_NATIVE_INT, whole_extent.size(),
																			whole_extent.data());

		const std::array<double, 3> origin{0.0, 0.0, 0.0};
		const std::array<double, 3> spacing{1.0, 1.0, 1.0};
		std::array<double, 9> direction{1.0, 0.0, 0.0, 0.0, 1.0,
																		0.0, 0.0, 0.0, 1.0};

		success &= write_vector_attribute(vtkhdf, "Origin", H5T_IEEE_F64LE,
																			H5T_NATIVE_DOUBLE, 3, origin.data());

		success &= write_vector_attribute(vtkhdf, "Spacing", H5T_IEEE_F64LE,
																			H5T_NATIVE_DOUBLE, 3, spacing.data());

		success &= write_vector_attribute(vtkhdf, "Direction", H5T_IEEE_F64LE,
																			H5T_NATIVE_DOUBLE, 9, direction.data());

		const std::uint64_t step_value = step;
		// Note: the VTK VTKHDF reader does not use these for animation;
		// only for post-processing.
		// Real time-driven animation requires the temporal VTKHDF layout (a /VTKHDF/Steps group).
		success &= write_scalar_attribute(vtkhdf, "SimulationStep", H5T_STD_U64LE,
																			H5T_NATIVE_UINT64, &step_value);
		success &= write_scalar_attribute(vtkhdf, "SimulationTime", H5T_IEEE_F64LE,
																			H5T_NATIVE_DOUBLE, &time);

		if (!success) {
			close();
			throw std::runtime_error("Failed to write VTKHDF metadata attributes");
		}
	}

	void close() noexcept {
		if (point_data != H5I_INVALID_HID) {
			H5Gclose(point_data);
			point_data = H5I_INVALID_HID;
		}

		if (vtkhdf != H5I_INVALID_HID) {
			H5Gclose(vtkhdf);
			vtkhdf = H5I_INVALID_HID;
		}

		if (file != H5I_INVALID_HID) {
			H5Fclose(file);
			file = H5I_INVALID_HID;
		}
	}

	~Hdf5VTKWriter() { close(); }

	// no copy allowed
	Hdf5VTKWriter(const Hdf5VTKWriter&) = delete;
	Hdf5VTKWriter& operator=(const Hdf5VTKWriter&) = delete;

	bool hdf5_error(const char* operation) const;

	bool write_string_attribute(hid_t object, std::string_view name,
															const std::string& value);

	bool write_vector_attribute(hid_t object, std::string_view name,
															hid_t file_type, hid_t memory_type,
															hsize_t length, const void* data);

	bool write_scalar_attribute(hid_t object, std::string_view name,
															hid_t file_type, hid_t memory_type,
															const void* data);

	bool write_scalar_double_dataset(std::string_view name,
																	 const std::vector<double>& data);

	bool write_vector_double_dataset(
		std::string_view name, const std::vector<std::array<double, 3>>& data);

	bool write_double_dataset(std::string_view name, int rank,
														const hsize_t* dimensions, const double* data,
														std::size_t data_size,
														const char* vtk_attribute = nullptr);

 private:
	hsize_t nx = 1;
	hsize_t ny = 1;
	hsize_t nz = 1;
	std::array<hsize_t, Dim> scalar_shape;
	std::array<hsize_t, Dim + 1> vector_shape;
	static constexpr std::size_t vtk_vector_dim = 3;

	hid_t point_data{H5I_INVALID_HID};
	hid_t vtkhdf{H5I_INVALID_HID};
	hid_t file{H5I_INVALID_HID};

	// struct Hdf5CompressionOptions {
	// 	int deflate_level = 1;	// 0 = disabled, 1-9 = enabled
	// 	bool shuffle = true;
	// 	std::size_t target_chunk_bytes = 1 << 20;	 // 1 MiB
	// }/model;
};

template <std::size_t Dim, typename Config_t>
bool Hdf5VTKWriter<Dim, Config_t>::hdf5_error(const char* operation) const {

	std::cerr << "HDF5 error in " << operation << '\n';
	H5Eprint2(H5E_DEFAULT, stderr);
	return false;
}

template <std::size_t Dim, typename Config_t>
bool Hdf5VTKWriter<Dim, Config_t>::write_vector_attribute(
	hid_t object, std::string_view name, hid_t file_type, hid_t memory_type,
	hsize_t length, const void* data) {
	hid_t space = H5Screate_simple(1, &length, nullptr);
	if (space < 0) {
		return hdf5_error("H5Screate_simple(attribute)");
	}

	const std::string name_string{name};
	hid_t attribute = H5Acreate2(object, name_string.c_str(), file_type, space,
															 H5P_DEFAULT, H5P_DEFAULT);
	if (attribute < 0) {
		H5Sclose(space);
		return hdf5_error("H5Acreate2");
	}

	const herr_t status = H5Awrite(attribute, memory_type, data);
	if (status < 0) {
		H5Aclose(attribute);
		H5Sclose(space);
		return hdf5_error("H5Awrite");
	}

	H5Aclose(attribute);
	H5Sclose(space);
	return true;
}

template <std::size_t Dim, typename Config_t>
bool Hdf5VTKWriter<Dim, Config_t>::write_scalar_attribute(hid_t object,
																													std::string_view name,
																													hid_t file_type,
																													hid_t memory_type,
																													const void* data) {
	hid_t space = H5Screate(H5S_SCALAR);
	if (space < 0) {
		return hdf5_error("H5Screate(scalar attribute)");
	}

	const std::string name_string{name};
	hid_t attribute = H5Acreate2(object, name_string.c_str(), file_type, space,
															 H5P_DEFAULT, H5P_DEFAULT);
	if (attribute < 0) {
		H5Sclose(space);
		return hdf5_error("H5Acreate2(scalar)");
	}

	const herr_t status = H5Awrite(attribute, memory_type, data);
	if (status < 0) {
		H5Aclose(attribute);
		H5Sclose(space);
		return hdf5_error("H5Awrite(scalar)");
	}

	H5Aclose(attribute);
	H5Sclose(space);
	return true;
}

template <std::size_t Dim, typename Config_t>
bool Hdf5VTKWriter<Dim, Config_t>::write_string_attribute(
	hid_t object, std::string_view name, const std::string& value) {
	hid_t space = H5Screate(H5S_SCALAR);
	//// hid_t object = vtkhdf;
	if (space < 0) {
		return hdf5_error("H5Screate(string attribute)");
	}

	hid_t string_type = H5Tcopy(H5T_C_S1);
	if (string_type < 0) {
		H5Sclose(space);
		return hdf5_error("H5Tcopy");
	}

	H5Tset_size(string_type, value.size() + 1);
	H5Tset_cset(string_type, H5T_CSET_UTF8);

	const std::string name_string{name};
	hid_t attribute = H5Acreate2(object, name_string.c_str(), string_type, space,
															 H5P_DEFAULT, H5P_DEFAULT);

	if (attribute < 0) {
		H5Tclose(string_type);
		H5Sclose(space);
		return hdf5_error("H5Acreate2(string)");
	}

	const herr_t status = H5Awrite(attribute, string_type, value.data());

	if (status < 0) {
		H5Aclose(attribute);
		H5Tclose(string_type);
		H5Sclose(space);
		return hdf5_error("H5Awrite(string)");
	}

	H5Aclose(attribute);
	H5Tclose(string_type);
	H5Sclose(space);
	return true;
}

template <std::size_t Dim, typename Config_t>
bool Hdf5VTKWriter<Dim, Config_t>::write_scalar_double_dataset(
	std::string_view name, const std::vector<double>& data) {
	// int rank = static_cast<int>(Dim);
	// const char* vtk_attribute = "Scalars";
	return this->write_double_dataset(name, static_cast<int>(Dim),
																		scalar_shape.data(), data.data(),
																		data.size(), "Scalars");
}

template <std::size_t Dim, typename Config_t>
bool Hdf5VTKWriter<Dim, Config_t>::write_vector_double_dataset(
	std::string_view name, const std::vector<std::array<double, 3>>& data) {
	// if (data.empty()) {
	// 	return hdf5_error("empty vector field");
	// }
	// int rank = static_cast<int>(Dim+1);
	// const char* vtk_attribute = "Vectors";
	return this->write_double_dataset(name, static_cast<int>(Dim + 1),
																		vector_shape.data(), data.front().data(),
																		data.size() * 3, "Vectors");
}

template <std::size_t Dim, typename Config_t>
bool Hdf5VTKWriter<Dim, Config_t>::write_double_dataset(
	std::string_view name, int rank, const hsize_t* dimensions,
	const double* data, std::size_t data_size, const char* vtk_attribute) {

	hid_t parent = point_data;
	std::size_t expected_size = 1;
	for (int d = 0; d < rank; ++d) {
		expected_size *= static_cast<std::size_t>(dimensions[d]);
	}

	if (data_size != expected_size) {
		return hdf5_error("dataset size mismatch");
	}

	hid_t dataspace = H5Screate_simple(rank, dimensions, nullptr);

	if (dataspace < 0) {
		return hdf5_error("H5Screate_simple(dataset)");
	}

	const std::string name_string{name};
	// VTKHDF stores portable IEEE-754 little-endian doubles.
	hid_t dataset = H5Dcreate2(parent, name_string.c_str(), H5T_IEEE_F64LE,
														 dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	if (dataset < 0) {
		H5Sclose(dataspace);
		return hdf5_error("H5Dcreate2");
	}

	// The memory buffer is native C++ double.
	bool success = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
													H5P_DEFAULT, data) >= 0;

	if (!success) {
		hdf5_error("H5Dwrite");
	}

	// This identifies active VTK arrays such as Scalars or Vectors.
	if (success && vtk_attribute != nullptr) {
		success = write_string_attribute(dataset, "Attribute", vtk_attribute);
	}

	H5Dclose(dataset);
	H5Sclose(dataspace);

	return success;
}

#endif