#ifndef WRITER_H
#define WRITER_H

#include <concepts>

template <typename Writer>
concept FieldWriter = requires(
	Writer& writer, std::string_view name, const std::vector<double>& scalar,
	const std::vector<std::array<double, 3>>& vector) {
	{ writer.write_scalar_double_dataset(name, scalar) } -> std::same_as<bool>;

	{ writer.write_vector_double_dataset(name, vector) } -> std::same_as<bool>;
};


#endif