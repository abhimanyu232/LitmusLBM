cmake_minimum_required(VERSION 3.28)

set(PROJECT_NAME LitmusLBM)
project(${PROJECT_NAME} VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -march=native -mtune=native")

### Clang Tidy
find_program(CLANG_TIDY NAMES "clang-tidy-20")

# Get GCC toolchain path from the C++ compiler (handles Spack, system, etc.)
get_filename_component(GCC_BIN_DIR "${CMAKE_CXX_COMPILER}" DIRECTORY)
get_filename_component(GCC_TOOLCHAIN "${GCC_BIN_DIR}" DIRECTORY)

if(CLANG_TIDY)
  message(STATUS "Clang-tidy found: ${CLANG_TIDY}")
  set(CMAKE_CXX_CLANG_TIDY
        ${CLANG_TIDY}
        "-config-file=${CMAKE_CURRENT_SOURCE_DIR}/.clang-tidy"
        "--extra-arg=-Wno-unknown-warning-option"
				"--extra-arg=--gcc-toolchain=${GCC_TOOLCHAIN}"
        "--extra-arg=-stdlib=libstdc++" # ensures that clang-tidy uses the right stdlib, else shows errors with gcc. does not affect compilation
				"--extra-arg=-std=c++23" 
				"--extra-arg=-fopenmp=libomp"
				"--extra-arg=-I/usr/lib/llvm-18/lib/clang/18/include" # ensures that clang-tidy uses the right include dir, else cant find omp.h. does not affect compilation
    ) 
else() ##-Wno-unknown-warning-option,-stdlib=libc++,-std=c++23
  message(WARNING "Clang-tidy not found. Static analysis disabled.")
endif()

### SET BUILD TYPE
## "Debug: {-g}" 
## "Release: {-O3 -DNDEBUG} " 
## "MinSizeRel: {-Os}" 
## "RelWithDebInfo: { -O2 -g -DNDEBUG}") [DEFAULT]
set(default_build_type "RelWithDebInfo")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

find_package(OpenMP REQUIRED)
###### hypothetical dependencies ######
# find_package(Boost 1.36.0 COMPONENTS filesystem system REQUIRED)
# target_link_libraries(example PUBLIC
#   ${Boost_LIBRARIES}
#   # here you can add any library dependencies
# )

### Define include directories
set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
### Compile header files into a library
## include_directories(${COMMON_INCLUDES})
### INSTEAD of the above, use interface library to include directories
add_library(LLBM INTERFACE) 
target_include_directories(LLBM INTERFACE ${COMMON_INCLUDES})
target_link_libraries(LLBM INTERFACE OpenMP::OpenMP_CXX)
target_compile_features(LLBM INTERFACE cxx_std_23)
set_target_properties(LLBM PROPERTIES CXX_CLANG_TIDY "${CMAKE_CXX_CLANG_TIDY}")

## IF I needed an external library, say header only Boost 
# target_link_libraries(lm_matlib INTERFACE Boost::Boost)


### Define Targets
add_subdirectory(src)



# set(TESTING OFF CACHE BOOL "True if tests shall be built")
# if (TESTING)
#   enable_testing()
#   add_subdirectory(tests)
# endif() 