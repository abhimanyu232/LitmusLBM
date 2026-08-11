# case: DOUBLE PERIODIC SHEAR LAYER
set(CASE_SETUP_NAME dbl_periodic_shear_layer_2D)
add_executable(${CASE_SETUP_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/dbl_periodic_shear_layer_2D.cpp)
target_link_libraries(${CASE_SETUP_NAME} PRIVATE LLBM)
install(
    TARGETS ${CASE_SETUP_NAME}
    RUNTIME DESTINATION .
)

set_target_properties(${CASE_SETUP_NAME} PROPERTIES CXX_CLANG_TIDY "${CMAKE_CXX_CLANG_TIDY}")
