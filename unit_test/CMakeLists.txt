set(UNIT_TEST_SOURCE

        rsa3d/MatrixTest.cpp
        rsa3d/PackingTest.cpp
        rsa3d/VectorTest.cpp)

add_library(unit_test_runner lib/catch.hpp lib/fakeit.hpp Main.cpp)
target_include_directories(unit_test_runner PUBLIC lib)

function(register_targets_for_dimensions SD AD)
    set(DIM_SUFFIX ".${SD}.${AD}")
    set(UNIT_TEST_EXEC "unit_test${DIM_SUFFIX}")
    set(RSA_LIB "rsa_main${DIM_SUFFIX}")

    add_executable(${UNIT_TEST_EXEC} ${UNIT_TEST_SOURCE})
    target_compile_definitions(${UNIT_TEST_EXEC} PUBLIC RSA_SPATIAL_DIMENSION=${SD} RSA_ANGULAR_DIMENSION=${AD})
    target_link_libraries(${UNIT_TEST_EXEC} unit_test_runner ${RSA_LIB})
endfunction()

register_targets_for_dimensions(3 0)    # 3D spheres, cuboids, platonic solids and so on
register_targets_for_dimensions(2 0)    # 2D disks
register_targets_for_dimensions(2 1)    # AnisotropicShape2D's