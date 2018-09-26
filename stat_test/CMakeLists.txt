set(STAT_TEST_SOURCE
        utility/Distrubution.h
        utility/IndependentPairFactory.cpp
        utility/InfoLooper.h
        utility/ParallelPairFactory.cpp
        utility/MockBC.h
        utility/ParalellInfoLooper.h
        utility/Quantity.cpp
        utility/ShapeGenerators.h
        utility/ShapePairFactory.h
        utility/UniformBallDistribution.cpp
        utility/UniformBoxDistribution.cpp
        AnisotropicShape2DExclusionDrawer.cpp
        AnisotropicShape2DExclusionTest.cpp
        Main.cpp
        OrderParamTest.cpp
        PackingOverlapsTest.cpp
        ShapeBCTest.cpp
        ShapeOverlapTest.cpp
        ShapePointInsideTest.cpp
        ShapeStaticSizesTest.cpp
        CuboidSpeedTest.cpp
        RectangleCase.cpp
        VectorSpeedTest.cpp)

function(register_targets_for_dimensions SD AD)
    set(DIM_SUFFIX ".${SD}.${AD}")
    set(STAT_TEST_EXEC "stat_test${DIM_SUFFIX}")
    set(RSA_LIB "rsa_main${DIM_SUFFIX}")

    add_executable(${STAT_TEST_EXEC} ${STAT_TEST_SOURCE})
    target_compile_definitions(${STAT_TEST_EXEC} PUBLIC RSA_SPATIAL_DIMENSION=${SD} RSA_ANGULAR_DIMENSION=${AD})
    target_link_libraries(${STAT_TEST_EXEC} ${RSA_LIB})
endfunction()

register_targets_for_dimensions(3 0)    # 3D spheres, cuboids, platonic solids and so on
register_targets_for_dimensions(2 0)    # 2D disks
register_targets_for_dimensions(2 1)    # AnisotropicShape2D's