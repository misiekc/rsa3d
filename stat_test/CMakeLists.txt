set(STAT_TEST_COMMON_SOURCE
        utils/Distrubution.h
        utils/IndependentPairFactory.cpp
        utils/InfoLooper.h
        utils/ParallelPairFactory.cpp
        utils/ParalellInfoLooper.h
        utils/ShapeGenerators.h
        utils/ShapePairFactory.h
        utils/TestExitCodes.h
        utils/UniformBallDistribution.cpp
        utils/UniformBoxDistribution.cpp
        AnisotropicShape2DExclusionDrawer.cpp
        AnisotropicShape2DExclusionTest.cpp
        Main.cpp
        PackingOverlapsTest.cpp
        ShapeBCTest.cpp
        ShapeOverlapTest.cpp
        ShapePointInsideTest.cpp
        ShapeStaticSizesTest.cpp
        ShapeUltimateTest.cpp
        VectorSpeedTest.cpp)

set(STAT_TEST_SOURCE.3.0
        OrderParamTest.cpp
        CuboidSpeedTest.cpp)

function(register_targets_for_dimensions SD AD)
    set(DIM_SUFFIX ".${SD}.${AD}")
    set(STAT_TEST_EXEC "stat_test${DIM_SUFFIX}")
    set(RSA_LIB "rsa_main${DIM_SUFFIX}")
    set(SHAPES_LIB "shapes${DIM_SUFFIX}")
    set(ADDITIONAL_SOURCE "STAT_TEST_SOURCE${DIM_SUFFIX}")

    if(DEFINED ${ADDITIONAL_SOURCE})
        add_executable(${STAT_TEST_EXEC} ${STAT_TEST_COMMON_SOURCE} ${${ADDITIONAL_SOURCE}})
    else()
        add_executable(${STAT_TEST_EXEC} ${STAT_TEST_COMMON_SOURCE})
    endif()

    target_compile_definitions(${STAT_TEST_EXEC} PUBLIC RSA_SPATIAL_DIMENSION=${SD} RSA_ANGULAR_DIMENSION=${AD})
    target_link_libraries(${STAT_TEST_EXEC} ${SHAPES_LIB} ${RSA_LIB})
endfunction()

register_targets_for_dimensions(3 0)    # 3D spheres, cuboids, platonic solids and so on
register_targets_for_dimensions(2 0)    # 2D disks
register_targets_for_dimensions(2 1)    # AnisotropicShape2D's
register_targets_for_dimensions(1 1)