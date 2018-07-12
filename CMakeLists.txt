cmake_minimum_required(VERSION 3.8)
project(rsa3d)

set(CMAKE_CXX_STANDARD 11)

set(MAIN_SOURCE
        rsa3d/analizator/Analyzer.cpp
        rsa3d/shape/shapes/cuboid/Cuboid.cpp
        rsa3d/shape/shapes/cuboid/CuboidOverlapStrategy.h
        rsa3d/shape/shapes/cuboid/MineOverlap.cpp
        rsa3d/shape/shapes/cuboid/OptimizedSATOverlap.cpp
        rsa3d/shape/shapes/cuboid/SATOverlap.cpp
        rsa3d/shape/shapes/cuboid/TriTriOverlap.cpp
        rsa3d/shape/shapes/regular_solid/Cuboctahedron.cpp
        rsa3d/shape/shapes/regular_solid/Dodecahedron.cpp
        rsa3d/shape/shapes/regular_solid/Icosahedron.cpp
        rsa3d/shape/shapes/regular_solid/Octahedron.cpp
        rsa3d/shape/shapes/regular_solid/RegularSolid.h
        rsa3d/shape/shapes/regular_solid/RegularSolid.tpp
        rsa3d/shape/shapes/regular_solid/Rhombicuboctahedron.cpp
        rsa3d/shape/shapes/regular_solid/SATOverlap.h
        rsa3d/shape/shapes/regular_solid/SnubCube.cpp
        rsa3d/shape/shapes/regular_solid/Tetrahedron.cpp
        rsa3d/shape/shapes/regular_solid/TriTriOverlap.h
        rsa3d/shape/shapes/regular_solid/TruncatedCube.cpp
        rsa3d/shape/shapes/regular_solid/TruncatedCuboctahedron.cpp
        rsa3d/shape/shapes/regular_solid/TruncatedOctahedron.cpp
        rsa3d/shape/shapes/polygon/HBPolygon.cpp
        rsa3d/shape/shapes/polygon/SBPolygon.cpp
        rsa3d/shape/shapes/Ellipse.cpp
        rsa3d/shape/shapes/OrientedCuboid.h
        rsa3d/shape/shapes/OrientedCuboid.tpp
        rsa3d/shape/shapes/polygon/Polygon.cpp
        rsa3d/shape/shapes/Rectangle.cpp
        rsa3d/shape/shapes/Sphere.h
        rsa3d/shape/shapes/Sphere.tpp
        rsa3d/shape/shapes/SpheroCylinder2D.cpp
        rsa3d/shape/AnisotropicShape2D.cpp
        rsa3d/shape/ConvexShape.h
        rsa3d/shape/ConvexShape.tpp
        rsa3d/shape/OrderCalculable.h
        rsa3d/shape/OverlapStrategy.h
        rsa3d/shape/OverlapStrategyShape.h
        rsa3d/shape/Shape.h
        rsa3d/shape/Shape.tpp
        rsa3d/shape/ShapeFactory.cpp
        rsa3d/surfaces/NBoxFBC.cpp
        rsa3d/surfaces/NBoxPBC.cpp
        rsa3d/BoundaryConditions.h
        rsa3d/Config.cpp
        rsa3d/Intersection.cpp
        rsa3d/Matrix.h
        rsa3d/NeighbourGrid.h
        rsa3d/Packing.cpp
        rsa3d/PackingGenerator.cpp
        rsa3d/Parameters.cpp
        rsa3d/Positioned.h
        rsa3d/Positioned.tpp
        rsa3d/RND.cpp
        rsa3d/Surface.cpp
        rsa3d/Timer.cpp
        rsa3d/Utils.cpp
        rsa3d/Vector.h
        rsa3d/Vector.tpp
        rsa3d/Voxel.h
        rsa3d/Voxel.cpp
        rsa3d/VoxelList.h
        rsa3d/VoxelList.cpp)

set(STATISTICS_SOURCE
        statistics/ASFRegression.cpp
        statistics/LinearRegression.cpp
        statistics/LogPlot.cpp
        statistics/Plot.cpp
        statistics/PowerRegression.cpp)

set(TEST_SOURCE
        test/utility/BallFactory.cpp
        test/utility/BoxFactory.cpp
        test/utility/InfoLooper.h
        test/utility/IsotropicFactory.h
        test/utility/MockBC.h
        test/utility/ParalellInfoLooper.h
        test/utility/Quantity.cpp
        test/utility/ShapeGenerators.h
        test/utility/ShapePairFactory.h
        test/AnisotropicShape2DExclusionDrawer.cpp
        test/AnisotropicShape2DExclusionTest.cpp
        test/PackingOverlapsTest.cpp
        test/ShapeOverlapTest.cpp
        test/ShapePointInsideTest.cpp
        test/ShapeStaticSizesTest.cpp
        test/CuboidSpeedTest.cpp
        test/RectangleCase.cpp
        test/VectorSpeedTest.cpp)

# Generates targets: rsa common code library, rsa test library, rsa exec, rsa test exec for given spatial and angular
# dimenstions. It is done so multiple binaries can coexist without recompiling every time parameters are changed.
function(register_targets_for_dimensions SD AD)
    set(DIM_SUFFIX ".${SD}.${AD}")

    set(RSA_LIB "librsa${DIM_SUFFIX}")
    set(TEST_LIB "libtest${DIM_SUFFIX}")
    set(RSA_EXEC "rsa${DIM_SUFFIX}")
    set(RSA_TEST_EXEC "rsa_test${DIM_SUFFIX}")

    add_library(${RSA_LIB} ${MAIN_SOURCE})
    target_compile_definitions(${RSA_LIB} PUBLIC RSA_SPATIAL_DIMENSION=${SD} RSA_ANGULAR_DIMENSION=${AD})

    add_library(${TEST_LIB} ${TEST_SOURCE})
    target_compile_definitions(${TEST_LIB} PUBLIC RSA_SPATIAL_DIMENSION=${SD} RSA_ANGULAR_DIMENSION=${AD})

    add_executable(${RSA_EXEC} rsa3d/Main.cpp)
    target_include_directories(${RSA_EXEC} PUBLIC statistics)
    target_link_libraries(${RSA_EXEC} ${RSA_LIB} libstatistics)

    add_executable(${RSA_TEST_EXEC} test/Main.cpp)
    target_link_libraries(${RSA_TEST_EXEC} ${TEST_LIB} ${RSA_LIB})
endfunction(register_targets_for_dimensions)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wall -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")   # -O3 in relase, in debug no optimisations
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp -static")    # link dependencies statically

add_library(libstatistics ${STATISTICS_SOURCE})

# Here come all dimensions for which binaries should be generated
register_targets_for_dimensions(3 0)    # 3D spheres, cuboids, platonic solids and so on
register_targets_for_dimensions(2 0)    # 2D disks
register_targets_for_dimensions(2 1)    # AnisotropicShape2D's