cmake_minimum_required(VERSION 3.8)
project(rsa3d)

set(CMAKE_CXX_STANDARD 11)

set(MAIN_SOURCE
        rsa3d/analizator/Analyzer.cpp
        rsa3d/shapes/cube_strategies/OptimizedSATOverlap.cpp
        rsa3d/shapes/cube_strategies/CuboidOverlapStrategy.h
        rsa3d/shapes/cube_strategies/MineOverlap.cpp
        rsa3d/shapes/cube_strategies/SATOverlap.cpp
        rsa3d/shapes/cube_strategies/TriTriOverlap.cpp
        rsa3d/shapes/platonic_solid/PlatonicSolid.h
        rsa3d/shapes/platonic_solid/PlatonicSolid.tpp
        rsa3d/shapes/platonic_solid/RegularTetrahedron.cpp
        rsa3d/shapes/platonic_solid/RegularOctahedron.cpp
        rsa3d/shapes/platonic_solid/RegularDodecahedron.cpp
        rsa3d/shapes/platonic_solid/RegularIcosahedron.cpp
        rsa3d/shapes/polygons/Segment2D.cpp
        rsa3d/shapes/Cuboid.cpp
        rsa3d/shapes/Ellipse.cpp
        rsa3d/shapes/OrientedCuboid.h
        rsa3d/shapes/polygons/Polygon.cpp
        rsa3d/shapes/Rectangle.cpp
        rsa3d/shapes/Sphere.h
        rsa3d/shapes/SpheroCylinder2D.cpp
        rsa3d/surfaces/NBoxFBC.cpp
        rsa3d/surfaces/NBoxPBC.cpp
        rsa3d/AnisotropicShape2D.cpp
        rsa3d/BoundaryConditions.cpp
        rsa3d/Config.cpp
        rsa3d/Intersection.cpp
        rsa3d/Matrix.h
        rsa3d/NeighbourGrid.h
        rsa3d/OrientedFace.cpp
        rsa3d/PackingGenerator.h
        rsa3d/Parameters.cpp
        rsa3d/Positioned.h
        rsa3d/RND.cpp
        rsa3d/Shape.h
        rsa3d/ShapeFactory.cpp
        rsa3d/Surface.cpp
        rsa3d/Timer.cpp
        rsa3d/Utils.cpp
        rsa3d/Vector.h
        rsa3d/Voxel.h
        rsa3d/VoxelList.h rsa3d/shapes/OverlapStrategyShape.h rsa3d/shapes/OverlapStrategy.h rsa3d/shapes/platonic_solid/SATOverlap.tpp rsa3d/shapes/platonic_solid/SATOverlap.h rsa3d/shapes/platonic_solid/TriTriOverlap.tpp rsa3d/shapes/platonic_solid/TriTriOverlap.h)

set(STATISTICS_SOURCE
        statistics/ASFRegression.cpp
        statistics/LinearRegression.cpp
        statistics/LogPlot.cpp
        statistics/Plot.cpp
        statistics/PowerRegression.cpp)

set(TEST_SOURCE
        test/utility/BallFactory.cpp
        test/utility/BoxFactory.cpp
        test/utility/ShapePairFactory.h
        test/utility/MockBC.h
        test/utility/Quantity.cpp
        test/AnisotropicShape2DExclusionDrawer.cpp
        test/AnisotropicShape2DExclusionTest.cpp
        test/ShapeOverlapTest.cpp
        test/ShapePointInsideTest.cpp
        test/CuboidSpeedTest.cpp
        test/RectangleCase.cpp
        test/VectorSpeedTest.cpp)

include_directories(statistics)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wall -pedantic")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
#set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -fopenmp")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")

add_library(rsa_src ${MAIN_SOURCE})
add_library(tests ${TEST_SOURCE})
add_library(statistics ${STATISTICS_SOURCE})

add_executable(rsa3d rsa3d/Main.cpp)
target_link_libraries(rsa3d rsa_src statistics)

add_executable(rsa_test test/Main.cpp)
target_link_libraries(rsa_test tests rsa_src)