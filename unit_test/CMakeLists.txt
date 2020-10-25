add_library(unit_test_runner lib/catch.hpp lib/fakeit.hpp Main.cpp matchers/VectorApproxMatcher.h)
target_include_directories(unit_test_runner PUBLIC lib)

add_executable(unit_test.2.0
        rsa3d/ConfigTest.cpp
        rsa3d/MatrixTest.cpp
        rsa3d/PackingTest.cpp
        rsa3d/ParametersTest.cpp
        rsa3d/ProgramArgumentsTest.cpp
        rsa3d/VectorTest.cpp
        rsa3d/utils/QuantityTest.cpp rsa3d/surface_functions/VirtualSineSurfaceFunctionTest.cpp)
target_compile_definitions(unit_test.2.0 PUBLIC RSA_SPATIAL_DIMENSION=2 RSA_ANGULAR_DIMENSION=0)
target_link_libraries(unit_test.2.0  unit_test_runner rsa_main.2.0 shapes.2.0 rsa_main.2.0)

add_executable(unit_test.2.1
        rsa3d/shape/shapes/polydisk/PolydiskTest.cpp
        rsa3d/shape/shapes/polygon/PolygonTest.cpp
        rsa3d/shape/shapes/polygon/TriangleTest.cpp rsa3d/shape/shapes/polydisk/KmerTest.cpp rsa3d/shape/shapes/RegularDiskopolygonTest.cpp)
target_compile_definitions(unit_test.2.1 PUBLIC RSA_SPATIAL_DIMENSION=2 RSA_ANGULAR_DIMENSION=1)
target_link_libraries(unit_test.2.1 unit_test_runner rsa_main.2.1 shapes.2.1 rsa_main.2.1)