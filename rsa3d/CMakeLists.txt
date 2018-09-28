set(RSA_SOURCE
        analizator/Analyzer.cpp
        analizator/ExclusionZoneVisualizer.cpp
        shape/shapes/cuboid/Cuboid.cpp
        shape/shapes/cuboid/CuboidOverlapStrategy.h
        shape/shapes/cuboid/MineOverlapCB.cpp
        shape/shapes/cuboid/OptimizedSATOverlapCB.cpp
        shape/shapes/cuboid/SATOverlapCB.cpp
        shape/shapes/cuboid/TriTriOverlapCB.cpp
        shape/shapes/regular_solid/Cuboctahedron.cpp
        shape/shapes/regular_solid/CubeToTetrahedron.cpp
        shape/shapes/regular_solid/Dodecahedron.cpp
        shape/shapes/regular_solid/Icosahedron.cpp
        shape/shapes/regular_solid/Icosidodecahedron.cpp
        shape/shapes/regular_solid/Octahedron.cpp
        shape/shapes/regular_solid/RegularSolid.h
        shape/shapes/regular_solid/RegularSolid.tpp
        shape/shapes/regular_solid/RegularSolidBase.cpp
        shape/shapes/regular_solid/Rhombicosidodecahedron.cpp
        shape/shapes/regular_solid/Rhombicuboctahedron.cpp
        shape/shapes/regular_solid/SATOverlapRS.cpp
        shape/shapes/regular_solid/SnubCube.cpp
        shape/shapes/regular_solid/SnubDodecahedron.cpp
        shape/shapes/regular_solid/Tetrahedron.cpp
        shape/shapes/regular_solid/TriTriOverlapRS.cpp
        shape/shapes/regular_solid/TruncatedCube.cpp
        shape/shapes/regular_solid/TruncatedCuboctahedron.cpp
        shape/shapes/regular_solid/TruncatedDodecahedron.cpp
        shape/shapes/regular_solid/TruncatedIcosahedron.cpp
        shape/shapes/regular_solid/TruncatedIcosidodecahedron.cpp
        shape/shapes/regular_solid/TruncatedOctahedron.cpp
        shape/shapes/regular_solid/TruncatedTetrahedron.cpp
        shape/shapes/regular_solid/UnoptimizedSATOverlapRS.cpp
        shape/shapes/polygon/HBPolygon.cpp
        shape/shapes/polygon/SBPolygon.cpp
        shape/shapes/Ellipse.cpp
        shape/shapes/Ellipsoid.cpp
        shape/shapes/OrientedCuboid.h
        shape/shapes/OrientedCuboid.tpp
        shape/shapes/Polydisk.cpp
        shape/shapes/polygon/Polygon.cpp
        shape/shapes/Rectangle.cpp
        shape/shapes/Sphere.h
        shape/shapes/Sphere.tpp
        shape/shapes/SpheroCylinder2D.cpp
        shape/AnisotropicShape2D.cpp
        shape/ConvexShape.h
        shape/ConvexShape.tpp
        shape/OrderCalculable.h
        shape/OrderParameters.cpp
        shape/OverlapStrategy.h
        shape/OverlapStrategyShape.h
        shape/Shape.h
        shape/Shape.tpp
        shape/ShapeFactory.cpp
        surfaces/NBoxFBC.cpp
        surfaces/NBoxPBC.cpp
        BoundaryConditions.h
        Config.cpp
        Intersection.cpp
        Matrix.h
        Matrix.tpp
        NeighbourGrid.h
        Packing.cpp
        PackingGenerator.cpp
        Parameters.cpp
        Positioned.h
        Positioned.tpp
        RND.cpp
        Surface.cpp
        ThreadLocalRND.cpp
        Timer.cpp
        Utils.cpp
        Vector.h
        Vector.tpp
        Voxel.h
        Voxel.cpp
        VoxelList.h
        VoxelList.cpp)

function(register_targets_for_dimensions SD AD)
    set(DIM_SUFFIX ".${SD}.${AD}")
    set(RSA_LIB "rsa_main${DIM_SUFFIX}")

    add_library(${RSA_LIB} ${RSA_SOURCE})
    target_compile_definitions(${RSA_LIB} PUBLIC RSA_SPATIAL_DIMENSION=${SD} RSA_ANGULAR_DIMENSION=${AD})

    set(RSA_EXEC "rsa${DIM_SUFFIX}")
    add_executable(${RSA_EXEC} Main.cpp)
    target_include_directories(${RSA_EXEC} PUBLIC statistics)
    target_link_libraries(${RSA_EXEC} ${RSA_LIB} statistics)
endfunction()

register_targets_for_dimensions(3 0)    # 3D spheres, cuboids, platonic solids and so on
register_targets_for_dimensions(2 0)    # 2D disks
register_targets_for_dimensions(2 1)    # AnisotropicShape2D's