set(RSA_SOURCE
        analizator/Analyzer.cpp
        analizator/DomainAnalyzer.cpp
        analizator/ExclusionZoneVisualizer.cpp
        geometry/xenocollide/BodyBuilder.cpp
        geometry/xenocollide/Collide.cpp
        geometry/xenocollide/CollideGeometry.cpp
        geometry/xenocollide/MapPtr.cpp
        geometry/xenocollide/Quat.cpp
        geometry/Geometry.cpp
        geometry/Matrix.h
        geometry/Matrix.tpp
        geometry/Vector.h
        geometry/Vector.tpp
        modes/AccuracySimulation.cpp
        modes/AnalyzeMode.cpp
        modes/DomainsMode.cpp
        modes/BCExpandMode.cpp
        modes/BoundariesSimulation.cpp
        modes/ContactFunctionMode.cpp
        modes/DatFileGenerationMode.cpp
        modes/DebugMode.cpp
        modes/DefaultSimulation.cpp
        modes/DensitySimulation.cpp
        modes/ExclusionZonesMode.cpp
        modes/PackingTestMode.cpp
        modes/PovrayMode.cpp
        modes/Simulation.cpp
        modes/Saturation.cpp
        modes/WolframMode.cpp
        modes/SurfaceAreaMode.cpp
        surface_functions/FlatSurfaceFunction.h
        surface_functions/SineSurfaceFunction.cpp
        surface_functions/VirtualSineSurfaceFunction.cpp
        utils/Assertions.h
        utils/OMPMacros.h
        utils/Quantity.cpp
        utils/Timer.cpp
        utils/Utils.cpp
        BoundaryConditions.h
        Config.cpp
        CurvedSurface.cpp
        CurvedSurfaceVoxelList.cpp
        boundary_conditions/FreeBC.h
        boundary_conditions/PeriodicBC.h
        NeighbourGrid.h
        OrientedCuboidVoxelList.cpp
        Packing.cpp
        PackingGenerator.cpp
        Parameters.cpp
        Positioned.h
        Positioned.tpp
        ProgramArguments.cpp
        ProgramMode.h
        RND.cpp
        Surface.cpp
        SurfaceFunction.h
        ThreadLocalRND.cpp
        VariableSizeVoxelList.cpp
        Voxel.cpp
        VoxelList.cpp)

add_subdirectory(shape)

function(register_targets_for_dimensions SD AD)
    set(DIM_SUFFIX ".${SD}.${AD}")
    set(RSA_LIB "rsa_main${DIM_SUFFIX}")
    set(SHAPES_LIB "shapes${DIM_SUFFIX}")
    set(RSA_EXEC "rsa${DIM_SUFFIX}")

    add_library(${RSA_LIB} ${RSA_SOURCE})
    target_compile_definitions(${RSA_LIB} PUBLIC RSA_SPATIAL_DIMENSION=${SD} RSA_ANGULAR_DIMENSION=${AD})

    add_executable(${RSA_EXEC} Main.cpp)
    target_include_directories(${RSA_EXEC} PUBLIC statistics)
    target_link_libraries(${RSA_EXEC} ${SHAPES_LIB} ${RSA_LIB} statistics)
endfunction()

register_targets_for_dimensions(3 0)    # 3D spheres, cuboids, platonic solids and so on
register_targets_for_dimensions(2 0)    # 2D disks
register_targets_for_dimensions(2 1)    # AnisotropicShape2D's
register_targets_for_dimensions(1 1)    # paris car parking problem