################################################
# EXECUTABLES:                                 #
#     rsa, rsa_stat_test, rsa_unit_test        #
#                                              #
# DEFAULT TARGET:                              #
#     rsa                                      #
#                                              #
# DOCUMENTATION:                               #
#     make doc                                 #
#                                              #
# CLEANING ALL BUILD FILES:                    #
#     make clean                               #
#                                              #
################################################


####################
# PROJECT SETTINGS #
####################

# RSA dimensions
RSA_SPATIAL_DIMENSION = 2
RSA_ANGULAR_DIMENSION = 1

# Compiler flags
CFLAGS = -Wall -pedantic -std=c++17 -I"$(CURDIR)/statistics" -I"$(CURDIR)/unit_test/lib" -O3 -fopenmp -DRSA_SPATIAL_DIMENSION=$(RSA_SPATIAL_DIMENSION) -DRSA_ANGULAR_DIMENSION=$(RSA_ANGULAR_DIMENSION)

# Linker flags
LFLAGS = -fopenmp

# Executables name
SUFFIX = .$(RSA_SPATIAL_DIMENSION).$(RSA_ANGULAR_DIMENSION)
EXEC = rsa$(SUFFIX)
STAT_TEST_EXEC = stat_test$(SUFFIX)
UNIT_TEST_EXEC = unit_test$(SUFFIX)

.DEFAULT_GOAL := $(EXEC)

# Statistics library name
LIBSTAT = libstat.a

# All packages in the project
PACKAGES = rsa3d/analizator/ \
           rsa3d/geometry/ \
           rsa3d/shape/ \
           rsa3d/shape/shapes/ \
           rsa3d/shape/shapes/cuboid/ \
           rsa3d/shape/shapes/polygon/ \
           rsa3d/shape/shapes/regular_solid/ \
           rsa3d/shape/shapes/spherocylinder/ \
           rsa3d/surfaces/ \
           rsa3d/utils/ \
           statistics/ \
           stat_test/ \
           stat_test/utils/ \
           unit_test/ \
           unit_test/rsa3d

# Source files list (without extensions) except of the Main.cpp file
# These objects will be linked with both main executable and test
# executables
OBJS_COMMON = rsa3d/Config \
       rsa3d/Packing \
       rsa3d/PackingGenerator \
       rsa3d/Parameters \
       rsa3d/ProgramArguments \
       rsa3d/RND \
       rsa3d/Surface \
       rsa3d/ThreadLocalRND \
       rsa3d/Voxel \
       rsa3d/VoxelList \
       rsa3d/geometry/Geometry \
       rsa3d/shape/shapes/cuboid/Cuboid \
       rsa3d/shape/shapes/cuboid/MineOverlapCB \
       rsa3d/shape/shapes/cuboid/OptimizedSATOverlapCB \
       rsa3d/shape/shapes/cuboid/SATOverlapCB \
       rsa3d/shape/shapes/cuboid/TriTriOverlapCB \
       rsa3d/shape/shapes/regular_solid/Cuboctahedron \
       rsa3d/shape/shapes/regular_solid/CubeToTetrahedron \
       rsa3d/shape/shapes/regular_solid/Dodecahedron \
       rsa3d/shape/shapes/regular_solid/Icosahedron \
       rsa3d/shape/shapes/regular_solid/Icosidodecahedron \
       rsa3d/shape/shapes/regular_solid/Octahedron \
       rsa3d/shape/shapes/regular_solid/RegularSolidBase \
       rsa3d/shape/shapes/regular_solid/Rhombicosidodecahedron \
       rsa3d/shape/shapes/regular_solid/Rhombicuboctahedron \
       rsa3d/shape/shapes/regular_solid/SATOverlapRS \
       rsa3d/shape/shapes/regular_solid/SnubCube \
       rsa3d/shape/shapes/regular_solid/SnubDodecahedron \
       rsa3d/shape/shapes/regular_solid/Tetrahedron \
       rsa3d/shape/shapes/regular_solid/TriTriOverlapRS \
       rsa3d/shape/shapes/regular_solid/TruncatedCube \
       rsa3d/shape/shapes/regular_solid/TruncatedCuboctahedron \
       rsa3d/shape/shapes/regular_solid/TruncatedDodecahedron \
       rsa3d/shape/shapes/regular_solid/TruncatedIcosahedron \
       rsa3d/shape/shapes/regular_solid/TruncatedIcosidodecahedron \
       rsa3d/shape/shapes/regular_solid/TruncatedOctahedron \
       rsa3d/shape/shapes/regular_solid/TruncatedTetrahedron \
       rsa3d/shape/shapes/regular_solid/UnoptimizedSATOverlapRS \
       rsa3d/shape/shapes/spherocylinder/SpheroCylinder2D \
       rsa3d/shape/shapes/spherocylinder/Stolen2DOverlapSC \
       rsa3d/shape/shapes/polygon/HBPolygon \
       rsa3d/shape/shapes/polygon/SBPolygon \
       rsa3d/shape/shapes/polygon/Polygon \
       rsa3d/shape/shapes/polygon/Triangle \
       rsa3d/shape/shapes/Ellipse \
       rsa3d/shape/shapes/Ellipsoid \
       rsa3d/shape/shapes/Polydisk \
       rsa3d/shape/shapes/Rectangle \
       rsa3d/shape/AnisotropicShape2D \
       rsa3d/shape/OrderParameters \
       rsa3d/shape/ShapeFactory \
       rsa3d/surfaces/NBoxFBC \
       rsa3d/surfaces/NBoxPBC \
       rsa3d/utils/Quantity \
       rsa3d/utils/Timer \
       rsa3d/utils/Utils

# Objects to be linked only with main executable
OBJS_MAIN = rsa3d/Main \
       		rsa3d/analizator/Analyzer \
       		rsa3d/analizator/ExclusionZoneVisualizer \

# Objects to be linked only with statistical tests executable
OBJS_STAT_TEST = stat_test/utils/IndependentPairFactory \
            stat_test/utils/ParallelPairFactory \
            stat_test/utils/UniformBallDistribution \
            stat_test/utils/UniformBoxDistribution\
       		stat_test/AnisotropicShape2DExclusionDrawer \
       		stat_test/AnisotropicShape2DExclusionTest \
       		stat_test/PackingOverlapsTest \
       		stat_test/CuboidSpeedTest \
       		stat_test/Main \
       		stat_test/OrderParamTest \
       		stat_test/ShapeBCTest \
       		stat_test/ShapeOverlapTest \
       		stat_test/ShapePointInsideTest \
       		stat_test/ShapeStaticSizesTest \
       		stat_test/ShapeUltimateTest \
       		stat_test/VectorSpeedTest

# Objects to be linked only with unit tests executable
OBJS_UNIT_TEST = unit_test/rsa3d/MatrixTest \
			unit_test/rsa3d/VectorTest \
            unit_test/rsa3d/PackingTest \
            unit_test/Main \

# Source files list (without extensions) for statistics lib
OBJS_STAT = statistics/ArrayFunction \
            statistics/ASFRegression \
            statistics/LinearRegression \
            statistics/LogPlot \
            statistics/Plot \
            statistics/PowerRegression

# build folder
BUILD = make-build

# OBJ files folder
OBJDIR = $(BUILD)/obj

# D files folder
DEPSDIR = $(BUILD)/deps

# Documentation dir
DOCDIR = doc

# Map all objects from source file lists to dependency files
OBJS_COMMON_D = $(OBJS_COMMON:%=$(OBJDIR)/%.o)
OBJS_MAIN_D = $(OBJS_MAIN:%=$(OBJDIR)/%.o)
OBJS_STAT_TEST_D = $(OBJS_STAT_TEST:%=$(OBJDIR)/%.o)
OBJS_UNIT_TEST_D = $(OBJS_UNIT_TEST:%=$(OBJDIR)/%.o)
OBJS_STAT_D = $(OBJS_STAT:%=$(OBJDIR)/%.o)

# creating folders
define make_dirs
@mkdir -p $(OBJDIR)
@mkdir -p $(DEPSDIR)
@for dir in $(PACKAGES); \
do \
mkdir -p $(OBJDIR)/$$dir && mkdir -p $(DEPSDIR)/$$dir; \
done
endef

####################
# COMPILATION      #
####################

# map all objects from sources lists to dependency folders and
# include dependencies
-include $(OBJS_COMMON:%=$(DEPSDIR)/%.d)
-include $(OBJS_MAIN:%=$(DEPSDIR)/%.d)
-include $(OBJS_STAT_TEST:%=$(DEPSDIR)/%.d)
-include $(OBJS_UNIT_TEST:%=$(DEPSDIR)/%.d)
-include $(OBJS_STAT:%=$(DEPSDIR)/%.d)

# compilation, generating dependency files
# make sure all necessary diractories are created before compilation
$(OBJDIR)/%.o: %.cpp
	$(make_dirs)
	@echo Compiling $*.cpp ...
	@g++ $(CFLAGS) -c $*.cpp -o $@ -MMD -MF $(DEPSDIR)/$*.d -MT '$@'

####################
# TARGETS          #
####################

# main executable
$(EXEC): $(OBJS_COMMON_D) $(OBJS_MAIN_D) $(LIBSTAT)
	@echo 'Linking binary $@'
	@g++ -o $@ $^ $(LFLAGS)

# statistical test executable
$(STAT_TEST_EXEC): $(OBJS_COMMON_D) $(OBJS_STAT_TEST_D)
	@echo 'Linking binary $@'
	@g++ -o $@ $^ $(LFLAGS)

# unit test executable
$(UNIT_TEST_EXEC): $(OBJS_COMMON_D) $(OBJS_UNIT_TEST_D)
	@echo 'Linking binary $@'
	@g++ -o $@ $^ $(LFLAGS)

# statistics library
$(LIBSTAT): $(OBJS_STAT_D)
	@echo 'Linking library $@'
	@ar crf $@ $^

.PHONY: all clean doc

# all executables
all: $(EXEC) $(STAT_TEST_EXEC) $(UNIT_TEST_EXEC)

# cleaning compilation results
clean:
	rm -rf $(EXEC) $(STAT_TEST_EXEC) $(UNIT_TEST_EXEC) $(LIBSTAT) $(BUILD) $(DOCDIR)
	
# documentation
doc:
	@mkdir -p $(DOCDIR)
	@doxygen
