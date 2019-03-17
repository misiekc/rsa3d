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
RSA_ANGULAR_DIMENSION = 0

# Compiler flags
CFLAGS = -Wall -pedantic -std=c++17 -I"$(CURDIR)/statistics" -I"$(CURDIR)/unit_test/lib" -O3 -fopenmp -DRSA_SPATIAL_DIMENSION=$(RSA_SPATIAL_DIMENSION) -DRSA_ANGULAR_DIMENSION=$(RSA_ANGULAR_DIMENSION)

# Linker flags
LFLAGS = -fopenmp

# Executable name
EXEC = rsa
.DEFAULT_GOAL := $(EXEC)

# Statistics library name
LIBSTAT = libstat.a

# All packages in the project
PACKAGES = rsa3d/analizator/ \
           rsa3d/shape/ \
           rsa3d/shape/shapes/ \
           rsa3d/shape/shapes/cuboid/ \
           rsa3d/shape/shapes/polygon/ \
           rsa3d/shape/shapes/regular_solid/ \
           rsa3d/shape/shapes/spherocylinder/ \
           rsa3d/surfaces/ \
           statistics/ \
           stat_test/ \
           stat_test/utility/ \
           unit_test/ \
           unit_test/rsa3d

# Source files list (without extensions) except of the Main.cpp file
# These objects will be linked with both main executable and test
# executables
OBJS_COMMON = rsa3d/Config \
       rsa3d/Intersection \
       rsa3d/Packing \
       rsa3d/PackingGenerator \
       rsa3d/Parameters \
       rsa3d/RND \
       rsa3d/Surface \
       rsa3d/ThreadLocalRND \
       rsa3d/Timer \
       rsa3d/Utils \
       rsa3d/Voxel \
       rsa3d/VoxelList \
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
       rsa3d/shape/shapes/Ellipse \
       rsa3d/shape/shapes/Ellipsoid \
       rsa3d/shape/shapes/Polydisk \
       rsa3d/shape/shapes/Rectangle \
       rsa3d/shape/AnisotropicShape2D \
       rsa3d/shape/OrderParameters \
       rsa3d/shape/ShapeFactory \
       rsa3d/surfaces/NBoxFBC \
       rsa3d/surfaces/NBoxPBC

# Objects to be linked only with main executable
OBJS_MAIN = rsa3d/Main \
       		rsa3d/analizator/Analyzer \
       		rsa3d/analizator/ExclusionZoneVisualizer \

# Objects to be linked only with statistical tests executable
OBJS_STAT_TEST = stat_test/utility/IndependentPairFactory \
            stat_test/utility/ParallelPairFactory \
            stat_test/utility/Quantity \
            stat_test/utility/UniformBallDistribution \
            stat_test/utility/UniformBoxDistribution\
       		stat_test/utility/Quantity \
       		stat_test/AnisotropicShape2DExclusionDrawer \
       		stat_test/AnisotropicShape2DExclusionTest \
       		stat_test/PackingOverlapsTest \
       		stat_test/CuboidSpeedTest \
       		stat_test/Main \
       		stat_test/OrderParamTest \
       		stat_test/RectangleCase \
       		stat_test/ShapeBCTest \
       		stat_test/ShapeOverlapTest \
       		stat_test/ShapePointInsideTest \
       		stat_test/ShapeStaticSizesTest \
       		stat_test/VectorSpeedTest

# Objects to be linked only with unit tests executable
OBJS_UNIT_TEST = unit_test/rsa3d/MatrixTest \
			unit_test/rsa3d/VectorTest \
            unit_test/rsa3d/PackingTest \
            unit_test/Main \

# Source files list (withous extensions) for statistics lib
OBJS_STAT = statistics/ASFRegression \
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
$(EXEC)_stat_test: $(OBJS_COMMON_D) $(OBJS_STAT_TEST_D)
	@echo 'Linking binary $@'
	@g++ -o $@ $^ $(LFLAGS)

# unit test executable
$(EXEC)_unit_test: $(OBJS_COMMON_D) $(OBJS_UNIT_TEST_D)
	@echo 'Linking binary $@'
	@g++ -o $@ $^ $(LFLAGS)

# statistics library
$(LIBSTAT): $(OBJS_STAT_D)
	@echo 'Linking library $@'
	@ar crf $@ $^

.PHONY: all clean doc

# all executables
all: $(EXEC) $(EXEC)_stat_test $(EXEC)_unit_test

# cleaning compilation results
clean:
	rm -rf $(EXEC) $(EXEC)_stat_test $(EXEC)_unit_test $(LIBSTAT) $(BUILD) $(DOCDIR)
	
# documentation
doc:
	@mkdir -p $(DOCDIR)
	@doxygen
