################################################
# COMPILATION AND CONSOLIDATION OF EXECUTABLE: #
#     make rsa                                 #
# IT'S THE DEFAULT TARGET, SO ONE CAN USE      #
#     make                                     #
#                                              #
# C & C OF CUSTOM TEST EXECUTABLE:             #
#     make rsa_test                            #
#                                              #
# MAKING ALL BINARIES:                         #
#     make all                                 #
#                                              #
# CLEANING ALL BUILD FILES (make-build EXEC):  #
#     make clean                               #
#                                              #
################################################


#####################
# COMPILER SETTINGS #
#####################

# Compiler to use
CC = g++

# Compiler flags
CFLAGS = -Wall -pedantic -std=c++11 -I"$(CURDIR)/statistics" -O3 -fopenmp -DRSA_SPATIAL_DIMENSION=2 -DRSA_ANGULAR_DIMENSION=1

# Linker flags
LFLAGS = -fopenmp

####################
# PROJECT SETTINGS #
####################

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
           rsa3d/shape/shapes/regular_solid/ \
           rsa3d/shape/shapes/polygon/ \
           rsa3d/surfaces/ \
           statistics/ \
           test/ \
           test/utility/

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
       rsa3d/shape/shapes/cuboid/MineOverlap \
       rsa3d/shape/shapes/cuboid/OptimizedSATOverlap \
       rsa3d/shape/shapes/cuboid/SATOverlap \
       rsa3d/shape/shapes/cuboid/TriTriOverlap \
       rsa3d/shape/shapes/regular_solid/Cuboctahedron \
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
       rsa3d/shape/shapes/polygon/HBPolygon \
       rsa3d/shape/shapes/polygon/SBPolygon \
       rsa3d/shape/shapes/polygon/Polygon \
       rsa3d/shape/shapes/Ellipse \
       rsa3d/shape/shapes/Ellipsoid \
       rsa3d/shape/shapes/Rectangle \
       rsa3d/shape/shapes/SpheroCylinder2D \
       rsa3d/shape/AnisotropicShape2D \
       rsa3d/shape/OrderParameters \
       rsa3d/shape/ShapeFactory \
       rsa3d/surfaces/NBoxFBC \
       rsa3d/surfaces/NBoxPBC

# Objects to be linked only with main executable
OBJS_MAIN = rsa3d/Main \
       		rsa3d/analizator/Analyzer \

# Objects to be linked only with custom tests executable
OBJS_TEST = test/utility/BallFactory \
       		test/utility/BoxFactory \
       		test/utility/Quantity \
       		test/AnisotropicShape2DExclusionDrawer \
       		test/AnisotropicShape2DExclusionTest \
       		test/PackingOverlapsTest \
       		test/CuboidSpeedTest \
       		test/Main \
       		test/OrderParamTest \
       		test/RectangleCase \
       		test/ShapeOverlapTest \
       		test/ShapePointInsideTest \
       		test/ShapeStaticSizesTest \
       		test/VectorSpeedTest

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
OBJS_TEST_D = $(OBJS_TEST:%=$(OBJDIR)/%.o)
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
-include $(OBJS_TEST:%=$(DEPSDIR)/%.d)
-include $(OBJS_STAT:%=$(DEPSDIR)/%.d)

# *.c and *.cpp compilation, generating dependency files
# make sure all necessary diractories are created before compilation
$(OBJDIR)/%.o: %.c
	$(make_dirs)
	$(CC) $(CFLAGS) -c $*.c -o $@ -MMD -MF $(DEPSDIR)/$*.d -MT '$@'

$(OBJDIR)/%.o: %.cpp
	$(make_dirs)
	$(CC) $(CFLAGS) -c $*.cpp -o $@ -MMD -MF $(DEPSDIR)/$*.d -MT '$@'

####################
# TARGETS          #
####################

# main executable
$(EXEC): $(OBJS_COMMON_D) $(OBJS_MAIN_D) $(LIBSTAT)
	@echo '## LINKING $(EXEC)'
	$(CC) -o $@ $^ $(LFLAGS)

# custom test executable
$(EXEC)_test: $(OBJS_COMMON_D) $(OBJS_TEST_D)
	@echo '## LINKING $(EXEC)_test'
	$(CC) -o $@ $^ $(LFLAGS)

# statistics library
$(LIBSTAT): $(OBJS_STAT_D)
	@echo '## MAKING $(LIBSTAT)'
	ar crf $@ $^

.PHONY: all clean doc

# all executables
all: $(EXEC) $(EXEC)_test

# cleaning compilation results
clean:
	rm -rf $(EXEC) $(EXEC)_test $(LIBSTAT) $(BUILD) $(DOCDIR)
	
# documentation
doc:
	@mkdir -p $(DOCDIR)
	@doxygen
