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
CFLAGS = -Wall -pedantic -std=c++14 -I"$(CURDIR)/statistics" -O3 -fopenmp

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
           rsa3d/shapes/ \
           rsa3d/shapes/cube_strategies/ \
           rsa3d/surfaces/ \
           statistics/ \
           test/ \
           test/utility/

# Source files list (without extensions) except of the Main.cpp file
# These objects will be linked with both main executable and test
# executables
OBJS_COMMON = rsa3d/AnisotropicShape2D \
       rsa3d/BoundaryConditions \
       rsa3d/Config \
       rsa3d/Intersection \
       rsa3d/OrientedFace \
       rsa3d/Parameters \
       rsa3d/RND \
       rsa3d/ShapeFactory \
       rsa3d/Surface \
       rsa3d/Timer \
       rsa3d/Utils \
       rsa3d/shapes/cube_strategies/MineOverlap \
       rsa3d/shapes/cube_strategies/OptimizedSATOverlap \
       rsa3d/shapes/cube_strategies/SATOverlap \
       rsa3d/shapes/cube_strategies/TriTriOverlap \
       rsa3d/shapes/Cuboid \
       rsa3d/shapes/Ellipse \
       rsa3d/shapes/Rectangle \
       rsa3d/shapes/SpheroCylinder2D \
       rsa3d/surfaces/NBoxFBC \
       rsa3d/surfaces/NBoxPBC \

# Objects to be linked only with main executable
OBJS_MAIN = rsa3d/Main \
       		rsa3d/analizator/Analyzer \

# Objects to be linked only with custom tests executable
OBJS_TEST = test/utility/BallFactory \
       		test/utility/BoxFactory \
       		test/utility/Quantity \
       		test/AnisotropicShape2DExclusionTest \
       		test/AnisotropicShape2DExclusionDrawer \
       		test/CuboidIntTest \
       		test/CuboidPointInsideTest \
       		test/CuboidSpeedTest \
       		test/Main \
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
	@mkdir $(DOCDIR)
	@doxygen
