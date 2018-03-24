#############################################
#         GENERAL PURPOSE MAKEFILE          #
#############################################
# COMPILATION AND CONSOLIDATION:            #
#     make                                  #
#                                           #
# CLEANING HELPER FILES (OBJ & D):          #
#     make clean                            #
#                                           #
# CLEANING ALL RESULT OF COMPILATION        #
#     make clean_all                        #
#############################################


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

# Libraries
LIBSTAT = libstat.a

# Source files list (without extensions)
OBJS = rsa3d/AnisotropicShape2D \
       rsa3d/BoundaryConditions \
       rsa3d/Config \
       rsa3d/DynamicMatrix \
       rsa3d/DynamicVector \
       rsa3d/Intersection \
       rsa3d/OrientedFace \
       rsa3d/Parameters \
       rsa3d/RND \
       rsa3d/ShapeFactory \
       rsa3d/Surface \
       rsa3d/Utils \
       rsa3d/shapes/cube_strategies/MineOverlap \
       rsa3d/shapes/cube_strategies/OptimizedSATOverlap \
       rsa3d/shapes/cube_strategies/SATOverlap \
       rsa3d/shapes/cube_strategies/TriTriOverlap \
       rsa3d/shapes/Cuboid \
       rsa3d/shapes/Ellipse \
       rsa3d/shapes/SpheroCylinder2D \
       rsa3d/surfaces/NBoxFBC \
       rsa3d/surfaces/NBoxPBC \

OBJS_MAIN = rsa3d/Main \
       		rsa3d/analizator/Analyzer \

OBJS_TEST = test/utility/BallFactory \
       		test/utility/BoxFactory \
       		test/utility/Quantity \
       		test/utility/Timer \
       		test/AnisotropicShape2DExclusionTest \
       		test/CuboidIntTest \
       		test/CuboidPointInsideTest \
       		test/CuboidSpeedTest \
       		test/Main \
       		test/TriangleIntTest \
       		test/VectorSpeedTest
       
# Source files list (withous extensions) for statistics lib
OBJS_STAT = statistics/ASFRegression \
            statistics/LinearRegression \
            statistics/LogPlot \
            statistics/Plot \
            statistics/PowerRegression

# A list of all subfolders in the project (all dirs in path
# are created automaticly, they don't have to be listed explicitly)
SUBDIRS = rsa3d/analizator/ \
          rsa3d/shapes/ \
          rsa3d/shapes/cube_strategies/ \
          rsa3d/surfaces/ \
          statistics/ \
          test/ \
          test/utility/

# build folder
BUILD = make-build

# OBJ files folder
OBJDIR = $(BUILD)/obj

# D files folder
DEPSDIR = $(BUILD)/deps

###################################
# MAKEFILE INTERNAL CONFIGURATION #
###################################

OBJSD = $(OBJS:%=$(OBJDIR)/%.o)
OBJS_MAIND = $(OBJS_MAIN:%=$(OBJDIR)/%.o)
OBJS_TESTD = $(OBJS_TEST:%=$(OBJDIR)/%.o)
OBJS_STATD = $(OBJS_STAT:%=$(OBJDIR)/%.o)
DEPS = $(OBJS:%=$(DEPSDIR)/%.d)
DEPS_MAIN = $(OBJS_MAIN:%=$(DEPSDIR)/%.d)
DEPS_TEST = $(OBJS_TEST:%=$(DEPSDIR)/%.d)
DEPS_STAT = $(OBJS_STAT:%=$(DEPSDIR)/%.d)

# creating folders
#-----------------------------------
define make_dirs
@mkdir -p $(OBJDIR)
@mkdir -p $(DEPSDIR)
@for dir in $(SUBDIRS); \
do \
mkdir -p $(OBJDIR)/$$dir && mkdir -p $(DEPSDIR)/$$dir; \
done
endef

####################
# TARGETS          #
####################

# executable linking
#-----------------------------------------------------------------
$(EXEC): $(OBJSD) $(OBJS_MAIND) $(LIBSTAT)
	@echo '## LINKING $(EXEC)'
	$(CC) -o $@ $^ $(LFLAGS)

$(EXEC)_test: $(OBJSD) $(OBJS_TESTD)
	@echo '## LINKING $(EXEC)_test'
	$(CC) -o $@ $^ $(LFLAGS)
	
$(LIBSTAT): $(OBJS_STATD)
	@echo '## MAKING $(LIBSTAT)'
	ar crf $@ $^

all: $(EXEC) $(EXEC)_test


# including dependencies - after $(EXEC) so make with no targets
# refer to executable linking
#----------------------------------------------------------------
-include $(DEPS)
-include $(DEPS_MAIN)
-include $(DEPS_TEST)
-include $(DEPS_STAT)

# *.c and *.cpp compilation, generating dependency files
#----------------------------------------------------------------
$(OBJDIR)/%.o: %.c
	$(make_dirs)
	$(CC) $(CFLAGS) -c $*.c -o $@ -MMD -MF $(DEPSDIR)/$*.d -MT '$@'

$(OBJDIR)/%.o: %.cpp
	$(make_dirs)
	$(CC) $(CFLAGS) -c $*.cpp -o $@ -MMD -MF $(DEPSDIR)/$*.d -MT '$@'

# cleaning compilation results
#---------------------------------------------------------------
clean:
	rm -rf $(EXEC) $(EXEC)_test $(LIBSTAT) $(BUILD)
