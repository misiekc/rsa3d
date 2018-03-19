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
       rsa3d/Main \
       rsa3d/OrientedFace \
       rsa3d/Parameters \
       rsa3d/RND \
       rsa3d/ShapeFactory \
       rsa3d/Surface \
       rsa3d/Utils \
       rsa3d/analizator/Analyzer \
       rsa3d/shapes/cube_strategies/MineOverlap \
       rsa3d/shapes/cube_strategies/OptimizedSATOverlap \
       rsa3d/shapes/cube_strategies/SATOverlap \
       rsa3d/shapes/cube_strategies/TriTriOverlap \
       rsa3d/shapes/Cuboid \
       rsa3d/shapes/Ellipse \
       rsa3d/shapes/SpheroCylinder2D \
       rsa3d/surfaces/NBoxFBC \
       rsa3d/surfaces/NBoxPBC \
       test/utility/BallFactory \
       test/utility/BoxFactory \
       test/utility/Quantity \
       test/utility/Timer \
       test/AnisotropicShape2DExclusionTest \
       test/CuboidIntTest \
       test/CuboidPointInsideTest \
       test/CuboidSpeedTest \
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
OBJS_STATD = $(OBJS_STAT:%=$(OBJDIR)/%.o)
DEPS = $(OBJS:%=$(DEPSDIR)/%.d)

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
$(EXEC): $(OBJSD) $(LIBSTAT)
	@echo '## LINKING'
	$(CC) -o $@ $^ $(LFLAGS)
	
$(LIBSTAT): $(OBJS_STATD)
	@echo '## MAKING $(LIBSTAT)'
	ar crf $@ $^

# including dependencies - after $(EXEC) so make with no targets
# refer to executable linking
#----------------------------------------------------------------
-include $(DEPS)

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
clean_all:
	rm -rf $(EXEC) $(LIBSTAT) $(BUILD)

clean:
	rm -rf $(BUILD)

