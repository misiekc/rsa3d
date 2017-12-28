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
CFLAGS = -Wall -pedantic -std=c++14 -I"$(CURDIR)/statistics" -O3

# Linker flags 
LFLAGS =

####################
# PROJECT SETTINGS #
####################

# Executable name
EXEC = rsa

# Libraries
LIBSTAT = libstat.a

# Source files list (withous extensions)
OBJS = rsa3d/BoundaryConditions \
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
       rsa3d/shapes/Cuboid \
       rsa3d/shapes/cube_strategies/MineOverlap \
       rsa3d/shapes/cube_strategies/OptimizedSATOverlap \
       rsa3d/shapes/cube_strategies/SATOverlap \
       rsa3d/shapes/cube_strategies/TriTriOverlap \
       rsa3d/surfaces/NBoxFBC \
       rsa3d/surfaces/NBoxPBC \
       rsa3d/tests/utility/BallFactory \
       rsa3d/tests/utility/BoxFactory \
       rsa3d/tests/utility/Quantity \
       rsa3d/tests/utility/Timer \
       rsa3d/tests/CuboidIntTest \
       rsa3d/tests/CuboidPointInsideTest \
       rsa3d/tests/CuboidSpeedTest \
       rsa3d/tests/TriangleIntTest
       
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
          rsa3d/shapes/cube_strategies \
          rsa3d/surfaces/ \
          rsa3d/tests/ \
          rsa3d/tests/utility \
          statistics/

# OBJ files folder
OBJDIR = build/obj

# D files folder
DEPSDIR = build/deps

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
	rm -rf $(EXEC) $(LIBSTAT) build

clean:
	rm -rf build

