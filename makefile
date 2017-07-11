#############################################
#    PLIK MAKEFILE OGÓLNEGO PRZEZNACZENIA   #
#############################################
# KOMPILACJA I LINKOWANIE PROGRAMU:         #
#     make                                  #
#                                           #
# CZYSZCZENIE PLIÓW POMOCNICZYCH (OBJ i D): #
#     make clean                            #
#                                           #
# CZYSZCZENIE WSZYSTKICH WYNIKÓW KOMPILACJI #
#     make clean_all                        #
#############################################


##########################
# USTAWIENIA KOMPILATORA #
##########################

# Kompilator
CC = g++

# Flagi kompilacji
CFLAGS = -Wall -pedantic -std=c++11 -I./ -g

# Flagi linkowania
LFLAGS =

#######################
# USTAWIENIA PROJEKTU #
#######################

# Nazwa pliku wykonywalnego
EXEC = rsa

# Lista plików źródłowych (bez rozszerzenia)
OBJS = rsa3d/BoundaryConditions \
       rsa3d/Intersection \
       rsa3d/Main \
       rsa3d/Matrix \
       rsa3d/NeighbourGrid \
       rsa3d/PackingGenerator \
       rsa3d/Parameters \
       rsa3d/Positioned \
       rsa3d/RND \
       rsa3d/Shape \
       rsa3d/ShapeFactory \
       rsa3d/Surface \
       rsa3d/Utils \
       rsa3d/Vector \
       rsa3d/Voxel \
       rsa3d/VoxelList \
       rsa3d/analizator/Analyzer \
       rsa3d/shapes/Sphere \
       rsa3d/shapes/Cuboid \
       rsa3d/surfaces/NBoxPBC \
       rsa3d/tests/CuboidIntTest \
       rsa3d/tests/TriangleIntTest \
       statistics/ASFRegression \
       statistics/LinearRegression \
       statistics/LogPlot \
       statistics/Plot \
       statistics/PowerRegression

# Lista wszystkich podfolderów w projekcie (foldery w ścieżce tworzone
# automatycznie)
SUBDIRS = rsa3d/analizator/ \
          rsa3d/shapes/ \
          rsa3d/surfaces/ \
          rsa3d/tests/ \
          statistics/

# Folder na pliki OBJ
OBJDIR = build/obj

# Folder na pliki D
DEPSDIR = build/deps

####################################
# WEWNĘTRZNA KONFIGURACJA MAKEFILE #
####################################

OBJSD = $(OBJS:%=$(OBJDIR)/%.o)
DEPS = $(OBJS:%=$(DEPSDIR)/%.d)

# tworzenie potrzebnych folderów
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
# CELE             #
####################

# linkowanie programu
#-----------------------------------------------------------------
$(EXEC): $(OBJSD)
	@echo '## LINKING'
	$(CC) -o $@ $^ $(LFLAGS)	

# załączenie zależności - po $(EXEC), aby puste make odnosiło się
# do pliku wykonywalnego
#----------------------------------------------------------------
-include $(DEPS)

# kompilacja plików *.c i *.cpp, tworzenie plików zależności
#----------------------------------------------------------------
$(OBJDIR)/%.o: %.c
	$(make_dirs)
	$(CC) $(CFLAGS) -c $*.c -o $@ -MMD -MF $(DEPSDIR)/$*.d -MT '$@'

$(OBJDIR)/%.o: %.cpp
	$(make_dirs)
	$(CC) $(CFLAGS) -c $*.cpp -o $@ -MMD -MF $(DEPSDIR)/$*.d -MT '$@'

# czyszczenie wyników kompilacji
#---------------------------------------------------------------
clean_all:
	rm -rf $(EXEC) build

clean:
	rm -rf build

