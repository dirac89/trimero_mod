######################################################################################################
#                TYPE : makefile for gsl
#
#    COMMENTS  :  For compilation and linking of a c++ code of the file "gsl_makefile_in_c++"
#
#     COMMENTS :  For more details about constructing makefiles take a look for example at
#                            http://stackoverflow.com/questions/2481269/how-to-make-simple-c-makefile
#
################################################################################

# source files, objects, headers and libraries _________________________________
SOURCES      = $(wildcard src/*.cpp)
OBJECTS      = $(addprefix build/obj/, $(SOURCES:.cpp=.o))
INCLUDE_PATH = -Isrc/ -I/usr/software/gsl-2.4
LIBS         = -lgsl -lgslcblas 
LIB_PATH     = -L/usr/software/gsl-2.4/

# compiler______________________________________________________________________
COMP         = g++
COMP_FLAGS   = $(INCLUDE_PATH) -O3

# linker________________________________________________________________________
LINK         = g++
LINK_FLAGS   = $(LIB_PATH) -static

# executable____________________________________________________________________
all: build/bin/trimer_exe

create_dirs: 
	mkdir -p build/obj/src/ build/bin/ 
	
build/bin/trimer_exe: create_dirs $(OBJECTS)
	$(LINK) $(OBJECTS) $(LIBS) $(LINK_FLAGS) -o $@

# rule for source to object conversion__________________________________________
build/obj/%.o: %.cpp
	$(COMP) $(COMP_FLAGS) -c $< -o $@

# remove objects and executable_________________________________________________
clean:
	rm -rf build/

# tool that generates dependencies______________________________________________
#dependencies:
#	makedepend $(COMP_FLAGS) $(SOURCES)

# ADD ALSO THE HEADER FILES USED BY EACH OBJECTIVE FILE ... for example
#main.o: include/file1.h
#src/file1.o: include/file1.h /opt/gsl/include/gsl/gsl_complex.h
#src/file1.o: /opt/gsl/include/gsl/gsl_complex_math.h
#src/file2.o: include/file2.h include/file1.h
#src/file3.o: include/file3.h include/file2.h include/file1.h
