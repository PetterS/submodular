# Petter Strandmark 2011

# You might have to change these to g++-4.5 or something similar if you have several versions installed
COMPILER = g++ -c
LINKER = g++

# Temporary directory to store object files; needs to exist
OBJDIR = ~/obj/submodularcygwin/

# Location of Clp if not installed where gcc can find it
CLPLIBDIR = ~/Programming/cygwin/coin-Clp/lib/
CLPINCLUDEDIR = ~/Programming/cygwin/coin-Clp/include/

# Names of executables
TARGET  = bin/submodular_gcc.exe

##################################################################################

SUBMOBJ = $(OBJDIR)PseudoBoolean.o $(OBJDIR)PseudoBoolean_create_g.o $(OBJDIR)PseudoBoolean_heuristic.o $(OBJDIR)PseudoBoolean_minimize.o $(OBJDIR)PseudoBoolean_reduce.o $(OBJDIR)PseudoBoolean_complete.o
QPBOOBJ = $(OBJDIR)QPBO.o $(OBJDIR)QPBO_extra.o $(OBJDIR)QPBO_maxflow.o $(OBJDIR)QPBO_postprocessing.o 
MAXFLOWOBJ = $(OBJDIR)graph.o $(OBJDIR)maxflow.o 
LIBDIR = $(OBJDIR)
INCLUDE = -I source/library -I thirdparty/Petter -I thirdparty/HOCR -I thirdparty/QPBO -I thirdparty/maxflow-v3.01.src
OPTIONS = -std=c++0x




all : $(TARGET)

$(TARGET): $(OBJDIR)main_program.o $(OBJDIR)submodular_tests.o $(OBJDIR)libpetter.a $(QPBOOBJ) $(MAXFLOWOBJ)
	$(LINKER) -L $(LIBDIR) -L $(CLPLIBDIR) $(OBJDIR)main_program.o  $(OBJDIR)submodular_tests.o $(QPBOOBJ) $(MAXFLOWOBJ) -o $(TARGET) -lpetter -lClp -lCoinUtils 

clean:
	rm -f $(TARGET)
	rm -f $(OBJDIR)*.o
	rm -f $(OBJDIR)*.la
	
test: $(TARGET)
	./$(TARGET)
	
# Main program; compile with -Wall -Werror -pedantic-errors
$(OBJDIR)main_program.o: source/main_program.cpp
	$(COMPILER) $(OPTIONS) -Wall -Werror -pedantic-errors $(INCLUDE) source/main_program.cpp -o $(OBJDIR)main_program.o
$(OBJDIR)submodular_tests.o: source/submodular_tests.cpp
	$(COMPILER) $(OPTIONS) -Wall -Werror -pedantic-errors $(INCLUDE) source/submodular_tests.cpp -o $(OBJDIR)submodular_tests.o
# Main program; -pedantic-errors not possible since it includes QPBO.h
$(OBJDIR)image_denoising.o: source/image_denoising.cpp
	$(COMPILER) $(OPTIONS) -Wall -Werror $(INCLUDE) source/image_denoising.cpp -o $(OBJDIR)image_denoising.o 
    
$(LIBDIR)libpetter.a : $(SUBMOBJ) $(OBJDIR)Petter-Color.o
	ar rcs $(LIBDIR)libpetter.a $(OBJDIR)Petter-Color.o $(SUBMOBJ)

# Library; compile with -Wall -Werror -pedantic-errors
$(OBJDIR)Petter-Color.o: thirdparty/Petter/Petter-Color.cc thirdparty/Petter/Petter-Color.h
	$(COMPILER) $(OPTIONS) -Wall -Werror -pedantic-errors $(INCLUDE) thirdparty/Petter/Petter-Color.cc -o $(OBJDIR)Petter-Color.o	

# Library; compile with -Wall -Werror
# -pedantic-errors not possible since it uses QPBO
$(OBJDIR)%.o: source/library/%.cpp source/library/PseudoBoolean.h
	$(COMPILER) $(OPTIONS) -Wall -Werror $(INCLUDE) -I $(CLPINCLUDEDIR)  $< -o $@

# Will generate warnings in g++ 4.5
$(OBJDIR)%.o: thirdparty/QPBO/%.cpp 
	$(COMPILER) $(OPTIONS) $(INCLUDE) $< -o $@

# Will generate warnings in g++ 4.5
$(OBJDIR)%.o: thirdparty/maxflow-v3.01.src/%.cpp 
	$(COMPILER) $(OPTIONS) $(INCLUDE) $< -o $@

	
