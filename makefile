# Source Directory
#SRCDIR = src/

# Put compiler command in FC

FC=mpiifort

# Put compile time options in FFLAGS


#FFLAGS = -g -fbacktrace -fno-align-commons -fbounds-check 
#FFLAGS = -g -fdefault-real-8 -fbacktrace -fno-align-commons -fbounds-check#FFLAGS = -fdefault-real-8   

FFLAGS = -Ofast -align -real-size 64 -mcmodel=large
#FFLAGS = -g -O0 -traceback -check all -debug all -align -real-size 64 -mcmodel=large

OBJ = main.o \
      hatnderivatives6thOrder.o \
      derivatives6thOrder.o \
      derivatives6thOrderPBC.o \
      energy.o \
      evolveeuler.o \
      fluxesNumRelSO3.o \
      covariantDerivsSO3.o \
      fieldStrengthsSO3.o \
      icMMbarTwistedNumerical.o \
      ghostsend.o \
      ghostrecvl.o \
      ghostrecvu.o \
      printenergy.o \
      evolveloop.o \
      processCoord.o \
      icgamma.o	\
      boundaryIndices.o \
####################################
#ghost cells seem to be phased out
#      ghostcellsX.o \
      ghostcellsY.o \
      ghostcellsZ.o \
      hatnghostcellsX.o \
      hatnghostcellsY.o \
      hatnghostcellsZ.o \
#####################################
#Below are evolution subroutines
#      emain.o \
      eevolveeuler.o \
      eghostcellsX.o \
      eghostcellsY.o \
      eghostcellsZ.o \
      eleapforward.o \
      eaverageforhalfstep.o \
      echangename.o \
      ederivatives6thOrder.o \
      ederivatives6thOrderPBC.o \
      ecovariantDerivs.o \
      eenergy.o \
      eprintenergy.o \
      effluxesNumRel.o \
      efieldStrengths.o \
      efluxesNumRel.o \
      einitialconditions.o \

# Program name
PROGRAM = cnSO3.out

# Default make command is all
all: $(PROGRAM) 

# Create the program executable
$(PROGRAM): $(OBJ) 
	$(FC) -o  $@ $^ 
	rm *.o

# Object files (.o) are generated when .f files are compiled
%.o: %.f
	$(FC) $(FFLAGS) -c $*.f

# Clean all the junk
clean:
	rm: -f *.o *.dvi *.aux *.log core.* *~

# In this file, FC, FFLAGS, and PROGRAM are macros (variables) that can be used 
# throughout the file. If we just type make in the terminal, the command all in 
# this file be executed. $(PROGRAM) in all command is it's dependency and this 
# will cause $(PROGRAM) command to execute which, in turn, depends on the 
# object files $(OBJ). $(PROGRAM) commands is then executed. But $(OBJ) again 
# have dependencies on %.f files to be generated. So then the %.o command is 
# executed. Flag -c says to generate the object file, the -o $@ says to put 
# output of the compilation in the file named on the left side of the :, $< is 
# the first term in the dependency list. $@ and $^ are left and right sides of 
# :, respectively.
