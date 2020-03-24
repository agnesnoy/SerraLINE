# Set up the flags you want in FFLAGS and compiler in FC

# flags for GNU 
#FFLAGS = -O0
FC = gfortran
FFLAGS = -g -Wall -O0 -fcheck=all #-std=f2003

# Suffixes rules to control how .o and .mod files are treated
.SUFFIXES:
.SUFFIXES: .f90 .o .mod 
.f90.o:
	$(FC) $(FFLAGS)  -c $< 
.f90.mod:
	$(FC) $(FFLAGS)  -c $< 

# Objects (OBJ) to be used

OBJ  = SerraLINE.f90 parms.o io_mod.o functions_mod.o
OBJ2 = Extract.f90 parms.o io_mod.o

#Compilation
all:	SerraLINE Extract

SerraLINE:	$(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $@

Extract:	$(OBJ2)
	$(FC) $(FFLAGS) $(OBJ2) -o $@

parms.o parms.mod: parms.f90
	$(FC) $(FFLAGS) -c parms.f90

#Instruction to clean mods and executables
clean:
	rm *.o *.mod SerraLINE Extract

