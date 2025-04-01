# Compiler
FC = gfortran

# Compiler flags
FFLAGS = -g -cpp

# Find all .F files and corresponding .o files
SRCS = $(wildcard *.F)
OBJS = $(filter-out SGPDR.o, $(SRCS:.F=.o))

# Default target
SGPDR: SGPDR.o $(OBJS)
	$(FC) $(FFLAGS) -o $@ SGPDR.o $(OBJS)

# Compile all .F files into .o files
%.o: %.F
	$(FC) $(FFLAGS) -c $<

# Clean up
clean:
	rm -f *.o SGPDR
