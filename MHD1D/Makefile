# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -g -std=c++14 -pedantic -O3

# Target executable
EXEC = run

# Source files 
SRCS = euler1D.C MHD1DSolver.C method.C EoS.C material.C

# Object files (replace .C with .o)
OBJS = $(SRCS:.C=.o)

# Header files
DEPS = euler1D.H MHD1DSolver.H method.H EoS.H material.H

# Default target
all: $(EXEC)

# Link object files to create the executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(EXEC) 

# Compile source files into object files
%.o: %.C $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule to remove generated files
.PHONY: clean
clean:
	rm -f $(OBJS) $(EXEC)
