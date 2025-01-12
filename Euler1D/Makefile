# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -g -std=c++14 -pedantic -O3

# Target executable
EXEC = run

# Source files
SRCS = euler1DSolver.C EoS.C euler1D.C method.C material.C

# Object files (replace .C with .o)
OBJS = $(SRCS:.C=.o)

# Header files
DEPS = euler1DSolver.H EoS.H euler1D.H method.H material.H

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
