# N-Body Gravitational Simulation Makefile
# Author: Makar Ulesov

# Compiler settings
CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic -O3
DEBUG_FLAGS := -g -O0 -DDEBUG

# Target executable
TARGET := nbody
SRCS := main.cpp
OBJS := $(SRCS:.cpp=.o)

# Default target
all: $(TARGET)

# Link
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Debug build
debug: CXXFLAGS += $(DEBUG_FLAGS)
debug: clean $(TARGET)

# Run simulation
run: $(TARGET)
	./$(TARGET)

# Generate visualization (requires Python with matplotlib and Pillow)
plot: simulation_output.csv
	python3 plot.py

# Full pipeline: build, run, and visualize
simulate: run plot

# Clean build artifacts
clean:
	rm -f $(TARGET) $(OBJS) simulation_output.csv

# Clean everything including generated visualizations
distclean: clean
	rm -f orbit_animation.gif

# Check dependencies for Python visualization
check-deps:
	@echo "Checking Python dependencies..."
	@python3 -c "import matplotlib" 2>/dev/null || echo "Missing: matplotlib (pip install matplotlib)"
	@python3 -c "import PIL" 2>/dev/null || echo "Missing: Pillow (pip install Pillow)"
	@python3 -c "import pandas" 2>/dev/null || echo "Missing: pandas (pip install pandas)"
	@echo "Done."

.PHONY: all debug run plot simulate clean distclean check-deps
