CXX = g++
CXXFLAGS_OPT = -O3 -g -Wall -std=c++20 -fopenmp
CXXFLAGS_DEBUG = -O2 -g -Wall -std=c++20 -fopenmp
TARGET = lbm.x
SRC = src/lbm.cpp

# O0 -> O1 : approx 1100% improvement or more. 
# O1 -> O2 : approx 20-25% improvement or more. 
# O2 -> O3 : approx 30% improvement or more. 
# -march=native -flto -funroll-loops -> currently no difference on top of O3

.PHONY: all debug clean run visualize simulate

all: $(TARGET)
	
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS_OPT) -o $(TARGET) $(SRC)

debug: $(SRC)
	$(CXX) $(CXXFLAGS_DEBUG) -o $(TARGET) $(SRC)

run: $(TARGET)
	./$(TARGET)

visualize:
	python3 visualize_lbm.py

clean:
	rm -f $(TARGET)

# Run everything: compile, simulate, and visualize
simulate: clean run visualize