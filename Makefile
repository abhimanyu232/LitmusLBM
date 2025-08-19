CXX = g++
CXXFLAGS = -O3 -g -Wall -std=c++20 -fopenmp
TARGET = lbm.x
SRC = lbm_ai.cpp

# O0 -> O1 : approx 1100% improvement or more. 
# O1 -> O2 : approx 20-25% improvement or more. 
# O2 -> O3 : approx 30% improvement or more. 
# -march=native -flto -funroll-loops -> currently no difference on top of O3

.PHONY: all clean run visualize simulate

all: $(TARGET)
	
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

run: $(TARGET)
	./$(TARGET)

visualize:
	python3 visualize_lbm.py

clean:
	rm -f $(TARGET)
	rm -rf lbm_ai/velocity_*.txt lbm_ai/vorticity_*.txt
	rm -rf plots/velocity_field_*.png  plots/vorticity_field_*.png plots/combined_fields_*.png 
	rm plot/*.mp4

# Run everything: compile, simulate, and visualize
simulate: clean run visualize