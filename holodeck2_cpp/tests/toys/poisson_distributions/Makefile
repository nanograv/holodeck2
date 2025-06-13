# Makefile to build the Poisson distribution toy example
# This Makefile is used to compile the Poisson distribution toy example in C++
# It includes the necessary flags and libraries for the compilation.
# Run `make` to build the project, `make clean` to remove object files and the executable.
# Run `make run` to execute the compiled program.

CXX = g++
# CXXFLAGS = -std=c++23 -Wall -Wextra -O2
# LDFLAGS = -lm
CXXFLAGS =  -std=c++23 \
			-Wall \
			-I/usr/local/include \
			-I/opt/homebrew/Cellar/boost/1.88.0/include

LDFLAGS =   -lm \
			# -L/usr/local/lib \
			# -L/opt/homebrew/Cellar/boost/1.88.0/lib

TARGET = poisson_distribution
SOURCES = poisson_distribution.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

run: $(TARGET)
	./$(TARGET)

.PHONY: all clean
