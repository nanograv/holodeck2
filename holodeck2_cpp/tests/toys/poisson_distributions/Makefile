# Makefile to build the Poisson distribution toy example
# This Makefile is used to compile the Poisson distribution toy example in C++
# It includes the necessary flags and libraries for the compilation.
# Run `make` to build the project, `make clean` to remove object files and the executable.
# Run `make run` to execute the compiled program.

CXX = g++

CXXFLAGS = -std=c++23 \
           -Wall \
           -I/usr/local/include \
           -I/opt/homebrew/Cellar/boost/1.88.0/include

LDFLAGS = -lm \
          # -L/usr/local/lib \
          # -L/opt/homebrew/Cellar/boost/1.88.0/lib

BUILDDIR = build
TARGETNAME = poisson_distribution
TARGET = $(BUILDDIR)/$(TARGETNAME)

SOURCES = poisson_distribution.cpp
OBJECTS = $(SOURCES:.cpp=.o)
BUILD_OBJECTS = $(addprefix $(BUILDDIR)/, $(OBJECTS))

all: $(TARGET)

# Rule to link final binary
$(TARGET): $(BUILD_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Rule to compile object files in build dir
$(BUILDDIR)/%.o: %.cpp
	@mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR)

run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run