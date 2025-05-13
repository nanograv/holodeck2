# ==============================================================================
#
#
# -- General:
# configure (make build files) :  cmake -S . -B build
#     build (w/ external tool) :  cmake --build build
#
# -- Presets:
# cmake --preset default
# cmake --build --preset default  [--clean-first]
#
# ==============================================================================

cmake_minimum_required(VERSION 3.16)
project(custom LANGUAGES CXX)

# Force x86_64 if you're using x86_64 Python (like Anaconda)
set(CMAKE_OSX_ARCHITECTURES "x86_64")
set(CMAKE_CXX_STANDARD 17)

# Get Python include path and extension suffix
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

add_library(custom MODULE src/custom_pymodule.cpp)

# Determine appropriate file-suffix for shared-library
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
    OUTPUT_VARIABLE PYTHON_EXT_SUFFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "PYTHON_EXT_SUFFIX = ${PYTHON_EXT_SUFFIX}")

set_target_properties(custom PROPERTIES
    PREFIX ""                               # remove (default) 'lib' prefix
    SUFFIX "${PYTHON_EXT_SUFFIX}"           # add appropriate suffix
    MACOSX_BUNDLE OFF                       # We need this to be a shared-library
    LINK_FLAGS "-undefined dynamic_lookup"  # avoid conflicts with libpython.dylib
)

# Include Python headers
target_include_directories(custom PRIVATE ${Python3_INCLUDE_DIRS})
