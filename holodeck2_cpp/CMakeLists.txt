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
project(holodeck2_cpp LANGUAGES CXX)

# Force x86_64 if you're using x86_64 Python (like Anaconda)
set(CMAKE_OSX_ARCHITECTURES "x86_64")
set(CMAKE_CXX_STANDARD 20)

set(DATA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data")


add_compile_definitions(DEBUG)



# ==============================================================================
# ====    Find Dependencies    ====
# ==============================================================================

# Get Python include path and extension suffix
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)


# ==============================================================================
# ====    Setup Module    ====
# ==============================================================================

add_library(holodeck2_cpp MODULE src/pymodule_holodeck2_cpp.cpp)

target_compile_definitions(holodeck2_cpp PUBLIC PATH_DATA_DIR="${DATA_PATH}")

# Determine appropriate file-suffix for shared-library
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
    OUTPUT_VARIABLE PYTHON_EXT_SUFFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "PYTHON_EXT_SUFFIX = ${PYTHON_EXT_SUFFIX}")

set_target_properties(holodeck2_cpp PROPERTIES
    PREFIX ""                               # remove (default) 'lib' prefix
    SUFFIX "${PYTHON_EXT_SUFFIX}"           # add appropriate suffix
    MACOSX_BUNDLE OFF                       # We need this to be a shared-library
    LINK_FLAGS "-undefined dynamic_lookup"  # avoid conflicts with libpython.dylib
)

# Include Python headers
target_include_directories(holodeck2_cpp PRIVATE ${Python3_INCLUDE_DIRS})

# Make symlink to CWD
add_custom_command(TARGET holodeck2_cpp POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            $<TARGET_FILE:holodeck2_cpp>
            ${CMAKE_CURRENT_SOURCE_DIR}/$<TARGET_FILE_NAME:holodeck2_cpp>
    COMMENT "Symlinking built Python extension to source dir"
)


# ==============================================================================
# ====    Setup Executables    ====
# ==============================================================================

add_executable(holodeck2cpp
    src/main.cpp
    src/sam.cpp
    src/physics.cpp
    src/cosmology.cpp
)

target_compile_definitions(holodeck2cpp PUBLIC PATH_DATA_DIR="${DATA_PATH}")

# Make symlink to executable in CWD
add_custom_command(TARGET holodeck2cpp POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            $<TARGET_FILE:holodeck2cpp>
            ${CMAKE_CURRENT_SOURCE_DIR}/$<TARGET_FILE_NAME:holodeck2cpp>
    COMMENT "Symlinking built c++ executable to source dir"
)

# ---- Tests

add_executable(tester
    tests/tester.cpp
    src/cosmology.cpp
    src/physics.cpp
)

target_compile_definitions(tester PUBLIC PATH_DATA_DIR="${DATA_PATH}")
