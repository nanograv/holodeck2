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

# ---- Python

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# ---- HDF5
# By default, Cmake is finding the system HDF5 (from homebrew) which is compiled for arm64, which
# is not compatible with x86_64 (which we're forcing, to work with the current anaconda python).
# To fix this, we need to set the HDF5_ROOT variable to point to the anaconda hdf5 installation.

# Try to detect conda and use its prefix for HDF5
if(DEFINED ENV{CONDA_PREFIX})
    set(CONDA_PREFIX $ENV{CONDA_PREFIX})

    # Anaconda does not include a `/cmake/hdf5` directory by default to tell CMake where to find
    # files, so we need to set them manually:
    set(HDF5_INCLUDE_DIRS "$ENV{CONDA_PREFIX}/include")
    set(HDF5_LIBRARIES "$ENV{CONDA_PREFIX}/lib/libhdf5.dylib")

    # This would work if anaconda did include a cmake directory, but it doesn't:
    # list(PREPEND CMAKE_PREFIX_PATH "${CONDA_PREFIX}")
    # message(STATUS "Adding Conda prefix to cmake path: ${CONDA_PREFIX}")
    # message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")

else()
    message(FATAL_ERROR "No `CONDA_PREFIX` set and no alternative HDF5 config provided.")

endif()

# NOTE: The following line would override the HDF5 paths we manually set above.
# find_package(HDF5 REQUIRED COMPONENTS C)
message(STATUS "HDF5 include: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5 libs:    ${HDF5_LIBRARIES}")


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
target_include_directories(holodeck2_cpp PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(holodeck2_cpp PRIVATE ${HDF5_LIBRARIES})

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
    src/cosmology.cpp
    src/main.cpp
    src/physics.cpp
    src/sam.cpp
    src/utils.cpp
)

target_compile_definitions(holodeck2cpp PUBLIC PATH_DATA_DIR="${DATA_PATH}")

target_include_directories(holodeck2cpp PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(holodeck2cpp PRIVATE ${HDF5_LIBRARIES})


# Make symlink to executable in CWD
add_custom_command(TARGET holodeck2cpp POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            $<TARGET_FILE:holodeck2cpp>
            ${CMAKE_CURRENT_SOURCE_DIR}/$<TARGET_FILE_NAME:holodeck2cpp>
    COMMENT "Symlinking built c++ executable to source dir"
)

# ---- Tests

add_executable(test_holodeck
    tests/tester.cpp
    src/cosmology.cpp
    src/physics.cpp
    src/utils.cpp
)

target_compile_definitions(test_holodeck PUBLIC PATH_DATA_DIR="${DATA_PATH}")
target_include_directories(test_holodeck PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(test_holodeck PRIVATE ${HDF5_LIBRARIES})
