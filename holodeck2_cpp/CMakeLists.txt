# ==============================================================================
#
#
# -- General:
# configure (make build files) :  cmake -S . -B build
#     build (w/ external tool) :  cmake --build build
#
# -- Presets:
# cmake --preset default
# cmake --build --preset default  [--clean-first] [-DDEBUG_MODE=ON]
#
# ==============================================================================

cmake_minimum_required(VERSION 3.16)
project(holodeck2_cpp LANGUAGES CXX)

# Force x86_64 if you're using x86_64 Python (like Anaconda)
set(CMAKE_OSX_ARCHITECTURES "x86_64")
set(CMAKE_CXX_STANDARD 20)

set(DATA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data")


# Use -DDEBUG_MODE=ON to enable debug mode
# option(DEBUG_MODE "Enable debug mode macros" OFF)

# if(DEBUG_MODE)
#     add_compile_definitions(DEBUG)
# endif()

add_compile_definitions(DEBUG)



# ==============================================================================
# ====    Find Dependencies    ====
# ==============================================================================

# Get Python include path and extension suffix
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# find_package(Qt6 REQUIRED COMPONENTS
#     Widgets  # Covers QApplication, QWidget, QTimer
#     Core     # Covers QTime
#     Gui      # Covers QPainter
# )

# find_package(Boost REQUIRED COMPONENTS
#     program_options      # for parsing configuration input files
# )


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
)

# target_link_libraries(warp PRIVATE
#     # Qt6
#     Qt6::Widgets
#     Qt6::Core
#     Qt6::Gui
#     # Boost
#     Boost::program_options
# )


add_custom_command(TARGET holodeck2cpp POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            $<TARGET_FILE:holodeck2cpp>
            ${CMAKE_CURRENT_SOURCE_DIR}/$<TARGET_FILE_NAME:holodeck2cpp>
    COMMENT "Symlinking built c++ executable to source dir"
)

# ---- Tests

add_executable(tester
    tests/tester.cpp
)

# target_link_libraries(tester PRIVATE holodeck2_cpp)
target_compile_definitions(tester PUBLIC PATH_DATA_DIR="${DATA_PATH}")

# target_compile_definitions(tester PRIVATE
#     PATH_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data"
# )


# target_include_directories(test_cosmology PRIVATE include)
# target_link_libraries(test_cosmology PRIVATE cosmo_core_lib)  # replace with actual target



# Add directories where *header files* can be found for a particular target
# target_include_directories(warp PUBLIC include)

# ---- Variables
# set(MY_VARIABLE "value")
# Use variable as${MY_VARIABLE}.

# ---- Modify C code with CMake variables
# configure_file (
#     "${PROJECT_SOURCE_DIR}/include/My/Version.h.in"
#     "${PROJECT_BINARY_DIR}/include/My/Version.h"
# )

