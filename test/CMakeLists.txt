cmake_minimum_required(VERSION 3.14)

project(cellmapTests LANGUAGES CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)

# ---- Dependencies ----

if (PROJECT_IS_TOP_LEVEL)
    find_package(cellmap REQUIRED)
    enable_testing()
endif ()

# ---- Tests ----

add_executable(cellmap_test source/cellmap_test.cpp)
target_link_libraries(cellmap_test PRIVATE cellmap::cellmap)
target_compile_features(cellmap_test PRIVATE cxx_std_20)

add_test(NAME cellmap_test COMMAND cellmap_test)

# ---- End-of-file commands ----

add_folders(Test)
