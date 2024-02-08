cmake_minimum_required(VERSION 3.14)

project(cheesemapTests LANGUAGES CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)

# ---- Dependencies ----

if (PROJECT_IS_TOP_LEVEL)
    find_package(cheesemap REQUIRED)
    enable_testing()
endif ()

# ---- Tests ----

# Find package for GTest
find_package(GTest REQUIRED)

# Get the names of all the tests from the source files
file(GLOB_RECURSE test_sources CONFIGURE_DEPENDS source/*.cpp)
foreach (test_source IN LISTS test_sources)
    get_filename_component(test_name ${test_source} NAME_WE)
    list(APPEND TEST_NAMES ${test_name})
endforeach ()

# Add a test for each source file
foreach (test_name IN LISTS TEST_NAMES)
    add_executable(${test_name} source/${test_name}.cpp)
    target_include_directories(${test_name} PRIVATE "include/")
    target_link_libraries(${test_name} PRIVATE cheesemap::cheesemap)
    target_link_libraries(${test_name} PRIVATE GTest::GTest GTest::Main)
    target_compile_features(${test_name} PRIVATE cxx_std_20)
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach ()

# ---- End-of-file commands ----

add_folders(Test)
