cmake_minimum_required(VERSION 3.14)

project(cellmapExamples CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)

if (PROJECT_IS_TOP_LEVEL)
    find_package(cellmap REQUIRED)
endif ()

add_custom_target(run-examples)

function(add_example NAME)
    add_executable("${NAME}" "${NAME}.cpp")
    target_include_directories("${NAME}" PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
    target_link_libraries("${NAME}" PRIVATE cellmap::cellmap)
    target_compile_features("${NAME}" PRIVATE cxx_std_20)
    add_custom_target("run_${NAME}" COMMAND "${NAME}" VERBATIM)
    add_dependencies("run_${NAME}" "${NAME}")
    add_dependencies(run-examples "run_${NAME}")
endfunction()

add_example(example)

add_folders(Example)
