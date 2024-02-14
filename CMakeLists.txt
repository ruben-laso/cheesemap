cmake_minimum_required(VERSION 3.14)

include(cmake/prelude.cmake)

project(
        cheesemap
        VERSION 0.1.0
        DESCRIPTION "Data structure for point clouds"
        HOMEPAGE_URL "https://example.com/"
        LANGUAGES NONE
)

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# ---- Declare library ----

add_library(cheesemap_cheesemap INTERFACE)
add_library(cheesemap::cheesemap ALIAS cheesemap_cheesemap)

set_property(
        TARGET cheesemap_cheesemap PROPERTY
        EXPORT_NAME cheesemap
)

target_include_directories(
        cheesemap_cheesemap ${warning_guard}
        INTERFACE
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
)

target_compile_features(cheesemap_cheesemap INTERFACE cxx_std_20)

# ---- Dependencies ----

# Required: range-v3
find_package(range-v3 REQUIRED)
target_link_libraries(
        cheesemap_cheesemap ${warning_guard}
        INTERFACE
        range-v3::range-v3
)

# Required: armadillo
find_package(Armadillo REQUIRED)
target_link_libraries(
        cheesemap_cheesemap ${warning_guard}
        INTERFACE
        armadillo
)

# ---- Install rules ----

if (NOT CMAKE_SKIP_INSTALL_RULES)
    include(cmake/install-rules.cmake)
endif ()

# ---- Examples ----

if (PROJECT_IS_TOP_LEVEL)
    option(BUILD_EXAMPLES "Build examples tree." "${cheesemap_DEVELOPER_MODE}")
    if (BUILD_EXAMPLES)
        add_subdirectory(example)
    endif ()
endif ()

# ---- Developer mode ----

if (NOT cheesemap_DEVELOPER_MODE)
    return()
elseif (NOT PROJECT_IS_TOP_LEVEL)
    message(
            AUTHOR_WARNING
            "Developer mode is intended for developers of cheesemap"
    )
endif ()

include(cmake/dev-mode.cmake)
