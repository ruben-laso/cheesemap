cmake_minimum_required(VERSION 3.14)

include(cmake/prelude.cmake)

project(
    cellmap
    VERSION 0.1.0
    DESCRIPTION "Data structure for point clouds"
    HOMEPAGE_URL "https://example.com/"
    LANGUAGES NONE
)

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# ---- Declare library ----

add_library(cellmap_cellmap INTERFACE)
add_library(cellmap::cellmap ALIAS cellmap_cellmap)

set_property(
    TARGET cellmap_cellmap PROPERTY
    EXPORT_NAME cellmap
)

target_include_directories(
    cellmap_cellmap ${warning_guard}
    INTERFACE
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
)

target_compile_features(cellmap_cellmap INTERFACE cxx_std_20)

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES)
  include(cmake/install-rules.cmake)
endif()

# ---- Examples ----

if(PROJECT_IS_TOP_LEVEL)
  option(BUILD_EXAMPLES "Build examples tree." "${cellmap_DEVELOPER_MODE}")
  if(BUILD_EXAMPLES)
    add_subdirectory(example)
  endif()
endif()

# ---- Developer mode ----

if(NOT cellmap_DEVELOPER_MODE)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  message(
      AUTHOR_WARNING
      "Developer mode is intended for developers of cellmap"
  )
endif()

include(cmake/dev-mode.cmake)
