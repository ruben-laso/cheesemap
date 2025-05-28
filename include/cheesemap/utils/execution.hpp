#pragma once

#if defined(__cpp_lib_execution) && (__cpp_lib_execution >= 201902L)
#include <execution>
#else
#warning "<execution> and parallel algorithms (e.g., par_unseq) not supported by this compiler. Cheesemap will use sequential algorithms."
#endif
