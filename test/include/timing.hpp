#pragma once

#include <chrono>

namespace chs::timing
{
	static double last_seconds = 0;
	static auto   last_start   = std::chrono::high_resolution_clock::now();
	static auto   last_end     = std::chrono::high_resolution_clock::now();
} // namespace chs::timing

#ifndef TIME_IT
#define TIME_IT(code)                                                                                                  \
	chs::timing::last_start = std::chrono::high_resolution_clock::now();                                           \
	code;                                                                                                          \
	chs::timing::last_end     = std::chrono::high_resolution_clock::now();                                         \
	chs::timing::last_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(chs::timing::last_end -  \
	                                                                                      chs::timing::last_start) \
	                                    .count();
#endif