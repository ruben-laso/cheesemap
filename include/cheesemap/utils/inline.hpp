#pragma once

#ifndef CHSINLINE
#if defined(_MSC_VER)
#define CHSINLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define CHSINLINE inline __attribute__((always_inline))
#else
#define CHSINLINE inline
#endif
#endif
