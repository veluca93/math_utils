#pragma once
#include <assert.h>
#include <stdlib.h>

inline size_t NextPowerOfTwo(size_t n) {
  assert(n != 0);
  return n == 1 ? 1 : 1ULL << (64 - __builtin_clzll(n - 1));
}
