#pragma once
#include <assert.h>
#include <limits>
#include <stdlib.h>

#define CONSTEXPR_FAIL(...) __builtin_unreachable()

inline size_t NextPowerOfTwo(size_t n) {
  assert(n != 0);
  return n == 1 ? 1 : 1ULL << (64 - __builtin_clzll(n - 1));
}

using mod_t = long long;

constexpr bool IsPrime(size_t p) {
  for (size_t a = 2; a * a <= p; a += 2) {
    if (p % a == 0)
      return false;
  }
  return true;
}

template <size_t mod> constexpr mod_t ModFastExp(mod_t v, size_t exp) {
  static_assert(mod <= std::numeric_limits<size_t>::max() / mod,
                "Modulo is too big!");
  mod_t ret = 1;
  while (exp > 0) {
    if (exp % 2 == 1) {
      ret = (ret * v) % mod;
    }
    v = (v * v) % mod;
    exp >>= 1;
  }
  return ret;
}

template <size_t mod> constexpr mod_t ModInverse(mod_t a) {
  static_assert(IsPrime(mod), "Modulo must be prime!");
  static_assert(mod <= std::numeric_limits<size_t>::max() / mod,
                "Modulo is too big!");
  return ModFastExp<mod>(a, mod - 2);
}
