#include "util.h"
#include <vector>

namespace detail {

template <size_t mod> constexpr size_t ModPow2Order(mod_t root) {
  size_t order = 1;
  for (size_t _ = 0; _ < 30; _++, order <<= 1) {
    if (ModFastExp<mod>(root, order) == 1) {
      return order;
    }
  }
  CONSTEXPR_FAIL("Invalid root!");
}

template <size_t mod, mod_t root, bool invert> void ntt(std::vector<mod_t> *a) {
  size_t n = a->size();
  assert(__builtin_popcountll(n) == 1);

  constexpr size_t root_order = ModPow2Order<mod>(root);
  constexpr mod_t root_inv = ModInverse<mod>(root);
  assert(n <= root_order);

  for (size_t i = 1, j = 0; i < n; i++) {
    size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) {
      j ^= bit;
    }
    j ^= bit;

    if (i < j) {
      std::swap((*a)[i], (*a)[j]);
    }
  }

  for (size_t len = 2; len <= n; len <<= 1) {
    mod_t wlen = invert ? root_inv : root;
    for (size_t i = len; i < root_order; i <<= 1) {
      wlen = (mod_t)(1LL * wlen * wlen % mod);
    }

    for (size_t i = 0; i < n; i += len) {
      size_t w = 1;
      for (size_t j = 0; j < len / 2; j++) {
        mod_t u = (*a)[i + j];
        mod_t v = (mod_t)(1LL * (*a)[i + j + len / 2] * w % mod);
        (*a)[i + j] = u + v < (mod_t)mod ? u + v : u + v - mod;
        (*a)[i + j + len / 2] = u - v >= 0 ? u - v : u - v + mod;
        w = (mod_t)(1LL * w * wlen % mod);
      }
    }
  }

  if (invert) {
    mod_t n_inv = ModInverse<mod>(n);
    for (mod_t &x : *a)
      x = (mod_t)(1LL * x * n_inv % mod);
  }
}
} // namespace detail

template <size_t mod, mod_t root> void NTT(std::vector<mod_t> *a) {
  a->resize(NextPowerOfTwo(a->size()));
  detail::ntt<mod, root, /*inverse=*/false>(a);
}

template <size_t mod, mod_t root> void INTT(std::vector<mod_t> *a) {
  a->resize(NextPowerOfTwo(a->size()));
  detail::ntt<mod, root, /*inverse=*/true>(a);
}
