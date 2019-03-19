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
  return 0;
}

template <size_t mod> constexpr size_t BestOrder() {
  return 1ULL << __builtin_ctzll(mod - 1);
}

template <size_t mod> constexpr size_t FindRoot() {
  size_t best = BestOrder<mod>();
  for (size_t root = 1; root < mod; root++) {
    if (ModPow2Order<mod>(root) == best) {
      return root;
    }
  }
  CONSTEXPR_FAIL("Could not find any root!");
}

template <size_t mod, bool invert> void ntt(std::vector<mod_t> *a) {
  size_t n = a->size();
  assert(__builtin_popcountll(n) == 1);

  constexpr size_t root = FindRoot<mod>();
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

template <size_t mod> void NTT(std::vector<mod_t> *a) {
  a->resize(NextPowerOfTwo(a->size()));
  detail::ntt<mod, /*inverse=*/false>(a);
}

template <size_t mod> void INTT(std::vector<mod_t> *a) {
  a->resize(NextPowerOfTwo(a->size()));
  detail::ntt<mod, /*inverse=*/true>(a);
}
