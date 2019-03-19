#include "ntt.h"
#include <algorithm>

using IPoly = std::vector<mod_t>;

template <size_t mod> IPoly PolyMulS(const IPoly &a, const IPoly &b) {
  IPoly ret(a.size() + b.size() - 1);
  for (size_t i = 0; i < a.size(); i++) {
    for (size_t j = 0; j < b.size(); j++) {
      ret[i + j] += a[i] * b[j] % mod;
    }
  }
  for (size_t i = 0; i < ret.size(); i++) {
    ret[i] %= mod;
  }
  return ret;
}

template <size_t mod> IPoly PolyMulF(IPoly a, IPoly b) {
  a.resize(a.size() + b.size() - 1);
  size_t res_size = a.size();
  NTT<mod>(&a);
  b.resize(a.size());
  NTT<mod>(&b);
  for (size_t i = 0; i < a.size(); i++) {
    a[i] = a[i] * b[i] % mod;
  }
  INTT<mod>(&a);
  a.resize(res_size);
  return a;
}

template <size_t mod> mod_t PolyVal(const IPoly &a, mod_t x) {
  mod_t ret = a.back();
  for (size_t i = a.size() - 1; i > 0; i--) {
    ret = (x * ret % mod + a[i - 1] + mod) % mod;
  }
  return ret;
}

template <size_t mod> IPoly PolyAdd(const IPoly &a, const IPoly &b) {
  IPoly ret(std::max(a.size(), b.size()));
  for (size_t i = 0; i < a.size(); i++) {
    ret[i] += a[i];
  }
  for (size_t i = 0; i < b.size(); i++) {
    ret[i] += b[i];
  }
  for (mod_t &v : ret) {
    v %= mod;
  }
  while (ret.size() > 1 && ret.back() == 0) {
    ret.pop_back();
  }
  return ret;
}

template <size_t mod> IPoly PolySub(const IPoly &a, const IPoly &b) {
  IPoly ret(std::max(a.size(), b.size()));
  for (size_t i = 0; i < a.size(); i++) {
    ret[i] += a[i];
  }
  for (size_t i = 0; i < b.size(); i++) {
    ret[i] -= b[i];
  }
  for (mod_t &v : ret) {
    v = (v + mod) % mod;
  }
  while (ret.size() > 1 && ret.back() == 0) {
    ret.pop_back();
  }
  return ret;
}

template <size_t mod> IPoly PolyMul(IPoly a, IPoly b) {
  // TODO: tune constant
  const constexpr size_t kNttMulThres = 32;
  if (a.size() < kNttMulThres || b.size() < kNttMulThres) {
    return PolyMulS<mod>(std::move(a), std::move(b));
  } else {
    return PolyMulF<mod>(std::move(a), std::move(b));
  }
}

template <size_t mod> void PolyNeg(IPoly *a) {
  for (mod_t &v : *a) {
    v = mod - v;
  }
}

IPoly PolyShift(const IPoly &a, long long shift) {
  long long ns = a.size() + shift;
  if (ns <= 0)
    return {};
  IPoly answer(ns);
  for (size_t i = 0; i < a.size(); i++) {
    long long pos = i + shift;
    if (pos < 0)
      continue;
    answer[pos] = a[i];
  }
  return answer;
}

// Computes the inverse of a modulo x^d.
template <size_t mod> IPoly PolyRecp(const IPoly &a, size_t d) {
  assert(a[0] != 0);
  IPoly inverse{ModInverse<mod>(a[0])};
  for (size_t inv_modp = 1; inv_modp < d; inv_modp *= 2) {
    IPoly inva = PolyMul<mod>(inverse, a);
    if (inva.size() > 2 * inv_modp) {
      inva.resize(2 * inv_modp);
    }
    PolyNeg<mod>(&inva);
    inva[0] += 2;
    inverse = PolyMul<mod>(inva, inverse);
    inverse.resize(2 * inv_modp);
  }
  inverse.resize(d);
  return inverse;
}

template <size_t mod> IPoly PolyIdiv(IPoly a, IPoly b) {
  if (a.size() < b.size()) {
    return a;
  }
  size_t rdegree = a.size() - b.size() + 1;
  std::reverse(b.begin(), b.end());
  b = PolyRecp<mod>(b, rdegree);
  std::reverse(a.begin(), a.end());
  b = PolyMul<mod>(a, b);
  b.resize(rdegree);
  std::reverse(b.begin(), b.end());
  return b;
}

template <size_t mod> IPoly PolyMod(IPoly a, IPoly b) {
  return PolySub<mod>(a, PolyMul<mod>(b, PolyIdiv<mod>(a, b)));
}
