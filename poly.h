#pragma once
#include "util.h"
#include <assert.h>
#include <complex>
#include <vector>

using Poly = std::vector<double>;

Poly PolyMulS(const Poly &a, const Poly &b);
Poly PolyMulF(const Poly &a, const Poly &b);

double PolyVal(const Poly &a, double x);
Poly PolyAdd(const Poly &a, const Poly &b);
Poly PolySub(const Poly &a, const Poly &b);
Poly PolyMul(Poly a, double v);
void PolyNeg(Poly *a);
void PolyPrint(const Poly &a);

inline Poly PolyMul(const Poly &a, const Poly &b) {
  // TODO: possibly tune this constant.
  const constexpr size_t kFftMulThres = 32;
  if (a.size() < kFftMulThres || b.size() < kFftMulThres) {
    return PolyMulS(a, b);
  } else {
    return PolyMulF(a, b);
  }
}
