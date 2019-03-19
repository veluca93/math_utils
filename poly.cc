#include "poly.h"
#include "fft.h"

void PolyPrint(const Poly &a) {
  assert(a.size() > 0);
  fprintf(stderr, "%f", a[0]);
  for (size_t i = 1; i < a.size(); i++) {
    fprintf(stderr, "%+fx**%lu", a[i], i);
  }
  fprintf(stderr, "\n");
}

double PolyVal(const Poly &a, double x) {
  double ret = a.back();
  for (size_t i = a.size() - 1; i > 0; i--) {
    ret *= x;
    ret += a[i - 1];
  }
  return ret;
}

Poly PolyAdd(const Poly &a, const Poly &b) {
  Poly ret(std::max(a.size(), b.size()));
  for (size_t i = 0; i < a.size(); i++) {
    ret[i] += a[i];
  }
  for (size_t i = 0; i < b.size(); i++) {
    ret[i] += b[i];
  }
  while (ret.size() > 1 && std::abs(ret.back()) < 1e-12) {
    ret.pop_back();
  }
  return ret;
}

Poly PolySub(const Poly &a, const Poly &b) {
  Poly ret(std::max(a.size(), b.size()));
  for (size_t i = 0; i < a.size(); i++) {
    ret[i] += a[i];
  }
  for (size_t i = 0; i < b.size(); i++) {
    ret[i] -= b[i];
  }
  while (ret.size() > 1 && std::abs(ret.back()) < 1e-12) {
    ret.pop_back();
  }
  return ret;
}

Poly PolyMul(Poly a, double v) {
  for (double &c : a) {
    c *= v;
  }
  return a;
}

void PolyNeg(Poly *a) {
  for (double &v : *a) {
    v = -v;
  }
}

Poly PolyMulS(const Poly &a, const Poly &b) {
  Poly ret(a.size() + b.size() - 1);
  for (size_t i = 0; i < a.size(); i++) {
    for (size_t j = 0; j < b.size(); j++) {
      ret[i + j] += a[i] * b[j];
    }
  }
  return ret;
}

Poly PolyMulF(const Poly &a, const Poly &b) {
  auto ap = Cplx(a);
  ap.resize(a.size() + b.size() - 1);
  FFT(&ap);
  auto bp = Cplx(b);
  bp.resize(ap.size());
  FFT(&bp);
  for (size_t i = 0; i < ap.size(); i++) {
    ap[i] *= bp[i];
  }
  IFFT(&ap);
  Poly ret = Real(ap);
  ret.resize(a.size() + b.size() - 1);
  return ret;
}
