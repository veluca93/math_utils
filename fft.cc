#include "fft.h"
#include "util.h"
#include <math.h>
#include <stdio.h>

std::vector<double> Real(const ComplexVector &cplx) {
  std::vector<double> ret(cplx.size());
  for (size_t i = 0; i < cplx.size(); i++) {
    ret[i] = cplx[i].real;
  }
  return ret;
}

std::vector<double> Imag(const ComplexVector &cplx) {
  std::vector<double> ret(cplx.size());
  for (size_t i = 0; i < cplx.size(); i++) {
    ret[i] = cplx[i].imag;
  }
  return ret;
}

ComplexVector Cplx(const std::vector<double> &real) {
  ComplexVector ret(real.size());
  for (size_t i = 0; i < real.size(); i++) {
    ret[i] = Complex{real[i], 0};
  }
  return ret;
}

ComplexVector Cplx(const std::vector<double> &real,
                   const std::vector<double> &imag) {
  assert(real.size() == imag.size());
  ComplexVector ret(real.size());
  for (size_t i = 0; i < real.size(); i++) {
    ret[i] = Complex{real[i], imag[i]};
  }
  return ret;
}

namespace {
// Taken from https://cp-algorithms.com/algebra/fft.html
template <bool invert> void FFTImpl(ComplexVector *a) {
  const size_t n = a->size();

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

  std::vector<Complex> ws(n + 1);
  for (size_t len = 2; len <= n; len <<= 1) {
    double ang = 2 * M_PI / len * (invert ? -1 : 1);
    for (size_t i = 0; i < len / 2; i++) {
      ws[i].real = std::cos(i * ang);
      ws[i].imag = std::sin(i * ang);
    }
    for (size_t i = 0; i < n; i += len) {
      for (size_t j = 0; j < len / 2; j++) {
        Complex w = ws[j];
        auto u = (*a)[i + j];
        auto v = (*a)[i + j + len / 2] * w;
        (*a)[i + j] = u + v;
        (*a)[i + j + len / 2] = u - v;
      }
    }
  }

  if (invert) {
    double inv = 1.0 / n;
    for (auto &x : *a) {
      x *= inv;
    }
  }
}
} // namespace

void FFT(ComplexVector *vals) {
  vals->resize(NextPowerOfTwo(vals->size()));
  FFTImpl</*invert=*/false>(vals);
}
void IFFT(ComplexVector *freqs) {
  assert(__builtin_popcountll(freqs->size()) == 1);
  FFTImpl</*invert=*/true>(freqs);
}
