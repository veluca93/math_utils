#pragma once
#include <assert.h>
#include <complex>
#include <vector>

using ComplexVector = std::vector<std::complex<double>>;

std::vector<double> Real(const ComplexVector &cplx);
std::vector<double> Imag(const ComplexVector &cplx);
ComplexVector Cplx(const std::vector<double> &real);
ComplexVector Cplx(const std::vector<double> &real,
                   const std::vector<double> &imag);

void FFT(ComplexVector *vals);
void IFFT(ComplexVector *freqs);

inline size_t NextPowerOfTwo(size_t n) {
  assert(n != 0);
  return n == 1 ? 1 : 1ULL << (64 - __builtin_clzll(n - 1));
}
