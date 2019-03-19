#pragma once
#include <vector>

struct Complex {
  double real;
  double imag;
  Complex operator*(const Complex &other) const {
    return {real * other.real - imag * other.imag,
            real * other.imag + imag * other.real};
  }
  Complex &operator*=(double m) {
    real *= m;
    imag *= m;
    return *this;
  }
  Complex &operator*=(const Complex &other) {
    double r = real;
    double i = imag;
    real = r * other.real - i * other.imag;
    imag = r * other.imag + i * other.real;
    return *this;
  }
  Complex operator+(const Complex &other) const {
    return {real + other.real, imag + other.imag};
  }
  Complex operator-(const Complex &other) const {
    return {real - other.real, imag - other.imag};
  }
};

using ComplexVector = std::vector<Complex>;

std::vector<double> Real(const ComplexVector &cplx);
std::vector<double> Imag(const ComplexVector &cplx);
ComplexVector Cplx(const std::vector<double> &real);
ComplexVector Cplx(const std::vector<double> &real,
                   const std::vector<double> &imag);

void FFT(ComplexVector *vals);
void IFFT(ComplexVector *freqs);
