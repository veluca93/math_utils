#include "fft.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <random>

namespace {

TEST(FFTTest, CplxReal) {
  std::mt19937 rng(123);
  std::vector<double> real(1000);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  for (double &d : real) {
    d = dist(rng);
  }
  EXPECT_EQ(Real(Cplx(real)), real);
}

TEST(FFTTest, CplxRealImag) {
  std::mt19937 rng(123);
  std::vector<double> real(1000);
  std::vector<double> imag(1000);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  for (double &d : real) {
    d = dist(rng);
  }
  for (double &d : imag) {
    d = dist(rng);
  }
  EXPECT_EQ(Real(Cplx(real, imag)), real);
  EXPECT_EQ(Imag(Cplx(real, imag)), imag);
}

TEST(FFTTest, NextPowerOfTwo) {
  for (size_t i = 0; i < 20; i++) {
    for (size_t j = (1UL << i) + 1; j <= (2UL << i); j++) {
      EXPECT_EQ(NextPowerOfTwo(j), 2UL << i);
    }
  }
}

MATCHER(DoubleNearM, "") {
  return std::fabs(std::get<0>(arg) - std::get<1>(arg)) < 1e-5;
}

TEST(FFTTest, FFTIFFTInverse) {
  std::mt19937 rng(123);
  std::vector<double> real(10003);
  std::vector<double> imag(10003);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  for (double &d : real) {
    d = dist(rng);
  }
  for (double &d : imag) {
    d = dist(rng);
  }
  auto cplx = Cplx(real, imag);
  FFT(&cplx);
  IFFT(&cplx);
  cplx.resize(real.size());
  EXPECT_THAT(Real(cplx), ::testing::Pointwise(DoubleNearM(), real));
  EXPECT_THAT(Imag(cplx), ::testing::Pointwise(DoubleNearM(), imag));
}

} // namespace
