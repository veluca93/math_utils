#include "fft.h"
#include "util.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <random>

namespace {

using ::testing::DoubleNear;

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

class FFTTestP : public ::testing::TestWithParam<size_t> {};

TEST_P(FFTTestP, FFTIFFTInverse) {
  std::mt19937 rng(123);
  std::vector<double> real(GetParam());
  std::vector<double> imag(GetParam());
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
  EXPECT_THAT(Real(cplx), ::testing::Pointwise(DoubleNear(1e-5), real));
  EXPECT_THAT(Imag(cplx), ::testing::Pointwise(DoubleNear(1e-5), imag));
}

TEST_P(FFTTestP, FFTIFFTInverseBigCoefficients) {
  std::mt19937 rng(123);
  std::vector<double> real(GetParam());
  std::vector<double> imag(GetParam());
  std::iota(real.begin(), real.end(), 1.0);
  std::iota(imag.begin(), imag.end(), 1.0);
  auto cplx = Cplx(real, imag);
  FFT(&cplx);
  IFFT(&cplx);
  cplx.resize(real.size());
  EXPECT_THAT(Real(cplx), ::testing::Pointwise(DoubleNear(1e-4), real));
  EXPECT_THAT(Imag(cplx), ::testing::Pointwise(DoubleNear(1e-4), imag));
}

INSTANTIATE_TEST_CASE_P(UpToMillion, FFTTestP,
                        ::testing::Range<size_t>(100000, 1000001, 100000));

} // namespace
