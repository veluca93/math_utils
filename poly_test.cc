#include "poly.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <random>

namespace {

using ::testing::DoubleNear;

TEST(PolyTest, PolyVal) {
  EXPECT_THAT(PolyVal({1, -1, 2}, 5), DoubleNear(46, 1e-5));
  EXPECT_THAT(PolyVal({1, 0, 0, 0, 0, 0, 2}, 5), DoubleNear(31251, 1e-5));
}

class PolyLinearTestP : public ::testing::TestWithParam<size_t> {};
class PolyQuadraticTestP : public ::testing::TestWithParam<size_t> {};

TEST_P(PolyLinearTestP, PolyAdd) {
  std::mt19937 rng(123);
  Poly a(GetParam()), b(GetParam() / 2);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  for (double &d : a) {
    d = dist(rng);
  }
  for (double &d : b) {
    d = dist(rng);
  }
  double x = dist(rng);
  Poly sum = PolyAdd(a, b);
  EXPECT_THAT(PolyVal(a, x) + PolyVal(b, x), DoubleNear(PolyVal(sum, x), 1e-5));
}

TEST_P(PolyLinearTestP, PolySub) {
  std::mt19937 rng(123);
  Poly a(GetParam()), b(GetParam() / 2);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  for (double &d : a) {
    d = dist(rng);
  }
  for (double &d : b) {
    d = dist(rng);
  }
  double x = dist(rng);
  Poly diff = PolySub(a, b);
  EXPECT_THAT(PolyVal(a, x) - PolyVal(b, x),
              DoubleNear(PolyVal(diff, x), 1e-5));
}

TEST_P(PolyLinearTestP, PolySubSelf) {
  std::mt19937 rng(123);
  Poly a(GetParam());
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  for (double &d : a) {
    d = dist(rng);
  }
  Poly diff = PolySub(a, a);
  EXPECT_EQ(diff.size(), 1);
  EXPECT_NEAR(diff[0], 0, 1e-12);
}

TEST_P(PolyLinearTestP, PolyMulF) {
  std::mt19937 rng(123);
  Poly a(GetParam()), b(GetParam() / 2);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  for (double &d : a) {
    d = dist(rng);
  }
  for (double &d : b) {
    d = dist(rng);
  }
  double x = dist(rng);
  Poly prod = PolyMulF(a, b);
  EXPECT_THAT(std::abs(prod.back()), ::testing::Ge(1e-12));
  EXPECT_THAT(PolyVal(a, x) * PolyVal(b, x),
              DoubleNear(PolyVal(prod, x), 1e-5));
}

TEST_P(PolyQuadraticTestP, PolyMulS) {
  std::mt19937 rng(123);
  Poly a(GetParam()), b(GetParam() / 2);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  for (double &d : a) {
    d = dist(rng);
  }
  for (double &d : b) {
    d = dist(rng);
  }
  double x = dist(rng);
  Poly prod = PolyMulS(a, b);
  EXPECT_THAT(std::abs(prod.back()), ::testing::Ge(1e-12));
  EXPECT_THAT(PolyVal(a, x) * PolyVal(b, x),
              DoubleNear(PolyVal(prod, x), 1e-5));
}

TEST_P(PolyQuadraticTestP, PolyMulSvF) {
  std::mt19937 rng(123);
  size_t sz = GetParam();
  Poly a(sz), b(sz + 2);
  std::uniform_real_distribution<double> dist(1 << 14, 1 << 15);
  for (double &d : a) {
    d = dist(rng);
  }
  for (double &d : b) {
    d = dist(rng);
  }
  Poly prod_s = PolyMulS(a, b);
  Poly prod_f = PolyMulF(a, b);
  EXPECT_THAT(prod_f, testing::Pointwise(DoubleNear(5e-1), prod_s));
}

INSTANTIATE_TEST_CASE_P(UpToMillion, PolyLinearTestP,
                        ::testing::Range<size_t>(100000, 1000001, 100000));
INSTANTIATE_TEST_CASE_P(UpTo30Thousand, PolyQuadraticTestP,
                        ::testing::Range<size_t>(3000, 36001, 3000));

} // namespace
