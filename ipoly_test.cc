#include "ipoly.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <random>

namespace {

constexpr size_t mod = 998244353;
using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::Lt;
using ::testing::Ne;

TEST(IPolyTest, PolyVal) {
  EXPECT_THAT(PolyVal<mod>({1, -1, 2}, 5), Eq(46));
  EXPECT_THAT(PolyVal<mod>({1, 0, 0, 0, 0, 0, 2}, 5), Eq(31251));
}

class PolyLinearTestP : public ::testing::TestWithParam<size_t> {};
class PolyQuadraticTestP : public ::testing::TestWithParam<size_t> {};

TEST_P(PolyLinearTestP, PolyAdd) {
  std::mt19937 rng(123);
  IPoly a(GetParam()), b(GetParam() / 2);
  std::uniform_int_distribution<mod_t> dist(0, mod);
  for (mod_t &d : a) {
    d = dist(rng);
  }
  for (mod_t &d : b) {
    d = dist(rng);
  }
  double x = dist(rng);
  IPoly sum = PolyAdd<mod>(a, b);
  EXPECT_THAT((PolyVal<mod>(a, x) + PolyVal<mod>(b, x)) % mod,
              Eq(PolyVal<mod>(sum, x)));
}

TEST_P(PolyLinearTestP, PolySub) {
  std::mt19937 rng(123);
  IPoly a(GetParam()), b(GetParam() / 2);
  std::uniform_int_distribution<mod_t> dist(0, mod);
  for (mod_t &d : a) {
    d = dist(rng);
  }
  for (mod_t &d : b) {
    d = dist(rng);
  }
  double x = dist(rng);
  IPoly diff = PolySub<mod>(a, b);
  EXPECT_THAT((PolyVal<mod>(a, x) - PolyVal<mod>(b, x) + mod) % mod,
              Eq(PolyVal<mod>(diff, x)));
}

TEST_P(PolyLinearTestP, PolySubSelf) {
  std::mt19937 rng(123);
  IPoly a(GetParam());
  std::uniform_int_distribution<mod_t> dist(0, mod);
  for (mod_t &d : a) {
    d = dist(rng);
  }
  IPoly diff = PolySub<mod>(a, a);
  EXPECT_THAT(diff, ElementsAre(0));
}

TEST_P(PolyLinearTestP, PolyMulF) {
  std::mt19937 rng(123);
  IPoly a(GetParam()), b(GetParam() / 2);
  std::uniform_int_distribution<mod_t> dist(0, mod);
  for (mod_t &d : a) {
    d = dist(rng);
  }
  for (mod_t &d : b) {
    d = dist(rng);
  }
  double x = dist(rng);
  IPoly prod = PolyMulF<mod>(a, b);
  EXPECT_THAT(std::abs(prod.back()), Ne(0));
  EXPECT_THAT((PolyVal<mod>(a, x) * PolyVal<mod>(b, x)) % mod,
              Eq(PolyVal<mod>(prod, x)));
}

TEST_P(PolyQuadraticTestP, PolyMulS) {
  std::mt19937 rng(123);
  IPoly a(GetParam()), b(GetParam() / 2);
  std::uniform_int_distribution<mod_t> dist(0, mod);
  for (mod_t &d : a) {
    d = dist(rng);
  }
  for (mod_t &d : b) {
    d = dist(rng);
  }
  if (a.back() == 0)
    a.back() = 1;
  if (b.back() == 0)
    b.back() = 1;
  double x = dist(rng);
  IPoly prod = PolyMulS<mod>(a, b);
  EXPECT_THAT(prod.back(), Ne(0));
  EXPECT_THAT((PolyVal<mod>(a, x) * PolyVal<mod>(b, x)) % mod,
              Eq(PolyVal<mod>(prod, x)));
}

TEST_P(PolyQuadraticTestP, PolyMulSvF) {
  std::mt19937 rng(123);
  IPoly a(GetParam()), b(GetParam() / 2);
  std::uniform_int_distribution<mod_t> dist(0, mod);
  for (mod_t &d : a) {
    d = dist(rng);
  }
  for (mod_t &d : b) {
    d = dist(rng);
  }
  IPoly prod_s = PolyMulS<mod>(a, b);
  IPoly prod_f = PolyMulF<mod>(a, b);
  EXPECT_THAT(prod_f, Eq(prod_s));
}

TEST_P(PolyLinearTestP, PolyRecp) {
  size_t sz = GetParam();
  IPoly a{1, mod - 2, 1};
  size_t modp = sz;
  auto recp = PolyRecp<mod>(a, modp);
  auto recp_times_a = PolyMul<mod>(recp, a);
  recp_times_a.resize(modp);
  IPoly expected(modp);
  expected[0] = 1;
  EXPECT_THAT(recp_times_a, Eq(expected));
}

TEST_P(PolyLinearTestP, PolyMod) {
  std::mt19937 rng(123);
  IPoly a(GetParam()), b(GetParam() * 3 / 4);
  std::uniform_int_distribution<mod_t> dist(0, mod);
  for (mod_t &d : a) {
    d = dist(rng);
  }
  for (mod_t &d : b) {
    d = dist(rng);
  }
  IPoly md = PolyMod<mod>(a, b);
  EXPECT_THAT(md.size(), Lt(b.size()));
}

INSTANTIATE_TEST_CASE_P(UpToMillion, PolyLinearTestP,
                        ::testing::Range<size_t>(100000, 1000001, 100000));
INSTANTIATE_TEST_CASE_P(UpTo30Thousand, PolyQuadraticTestP,
                        ::testing::Range<size_t>(3000, 36001, 3000));

} // namespace
