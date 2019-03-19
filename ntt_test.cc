#include "ntt.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <random>

constexpr size_t mod = 998244353;
constexpr mod_t root = 31;

namespace {

class NTTTestP : public ::testing::TestWithParam<size_t> {};

TEST_P(NTTTestP, NTTINTTInverse) {
  std::mt19937 rng(123);
  std::vector<mod_t> v(GetParam());
  std::uniform_int_distribution<mod_t> dist(0, mod);
  for (mod_t &d : v) {
    d = dist(rng);
  }
  auto v_copy = v;
  NTT<mod, root>(&v_copy);
  INTT<mod, root>(&v_copy);
  v_copy.resize(v.size());
  EXPECT_THAT(v_copy, ::testing::Eq(v));
}

INSTANTIATE_TEST_CASE_P(UpToMillion, NTTTestP,
                        ::testing::Range<size_t>(100000, 1000001, 100000));
} // namespace
