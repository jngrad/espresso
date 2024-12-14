/*
 * Copyright (C) 2024 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_MODULE floating-point exceptions test
#define BOOST_TEST_DYN_LINK

#include <memory>
#include <mutex>
#include <optional>
#include <utility>

class fe_trap {
  struct global_state_params {
    std::weak_ptr<fe_trap> observer;
    std::mutex mutex;
  };
  static global_state_params global_state;

  struct scoped_instance {
    explicit scoped_instance(std::shared_ptr<fe_trap> ptr)
        : m_resource{std::move(ptr)} {}
    scoped_instance(scoped_instance const &) = delete;
    scoped_instance(scoped_instance &&) noexcept = default;
    scoped_instance &operator=(scoped_instance const &) = delete;
    scoped_instance &operator=(scoped_instance &&) noexcept = default;
    bool is_unique() const { return m_resource->is_unique(); }
    int get_flags() const { return m_resource->get_flags(); }

  private:
    std::shared_ptr<fe_trap> m_resource;
  };

  struct deleter {
    void operator()(fe_trap *ptr) { delete ptr; }
  };
  friend deleter;

  int m_flags;
  bool m_unique;

  fe_trap(std::optional<int> excepts, bool unique);
  ~fe_trap();

  static int parse_excepts(std::optional<int> excepts);

public:
  fe_trap(fe_trap const &) = delete;
  fe_trap(fe_trap &&) noexcept = delete;
  fe_trap &operator=(fe_trap const &) = delete;
  fe_trap &operator=(fe_trap &&) noexcept = delete;
  int get_flags() const { return m_flags; }
  bool is_unique() const { return m_unique; }
  static scoped_instance make_unique_scoped(std::optional<int> excepts = std::nullopt);
  static scoped_instance make_shared_scoped(std::optional<int> excepts = std::nullopt);
};


#include <cassert>
#include <cfenv>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <thread>

fe_trap::global_state_params fe_trap::global_state{{}, {}};

fe_trap::fe_trap(std::optional<int> excepts, bool unique) {
#if defined(__STDC_IEC_559__) and defined(__GLIBC__)
#define ESPRESSO_FENV_DEFINED_CTOR
  m_flags = parse_excepts(excepts);
  [[maybe_unused]] int status = feenableexcept(m_flags);
  assert(status == 0);
#endif // defined(__STDC_IEC_559__) and defined(__GLIBC__)
#if defined(__arm64__) and defined(__APPLE__)
#define ESPRESSO_FENV_DEFINED_CTOR
    using fpcr_t = decltype(std::fenv_t::__fpcr);
   m_flags = parse_excepts(excepts);
printf("%i -> %i\n", __fpcr_trap_invalid, FE_INVALID<<8);
printf("%i -> %i\n", __fpcr_trap_divbyzero, FE_DIVBYZERO<<8);
printf("%i -> %i\n", __fpcr_trap_overflow, FE_OVERFLOW<<8);
printf("%i -> %i\n", __fpcr_trap_underflow, FE_UNDERFLOW<<8);
  std::fenv_t env;
  [[maybe_unused]] auto const ret1 = std::fegetenv(&env);
  assert(ret1 == 0u);
  env.__fpcr |= static_cast<fpcr_t>(m_flags);
  [[maybe_unused]] auto const ret2 =  std::fesetenv(&env);
  assert(ret2 == 0u);
#endif // defined(__arm64__) and defined(__APPLE__)
  m_unique = unique;

  // sentinel if no implementation-defined solution was found
#if defined(ESPRESSO_FENV_DEFINED_CTOR)
#undef ESPRESSO_FENV_DEFINED_CTOR
#else
//#error "FE not supported"
#endif
}

fe_trap::~fe_trap() {
#if defined(__STDC_IEC_559__) and defined(__GLIBC__)
#define ESPRESSO_FENV_DEFINED_DTOR
  [[maybe_unused]] int status = fedisableexcept(m_flags);
  assert(status == 0 or status == m_flags);
#endif // defined(__STDC_IEC_559__) and defined(__GLIBC__)
#if defined(__arm64__) and defined(__APPLE__)
#define ESPRESSO_FENV_DEFINED_DTOR
    using fpcr_t = decltype(std::fenv_t::__fpcr);
    std::fenv_t env;
    [[maybe_unused]] auto const ret1 = std::fegetenv(&env);
    assert(ret1 == 0u);
printf("disable: %llu\n",env.__fpcr);
    assert((env.__fpcr & static_cast<fpcr_t>(m_flags)) == static_cast<fpcr_t>(m_flags));
    env.__fpcr &= static_cast<fpcr_t>(~m_flags);
    [[maybe_unused]] auto const ret2 =  std::fesetenv(&env);
    assert(ret2 == 0u);
#endif // defined(__arm64__) and defined(__APPLE__)

  // sentinel if no implementation-defined solution was found
#if defined(ESPRESSO_FENV_DEFINED_DTOR)
#undef ESPRESSO_FENV_DEFINED_DTOR
#else
//#error "FE not supported"
#endif
}

int fe_trap::parse_excepts(std::optional<int> excepts) {
  auto const fallback = FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW;
  int retval = excepts ? *excepts : fallback;
#if defined(__arm64__) and defined(__APPLE__)
  retval <<= 8;
#endif
  return retval;
}

fe_trap::scoped_instance fe_trap::make_unique_scoped(std::optional<int> excepts) {
  std::lock_guard<std::mutex> lock(fe_trap::global_state.mutex);
  if (fe_trap::global_state.observer.lock()) {
    throw std::runtime_error("Cannot create more than 1 instance of fe_trap");
  }
  auto raw_ptr = new fe_trap(excepts, true);
  auto watched = std::shared_ptr<fe_trap>(raw_ptr, deleter{});
  fe_trap::global_state.observer = watched;
  return fe_trap::scoped_instance(watched);
}

fe_trap::scoped_instance fe_trap::make_shared_scoped(std::optional<int> excepts) {
  std::lock_guard<std::mutex> lock(fe_trap::global_state.mutex);
  if (auto watched = fe_trap::global_state.observer.lock()) {
    if (watched->is_unique()) {
      throw std::runtime_error("Cannot create more than 1 instance of fe_trap");
    }
    if (watched->get_flags() != parse_excepts(excepts)) {
      throw std::invalid_argument("Cannot mix different exceptions with fe_trap");
    }
    return fe_trap::scoped_instance(watched);
  }
  auto raw_ptr = new fe_trap(excepts, false);
  auto watched = std::shared_ptr<fe_trap>(raw_ptr, deleter{});
  fe_trap::global_state.observer = watched;
  return fe_trap::scoped_instance(watched);
}



#include <cstdio>
#include <cfenv>

#if defined(__STDC_IEC_559__) and defined(__GLIBC__)
#define ESPRESSO_FPE_IS_SUPPORTED
#elif defined(__arm64__) and defined(__APPLE__)
#define ESPRESSO_FPE_IS_SUPPORTED
#endif

#if defined(ESPRESSO_FPE_IS_SUPPORTED) and not defined(__APPLE__)
#include <boost/test/unit_test.hpp>
#include <cassert>
#include <cmath>
#include <csetjmp>
#include <csignal>
#include <initializer_list>

static volatile std::sig_atomic_t last_signal_status = 0;
static volatile std::sig_atomic_t last_signal_code = 0;
static std::jmp_buf jmp_env;

void fpe_signal_handler(int signal) {
  ::last_signal_status = signal;
  siglongjmp(::jmp_env, 1);
}

static void fpe_signal_handler(int signum, siginfo_t *sip, void *) {
  ::last_signal_status = signum;
  ::last_signal_code = sip->si_code;
  siglongjmp(::jmp_env, 1);
}

BOOST_AUTO_TEST_CASE(trap_by_signal) {
  double volatile denominator = 0.;
  double volatile value = 1.;
  value = 0. / denominator;
  BOOST_REQUIRE(std::isnan(value));
  BOOST_REQUIRE_EQUAL(::last_signal_status, 0);
  BOOST_REQUIRE_EQUAL(::last_signal_code, 0);

#if defined(SIGFPE) and (SIGFPE == SIGILL)
  auto constexpr ref_signal_status = SIGILL;
  auto constexpr ref_signal_code = ILL_ILLTRP;
#else
  auto constexpr ref_signal_status = SIGFPE;
  auto constexpr ref_signal_code = FPE_FLTINV;
#endif

    struct sigaction act;
    act.sa_sigaction = fpe_signal_handler;
    sigemptyset(&act.sa_mask);
    act.sa_flags = SA_SIGINFO;
    sigaction(ref_signal_status, &act, nullptr);

    {
      auto const trap = fe_trap::make_unique_scoped();
      BOOST_REQUIRE(trap.is_unique());
      value = 0.;
      while (sigsetjmp(::jmp_env, 1) == 0) {
        value = 2.;
        value = 0. / denominator;
      }
      BOOST_CHECK_EQUAL(::last_signal_status, ref_signal_status);
      BOOST_CHECK_EQUAL(::last_signal_code, ref_signal_code);
      BOOST_REQUIRE(not std::isnan(value));
      BOOST_REQUIRE_EQUAL(value, 2.);
      ::last_signal_status = 0;
    }
    {
      auto const trap = fe_trap::make_shared_scoped();
      BOOST_REQUIRE(not trap.is_unique());
      value = 0.;
      while (sigsetjmp(::jmp_env, 1) == 0) {
        value = 2.;
        value = 0. / denominator;
      }
      BOOST_CHECK_EQUAL(::last_signal_status, ref_signal_status);
      BOOST_CHECK_EQUAL(::last_signal_code, ref_signal_code);
      BOOST_REQUIRE(not std::isnan(value));
      BOOST_REQUIRE_EQUAL(value, 2.);
      ::last_signal_status = 0;
    }
    {
      auto const trap1 = fe_trap::make_shared_scoped();
      {
          auto const trap2 = fe_trap::make_shared_scoped();
          BOOST_REQUIRE_EQUAL(trap1.get_flags(), trap2.get_flags());
          value = 0.;
          while (sigsetjmp(::jmp_env, 1) == 0) {
            value = 2.;
            value = 0. / denominator;
          }
          BOOST_CHECK_EQUAL(::last_signal_status, ref_signal_status);
          BOOST_CHECK_EQUAL(::last_signal_code, ref_signal_code);
          BOOST_REQUIRE(not std::isnan(value));
          BOOST_REQUIRE_EQUAL(value, 2.);
          ::last_signal_status = 0;
      }
    }
  value = 1. / denominator;
}
#endif

#if defined(ESPRESSO_FPE_IS_SUPPORTED)
BOOST_AUTO_TEST_CASE(exceptions) {
  {
    auto const trap = fe_trap::make_unique_scoped();
    BOOST_REQUIRE_THROW(fe_trap::make_unique_scoped(), std::runtime_error);
    BOOST_REQUIRE_THROW(fe_trap::make_shared_scoped(), std::runtime_error);
  }
  {
    auto const trap = fe_trap::make_shared_scoped();
    BOOST_REQUIRE_THROW(fe_trap::make_unique_scoped(), std::runtime_error);
    for (int other_excepts : {FE_INEXACT, FE_ALL_EXCEPT}) {
      if (trap.get_flags() != other_excepts) {
        BOOST_REQUIRE_THROW(fe_trap::make_shared_scoped(FE_ALL_EXCEPT),
                            std::invalid_argument);
      }
    }
  }
}
#else
#include <cstdio>
#include <cfenv>
int main() {
#if defined(__STDC_IEC_559__)
puts("defined(__STDC_IEC_559__)");
#else
puts("!defined(__STDC_IEC_559__)");
#endif
#if defined(__GLIBC__)
puts("defined(__GLIBC__)");
#else
puts("!defined(__GLIBC__)");
#endif
#if defined(__FAST_MATH__)
puts("defined(__FAST_MATH__)");
puts("!defined(__FAST_MATH__)");
#endif
#if defined(__arm64__)
puts("defined(__arm64__)");
#else
puts("!defined(__arm64__)");
#endif
#if defined(__APPLE__)
puts("defined(__APPLE__)");
#else
puts("!defined(__APPLE__)");
#endif
#if defined(HAVE_FPCR)
puts("defined(HAVE_FPCR)");
#else
puts("!defined(HAVE_FPCR)");
#endif
#if defined(HAVE_FEENABLEEXCEPT)
puts("defined(HAVE_FEENABLEEXCEPT)");
#else
puts("!defined(HAVE_FEENABLEEXCEPT)");
#endif
return 2;
}
#endif // defined(ESPRESSO_FPE_IS_SUPPORTED)
