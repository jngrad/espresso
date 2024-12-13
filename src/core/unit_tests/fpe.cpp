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
#if defined(__STDC_IEC_559__)
#if defined(__GLIBC__)
#define ESPRESSO_FENV_DEFINED_CTOR
  m_flags = parse_excepts(excepts);
  [[maybe_unused]] int status = feenableexcept(m_flags);
  assert(status == 0);
#endif // __GLIBC__
#endif // __STDC_IEC_559__
#if defined(__arm64__) and defined(__APPLE__)
#define ESPRESSO_FENV_DEFINED_CTOR
{
   m_flags = parse_excepts(excepts);
printf("%i -> %i\n", __fpcr_trap_invalid, FE_INVALID<<8);
printf("%i -> %i\n", __fpcr_trap_divbyzero, FE_DIVBYZERO<<8);
printf("%i -> %i\n", __fpcr_trap_overflow, FE_OVERFLOW<<8);
printf("%i -> %i\n", __fpcr_trap_underflow, FE_UNDERFLOW<<8);
//   m_flags = __fpcr_trap_divbyzero | __fpcr_trap_invalid | __fpcr_trap_overflow | __fpcr_trap_underflow ;
    std::fenv_t env;
    auto ret1 = std::fegetenv(&env);
assert(ret1 == 0);
printf("%i : %i  \n", FE_INVALID, __fpcr_trap_invalid);
 env.__fpcr |= m_flags;
  auto ret2 =  std::fesetenv(&env);
assert(ret2 == 0);
}
#endif
  m_unique = unique;

  // sentinel if no implementation-defined solution was found
#if defined(ESPRESSO_FENV_DEFINED_CTOR)
#undef ESPRESSO_FENV_DEFINED_CTOR
#else
#endif
}

fe_trap::~fe_trap() {
#if defined(__STDC_IEC_559__)
#if defined(__GLIBC__)
#define ESPRESSO_FENV_DEFINED_DTOR
  [[maybe_unused]] int status = fedisableexcept(m_flags);
  assert(status == 0 or status == m_flags);
#endif // __GLIBC__
#endif // __STDC_IEC_559__
#if defined(__arm64__) and defined(__APPLE__)
#define ESPRESSO_FENV_DEFINED_DTOR
{

    std::fenv_t env;
    auto ret1 = std::fegetenv(&env);
assert(ret1 == 0);
printf("disable: %i\n",env.__fpcr);
assert((env.__fpcr & m_flags) == m_flags);
 env.__fpcr &= ~m_flags;
  auto ret2 =  std::fesetenv(&env);
assert(ret2 == 0);
}
#endif

  // sentinel if no implementation-defined solution was found
#if defined(ESPRESSO_FENV_DEFINED_DTOR)
#undef ESPRESSO_FENV_DEFINED_DTOR
#else
#endif
}

int fe_trap::parse_excepts(std::optional<int> excepts) {
  int retval = 0;
#if defined(__STDC_IEC_559__)
#if defined(__GLIBC__)
#define ESPRESSO_FENV_DEFINED_CTOR
  retval = excepts ? *excepts : (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
#endif // __GLIBC__
#endif // __STDC_IEC_559__
  retval = excepts ? *excepts : (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW)<<8; // __fpcr_trap_divbyzero | __fpcr_trap_invalid | __fpcr_trap_overflow | __fpcr_trap_underflow ;
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
printf("%i != %i \n", watched->get_flags(), parse_excepts(excepts));
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

#if defined(__arm64__) and defined(__APPLE__)
//#if defined(__STDC_IEC_559__) and defined(__GLIBC__) and !defined(__FAST_MATH__)
//#include <boost/test/unit_test.hpp>
#define BOOST_CHECK_EQUAL(x, y) if (x != y) printf("%i != %i\n", x, y);
#define BOOST_REQUIRE_EQUAL(x, y) if (x != y) printf("%i != %i\n", x, y);
#define BOOST_CHECK(x) if (not x) printf("%i\n", (int)(x));
#define BOOST_REQUIRE(x) if (not x) printf("%i\n", (int)(x));
#include <cassert>
#include <cmath>
#include <csetjmp>
#include <csignal>
#include <initializer_list>

#ifdef SIGFPE
#undef SIGFPE
#define SIGFPE SIGILL
#endif

static volatile std::sig_atomic_t last_signal_status;
static std::jmp_buf jmp_env;

void sigfpe_handler(int signal) {
  ::last_signal_status = signal;
  siglongjmp(::jmp_env, 1);
}

//BOOST_AUTO_TEST_CASE(trap_by_signal) {
int main()  {
  double volatile denominator = 0.;
  double volatile value = 1.;
  printf("SIGILL: %i\n", SIGILL);
  printf("SIGFPE: %i\n", SIGFPE);
  value = 0. / denominator;
  BOOST_REQUIRE(std::isnan(value));
  BOOST_REQUIRE_EQUAL(::last_signal_status, 0);
  [[maybe_unused]] auto sig_ret = std::signal(SIGFPE, sigfpe_handler);
  assert(sig_ret != SIG_ERR);
  //std::signal(SIGILL, sigfpe_handler);
    {
      auto const trap = fe_trap::make_unique_scoped();
      BOOST_REQUIRE(trap.is_unique());
      value = 0.;
      while (sigsetjmp(::jmp_env, 1) == 0) {
        value = 2.;
        value = 0. / denominator;
      }
      BOOST_CHECK_EQUAL(::last_signal_status, SIGFPE);
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
      BOOST_CHECK_EQUAL(::last_signal_status, SIGFPE);
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
          BOOST_CHECK_EQUAL(::last_signal_status, SIGFPE);
          BOOST_REQUIRE(not std::isnan(value));
          BOOST_REQUIRE_EQUAL(value, 2.);
          ::last_signal_status = 0;
      }
    }
  std::signal(SIGFPE, SIG_DFL);
  std::signal(SIGILL, SIG_DFL);
  value = 1. / denominator;

{
    auto const trap = fe_trap::make_unique_scoped();
//    fe_trap::make_unique_scoped();
//BOOST_REQUIRE_THROW(fe_trap::make_unique_scoped(), std::runtime_error);
//fe_trap::make_shared_scoped();
  //  BOOST_REQUIRE_THROW(fe_trap::make_shared_scoped(), std::runtime_error);
  }
  {
    auto const trap = fe_trap::make_shared_scoped();
    //BOOST_REQUIRE_THROW(fe_trap::make_unique_scoped(), std::runtime_error);
    for (int other_excepts : {FE_INEXACT, FE_ALL_EXCEPT}) {
      if (trap.get_flags() != other_excepts) {
//fe_trap::make_shared_scoped(FE_ALL_EXCEPT);
       // BOOST_REQUIRE_THROW(fe_trap::make_shared_scoped(FE_ALL_EXCEPT),
         //                   std::invalid_argument);
      }
    }
  }

puts("out");
}
/*
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
*/
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
#else
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
#endif
