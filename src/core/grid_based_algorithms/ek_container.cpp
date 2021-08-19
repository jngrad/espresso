#include "ek_container.hpp"

namespace EK {
EKContainer<EKinWalberlaBase<double>> ek_container;

double get_tau() { return ek_container.get_tau(); }

void propagate() {
  // first calculate the charge for the potential, for that get all the
  // field-ids from the ekspecies pass the potential-field-id to the
  // flux-kernels of the eks for this the integrate function has to be split
  // with a public interface to diffusive and advective-flux this should also
  // allow the back-coupling to the LB with a field-id

  if (ek_container.empty()) {
    return;
  }

  ek_container.reset_charge();
  std::for_each(ek_container.begin(), ek_container.end(), [](auto const &ek) {
    ek_container.add_charge(ek->get_density_id(), ek->get_valency());
  });

  ek_container.solve_poisson();

  std::for_each(ek_container.begin(), ek_container.end(),
                [](auto const &ek) { ek->integrate(); });
}
} // namespace EK