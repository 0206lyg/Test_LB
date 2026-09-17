// SPDX-License-Identifier: GPL-2.0-or-later
// Cross viscosity law for OpenLB 1.9r0. All model inputs here are lattice units.
#ifndef SLURRY_CMC_CMC_CROSS_MODEL_H
#define SLURRY_CMC_CMC_CROSS_MODEL_H

#include <olb.h>

// SLURRY SCOPE BEGIN
namespace slurry { namespace cmc {

namespace cmc {
struct NU_ZERO : public olb::descriptors::FIELD_BASE<1> {
  template<typename T, typename D, typename F>
  static constexpr auto isValid(olb::FieldD<T,D,F> v) { return v > 0; }
};
struct NU_INF : public olb::descriptors::FIELD_BASE<1> {
  template<typename T, typename D, typename F>
  static constexpr auto isValid(olb::FieldD<T,D,F> v) { return v >= 0; }
};
struct LAMBDA : public olb::descriptors::FIELD_BASE<1> {
  template<typename T, typename D, typename F>
  static constexpr auto isValid(olb::FieldD<T,D,F> v) { return v > 0; }
};
struct M : public olb::descriptors::FIELD_BASE<1> {
  template<typename T, typename D, typename F>
  static constexpr auto isValid(olb::FieldD<T,D,F> v) { return v > 0; }
};

struct CrossModel {
  using parameters = olb::meta::list<NU_ZERO,NU_INF,LAMBDA,M>;
  static std::string getName() { return "CMC::Cross"; }
  template<typename T, typename PARAMETERS>
  static T computeViscosity(PARAMETERS& p, T gamma) any_platform {
    const T nu0 = p.template get<NU_ZERO>();
    const T nuInf = p.template get<NU_INF>();
    return nuInf + (nu0-nuInf) /
      (T{1} + olb::util::pow(p.template get<LAMBDA>() * gamma,
                           p.template get<M>()));
  }
};

template<typename T, typename D, typename MOMENTA=olb::momenta::BulkTuple>
using CrossBGKdynamics = olb::dynamics::Tuple<
  T,D,MOMENTA,olb::equilibria::SecondOrder,
  olb::visco::OmegaFromCell<olb::collision::BGK,CrossModel>>;
} // namespace cmc

} } // SLURRY SCOPE END
#endif
