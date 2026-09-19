// Standalone conservative local-adhesion checks. No OpenLB or PETSc required.
/* Run from the repository root:
   c++ -std=c++17 -O2 -Wall -Wextra -pedantic \
     -I olb-1.9r0/src/slurry/gr_re2 \
     tests/graphite_adhesion/adhesion_tests.cpp -o /tmp/adhesion_tests
   /tmp/adhesion_tests
*/
#include "reSquaredPotential.h"
#include "roughContact.h"

#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace g = slurry::gr_re2::graphite;
namespace {
constexpr double pi = 3.1415926535897932384626433832795;

void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double absolute, double relative,
          const std::string& message) {
  if (!std::isfinite(actual) ||
      std::abs(actual - expected) > absolute + relative * std::abs(expected)) {
    std::ostringstream out;
    out << std::setprecision(17) << message << ": actual=" << actual
        << ", expected=" << expected;
    throw std::runtime_error(out.str());
  }
}

void nearVector(const g::Vec3& actual, const g::Vec3& expected,
                double absolute, double relative, const std::string& message) {
  for (int axis = 0; axis < 3; ++axis)
    near(actual[axis], expected[axis], absolute, relative,
         message + " component " + std::to_string(axis));
}

void domainError(const std::function<void()>& operation, const std::string& message) {
  try { operation(); }
  catch (const std::domain_error&) { return; }
  throw std::runtime_error(message);
}

g::PairParameters strongParameters() {
  g::PairParameters p;
  p.hamaker = 9.9e-20;
  p.sigma = 4.197e-10;
  p.roughnessGap = 2.e-9;
  p.localGap = 3.e-10;
  p.localGapFraction = 1.;
  p.localSwitchExcessGap = 2.e-9;
  p.localCutoffExcessGap = 10.e-9;
  return p;
}

g::Vec3 support(const g::Body& body, const g::Vec3& normal) {
  g::Vec3 squared{};
  for (int k = 0; k < 3; ++k) squared[k] = body.axes[k] * body.axes[k];
  const auto qn = g::mul(g::rotatedDiagonal(body.rotation, squared), normal);
  return g::scale(qn, 1. / std::sqrt(g::dot(normal, qn)));
}

using Bodies = std::pair<g::Body, g::Body>;

Bodies atGap(double gap, const g::Mat3& ri = g::identity3,
             const g::Mat3& rj = g::identity3,
             const g::Vec3& normal = {0., 0., 1.}) {
  Bodies pair;
  pair.first.rotation = ri;
  pair.second.rotation = rj;
  pair.first.position = {2.e-6, -3.e-6, 1.e-6};
  const auto n = g::normalized(normal);
  pair.second.position = g::add(pair.first.position,
      g::add(g::add(support(pair.first, n), support(pair.second, n)),
             g::scale(n, gap)));
  return pair;
}

Bodies genericPair(double gap) {
  return atGap(gap, g::rotationIncrement({.19, -.28, .08}),
               g::rotationIncrement({-.11, .21, -.15}), {.18, -.12, 1.});
}

g::PairResult evaluate(const Bodies& pair, const g::PairParameters& p) {
  return g::evaluatePair(pair.first, pair.second, p);
}

void legacyAndFarField() {
  // Golden values captured from the pre-local-adhesion production evaluator,
  // for two tilted graphite ellipsoids. These check orientation derivatives as
  // well as the old force scale; there is no second implementation in the test.
  g::PairParameters legacy;
  require(legacy.localGapFraction == 0., "local adhesion must be opt-in outside pure_gr");
  const auto pair = genericPair(2.e-9);
  const auto actual = evaluate(pair, legacy);
  near(actual.ua, -1.0084121104705865e-18, 1.e-29, 2.e-10, "legacy attraction energy");
  near(actual.ur, 5.5434855511554439e-20, 1.e-30, 2.e-10, "legacy repulsion energy");
  nearVector(actual.forceI,
      {5.5735062538339617e-11, -3.739650819127399e-11, 3.0693649415272606e-10},
      1.e-21, 2.e-10, "legacy force");
  nearVector(actual.torqueI,
      {8.2303332401790989e-17, -4.3066191738222818e-16, -6.6782724534318915e-17},
      1.e-27, 2.e-10, "legacy torque I");
  nearVector(actual.torqueJ,
      {-4.1832054712562659e-16, 7.9977988044265892e-17, 8.507180435212385e-17},
      1.e-27, 2.e-10, "legacy torque J");

  auto local = strongParameters();
  auto off = local;
  off.localGapFraction = 0.;
  for (double gap : {13.e-9, 50.e-9, 450.e-9, 501.e-9}) {
    const auto bodies = genericPair(gap);
    const auto a = evaluate(bodies, local), b = evaluate(bodies, off);
    require(a.energy == b.energy && a.forceI == b.forceI &&
            a.torqueI == b.torqueI && a.torqueJ == b.torqueJ,
            "local correction must leave the far field exactly unchanged");
  }
}

void contactForceScale() {
  const auto p = strongParameters();
  const auto edge = g::rotationIncrement({0., .5 * pi, 0.});
  const std::vector<Bodies> pairs{atGap(p.roughnessGap),
      atGap(p.roughnessGap, edge), atGap(p.roughnessGap, edge, edge)};
  const double force[] = {936.214833e-9, 25.950663e-9, 13.748022e-9};
  for (std::size_t i = 0; i < pairs.size(); ++i) {
    const auto result = evaluate(pairs[i], p);
    near(g::dot(result.forceI, result.normal), force[i], 0., 8.e-5,
         "strong contact force scale " + std::to_string(i));
    nearVector(g::add(result.forceAttractiveI, result.forceRepulsiveI),
               result.forceI, 0., 0., "branch force sum");
    near(result.ua + result.ur, result.energy, 0., 0., "branch energy sum");
  }
  // Independent near-contact continuum check for two aligned oblate faces:
  // R_eff=a^2/(2c), w(D)=A/(12 pi D^2)-A sigma^6/(360 pi D^8).
  // The finite-size RE2 result should approach 2 pi R_eff w(D), including the
  // repulsive branch. This tests physical normalization, not copied AD code.
  const auto& body = pairs[0].first;
  const double radius = body.axes[0] * body.axes[0] / (2. * body.axes[2]);
  const double d = p.localGap;
  const double work = p.hamaker / (12. * pi * d * d) *
      (1. - std::pow(p.sigma / d, 6) / 30.);
  near(g::dot(evaluate(pairs[0], p).forceI, {0., 0., 1.}),
       2. * pi * radius * work, 0., .002, "Derjaguin force normalization");
}

void gradientsAndAngularMomentum() {
  auto p = strongParameters();
  // Exercise full replacement, blended replacement, and the switched region.
  for (double fraction : {1., .37}) {
    p.localGapFraction = fraction;
    for (double gap : {2.e-9, 4.7e-9, 9.e-9}) {
      const auto bodies = genericPair(gap);
      const auto result = evaluate(bodies, p);
      constexpr double dx = 2.e-14, angle = 2.e-8;
      for (int axis = 0; axis < 3; ++axis) {
        auto plus = bodies, minus = bodies;
        plus.second.position[axis] += dx;
        minus.second.position[axis] -= dx;
        const double derivative = (evaluate(plus, p).energy - evaluate(minus, p).energy) / (2. * dx);
        near(result.forceI[axis], derivative, 5.e-15, 3.e-5,
             "force is the energy gradient");
        for (int body = 0; body < 2; ++body) {
          plus = bodies;
          minus = bodies;
          g::Vec3 rotation{};
          rotation[axis] = angle;
          auto& bp = body == 0 ? plus.first : plus.second;
          auto& bm = body == 0 ? minus.first : minus.second;
          bp.rotation = g::rotateLaboratory(bp.rotation, rotation);
          bm.rotation = g::rotateLaboratory(bm.rotation, g::scale(rotation, -1.));
          const double torque = -(evaluate(plus, p).energy - evaluate(minus, p).energy) / (2. * angle);
          near(body == 0 ? result.torqueI[axis] : result.torqueJ[axis], torque,
               5.e-21, 3.e-5, "torque is the rotation energy gradient");
        }
      }
      const auto moment = g::sub(g::add(result.torqueI, result.torqueJ),
          g::cross(g::sub(bodies.second.position, bodies.first.position), result.forceI));
      nearVector(moment, {}, 2.e-25 + 2.e-11 * g::norm(result.torqueI), 0.,
                 "internal pair conserves total angular momentum");
    }
  }
}

void releaseAndSwitches() {
  const auto p = strongParameters();
  // Follow separation at fixed orientation through both local switch ends.
  // A repulsive barrier introduced by the switch would block physical release.
  for (bool tilted : {false, true}) {
    double previous = -std::numeric_limits<double>::infinity();
    for (int k = 0; k <= 240; ++k) {
      const double gap = p.roughnessGap + k * .05e-9;
      const auto result = evaluate(tilted ? genericPair(gap) : atGap(gap), p);
      require(g::dot(result.forceI, result.normal) > 0.,
              "separation must remain attractive through the local switch");
      require(result.energy > previous, "release energy must rise monotonically");
      previous = result.energy;
    }
  }
  for (double edge : {p.roughnessGap + p.localSwitchExcessGap,
                      p.roughnessGap + p.localCutoffExcessGap}) {
    constexpr double delta = 1.e-15;
    const auto before = evaluate(genericPair(edge - delta), p);
    const auto after = evaluate(genericPair(edge + delta), p);
    near(before.energy, after.energy, 1.e-27, 5.e-6, "local switch energy continuity");
    nearVector(before.forceI, after.forceI, 1.e-18, 5.e-6, "local switch force continuity");
    nearVector(before.torqueI, after.torqueI, 1.e-24, 5.e-6, "local switch torque continuity");
  }
}

void parameterAndTrialDomains() {
  const auto good = strongParameters();
  const double floor = good.roughnessGap - good.localGap;
  near(g::minimumPairGap(good), floor, 0., 0., "local gap trial-domain floor");
  domainError([&] { evaluate(atGap(floor - 1.e-12), good); },
              "a trial with nonpositive local distance must be rejected");
  require(std::isfinite(evaluate(atGap(floor + 1.e-12), good).energy),
          "positive local trial distance must remain evaluable");
  auto off = good;
  off.localGapFraction = 0.;
  near(g::minimumPairGap(off), 0., 0., 0., "legacy trial-domain floor");
  require(std::isfinite(evaluate(atGap(1.e-9), off).energy),
          "disabled local adhesion must not restrict positive geometric gaps");
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const std::vector<std::function<void(g::PairParameters&)>> corruptions{
      [](auto& p) { p.localGapFraction = -1.; },
      [](auto& p) { p.localGapFraction = 1.01; },
      [=](auto& p) { p.localGapFraction = nan; },
      [](auto& p) { p.localGap = 0.; },
      [](auto& p) { p.localGap = 3.e-9; },
      [](auto& p) { p.localSwitchExcessGap = -1.e-9; },
      [](auto& p) { p.localCutoffExcessGap = p.localSwitchExcessGap; },
      [](auto& p) { p.switchGap = 5.e-9; },
      [=](auto& p) { p.hamaker = nan; },
      [=](auto& p) { p.sigma = std::numeric_limits<double>::infinity(); }};
  for (const auto& corrupt : corruptions) {
    auto p = good;
    corrupt(p);
    domainError([&] { g::validatePairParameters(p); }, "invalid parameters must fail validation");
  }
}

void contactCapacity() {
  const auto p = strongParameters();
  auto bodies = atGap(p.roughnessGap);
  const auto interaction = evaluate(bodies, p);
  const double adhesion = g::dot(interaction.forceI, interaction.normal);
  g::RoughContactSettings settings;
  settings.enabled = true;
  settings.friction = 1.;
  settings.tangentialStiffness = 80.;
  constexpr double dt = 1.e-3;
  bodies.first.velocity = {10.e-9 / dt, 0., 0.};
  const auto contact = [&](const g::RoughContactSettings& s) {
    return g::roughContact(bodies.first, bodies.second, interaction.normal,
        interaction.leverI, interaction.leverJ, adhesion, adhesion, dt, {}, s);
  };
  const auto stick = contact(settings);
  require(stick.active && !stick.sliding, "10 nm displacement stays below the strong sliding cap");
  near(g::norm(stick.tangentForce), 800.e-9, 1.e-18, 1.e-12, "80 N/m mobilizes 800 nN at 10 nm");
  near(stick.candidateState.rollingCap, settings.rollingLength * adhesion,
       1.e-27, 1.e-12, "rolling birth cap uses the corrected total adhesion");
  auto soft = settings;
  soft.tangentialStiffness = 9.;
  near(g::norm(contact(soft).tangentForce), 90.e-9, 1.e-18, 1.e-12,
       "legacy stiffness mobilizes 90 nN at 10 nm");
  bodies.first.velocity = {20.e-9 / dt, 0., 0.};
  const auto slide = contact(settings);
  require(slide.active && slide.sliding, "sliding yield does not detach the contact");
  near(g::norm(slide.tangentForce), settings.friction * adhesion,
       1.e-18, 1.e-12, "sliding force stops at mu N");
  require(slide.plasticSlipWork > 0., "sliding dissipates positive work");
  bodies.first.velocity = {};
  bodies.first.omega = {2. * settings.rollingYieldAngle / dt, 0., 0.};
  const auto rolling = contact(settings);
  require(rolling.rolling, "rotation beyond the yield angle activates plastic rolling");
  near(g::norm(rolling.rollingTorque), settings.rollingLength * adhesion,
       1.e-27, 1.e-12, "rolling torque stops at the corrected birth cap");
  require(rolling.plasticRollWork > 0., "rolling dissipates positive work");
}
} // namespace

int main() {
  const std::vector<std::pair<std::string, std::function<void()>>> checks{
      {"legacy and far-field preservation", legacyAndFarField},
      {"contact force and continuum scale", contactForceScale},
      {"energy gradients and angular momentum", gradientsAndAngularMomentum},
      {"release and local switch continuity", releaseAndSwitches},
      {"parameter validation and trial domains", parameterAndTrialDomains},
      {"sliding and rolling capacity", contactCapacity}};
  int failed = 0;
  for (const auto& check : checks) {
    try { check.second(); std::cout << "PASS " << check.first << '\n'; }
    catch (const std::exception& error) {
      ++failed;
      std::cerr << "FAIL " << check.first << ": " << error.what() << '\n';
    }
  }
  return failed ? 1 : 0;
}
