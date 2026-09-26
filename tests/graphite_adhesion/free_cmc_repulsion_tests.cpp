// Standalone conservative free-CMC interaction checks.
/* c++ -std=c++17 -O2 -Wall -Wextra -pedantic \
     -I olb-1.9r0/src/slurry/gr_re2 \
     tests/graphite_adhesion/free_cmc_repulsion_tests.cpp -o /tmp/free_cmc_repulsion_tests
   /tmp/free_cmc_repulsion_tests */
#include "reSquaredPotential.h"

#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace g = slurry::gr_re2::graphite;
namespace {
constexpr double pi = 3.1415926535897932384626433832795;
using Bodies = std::pair<g::Body, g::Body>;

void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}
void near(double actual, double expected, double absolute, double relative,
          const std::string& message) {
  if (!std::isfinite(actual) || !std::isfinite(expected) ||
      std::abs(actual - expected) > absolute + relative * std::abs(expected)) {
    std::ostringstream out;
    out << std::setprecision(17) << message << ": actual=" << actual
        << ", expected=" << expected;
    throw std::runtime_error(out.str());
  }
}
void nearVector(const g::Vec3& actual, const g::Vec3& expected,
                double absolute, double relative, const std::string& message) {
  for (int k = 0; k < 3; ++k)
    near(actual[k], expected[k], absolute, relative,
         message + " component " + std::to_string(k));
}
void rejects(const std::function<void()>& operation, const std::string& message) {
  try { operation(); }
  catch (const std::domain_error&) { return; }
  throw std::runtime_error(message);
}
g::PairParameters freeOnly() {
  g::PairParameters p;
  p.hamaker = 0.;
  p.sigma = .4197e-9;
  p.freeCmcRepulsionPressure = 52000.;
  return p;
}
g::Vec3 support(const g::Body& b, const g::Vec3& normal) {
  g::Vec3 squared{};
  for (int k = 0; k < 3; ++k) squared[k] = b.axes[k] * b.axes[k];
  const auto qn = g::mul(g::rotatedDiagonal(b.rotation, squared), normal);
  return g::scale(qn, 1. / std::sqrt(g::dot(normal, qn)));
}
Bodies atGap(double gap, bool generic = false) {
  Bodies b;
  g::Vec3 n{0., 0., 1.};
  if (generic) {
    b.first.axes = {1.65e-6, 1.33e-6, .20e-6};
    b.second.axes = {1.45e-6, 1.75e-6, .23e-6};
    b.first.rotation = g::rotationIncrement({.19, -.28, .08});
    b.second.rotation = g::rotationIncrement({-.11, .21, -.15});
    n = g::normalized({.18, -.12, 1.});
  }
  b.first.position = {2.e-6, -3.e-6, 1.e-6};
  b.second.position = g::add(b.first.position,
      g::add(g::add(support(b.first, n), support(b.second, n)),
             g::scale(n, gap)));
  return b;
}
g::PairResult evaluate(const Bodies& b, const g::PairParameters& p) {
  return g::evaluatePair(b.first, b.second, p);
}
Bodies perturb(Bodies b, int coordinate, double amount) {
  if (coordinate < 3) b.second.position[coordinate] += amount;
  else {
    g::Vec3 angle{};
    angle[(coordinate - 3) % 3] = amount;
    auto& body = coordinate < 6 ? b.first : b.second;
    body.rotation = g::rotateLaboratory(body.rotation, angle);
  }
  return b;
}
double gradient(const g::PairResult& result, int coordinate) {
  if (coordinate < 3) return result.forceI[coordinate];
  return -(coordinate < 6 ? result.torqueI[coordinate - 3]
                          : result.torqueJ[coordinate - 6]);
}
bool samePhysics(const g::PairResult& a, const g::PairResult& b) {
  return a.energy == b.energy && a.ua == b.ua && a.ur == b.ur &&
      a.forceI == b.forceI && a.torqueI == b.torqueI && a.torqueJ == b.torqueJ &&
      a.forceAttractiveI == b.forceAttractiveI && a.forceRepulsiveI == b.forceRepulsiveI &&
      a.torqueAttractiveI == b.torqueAttractiveI && a.torqueAttractiveJ == b.torqueAttractiveJ &&
      a.torqueRepulsiveI == b.torqueRepulsiveI && a.torqueRepulsiveJ == b.torqueRepulsiveJ;
}

void analyticFaceToFace() {
  const auto p = freeOnly();
  // Aligned identical oblates: lambda_D = a^2/c, independent of separation.
  const auto contact = atGap(p.roughnessGap);
  const double lambda = std::pow(contact.first.axes[0], 2) / contact.first.axes[2];
  for (double h : {.9e-9, 2.e-9, 2.67e-9, 10.e-9, 20.e-9, 25.e-9, 45.e-9}) {
    const auto r = evaluate(atGap(h), p);
    const double pressure = p.freeCmcRepulsionPressure *
        std::exp(-(h - p.roughnessGap) / p.freeCmcRepulsionLength);
    const double force = pi * lambda * pressure * p.freeCmcRepulsionLength;
    near(r.ur, force * p.freeCmcRepulsionLength, 1.e-30, 2.e-10,
         "analytic FF repulsive energy, including beyond curvature cutoff");
    nearVector(r.forceI, {0., 0., -force}, 1.e-21, 2.e-10,
               "analytic FF repulsive force");
    nearVector(r.torqueI, {}, 1.e-26, 0., "centered first torque vanishes");
    nearVector(r.torqueJ, {}, 1.e-26, 0., "centered second torque vanishes");
    require(r.ua == 0. && r.energy == r.ur, "free-CMC energy belongs to repulsive branch");
  }
  // Positive-gap Newton trials below h0 retain the exponential, not a clamp.
  require(evaluate(atGap(.9e-9), p).ur > evaluate(contact, p).ur,
          "trial energy increases below h0 without distance clamping");
}

void conservativeGeometry() {
  auto p = freeOnly();
  // A longer test length keeps the far-switch term measurable; the same law is
  // exercised at every gap, with generic offsets and unequal rotated bodies.
  p.freeCmcRepulsionLength = 60.e-9;
  for (double h : {1.7e-9, 2.e-9, 2.4e-9, 8.e-9, 19.e-9, 20.e-9, 35.e-9, 450.e-9}) {
    const auto b = atGap(h, true);
    const auto result = evaluate(b, p);
    for (int coordinate = 0; coordinate < 9; ++coordinate) {
      const double epsilon = coordinate < 3 ? 2.e-13 : 2.e-7;
      const double numerical =
          (evaluate(perturb(b, coordinate, epsilon), p).energy -
           evaluate(perturb(b, coordinate, -epsilon), p).energy) / (2. * epsilon);
      near(gradient(result, coordinate), numerical,
           coordinate < 3 ? 1.e-17 : 1.e-23, 5.e-5,
           "translation/rotation energy gradient " + std::to_string(coordinate));
    }
    const auto moment = g::sub(g::add(result.torqueI, result.torqueJ),
        g::cross(g::sub(b.second.position, b.first.position), result.forceI));
    nearVector(moment, {}, 2.e-24 + 2.e-11 * g::norm(result.torqueI), 0.,
               "internal angular momentum conservation");
    require(result.forceI == result.forceRepulsiveI &&
            result.torqueI == result.torqueRepulsiveI &&
            result.torqueJ == result.torqueRepulsiveJ,
            "force and torque repulsive-branch bookkeeping");
  }
}

void disabledAndAdhesionUnchanged() {
  g::PairParameters p;
  p.sigma = .4197e-9;
  p.surfaceAdhesion = true;
  p.adhesionWork = .007666863221834307;
  auto dormant = p;
  dormant.freeCmcRepulsionLength = 1.e-6;
  auto active = p;
  active.freeCmcRepulsionPressure = 52000.;
  for (double h : {1.7e-9, 2.e-9, 2.4e-9, 8.e-9, 19.e-9, 25.e-9, 450.e-9}) {
    const auto b = atGap(h, true);
    const auto before = evaluate(b, p), zero = evaluate(b, dormant), after = evaluate(b, active);
    require(samePhysics(before, zero), "zero pressure ignores decay length bitwise");
    require(before.ua == after.ua && before.forceAttractiveI == after.forceAttractiveI &&
            before.torqueAttractiveI == after.torqueAttractiveI &&
            before.torqueAttractiveJ == after.torqueAttractiveJ,
            "free pressure does not change adhesion or its derivatives");
  }
  // Zero pressure also leaves the original RE2 and legacy local-gap modes alone.
  for (bool legacy : {false, true}) {
    p.surfaceAdhesion = false;
    p.localGapFraction = legacy ? .1 : 0.;
    dormant = p;
    dormant.freeCmcRepulsionLength = 1.e-6;
    require(samePhysics(evaluate(atGap(2.2e-9, true), p),
                        evaluate(atGap(2.2e-9, true), dormant)),
            "disabled free-CMC term preserves legacy modes");
  }
}

void pressureLinearityAndObjectivity() {
  auto p = freeOnly(), twice = p;
  twice.freeCmcRepulsionPressure *= 2.;
  const auto b = atGap(11.e-9, true);
  const auto result = evaluate(b, p), doubled = evaluate(b, twice);
  near(doubled.energy, 2. * result.energy, 1.e-29, 2.e-14, "pressure scales energy linearly");
  nearVector(doubled.forceI, g::scale(result.forceI, 2.), 1.e-22, 2.e-14,
             "pressure scales force linearly");
  nearVector(doubled.torqueI, g::scale(result.torqueI, 2.), 1.e-28, 2.e-14,
             "pressure scales torque linearly");
  const auto rotation = g::rotationIncrement({.71, -.36, .43});
  Bodies rotated = b;
  for (auto* body : {&rotated.first, &rotated.second}) {
    body->position = g::add(g::mul(rotation, body->position), {4.e-6, -2.e-6, 3.e-6});
    body->rotation = g::multiply(rotation, body->rotation);
  }
  const auto transformed = evaluate(rotated, p);
  near(transformed.energy, result.energy, 1.e-28, 2.e-9, "rigid-motion invariant energy");
  nearVector(transformed.forceI, g::mul(rotation, result.forceI), 1.e-19, 2.e-8,
             "rigid-motion covariant force");
  nearVector(transformed.torqueI, g::mul(rotation, result.torqueI), 1.e-25, 2.e-8,
             "rigid-motion covariant torque");
  const auto swapped = evaluate({b.second, b.first}, p);
  near(swapped.energy, result.energy, 1.e-28, 2.e-10, "particle exchange energy");
  nearVector(swapped.forceI, g::scale(result.forceI, -1.), 1.e-19, 2.e-9,
             "particle exchange force");
  nearVector(swapped.torqueI, result.torqueJ, 1.e-25, 2.e-9,
             "particle exchange torque");
}

void farSwitch() {
  auto p = freeOnly();
  p.freeCmcRepulsionLength = 100.e-9;
  const auto b = atGap(p.roughnessGap);
  const double lambda = std::pow(b.first.axes[0], 2) / b.first.axes[2];
  for (double h : {399.e-9, 400.e-9, 401.e-9, 450.e-9, 499.e-9, 500.e-9, 501.e-9}) {
    const auto r = evaluate(atGap(h), p);
    double sw = 1., derivative = 0.;
    if (h >= p.cutoffGap) sw = 0.;
    else if (h > p.switchGap) {
      const double t = (h - p.switchGap) / (p.cutoffGap - p.switchGap);
      sw = 1. - 10.*t*t*t + 15.*t*t*t*t - 6.*t*t*t*t*t;
      derivative = -30.*t*t*(1.-t)*(1.-t)/(p.cutoffGap-p.switchGap);
    }
    const double unswitched = pi * lambda * p.freeCmcRepulsionPressure *
        std::pow(p.freeCmcRepulsionLength, 2) *
        std::exp(-(h-p.roughnessGap)/p.freeCmcRepulsionLength);
    near(r.energy, unswitched * sw, 2.e-29, 2.e-9, "far switch multiplies full energy");
    near(r.forceI[2], unswitched * (derivative-sw/p.freeCmcRepulsionLength),
         2.e-20, 2.e-9, "far-switch derivative contributes to force");
  }
  for (double h : {p.switchGap, p.cutoffGap}) {
    const auto before = evaluate(atGap(h-1.e-13, true), p);
    const auto after = evaluate(atGap(h+1.e-13, true), p);
    near(before.energy, after.energy, 1.e-25, 5.e-6, "far-boundary energy continuity");
    nearVector(before.forceI, after.forceI, 2.e-17, 5.e-5, "far-boundary force continuity");
    nearVector(before.torqueI, after.torqueI, 2.e-23, 5.e-5, "far-boundary torque continuity");
  }
}

void invalidInputs() {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const std::vector<std::function<void(g::PairParameters&)>> corruptions{
      [](auto& p) { p.freeCmcRepulsionPressure = -1.; },
      [=](auto& p) { p.freeCmcRepulsionPressure = nan; },
      [=](auto& p) { p.freeCmcRepulsionPressure = inf; },
      [](auto& p) { p.freeCmcRepulsionLength = 0.; },
      [](auto& p) { p.freeCmcRepulsionLength = -5.e-9; },
      [=](auto& p) { p.freeCmcRepulsionLength = nan; },
      [=](auto& p) { p.freeCmcRepulsionLength = inf; },
      [](auto& p) { p.roughnessGap = 0.; },
      [](auto& p) { p.roughnessGap = -1.e-9; },
      [=](auto& p) { p.roughnessGap = nan; },
      [](auto& p) { p.freeCmcRepulsionPressure = 1.e308; p.freeCmcRepulsionLength = 10.; },
      [](auto& p) { p.freeCmcRepulsionPressure = 1.e200; p.freeCmcRepulsionLength = 1.e100; },
      [](auto& p) { p.freeCmcRepulsionLength = std::numeric_limits<double>::denorm_min(); }};
  for (const auto& corrupt : corruptions) {
    auto p = freeOnly();
    corrupt(p);
    rejects([&] { g::validatePairParameters(p); }, "invalid free-CMC parameters rejected");
  }
  auto disabled = freeOnly();
  disabled.freeCmcRepulsionPressure = 0.;
  disabled.freeCmcRepulsionLength = 0.;
  rejects([&] { g::validatePairParameters(disabled); }, "length remains valid input when disabled");
}
} // namespace

int main() {
  const std::vector<std::pair<std::string, std::function<void()>>> checks{
      {"analytic FF force, extended range, and unclamped trials", analyticFaceToFace},
      {"all nine conservative derivatives and angular momentum", conservativeGeometry},
      {"zero pressure and unchanged adhesion", disabledAndAdhesionUnchanged},
      {"pressure linearity, covariance, and particle exchange", pressureLinearityAndObjectivity},
      {"far switch energy, derivatives, and continuity", farSwitch},
      {"invalid free-CMC input rejection", invalidInputs}};
  int failures = 0;
  for (const auto& check : checks) {
    try { check.second(); std::cout << "PASS " << check.first << '\n'; }
    catch (const std::exception& e) {
      ++failures;
      std::cerr << "FAIL " << check.first << ": " << e.what() << '\n';
    }
  }
  return failures ? 1 : 0;
}
