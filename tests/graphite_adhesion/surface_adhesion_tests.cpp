// Standalone checks of the curvature-resolved finite-range adhesion model.
/* c++ -std=c++17 -O2 -Wall -Wextra -pedantic \
     -I olb-1.9r0/src/slurry/gr_re2 \
     tests/graphite_adhesion/surface_adhesion_tests.cpp -o /tmp/surface_adhesion_tests
   /tmp/surface_adhesion_tests */
#include "reSquaredPotential.h"
#include "contactCurvature.h"
#include "nearFieldResistance.h"

#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
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
g::PairParameters parameters() {
  g::PairParameters p;
  p.surfaceAdhesion = true;
  p.hamaker = 9.9e-20;
  p.sigma = .4197e-9;
  p.roughnessGap = 2.e-9;
  p.adhesionWork = .0219;
  p.adhesionRange = .77e-9;
  p.curvatureSwitchGap = 5.e-9;
  p.curvatureCutoffGap = 20.e-9;
  return p;
}
// Independent planar integrated 12-6 LJ reference, evaluated in long double.
double planarWork(const g::PairParameters& p, double separation) {
  const long double h = separation;
  const long double ratio = static_cast<long double>(p.sigma) / h;
  return static_cast<double>(p.hamaker / (12.L * pi * h * h) *
                             (1.L - std::pow(ratio, 6) / 30.L));
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
// Solve the symmetric spheroid support geometry without calling the production
// closest-gap routine. The horizontal centre displacement is prescribed.
Bodies offsetPair(double offset, double gap) {
  const double a = 1.65e-6, c = .20e-6;
  double lo = 0., hi = .5 * pi;
  for (int iteration = 0; iteration < 100; ++iteration) {
    const double theta = .5 * (lo + hi), sn = std::sin(theta), cs = std::cos(theta);
    const double projected = std::sqrt(a * a * sn * sn + c * c * cs * cs);
    const double x = 2. * a * a * sn / projected + gap * sn;
    (x < offset ? lo : hi) = theta;
  }
  const double theta = .5 * (lo + hi), sn = std::sin(theta), cs = std::cos(theta);
  const double projected = std::sqrt(a * a * sn * sn + c * c * cs * cs);
  Bodies b;
  b.second.position = {offset, 0., 2. * c * c * cs / projected + gap * cs};
  return b;
}
g::ContactCurvatureResult curvature(const Bodies& b) {
  return g::contactCurvature(b.first, b.second,
                            g::closestEllipsoidGap(b.first, b.second));
}
double lubricationRootDet(const Bodies& b) {
  const auto n = g::closestEllipsoidGap(b.first, b.second).normal;
  const auto e1 = g::normalized(g::cross(n, std::abs(n[0]) < .8
                                      ? g::Vec3{1., 0., 0.} : g::Vec3{0., 1., 0.}));
  const auto e2 = g::cross(n, e1);
  auto c = g::nfCurvature(b.first, n);
  const auto cj = g::nfCurvature(b.second, n);
  for (int k = 0; k < 9; ++k) c[k] += cj[k];
  return std::sqrt(g::dot(e1, g::mul(c, e1)) * g::dot(e2, g::mul(c, e2)) -
                   std::pow(g::dot(e1, g::mul(c, e2)), 2));
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

void normalizationAndRange() {
  const auto p = parameters();
  const auto b = atGap(p.roughnessGap);
  const double radius = b.first.axes[0] * b.first.axes[0] / (2. * b.first.axes[2]);
  const double geometricFactor = 2. * pi * radius;
  const double excess = p.adhesionWork - planarWork(p, p.roughnessGap);
  near(g::dot(evaluate(b, p).forceI, {0., 0., 1.}),
       geometricFactor * p.adhesionWork, 1.e-18, 2.e-11,
       "FF contact force = 2 pi R_eff W, including baseline repulsion");
  near(curvature(b).length, 2. * radius, 1.e-18, 1.e-11,
       "centered oblate analytic curvature");
  for (double openingFraction : {0., .25, .7, 1., 1.3}) {
    const double opening = openingFraction * p.adhesionRange;
    const auto opened = atGap(p.roughnessGap + opening);
    const double expectedWork = planarWork(p, p.roughnessGap + opening) +
        excess * std::max(0., 1. - openingFraction);
    near(g::dot(evaluate(opened, p).forceI, {0., 0., 1.}),
         geometricFactor * expectedWork, 1.e-18, 2.e-10,
         "opening force equals integrated finite-range surface traction");
  }
  // Adhesion range changes separation work, but not force at closed contact.
  auto shortRange = p;
  shortRange.adhesionRange = .31e-9;
  near(evaluate(b, shortRange).forceI[2], evaluate(b, p).forceI[2],
       1.e-18, 1.e-11, "range does not rescale contact force");
  near(evaluate(b, p).energy - evaluate(b, shortRange).energy,
       -.5 * geometricFactor * excess * (p.adhesionRange - shortRange.adhesionRange),
       1.e-27, 1.e-11, "range controls pair separation energy");
}

void allNineDerivatives() {
  const auto p = parameters();
  // Contact, open cohesive zone, free background, curvature blend, and far switch.
  for (double gap : {1.7e-9, 2.e-9, 2.4e-9, 4.e-9, 9.e-9, 15.e-9, 19.e-9, 43.e-9, 450.e-9}) {
    const auto b = atGap(gap, true);
    const auto result = evaluate(b, p);
    for (int coordinate = 0; coordinate < 9; ++coordinate) {
      const double epsilon = coordinate < 3 ? 2.e-14 : 2.e-8;
      const double derivative =
          (evaluate(perturb(b, coordinate, epsilon), p).energy -
           evaluate(perturb(b, coordinate, -epsilon), p).energy) / (2. * epsilon);
      near(gradient(result, coordinate), derivative,
           coordinate < 3 ? 2.e-16 : 2.e-22, 4.e-5,
           "all-nine conservative derivatives, coordinate " + std::to_string(coordinate));
    }
    nearVector(g::add(result.forceAttractiveI, result.forceRepulsiveI),
               result.forceI, 0., 0., "force branch bookkeeping");
    near(result.ua + result.ur, result.energy, 0., 0., "energy branch bookkeeping");
    const auto moment = g::sub(g::add(result.torqueI, result.torqueJ),
        g::cross(g::sub(b.second.position, b.first.position), result.forceI));
    nearVector(moment, {}, 2.e-25 + 2.e-11 * g::norm(result.torqueI), 0.,
               "internal angular momentum conservation");
  }
}

void curvatureDerivatives() {
  for (double gap : {2.e-9, 9.e-9, 35.e-9}) {
    const auto b = atGap(gap, true);
    const auto c = curvature(b);
    near(c.rootDet, lubricationRootDet(b), 1.e-7, 2.e-12,
         "adhesion and lubrication use the same tangent-plane curvature");
    for (int coordinate = 0; coordinate < 9; ++coordinate) {
      const double epsilon = coordinate < 3 ? 1.e-12 : 1.e-7;
      // Finite differencing the existing lubrication geometry is deliberately
      // separate from the production implicit derivative implementation.
      const double numerical =
          (2. / lubricationRootDet(perturb(b, coordinate, epsilon)) -
           2. / lubricationRootDet(perturb(b, coordinate, -epsilon))) / (2. * epsilon);
      near(c.lengthDerivative[coordinate], numerical,
           coordinate < 3 ? 2.e-8 : 2.e-14, 2.e-5,
           "implicit closest-normal curvature derivative " + std::to_string(coordinate));
    }
  }
}

void covarianceAndSwap() {
  const auto p = parameters();
  const auto rotation = g::rotationIncrement({.71, -.36, .43});
  for (double gap : {2.e-9, 2.4e-9, 31.e-9, 450.e-9}) {
    const auto b = atGap(gap, true);
    const auto original = evaluate(b, p);
    Bodies transformed = b;
    for (auto* body : {&transformed.first, &transformed.second}) {
      body->position = g::add(g::mul(rotation, body->position), {4.e-6, -2.e-6, 3.e-6});
      body->rotation = g::multiply(rotation, body->rotation);
    }
    const auto rotated = evaluate(transformed, p);
    near(rotated.energy, original.energy, 1.e-28, 2.e-9, "energy is objective");
    nearVector(rotated.forceI, g::mul(rotation, original.forceI),
               1.e-18, 2.e-8, "force rotates objectively");
    nearVector(rotated.torqueI, g::mul(rotation, original.torqueI),
               1.e-25, 2.e-8, "first torque rotates objectively");
    nearVector(rotated.torqueJ, g::mul(rotation, original.torqueJ),
               1.e-25, 2.e-8, "second torque rotates objectively");
    const auto swapped = evaluate({b.second, b.first}, p);
    near(swapped.energy, original.energy, 1.e-28, 2.e-10, "pair swap energy symmetry");
    nearVector(swapped.forceI, g::scale(original.forceI, -1.),
               1.e-18, 2.e-9, "pair swap force antisymmetry");
    nearVector(swapped.torqueI, original.torqueJ, 1.e-25, 2.e-9,
               "pair swap torque exchange");
  }
}

void switchesAndFarField() {
  const auto p = parameters();
  auto legacy = p;
  legacy.surfaceAdhesion = false;
  for (double gap : {50.e-9 + 1.e-14, 70.e-9, 410.e-9, 499.e-9, 501.e-9}) {
    const auto b = atGap(gap, true);
    const auto a = evaluate(b, p), old = evaluate(b, legacy);
    require(a.energy == old.energy && a.forceI == old.forceI &&
            a.torqueI == old.torqueI && a.torqueJ == old.torqueJ,
            "outside curvature cutoff the old RE2 law is bitwise unchanged");
  }
  for (double edge : {p.roughnessGap, p.roughnessGap + p.adhesionRange,
                      p.curvatureSwitchGap, p.curvatureCutoffGap,
                      p.switchGap, p.cutoffGap}) {
    constexpr double epsilon = 1.e-16;
    const auto before = evaluate(atGap(edge - epsilon, true), p);
    const auto after = evaluate(atGap(edge + epsilon, true), p);
    near(before.energy, after.energy, 1.e-25, 5.e-6, "junction energy continuity");
    nearVector(before.forceI, after.forceI, 2.e-13, 5.e-5, "junction force continuity");
    nearVector(before.torqueI, after.torqueI, 2.e-19, 5.e-5, "junction torque continuity");
  }
  // Opening along the actual contact normal exercises FF, offset FF, EF,
  // EE, and a general triaxial contact without changing the near curvature.
  // This catches an artificial repulsive barrier introduced by energy matching.
  std::vector<Bodies> startingPairs{atGap(p.roughnessGap),
      offsetPair(.25 * 1.65e-6, p.roughnessGap),
      offsetPair(.5 * 1.65e-6, p.roughnessGap),
      offsetPair(1.65e-6, p.roughnessGap), atGap(p.roughnessGap, true)};
  for (bool bothEdges : {false, true}) {
    Bodies b;
    b.first.rotation = g::rotationIncrement({0., .5 * pi, 0.});
    if (bothEdges) b.second.rotation = b.first.rotation;
    const g::Vec3 n{0., 0., 1.};
    b.second.position = g::add(g::add(support(b.first, n), support(b.second, n)),
                               g::scale(n, p.roughnessGap));
    startingPairs.push_back(b);
  }
  // Keep a modest deterministic subset of the independent orientation audit.
  std::mt19937_64 random(975350);
  std::normal_distribution<double> normal(0., 1.);
  std::uniform_real_distribution<double> angle(0., 2. * pi);
  const auto unitVector = [&] {
    return g::normalized({normal(random), normal(random), normal(random)});
  };
  for (int sample = 0; sample < 32; ++sample) {
    Bodies b;
    const auto firstAxis = unitVector();
    const double firstAngle = angle(random);
    const auto secondAxis = unitVector();
    const double secondAngle = angle(random);
    b.first.rotation = g::rotationIncrement(g::scale(firstAxis, firstAngle));
    b.second.rotation = g::rotationIncrement(g::scale(secondAxis, secondAngle));
    const auto n = unitVector();
    b.second.position = g::add(g::add(support(b.first, n), support(b.second, n)),
                               g::scale(n, p.roughnessGap));
    startingPairs.push_back(b);
  }
  for (const auto& initial : startingPairs) {
    const auto n = g::closestEllipsoidGap(initial.first, initial.second).normal;
    double previous = -std::numeric_limits<double>::infinity();
    for (int k = 0; k < 250; ++k) {
      auto b = initial;
      b.second.position = g::add(b.second.position, g::scale(n, k * .25e-9));
      const auto result = evaluate(b, p);
      require(g::dot(result.forceI, n) > 0., "separation remains attractive across matching");
      require(result.energy > previous, "separation has no artificial matching barrier");
      previous = result.energy;
    }
  }
}

void offsetGeometry() {
  auto p = parameters();
  // Exactly zero excess isolates the geometric change. Expected energies below
  // are analytic spheroid results, independent of the production helper.
  p.adhesionWork = g::surfaceBackgroundWork(p);
  auto legacy = p;
  legacy.surfaceAdhesion = false;
  const double a = 1.65e-6, c = .20e-6, gap = p.roughnessGap;
  const double fractions[] = {0., .05, .10, .25, .5, 1.};
  const double oldRatios[] = {1., .961, .860, .498463, .206871, .075295};
  const double alignedFiniteRatio = evaluate(atGap(gap), legacy).ua /
                                    (-p.hamaker * a * a / (12. * c * gap));
  for (int index = 0; index < 6; ++index) {
    const double offset = fractions[index] * a;
    double lo = 0., hi = .5 * pi;
    for (int iteration = 0; iteration < 100; ++iteration) {
      const double theta = .5 * (lo + hi), sn = std::sin(theta), cs = std::cos(theta);
      const double projected = std::sqrt(a * a * sn * sn + c * c * cs * cs);
      const double x = 2. * a * a * sn / projected + gap * sn;
      (x < offset ? lo : hi) = theta;
    }
    const double theta = .5 * (lo + hi), sn = std::sin(theta), cs = std::cos(theta);
    const double projectedSquared = a * a * sn * sn + c * c * cs * cs;
    const double analyticLength = a * a * c / projectedSquared;
    Bodies b;
    b.second.position = {offset, 0., 2. * c * c * cs / std::sqrt(projectedSquared) + gap * cs};
    const double expected = -p.hamaker * analyticLength / (12. * gap);
    const auto result = evaluate(b, p);
    near(result.gap, gap, 1.e-18, 0., "fixed-gap offset geometry");
    near(curvature(b).length, analyticLength, 1.e-18, 2.e-10,
         "offset curvature equals analytic spheroid curvature");
    near(result.ua, expected, 1.e-27, 2.e-9,
         "offset attraction uses local surface curvature without centerline attenuation");
    const double observedRatio = evaluate(b, legacy).ua / expected / alignedFiniteRatio;
    near(observedRatio, oldRatios[index], 0., .005,
         "legacy offset attenuation reproduces the reviewed discrepancy");
  }
}

void limitsAndDomains() {
  auto p = parameters();
  p.adhesionWork = g::surfaceBackgroundWork(p);
  const auto b = atGap(p.roughnessGap);
  const double radius = b.first.axes[0] * b.first.axes[0] / (2. * b.first.axes[2]);
  near(evaluate(b, p).forceI[2], 2. * pi * radius * planarWork(p, p.roughnessGap),
       1.e-18, 2.e-11, "zero excess leaves the baseline contact force");
  p.hamaker = 0.;
  near(evaluate(b, p).forceI[2], 2. * pi * radius * p.adhesionWork,
       1.e-18, 2.e-11, "zero dispersion leaves the declared surface adhesion");
  const auto zero = evaluate(atGap(p.roughnessGap + 2. * p.adhesionRange), p);
  require(zero.energy == 0. && g::norm(zero.forceI) == 0. && g::norm(zero.torqueI) == 0.,
          "zero dispersion gives zero interaction outside the cohesive range");
  p = parameters();
  near(g::minimumPairGap(p), 0., 0., 0., "new law has no shifted-gap singularity");
  require(std::isfinite(evaluate(atGap(.9e-9), p).energy),
          "positive-gap solver trials below h0 have a smooth energy extension");
  rejects([&] { evaluate(atGap(-.1e-9), p); }, "actual overlap is outside the pair domain");
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  const std::vector<std::function<void(g::PairParameters&)>> corruptions{
      [](auto& x) { x.localGapFraction = 1.; },
      [](auto& x) { x.roughnessGap = 0.; },
      [](auto& x) { x.adhesionWork = 0.; },
      [=](auto& x) { x.adhesionWork = nan; },
      [](auto& x) { x.adhesionRange = 0.; },
      [=](auto& x) { x.adhesionRange = nan; },
      [](auto& x) { x.curvatureSwitchGap = x.roughnessGap * .5; },
      [](auto& x) { x.curvatureCutoffGap = x.curvatureSwitchGap; },
      [](auto& x) { x.curvatureCutoffGap = x.cutoffGap; }};
  for (const auto& corrupt : corruptions) {
    auto invalid = p;
    corrupt(invalid);
    rejects([&] { g::validatePairParameters(invalid); },
            "incompatible or invalid new-mode parameters must be rejected");
  }
}
} // namespace

int main() {
  const std::vector<std::pair<std::string, std::function<void()>>> checks{
      {"surface-work normalization and finite range", normalizationAndRange},
      {"all nine energy derivatives and angular momentum", allNineDerivatives},
      {"curvature/lubrication agreement and implicit derivatives", curvatureDerivatives},
      {"rigid-motion covariance and particle exchange", covarianceAndSwap},
      {"junction continuity and unchanged far field", switchesAndFarField},
      {"analytic offset-flake attraction correction", offsetGeometry},
      {"zero-excess limits and invalid parameter rejection", limitsAndDomains}};
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
