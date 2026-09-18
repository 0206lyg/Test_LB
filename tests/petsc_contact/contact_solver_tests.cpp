// Standalone physical regression fixtures for the production particle solver.
// This translation unit deliberately does not include or link OpenLB.
#include "particleSubsteps.h"
#ifdef SLURRY_USE_PETSC
#include <petscsys.h>
#endif

#include <functional>
#include <cstdlib>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace g = slurry::gr_re2::graphite;
namespace {

constexpr double radius = 1.e-6;
constexpr double mass = 2.e-11;
constexpr double step = 1.e-6;
constexpr double load = 1.e-9;

void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double absolute, double relative,
          const std::string& message) {
  if (!std::isfinite(actual) ||
      std::abs(actual - expected) > absolute + relative * std::abs(expected)) {
    std::ostringstream out;
    out << std::setprecision(17) << message << ": actual=" << actual
        << ", expected=" << expected << ", allowed="
        << absolute + relative * std::abs(expected);
    throw std::runtime_error(out.str());
  }
}

void nearVector(const g::Vec3& actual, const g::Vec3& expected, double absolute,
                double relative, const std::string& message) {
  for (int axis = 0; axis < 3; ++axis)
    near(actual[axis], expected[axis], absolute, relative,
         message + " component " + std::to_string(axis));
}

bool identicalBody(const g::Body& a, const g::Body& b) {
  return a.position == b.position && a.velocity == b.velocity && a.omega == b.omega &&
      a.rotation == b.rotation && a.axes == b.axes && a.inertiaBody == b.inertiaBody &&
      a.mass == b.mass;
}

bool identicalContact(const g::RoughContactState& a, const g::RoughContactState& b) {
  return a.active == b.active && a.sliding == b.sliding && a.rolling == b.rolling &&
      a.normal == b.normal && a.elasticSlip == b.elasticSlip &&
      a.elasticRoll == b.elasticRoll && a.rollingCap == b.rollingCap &&
      a.rollingStiffness == b.rollingStiffness && a.normalLoad == b.normalLoad &&
      a.elasticEnergy == b.elasticEnergy && a.plasticSlipWork == b.plasticSlipWork &&
      a.plasticRollWork == b.plasticRollWork && a.releasedEnergy == b.releasedEnergy &&
      a.stepDuration == b.stepDuration && a.tangentForce == b.tangentForce &&
      a.rollingTorque == b.rollingTorque && a.leverI == b.leverI && a.leverJ == b.leverJ;
}

#ifdef SLURRY_USE_PETSC
struct TemporaryDiagnostics {
  std::filesystem::path directory;
  TemporaryDiagnostics() {
    const auto tag = std::chrono::steady_clock::now().time_since_epoch().count();
    directory = std::filesystem::temp_directory_path() /
        ("gr_petsc_fixture_" + std::to_string(tag));
    std::filesystem::create_directory(directory);
  }
  ~TemporaryDiagnostics() {
    std::error_code error;
    std::filesystem::remove_all(directory, error);
  }
  std::string prefix() const { return (directory / "particle_solver_rank0").string(); }
};
#endif

g::ParticleStepSettings settings() {
  g::ParticleStepSettings s;
#ifdef SLURRY_USE_PETSC
  s.solverBackend = "petsc";
  if (const char* prefix = std::getenv("CONTACT_TEST_DIAGNOSTICS"))
    s.solverDiagnosticsPrefix = prefix;
#endif
  s.box = {20.e-6, 20.e-6, 20.e-6};
  s.shearRate = 0.;
  // Isolate the rigid rough-contact law from attraction and lubrication.
  s.pair.hamaker = 0.;
  s.nearField.enabled = false;
  s.rough.enabled = true;
  s.rough.friction = 1.;
  s.maxSubsteps = 1;
  s.maxNewtonIterations = 80;
  s.maxKrylovIterations = 200;
  s.maxLineSearch = 24;
  s.relativeTolerance = 1.e-7;
  s.forceAbsoluteTolerance = 1.e-16;
  s.torqueAbsoluteTolerance = 1.e-22;
  s.contactGapTolerance = 1.e-13;
  return s;
}

g::Body sphere(const g::Vec3& position) {
  g::Body body;
  body.position = position;
  body.axes = {radius, radius, radius};
  body.mass = mass;
  const double inertia = .4 * mass * radius * radius;
  body.inertiaBody = {inertia, inertia, inertia};
  return body;
}

std::vector<g::Body> pairAtGap(double gap) {
  const double distance = 2. * radius + gap;
  return {sphere({8.e-6 - .5 * distance, 10.e-6, 10.e-6}),
          sphere({8.e-6 + .5 * distance, 10.e-6, 10.e-6})};
}

g::Vec3 momentum(const std::vector<g::Body>& bodies) {
  g::Vec3 total{};
  for (const auto& body : bodies)
    total = g::add(total, g::scale(body.velocity, body.mass));
  return total;
}

double kineticEnergy(const std::vector<g::Body>& bodies) {
  double total = 0.;
  for (const auto& body : bodies) {
    total += .5 * body.mass * g::dot(body.velocity, body.velocity);
    total += .5 * g::dot(body.omega, g::mul(
        g::rotatedDiagonal(body.rotation, body.inertiaBody), body.omega));
  }
  return total;
}

void accepted(const g::ParticleStepDiagnostics& d) {
  require(d.substeps == 1, "fixture must converge without timestep subdivision");
  require(std::isfinite(d.maxForceResidualRatio) && d.maxForceResidualRatio <= 1.,
          "accepted force balance must satisfy the configured tolerance");
  require(std::isfinite(d.maxTorqueResidualRatio) && d.maxTorqueResidualRatio <= 1.,
          "accepted torque balance must satisfy the configured tolerance");
  require(d.contactGapViolation <= settings().contactGapTolerance,
          "accepted state penetrates the rough-contact gap");
}

void freeFlight() {
  auto s = settings();
  s.rough.enabled = false;
  std::vector<g::Body> bodies{sphere({5.e-6, 6.e-6, 7.e-6})};
  bodies[0].velocity = {2.e-4, -3.e-4, 4.e-4};
  bodies[0].omega = {0., 0., 12.};
  const auto initial = bodies[0];
  const std::vector<g::Vec3> zero(1);
  const auto d = g::advanceParticles(bodies, zero, zero, step, 0., s);
  accepted(d);
  nearVector(bodies[0].position,
             g::add(initial.position, g::scale(initial.velocity, step)),
             1.e-17, 1.e-12, "free flight displacement");
  nearVector(bodies[0].velocity, initial.velocity, 1.e-14, 1.e-10,
             "free flight velocity");
  nearVector(bodies[0].omega, initial.omega, 1.e-10, 1.e-10,
             "free spherical rotation");
  const auto expected = g::rotationIncrement(g::scale(initial.omega, step));
  for (int k = 0; k < 9; ++k)
    near(bodies[0].rotation[k], expected[k], 1.e-12, 1.e-10,
         "free rotation matrix");
}

void forcedFreeMotion() {
  auto s = settings();
  s.rough.enabled = false;
  std::vector<g::Body> bodies{sphere({5.e-6, 6.e-6, 7.e-6})};
  bodies[0].velocity = {2.e-5, -1.e-5, 0.};
  const auto initial = bodies[0];
  const std::vector<g::Vec3> force{{load, -.5 * load, .25 * load}};
  const std::vector<g::Vec3> torque{{0., 0., 1.e-18}};
  const auto d = g::advanceParticles(bodies, force, torque, step, 0., s);
  accepted(d);
  const auto expectedVelocity = g::add(initial.velocity,
      g::scale(force[0], step / mass));
  nearVector(bodies[0].velocity, expectedVelocity, 1.e-11, 1.e-6,
             "backward Euler forced velocity");
  nearVector(bodies[0].position,
             g::add(initial.position, g::scale(expectedVelocity, step)),
             1.e-16, 1.e-10, "backward Euler forced displacement");
  near(bodies[0].omega[2], step * torque[0][2] / initial.inertiaBody[2],
       1.e-6, 1.e-6, "angular impulse");
}

void normalEquilibrium() {
  const auto s = settings();
  auto bodies = pairAtGap(s.rough.gap);
  const auto initial = bodies;
  const std::vector<g::Vec3> force{{load, 0., 0.}, {-load, 0., 0.}}, zero(2);
  std::vector<g::GapCache> cache;
  g::PersistentContactState history;
  for (int k = 0; k < 3; ++k) {
    accepted(g::advanceParticles(bodies, force, zero, step, k * step, s,
                                &cache, &history));
    require(history.size() == 1 && history[0].active,
            "compressive load must activate the pair");
    near(history[0].normalLoad, load, 3.e-15, 3.e-6,
         "normal reaction balances the external load");
    near(g::closestEllipsoidGap(bodies[0], bodies[1]).gap, s.rough.gap,
         s.contactGapTolerance, 0., "equilibrium rough-contact gap");
    for (int i = 0; i < 2; ++i)
      nearVector(bodies[i].position, initial[i].position, 3.e-14, 0.,
                 "equilibrium particle position");
    nearVector(momentum(bodies), {}, 1.e-23, 0.,
               "internal contact conserves pair linear momentum");
  }
}

void approachingContact() {
  const auto s = settings();
  constexpr double initialExtraGap = 1.e-10;
  auto bodies = pairAtGap(s.rough.gap + initialExtraGap);
  bodies[0].velocity = {1.e-4, 0., 0.};
  bodies[1].velocity = {-1.e-4, 0., 0.};
  const double initialEnergy = kineticEnergy(bodies);
  const std::vector<g::Vec3> zero(2);
  g::PersistentContactState history;
  accepted(g::advanceParticles(bodies, zero, zero, step, 0., s, nullptr, &history));
  require(history[0].active, "approaching pair must create a contact");
  near(g::closestEllipsoidGap(bodies[0], bodies[1]).gap, s.rough.gap,
       s.contactGapTolerance, 0., "new contact gap");
  const double expectedSpeed = initialExtraGap / (2. * step);
  near(bodies[0].velocity[0], expectedSpeed, 8.e-8, 1.e-6,
       "inelastic constrained step speed");
  near(bodies[1].velocity[0], -expectedSpeed, 8.e-8, 1.e-6,
       "opposing inelastic constrained step speed");
  near(history[0].normalLoad, mass * (1.e-4 - expectedSpeed) / step,
       2.e-12, 1.e-6, "normal contact impulse");
  nearVector(momentum(bodies), {}, 1.e-23, 0., "contact impulse momentum balance");
  require(kineticEnergy(bodies) <= initialEnergy * (1. + 1.e-8),
          "passive contact must not increase kinetic energy");
}

g::RoughContactState preloadedState(const g::ParticleStepSettings& s) {
  g::RoughContactState state;
  state.active = true;
  state.normal = {1., 0., 0.};
  state.normalLoad = load;
  state.elasticSlip = {0., 2.e-11, 0.};
  state.elasticRoll = {0., 0., 2.e-4};
  state.rollingCap = 1.e-16;
  state.rollingStiffness = state.rollingCap / s.rough.rollingYieldAngle;
  state.tangentForce = g::scale(state.elasticSlip, -s.rough.tangentialStiffness);
  state.rollingTorque = g::scale(state.elasticRoll, -state.rollingStiffness);
  state.elasticEnergy = .5 * s.rough.tangentialStiffness *
      g::dot(state.elasticSlip, state.elasticSlip) + .5 * state.rollingStiffness *
      g::dot(state.elasticRoll, state.elasticRoll);
  return state;
}

void separatingContactClearsHistory() {
  const auto s = settings();
  auto bodies = pairAtGap(s.rough.gap);
  const auto initial = bodies;
  const auto state = preloadedState(s);
  g::PersistentContactState history{state};
  const std::vector<g::Vec3> force{{-load, 0., 0.}, {load, 0., 0.}}, zero(2);
  accepted(g::advanceParticles(bodies, force, zero, step, 0., s, nullptr, &history));
  require(!history[0].active, "separating normal contact must release");
  near(history[0].normalLoad, 0., 1.e-20, 0., "released normal load");
  nearVector(history[0].elasticSlip, {}, 1.e-20, 0., "released sliding history");
  nearVector(history[0].elasticRoll, {}, 1.e-15, 0., "released rolling history");
  near(history[0].releasedEnergy, state.elasticEnergy, 1.e-29, 1.e-8,
       "released contact elastic energy is accounted once");
  for (int i = 0; i < 2; ++i) {
    const auto expectedVelocity = g::scale(force[i], step / mass);
    nearVector(bodies[i].velocity, expectedVelocity, 2.e-9, 1.e-6,
               "detached particle feels no stale friction");
    nearVector(bodies[i].position,
               g::add(initial[i].position, g::scale(expectedVelocity, step)),
               5.e-15, 0., "detached particle displacement");
  }
  require(g::closestEllipsoidGap(bodies[0], bodies[1]).gap > s.rough.gap,
          "separating bodies have a positive gap slack");
  accepted(g::advanceParticles(bodies, force, zero, step, step, s, nullptr, &history));
  near(history[0].releasedEnergy, 0., 1.e-29, 0.,
       "released elasticity is not counted a second time");
}

void slidingRollingHistoryPersistence() {
  const auto s = settings();
  auto bodies = pairAtGap(s.rough.gap);
  const auto initial = bodies;
  const auto state = preloadedState(s);
  g::PersistentContactState history{state};
  // Perturb only the initial normal-load guess so this also needs a solve.
  history[0].normalLoad = .8 * load;
  const double arm = radius + .5 * s.rough.gap;
  const double tangential = s.rough.tangentialStiffness * state.elasticSlip[1];
  const double rolling = state.rollingStiffness * state.elasticRoll[2];
  const std::vector<g::Vec3> force{{load, tangential, 0.}, {-load, -tangential, 0.}};
  const std::vector<g::Vec3> torque{{0., 0., arm * tangential + rolling},
                                     {0., 0., arm * tangential - rolling}};
  for (int k = 0; k < 3; ++k) {
    accepted(g::advanceParticles(bodies, force, torque, step, k * step, s,
                                nullptr, &history));
    require(history[0].active && !history[0].sliding && !history[0].rolling,
            "balanced elastic contact must retain its sticking branch");
    nearVector(history[0].elasticSlip, state.elasticSlip, 2.e-15, 1.e-5,
               "accepted step retains sliding spring history");
    nearVector(history[0].elasticRoll, state.elasticRoll, 2.e-8, 1.e-5,
               "accepted step retains rolling spring history");
    near(history[0].rollingCap, state.rollingCap, 1.e-25, 1.e-10,
         "rolling cap remains the value at contact birth");
    near(history[0].normalLoad, load, 3.e-15, 3.e-6,
         "preloaded normal equilibrium");
    near(history[0].elasticEnergy, state.elasticEnergy, 1.e-25, 3.e-4,
         "balanced contact retains stored elastic energy");
    for (int i = 0; i < 2; ++i)
      nearVector(bodies[i].position, initial[i].position, 3.e-14, 0.,
                 "balanced tangential and rolling loads do not move the bodies");
  }
}

void rigidRotationTransportsHistory() {
  const auto s = settings();
  auto bodies = pairAtGap(s.rough.gap);
  const auto state = preloadedState(s);
  const double angle = .2;
  const g::Vec3 omega{0., 0., angle / step};
  const auto rotation = g::rotationIncrement({0., 0., angle});
  const g::Vec3 centre{8.e-6, 10.e-6, 10.e-6};
  for (auto& body : bodies) {
    body.position = g::add(centre, g::mul(rotation, g::sub(body.position, centre)));
    body.rotation = rotation;
    body.velocity = g::cross(omega, g::sub(body.position, centre));
    body.omega = omega;
  }
  const auto normal = g::mul(rotation, g::Vec3{1., 0., 0.});
  const auto out = g::roughContact(bodies[0], bodies[1], normal,
      g::scale(normal, radius), g::scale(normal, -radius), load,
      0., step, state, s.rough, true);
  nearVector(out.candidateState.elasticSlip, g::mul(rotation, state.elasticSlip),
             2.e-20, 1.e-9, "rigid frame rotation transports tangential history");
  nearVector(out.candidateState.elasticRoll, g::mul(rotation, state.elasticRoll),
             1.e-15, 1.e-9, "rigid frame rotation transports rolling history");
  nearVector(out.tangentForce, g::mul(rotation, state.tangentForce),
             2.e-19, 1.e-9, "tangential force objectivity");
  near(out.elasticEnergy, state.elasticEnergy, 1.e-29, 1.e-8,
       "common rigid rotation cannot create elastic energy");
}

void integratedRollingResponse(bool yielding) {
  const auto s = settings();
  auto bodies = pairAtGap(s.rough.gap);
  auto state = preloadedState(s);
  state.elasticSlip = {};
  state.elasticRoll = {};
  state.tangentForce = {};
  state.rollingTorque = {};
  state.elasticEnergy = 0.;
  g::PersistentContactState history{state};
  const double initialOmega = yielding ? 1.e4 : 100.;
  bodies[0].omega = {0., 0., initialOmega};
  bodies[1].omega = {0., 0., -initialOmega};
  const double initialEnergy = kineticEnergy(bodies);
  const std::vector<g::Vec3> force{{load, 0., 0.}, {-load, 0., 0.}}, zero(2);
  accepted(g::advanceParticles(bodies, force, zero, step, 0., s, nullptr, &history));
  require(history[0].active, "opposing spins retain the loaded normal contact");
  require(history[0].rolling == yielding,
          "integrated rolling response selects the expected yield branch");
  const double inertia = bodies[0].inertiaBody[2];
  const double expectedOmega = yielding
      ? initialOmega - state.rollingCap * step / inertia
      : initialOmega / (1. + 2. * state.rollingStiffness * step * step / inertia);
  near(bodies[0].omega[2], expectedOmega, 2.e-5, 2.e-6,
       "implicit opposing-spin rolling response");
  near(bodies[1].omega[2], -expectedOmega, 2.e-5, 2.e-6,
       "opposite rolling angular impulse");
  const double expectedRoll = yielding ? s.rough.rollingYieldAngle
                                      : 2. * step * expectedOmega;
  near(history[0].elasticRoll[2], expectedRoll, 2.e-8, 2.e-6,
       "committed rolling spring follows the implicit motion");
  nearVector(history[0].elasticSlip, {}, 1.e-15, 0.,
             "opposing spins create rolling without artificial tangential slip");
  require(kineticEnergy(bodies) + history[0].elasticEnergy <=
              initialEnergy * (1. + 1.e-6),
          "passive rolling cannot create kinetic plus elastic energy");
  if (yielding)
    require(history[0].plasticRollWork > 0.,
            "integrated rolling yield must record plastic dissipation");
}

void slidingAndRollingYield() {
  const auto s = settings();
  auto bodies = pairAtGap(s.rough.gap);
  const auto state = preloadedState(s);
  bodies[0].velocity = {0., 1.e-2, 0.};
  bodies[0].omega = {0., 0., 5.e4};
  const auto out = g::roughContact(bodies[0], bodies[1], {1., 0., 0.},
      {radius, 0., 0.}, {-radius, 0., 0.}, load, 0., step, state, s.rough, true);
  require(out.sliding && out.rolling, "large slip and roll must activate both yield limits");
  near(g::norm(out.tangentForce), s.rough.friction * load, 1.e-20, 1.e-10,
       "Coulomb force cap");
  near(g::norm(out.rollingTorque), state.rollingCap, 1.e-27, 1.e-10,
       "rolling moment cap");
  require(out.plasticSlipWork > 0. && out.plasticRollWork > 0.,
          "both yielded contact modes must dissipate work");
  // Internal force and common-point moments have zero total angular resultant.
  const auto angular = g::add(g::add(out.torqueI, out.torqueJ),
      g::cross(g::sub(bodies[0].position, bodies[1].position), out.forceI));
  nearVector(angular, {}, 1.e-27, 0., "contact common-point angular conservation");
}

void threeBodyContactChain() {
  const auto s = settings();
  const double distance = 2. * radius + s.rough.gap;
  std::vector<g::Body> bodies{sphere({8.e-6 - distance, 10.e-6, 10.e-6}),
                            sphere({8.e-6, 10.e-6, 10.e-6}),
                            sphere({8.e-6 + distance, 10.e-6, 10.e-6})};
  const auto initial = bodies;
  const std::vector<g::Vec3> force{{load, 0., 0.}, {}, {-load, 0., 0.}}, zero(3);
  g::PersistentContactState history;
  for (int k = 0; k < 3; ++k) {
    const auto d = g::advanceParticles(bodies, force, zero, step, k * step, s,
                                     nullptr, &history);
    accepted(d);
    require(d.contacts == 2 && history.size() == 3,
            "three-body chain must support exactly two normal contacts");
    require(history[0].active && !history[1].active && history[2].active,
            "contact identity must follow ordered particle pairs");
    near(history[0].normalLoad, load, 5.e-15, 5.e-6, "left chain reaction");
    near(history[2].normalLoad, load, 5.e-15, 5.e-6, "right chain reaction");
    for (int i = 0; i < 3; ++i)
      nearVector(bodies[i].position, initial[i].position, 5.e-14, 0.,
                 "balanced multi-contact particle position");
    for (int i = 0; i < 2; ++i)
      near(g::closestEllipsoidGap(bodies[i], bodies[i+1]).gap, s.rough.gap,
           s.contactGapTolerance, 0., "multi-contact no penetration");
    nearVector(momentum(bodies), {}, 1.e-23, 0., "chain momentum conservation");
  }
}

void decreasingNormalLoadProjectsStoredSlip() {
  const auto s = settings();
  const auto bodies = pairAtGap(s.rough.gap);
  // A contact can leave stick without a velocity increment: reducing its
  // normal reaction shrinks the Coulomb disk around its stored elastic force.
  // Check both tangent directions and values immediately below/above yield.
  for (const double sign : {-1., 1.}) {
    for (const double offset : {-1.e-10, 0., 1.e-10}) {
      auto state = preloadedState(s);
      state.elasticSlip = {0., sign * (1. + offset) *
          s.rough.friction * load / s.rough.tangentialStiffness, 0.};
      state.tangentForce = g::scale(state.elasticSlip, -s.rough.tangentialStiffness);
      state.elasticEnergy = .5 * s.rough.tangentialStiffness *
          g::dot(state.elasticSlip, state.elasticSlip) + .5 * state.rollingStiffness *
          g::dot(state.elasticRoll, state.elasticRoll);
      const double reducedLoad = .7 * load;
      const auto out = g::roughContact(bodies[0], bodies[1], {1., 0., 0.},
          {radius, 0., 0.}, {-radius, 0., 0.}, reducedLoad, 0., step,
          state, s.rough, true);
      require(out.sliding, "normal unloading must yield the stored tangential spring");
      near(out.tangentForce[1], -sign * s.rough.friction * reducedLoad,
           1.e-22, 1.e-12, "unloading tracks the current Coulomb cap");
      near(g::norm(out.candidateState.elasticSlip),
           s.rough.friction * reducedLoad / s.rough.tangentialStiffness,
           1.e-23, 1.e-12, "unloading commits the projected elastic slip");
      require(out.plasticSlipWork > 0., "shrinking Coulomb disk records plastic slip");
      nearVector(out.candidateState.elasticRoll, state.elasticRoll,
                 1.e-16, 1.e-12, "normal unloading preserves stored elastic roll");
      near(out.candidateState.rollingCap, state.rollingCap, 1.e-28, 1.e-12,
           "normal unloading preserves the adhesive birth rolling cap");
      nearVector(out.rollingTorque, state.rollingTorque, 1.e-28, 1.e-12,
                 "normal unloading does not rescale adhesive rolling torque");
    }
  }
}

void coupledFrictionLoadReversal() {
  auto s = settings();
  s.maxSubsteps = 64;
  const double distance = 2. * radius + s.rough.gap;
  std::vector<g::Body> bodies{sphere({8.e-6 - distance, 10.e-6, 10.e-6}),
                            sphere({8.e-6, 10.e-6, 10.e-6}),
                            sphere({8.e-6 + distance, 10.e-6, 10.e-6})};
  g::PersistentContactState history(3);
  auto state = preloadedState(s);
  state.elasticSlip = {0., (1. - 1.e-12) * load / s.rough.tangentialStiffness, 0.};
  state.elasticRoll = {0., 0., .9999 * s.rough.rollingYieldAngle};
  state.tangentForce = g::scale(state.elasticSlip, -s.rough.tangentialStiffness);
  state.rollingTorque = g::scale(state.elasticRoll, -state.rollingStiffness);
  state.elasticEnergy = .5 * s.rough.tangentialStiffness *
      g::dot(state.elasticSlip, state.elasticSlip) + .5 * state.rollingStiffness *
      g::dot(state.elasticRoll, state.elasticRoll);
  history[0] = state;
  history[2] = state;
  const double arm = .5 * distance;
  const double normalSchedule[] = {1., .7, .4, 1.2, 1.2, .7, .7, 1.};
  const double shearSchedule[] = {1., 1., 1., -1., -1., -1., 1., 1.};
  bool sawSliding = false, sawRolling = false, sawPositive = false, sawNegative = false;
  double slipWork = 0., rollWork = 0.;
  std::vector<g::GapCache> cache;
  for (int k = 0; k < 8; ++k) {
    const double normal = normalSchedule[k] * load;
    const double shear = shearSchedule[k] * load;
    const double appliedRoll = k == 0 ? 0. : (k < 4 ? 2. : -2.) * state.rollingCap;
    const std::vector<g::Vec3> force{{normal, shear, 0.}, {}, {-normal, -shear, 0.}};
    const std::vector<g::Vec3> torque{{0., 0., arm * shear + appliedRoll},
        {0., 0., 2. * arm * shear}, {0., 0., arm * shear - appliedRoll}};
    const auto d = g::advanceParticles(bodies, force, torque, step, k * step,
                                      s, &cache, &history);
    require(d.maxForceResidualRatio <= 1. && d.maxTorqueResidualRatio <= 1. &&
            d.contactGapViolation <= s.contactGapTolerance,
            "coupled load reversal satisfies all physical acceptance criteria");
    require(history[0].active && !history[1].active && history[2].active,
            "coupled chain retains its two loaded contacts");
    for (const std::size_t p : {std::size_t{0}, std::size_t{2}}) {
      const auto& contact = history[p];
      require(contact.normalLoad >= 0., "coupled contact reaction remains compressive");
      require(g::norm(contact.tangentForce) <=
                  s.rough.friction * contact.normalLoad + s.forceAbsoluteTolerance,
              "coupled contact obeys its current Coulomb force cap");
      require(g::norm(contact.rollingTorque) <=
                  contact.rollingCap + s.torqueAbsoluteTolerance,
              "coupled contact obeys the adhesive rolling cap");
      near(contact.rollingCap, state.rollingCap, 1.e-27, 1.e-11,
           "normal and shear cycling preserve the rolling birth cap");
      sawSliding = sawSliding || contact.sliding;
      sawRolling = sawRolling || contact.rolling;
      sawPositive = sawPositive || contact.tangentForce[1] > .1 * load;
      sawNegative = sawNegative || contact.tangentForce[1] < -.1 * load;
      slipWork += contact.plasticSlipWork;
      rollWork += contact.plasticRollWork;
    }
    nearVector(momentum(bodies), {}, 3.e-20, 0.,
               "coupled contact cycling conserves total linear momentum");
    for (int i = 0; i < 2; ++i)
      near(g::closestEllipsoidGap(bodies[i], bodies[i + 1]).gap, s.rough.gap,
           s.contactGapTolerance, 0., "coupled chain no penetration");
  }
  require(sawSliding && slipWork > 0., "load cycling exercises dissipative sliding");
  require(sawPositive && sawNegative, "shear reversal changes the tangential force sign");
  require(sawRolling && rollWork > 0., "load cycling exercises dissipative rolling");
}

void failedStepDoesNotCommit() {
  auto s = settings();
#ifdef SLURRY_USE_PETSC
  TemporaryDiagnostics files;
  if (s.solverDiagnosticsPrefix.empty()) s.solverDiagnosticsPrefix = files.prefix();
#endif
  s.maxNewtonIterations = 1;
  s.relativeTolerance = 1.e-10;
  s.forceAbsoluteTolerance = 1.e-20;
  s.torqueAbsoluteTolerance = 1.e-26;
  s.contactGapTolerance = 1.e-16;
  auto bodies = pairAtGap(s.rough.gap);
  bodies[0].velocity = {1.e-3, 3.e-3, 0.};
  bodies[1].velocity = {-1.e-3, -3.e-3, 0.};
  bodies[0].omega = {10., 20., 30.};
  bodies[1].omega = {-20., 30., 10.};
  const auto initial = bodies;
  g::PersistentContactState history{preloadedState(s)};
  const auto originalHistory = history;
  std::vector<g::GapCache> cache(1);
  cache[0].normal = {1., 0., 0.};
  cache[0].valid = true;
  const auto originalCache = cache;
  const std::vector<g::Vec3> zero(2);
  bool threw = false;
  try {
    g::advanceParticles(bodies, zero, zero, step, 0., s, &cache, &history);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  require(threw, "one-iteration nonlinear fixture must report an unconverged step");
  for (std::size_t i = 0; i < bodies.size(); ++i) {
    require(identicalBody(bodies[i], initial[i]),
            "failed solve must preserve every caller-owned body value exactly");
  }
  require(history.size() == originalHistory.size(),
          "failed solve must preserve the caller's history size");
  require(identicalContact(history[0], originalHistory[0]),
          "failed solve must preserve every caller-owned contact value exactly");
  require(cache.size() == originalCache.size() &&
          cache[0].normal == originalCache[0].normal &&
          cache[0].valid == originalCache[0].valid,
          "failed solve must preserve the accepted gap cache");
#ifdef SLURRY_USE_PETSC
  for (const auto* suffix : {"_trace.csv", "_outcomes.csv", "_failure.dat"}) {
    const auto path = s.solverDiagnosticsPrefix + suffix;
    require(std::filesystem::is_regular_file(path) && std::filesystem::file_size(path) > 0,
            std::string("final failed outer step must save ") + suffix);
  }
  const auto replay = g::particle_detail::readParticleReplay(
      s.solverDiagnosticsPrefix + "_failure.dat");
  require(replay.bodies.size() == initial.size() && replay.contacts.size() == 1 &&
          replay.cache.size() == 1, "failure replay preserves state dimensions");
  for (std::size_t i = 0; i < initial.size(); ++i)
    require(identicalBody(replay.bodies[i], initial[i]),
            "failure replay roundtrips every body value exactly");
  require(identicalContact(replay.contacts[0], originalHistory[0]),
          "failure replay roundtrips every contact-history value exactly");
  require(replay.cache[0].normal == originalCache[0].normal &&
          replay.cache[0].valid == originalCache[0].valid,
          "failure replay roundtrips the accepted geometric cache");
  require(replay.force == zero && replay.torque == zero && replay.dt == step &&
          replay.outerDt == step && replay.time == 0. && replay.outerTime == 0. &&
          replay.count == 1 && replay.substep == 0,
          "failure replay roundtrips external loads and timestep identity");
  require(replay.settings.relativeTolerance == s.relativeTolerance &&
          replay.settings.forceAbsoluteTolerance == s.forceAbsoluteTolerance &&
          replay.settings.torqueAbsoluteTolerance == s.torqueAbsoluteTolerance &&
          replay.settings.contactGapTolerance == s.contactGapTolerance &&
          replay.settings.maxNewtonIterations == s.maxNewtonIterations &&
          replay.settings.maxSubsteps == s.maxSubsteps,
          "failure replay preserves the exact failed solve criteria");
#endif
}

#ifdef SLURRY_USE_PETSC
void recoveredRetryWritesNoDiagnostics() {
  auto s = settings();
  TemporaryDiagnostics files;
  s.solverDiagnosticsPrefix = files.prefix();
  s.maxNewtonIterations = 4;
  s.maxSubsteps = 32;
  auto bodies = pairAtGap(s.rough.gap);
  bodies[0].velocity = {1.e-3, 3.e-3, 0.};
  bodies[1].velocity = {-1.e-3, -3.e-3, 0.};
  bodies[0].omega = {10., 20., 30.};
  bodies[1].omega = {-20., 30., 10.};
  g::PersistentContactState history{preloadedState(s)};
  const std::vector<g::Vec3> zero(2);
  const auto d = g::advanceParticles(bodies, zero, zero, step, 0., s, nullptr, &history);
  require(d.substeps > 1 && d.substeps <= s.maxSubsteps,
          "recovery fixture must actually retry with smaller particle timesteps");
  require(d.maxForceResidualRatio <= 1. && d.maxTorqueResidualRatio <= 1. &&
          d.contactGapViolation <= s.contactGapTolerance,
          "recovered outer step must satisfy all physical criteria");
  require(std::filesystem::is_empty(files.directory),
          "recoverable failed attempts must not write trace, outcome, or replay files");
}

void savedFrictionBranchReplay(const std::string& input) {
  // The real failure state is supplied explicitly and remains outside the test
  // source tree. Run exactly its saved substep, retaining every physical input,
  // convergence tolerance, and iteration budget without further subdivision.
  auto replay = g::particle_detail::readParticleReplay(input);
  TemporaryDiagnostics files;
  replay.settings.solverBackend = "petsc";
  replay.settings.solverDiagnosticsPrefix = files.prefix();
  std::vector<g::Body> output;
  g::ParticleStepDiagnostics d;
  std::string error;
  const bool success = g::particle_detail::dispatchImplicitStep(
      replay.bodies, replay.force, replay.torque, replay.dt, replay.time,
      replay.settings, replay.cache, replay.contacts, output, d, error,
      replay.outerTime, replay.outerDt, replay.count, replay.substep);
  require(success, "saved friction-boundary substep must converge: " + error);
  require(std::isfinite(d.maxForceResidualRatio) && d.maxForceResidualRatio <= 1. &&
          std::isfinite(d.maxTorqueResidualRatio) && d.maxTorqueResidualRatio <= 1. &&
          std::isfinite(d.contactGapViolation) &&
          d.contactGapViolation <= replay.settings.contactGapTolerance,
          "saved substep must satisfy its original force, torque, and gap tolerances");
  require(d.frictionBranchCorrections > 0 &&
          d.frictionBranchAttempts >= d.frictionBranchCorrections,
          "saved fixture must exercise an accepted friction-branch correction");
  require(d.newtonIterations <= replay.settings.maxNewtonIterations,
          "friction-branch correction must remain within the original iteration budget");
  require(output.size() == replay.bodies.size(), "replay preserves particle count");
  for (const auto& body : output)
    require(g::finite(body.position) && g::finite(body.velocity) && g::finite(body.omega),
            "accepted replay state has finite positions and velocities");
  for (const auto& contact : replay.contacts) {
    if (!contact.active) continue;
    require(contact.normalLoad >= 0., "replayed contact reaction remains compressive");
    require(g::norm(contact.tangentForce) <= replay.settings.rough.friction *
                contact.normalLoad + replay.settings.forceAbsoluteTolerance,
            "replayed contact retains the unmodified Coulomb force cap");
    require(g::norm(contact.rollingTorque) <=
                contact.rollingCap + replay.settings.torqueAbsoluteTolerance,
            "replayed contact retains the unmodified rolling torque cap");
  }
  require(std::filesystem::is_empty(files.directory),
          "successful branch recovery must not emit failure files");
  std::cout << "REPLAY REGRESSION: particles=" << output.size()
            << "; Newton=" << d.newtonIterations
            << "; branch_attempts=" << d.frictionBranchAttempts
            << "; branch_corrections=" << d.frictionBranchCorrections
            << "; force_ratio=" << d.maxForceResidualRatio
            << "; torque_ratio=" << d.maxTorqueResidualRatio
            << "; gap_violation_m=" << d.contactGapViolation << '\n';
}
#endif

void productionOblateContact() {
  auto s = settings();
  s.pair.hamaker = 9.9e-20;
  s.nearField.enabled = true;
  s.nearField.viscosity = .000890;
  s.nearField.matchingGap = 50.e-9;
  s.relativeTolerance = .01;
  s.forceAbsoluteTolerance = 1.e-13;
  s.torqueAbsoluteTolerance = 1.65e-19;
  s.contactGapTolerance = 1.e-12;
  s.maxSubsteps = 256;
  s.maxNewtonIterations = 60;
  constexpr double a = 1.65e-6, c = .20e-6, actualDt = 1.1547e-6;
  const double actualMass = (4. * std::acos(-1.) / 3.) * a * a * c * 2400. * 824.62;
  std::vector<g::Body> bodies(2);
  for (auto& body : bodies) {
    body.axes = {a, a, c};
    body.mass = actualMass;
    body.inertiaBody = {actualMass * (a*a+c*c) / 5.,
                       actualMass * (a*a+c*c) / 5., 2.*actualMass*a*a / 5.};
  }
  bodies[0].position = {8.e-6, 10.e-6, 10.e-6 - c - .5*s.rough.gap};
  bodies[1].position = {8.e-6, 10.e-6, 10.e-6 + c + .5*s.rough.gap};
  bodies[0].velocity = {1.e-4, 0., 0.};
  bodies[1].velocity = {-1.e-4, 0., 0.};
  bodies[0].omega = {0., 20., 0.};
  bodies[1].omega = {0., -20., 0.};
  const auto pair = g::evaluatePair(bodies[0], bodies[1], s.pair);
  const double birthAttraction = g::dot(pair.forceI, pair.normal);
  require(birthAttraction > 0., "production RE2 contact has an adhesive birth force");
  const std::vector<g::Vec3> zero(2);
  std::vector<g::GapCache> cache;
  g::PersistentContactState history;
  double birthRollingCap = 0.;
  double momentumAllowance = 0.;
  for (int k = 0; k < 3; ++k) {
    const auto before = bodies;
    const auto d = g::advanceParticles(bodies, zero, zero, actualDt, k*actualDt, s,
                                     &cache, &history);
    require(d.substeps == 1, "production pair should solve one actual LB interval directly");
    require(d.maxForceResidualRatio <= 1. && d.maxTorqueResidualRatio <= 1. &&
            d.contactGapViolation <= s.contactGapTolerance,
            "production oblate step must satisfy the user's physical tolerances");
    require(d.contacts == 1 && history[0].active && history[0].normalLoad > 0.,
            "adhesive oblate pair retains a compressive contact reaction");
    require(history[0].rollingCap > 0. && g::norm(history[0].elasticRoll) > 0.,
            "new oblate contact creates and advances adhesive rolling history");
    require(d.lubricationDissipation > 0.,
            "water lubrication participates in the oblate contact dynamics");
    if (k == 0) birthRollingCap = history[0].rollingCap;
    near(history[0].rollingCap, birthRollingCap, 1.e-26, 1.e-10,
         "oblate rolling cap remains fixed at contact birth");
    near(birthRollingCap, s.rough.rollingLength * birthAttraction,
         2.e-18, 1.e-3, "oblate rolling cap comes from the RE2 birth attraction");
    near(g::closestEllipsoidGap(bodies[0], bodies[1]).gap, s.rough.gap,
         s.contactGapTolerance, 0., "production oblate no penetration");
    // If R=I-F and |R| <= atol + rtol*max(|I|,|F|), then
    // |R| <= (atol + rtol*|I|)/(1-rtol). Accumulate the corresponding
    // permitted impulse instead of demanding tighter conservation than the
    // user's force tolerance. Pair forces themselves are equal and opposite.
    for (std::size_t i = 0; i < bodies.size(); ++i) {
      double inertialForce = 0.;
      for (int axis = 0; axis < 3; ++axis)
        inertialForce = std::max(inertialForce, bodies[i].mass *
            std::abs(bodies[i].velocity[axis] - before[i].velocity[axis]) / actualDt);
      momentumAllowance += actualDt * (s.forceAbsoluteTolerance +
          s.relativeTolerance * inertialForce) / (1. - s.relativeTolerance);
    }
    nearVector(momentum(bodies), {}, momentumAllowance, 0.,
               "oblate pair internal forces conserve total momentum");
  }
}

} // namespace

int main(int argc, char** argv) {
  std::string selected, replayFixture;
  std::vector<char*> runtimeArguments{argv[0]};
  for (int i = 1; i < argc; ++i) {
    const std::string argument = argv[i];
    if (argument == "--case" || argument == "--replay-fixture") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << argument << '\n';
        return 1;
      }
      if (argument == "--case") selected = argv[++i];
      else replayFixture = argv[++i];
    }
    else runtimeArguments.push_back(argv[i]);
  }
#ifndef SLURRY_USE_PETSC
  if (!replayFixture.empty()) {
    std::cerr << "--replay-fixture requires a PETSc-enabled test binary\n";
    return 1;
  }
#endif
#ifdef SLURRY_USE_PETSC
  int runtimeCount = static_cast<int>(runtimeArguments.size());
  runtimeArguments.push_back(nullptr);
  char** runtimeArgv = runtimeArguments.data();
  if (PetscInitialize(&runtimeCount, &runtimeArgv, nullptr, nullptr)) return 2;
#endif
  std::vector<std::pair<std::string, std::function<void()>>> tests{
      {"free_flight", freeFlight},
      {"forced_free_motion", forcedFreeMotion},
      {"normal_equilibrium", normalEquilibrium},
      {"approaching_contact", approachingContact},
      {"separating_contact_clears_history", separatingContactClearsHistory},
      {"sliding_rolling_history_persistence", slidingRollingHistoryPersistence},
      {"integrated_elastic_rolling", [] { integratedRollingResponse(false); }},
      {"integrated_yielded_rolling", [] { integratedRollingResponse(true); }},
      {"rigid_rotation_transports_history", rigidRotationTransportsHistory},
      {"sliding_and_rolling_yield", slidingAndRollingYield},
      {"three_body_contact_chain", threeBodyContactChain},
      {"decreasing_normal_load_projects_stored_slip", decreasingNormalLoadProjectsStoredSlip},
      {"coupled_friction_load_reversal", coupledFrictionLoadReversal},
      {"failed_step_does_not_commit", failedStepDoesNotCommit},
#ifdef SLURRY_USE_PETSC
      {"recovered_retry_writes_no_diagnostics", recoveredRetryWritesNoDiagnostics},
#endif
      {"production_oblate_contact", productionOblateContact}};
#ifdef SLURRY_USE_PETSC
  if (!replayFixture.empty())
    tests.emplace_back("saved_friction_branch_replay",
                       [replayFixture] { savedFrictionBranchReplay(replayFixture); });
#endif
  int failed = 0, ran = 0;
  for (const auto& test : tests) {
    if (!selected.empty() && test.first != selected) continue;
    ++ran;
    try {
      test.second();
      std::cout << "PASS " << test.first << '\n';
    } catch (const std::exception& error) {
      ++failed;
      std::cerr << "FAIL " << test.first << ": " << error.what() << '\n';
    }
  }
  if (!ran) { std::cerr << "Unknown --case: " << selected << '\n'; failed = 1; }
  std::cout << "RESULT " << (ran - failed) << '/' << ran << " passed\n";
#ifdef SLURRY_USE_PETSC
  if (PetscFinalize()) return 2;
#endif
  return failed ? 1 : 0;
}
