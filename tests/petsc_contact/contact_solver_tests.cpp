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

struct ScopedPetscOptions {
  struct Saved { std::string name, value; bool present; };
  std::vector<Saved> saved;
  explicit ScopedPetscOptions(const std::vector<std::pair<std::string, std::string>>& options) {
    try {
      for (const auto& option : options) {
        char value[4096]{};
        PetscBool present = PETSC_FALSE;
        require(PetscOptionsGetString(nullptr, nullptr, option.first.c_str(), value,
                    sizeof(value), &present) == 0, "read existing PETSc test option");
        saved.push_back({option.first, value, present == PETSC_TRUE});
        require(PetscOptionsSetValue(nullptr, option.first.c_str(), option.second.c_str()) == 0,
                "set scoped PETSc test option");
      }
    } catch (...) { restore(); throw; }
  }
  ScopedPetscOptions(const ScopedPetscOptions&) = delete;
  ScopedPetscOptions& operator=(const ScopedPetscOptions&) = delete;
  void restore() {
    for (auto it = saved.rbegin(); it != saved.rend(); ++it) {
      if (it->present) PetscOptionsSetValue(nullptr, it->name.c_str(), it->value.c_str());
      else PetscOptionsClearValue(nullptr, it->name.c_str());
    }
    saved.clear();
  }
  ~ScopedPetscOptions() { restore(); }
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

g::RoughContactState zeroNormalRollingState(const g::ParticleStepSettings& s) {
  auto state = preloadedState(s);
  state.normalLoad = 0.;
  state.elasticSlip = {};
  state.tangentForce = {};
  state.elasticRoll = {0., 0., .25 * s.rough.rollingYieldAngle};
  state.rollingTorque = g::scale(state.elasticRoll, -state.rollingStiffness);
  state.elasticEnergy = .5 * state.rollingStiffness *
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

void zeroNormalLoadRetainsRollingUntilRelease() {
  const auto s = settings();
  auto bodies = pairAtGap(s.rough.gap);
  const auto state = zeroNormalRollingState(s);
  g::PersistentContactState history{state};
  const std::vector<g::Vec3> zero(2);
  const std::vector<g::Vec3> torque{g::scale(state.rollingTorque, -1.),
                                   state.rollingTorque};
  std::vector<g::GapCache> cache;
  for (int k = 0; k < 3; ++k) {
    const auto d = g::advanceParticles(bodies, zero, torque, step, k * step, s,
                                      &cache, &history);
    accepted(d);
    require(d.contactReleases == 0,
            "zero normal load alone cannot trigger a contact release update");
    require(history[0].active, "zero normal load must not erase a closed adhesive contact");
    near(history[0].normalLoad, 0., s.forceAbsoluteTolerance, 0.,
         "closed adhesive rolling equilibrium needs no normal reaction");
    nearVector(history[0].rollingTorque, state.rollingTorque, 2.e-22, 1.e-7,
               "finite adhesive rolling torque survives zero normal reaction");
    near(history[0].rollingCap, state.rollingCap, 1.e-27, 1.e-12,
         "zero normal load preserves the rolling cap at contact birth");
    near(history[0].releasedEnergy, 0., 1.e-29, 0.,
         "a closed zero-load contact does not release stored energy");
  }
  const auto beforeRelease = bodies;
  const double storedEnergy = history[0].elasticEnergy;
  // The free normal displacement is only twice the configured gap tolerance.
  // Losing the finite rolling torque at this transition must be solved in the
  // same step, rather than accepting a stale angular balance or erasing history
  // during an intermediate line-search trial.
  const double separatingForce = mass * s.contactGapTolerance / (step * step);
  const std::vector<g::Vec3> force{{-separatingForce, 0., 0.},
                                  {separatingForce, 0., 0.}};
  const auto released = g::advanceParticles(bodies, force, torque, step, 3. * step,
                                            s, &cache, &history);
  accepted(released);
  require(released.contactStateUpdates > 0 && released.contactReleases > 0,
          "release fixture exercises the normal-contact state update");
  require(!history[0].active, "accepted separation releases the adhesive contact");
  require(g::closestEllipsoidGap(bodies[0], bodies[1]).gap >
              s.rough.gap + s.contactGapTolerance,
          "release fixture actually crosses the normal-contact tolerance");
  nearVector(history[0].elasticSlip, {}, 1.e-24, 0., "release clears stored slip");
  nearVector(history[0].elasticRoll, {}, 1.e-18, 0., "release clears stored roll");
  nearVector(history[0].rollingTorque, {}, 1.e-28, 0., "release clears rolling torque");
  near(history[0].rollingCap, 0., 1.e-28, 0., "release clears the old birth cap");
  near(history[0].releasedEnergy, storedEnergy, 1.e-28, 1.e-8,
       "accepted release accounts for the stored adhesive energy exactly once");
  for (int i = 0; i < 2; ++i) {
    nearVector(bodies[i].velocity,
        g::add(beforeRelease[i].velocity, g::scale(force[i], step / mass)),
        2.e-11, 1.e-7, "released bodies follow the external normal impulse");
    nearVector(bodies[i].omega,
        g::add(beforeRelease[i].omega,
               g::scale(torque[i], step / beforeRelease[i].inertiaBody[2])),
        2.e-5, 1.e-6, "released rolling torque participates in the angular solve");
  }
  accepted(g::advanceParticles(bodies, force, torque, step, 4. * step, s,
                              &cache, &history));
  near(history[0].releasedEnergy, 0., 1.e-29, 0.,
       "a later open step cannot release the same elastic energy again");
}

void recontactCreatesFreshAdhesiveHistory() {
  const auto s = settings();
  const auto bodies = pairAtGap(s.rough.gap);
  const auto old = preloadedState(s);
  const auto released = g::roughContact(bodies[0], bodies[1], {1., 0., 0.},
      {radius, 0., 0.}, {-radius, 0., 0.}, 0., 0., step, old, s.rough, false);
  require(!released.candidateState.active, "release candidate is inactive");
  const auto stillClosed = g::roughContact(bodies[0], bodies[1], {1., 0., 0.},
      {radius, 0., 0.}, {-radius, 0., 0.}, load, 0., step, old, s.rough, true);
  nearVector(stillClosed.candidateState.elasticRoll, old.elasticRoll, 1.e-18, 1.e-12,
             "an uncommitted release trial cannot erase accepted rolling history");
  near(stillClosed.candidateState.rollingCap, old.rollingCap, 1.e-28, 1.e-12,
       "an uncommitted release trial cannot erase the accepted birth cap");
  constexpr double newBirthForce = 3.e-9;
  const auto reborn = g::roughContact(bodies[0], bodies[1], {1., 0., 0.},
      {radius, 0., 0.}, {-radius, 0., 0.}, load, newBirthForce, step,
      released.candidateState, s.rough, true);
  require(reborn.candidateState.active, "recontact creates a fresh active state");
  near(reborn.candidateState.rollingCap, s.rough.rollingLength * newBirthForce,
       1.e-28, 1.e-12, "recontact uses its new adhesive birth force");
  near(reborn.candidateState.rollingStiffness,
       s.rough.rollingLength * newBirthForce / s.rough.rollingYieldAngle,
       1.e-26, 1.e-12, "recontact creates its own rolling stiffness");
  nearVector(reborn.candidateState.elasticSlip, {}, 1.e-24, 0.,
             "recontact cannot inherit a released sliding spring");
  nearVector(reborn.candidateState.elasticRoll, {}, 1.e-18, 0.,
             "recontact cannot inherit a released rolling spring");
  nearVector(reborn.rollingTorque, {}, 1.e-28, 0.,
             "recontact at rest has no stale rolling torque");
  near(reborn.releasedEnergy, 0., 1.e-29, 0.,
       "recontact cannot count an earlier release twice");
}

void releasedAdhesiveContactStaysOpenInsideGapTolerance() {
  const auto s = settings();
  auto ellipsoid = sphere({});
  ellipsoid.axes = {radius, .5 * radius, .5 * radius};
  ellipsoid.rotation = g::rotationIncrement({0., 0., .25 * std::acos(-1.)});
  ellipsoid.inertiaBody = {.1 * mass * radius * radius,
                       .25 * mass * radius * radius,
                       .25 * mass * radius * radius};
  const g::Vec3 normal{1., 0., 0.}, centre{8.e-6, 10.e-6, 10.e-6};
  const auto shape = g::rotatedDiagonal(ellipsoid.rotation,
      {radius * radius, .25 * radius * radius, .25 * radius * radius});
  const auto shapeNormal = g::mul(shape, normal);
  const auto support = g::scale(shapeNormal, 1. / std::sqrt(g::dot(normal, shapeNormal)));
  ellipsoid.position = g::sub(g::sub(centre, support), g::scale(normal, .5 * s.rough.gap));
  std::vector<g::Body> bodies{ellipsoid,
      sphere(g::add(centre, g::scale(normal, radius + .5 * s.rough.gap)))};
  const auto initial = bodies;
  near(g::closestEllipsoidGap(bodies[0], bodies[1]).gap, s.rough.gap,
       1.e-18, 0., "aspherical chatter fixture starts at the rough-contact surface");

  // A tilted ellipsoid couples an applied rolling torque to its normal gap.
  // With the old spring retained, the normal loading alone opens by 1.25*tol.
  // Releasing that spring lets the external torque rotate the ellipsoid, so
  // the gap then returns to approximately +0.75*tol without compressive load.
  // This is an admissible open state, not a reason to resurrect the old spring.
  const double supportDerivative = -.75 * radius * radius /
      (2. * std::sqrt(.625 * radius * radius));
  const double releasedAngle = .5 * s.contactGapTolerance / supportDerivative;
  const double rollingMagnitude = -releasedAngle * ellipsoid.inertiaBody[2] / (step * step);
  auto state = zeroNormalRollingState(s);
  state.normal = normal;
  state.elasticRoll = {0., 0., -.25 * s.rough.rollingYieldAngle};
  state.rollingStiffness = rollingMagnitude / (.25 * s.rough.rollingYieldAngle);
  state.rollingCap = state.rollingStiffness * s.rough.rollingYieldAngle;
  state.rollingTorque = g::scale(state.elasticRoll, -state.rollingStiffness);
  state.elasticEnergy = .5 * state.rollingStiffness * g::dot(state.elasticRoll, state.elasticRoll);
  g::PersistentContactState history{state};
  std::vector<g::GapCache> cache(1);
  cache[0].valid = true;
  cache[0].normal = normal;
  const double openingForce = .625 * mass * s.contactGapTolerance / (step * step);
  const std::vector<g::Vec3> force{{-openingForce, 0., 0.}, {openingForce, 0., 0.}};
  const std::vector<g::Vec3> torque{g::scale(state.rollingTorque, -1.), state.rollingTorque};
  const auto d = g::advanceParticles(bodies, force, torque, step, 0., s, &cache, &history);
  accepted(d);
  require(d.contactReleases == 1 && d.contactActivations == 0,
          "an unloaded released contact must not chatter back into its old state");
  require(!history[0].active, "released contact stays open within the existing gap tolerance");
  const double slack = g::closestEllipsoidGap(bodies[0], bodies[1]).gap - s.rough.gap;
  require(slack > .5 * s.contactGapTolerance && slack < .95 * s.contactGapTolerance,
          "rolling-coupled endpoint lies inside the positive gap tolerance band");
  near(history[0].normalLoad, 0., 0., 0., "open chatter endpoint has zero normal reaction");
  nearVector(history[0].elasticSlip, {}, 0., 0., "released slip cannot be resurrected");
  nearVector(history[0].elasticRoll, {}, 0., 0., "released rolling spring cannot be resurrected");
  near(history[0].rollingCap, 0., 0., 0., "released birth cap cannot be resurrected");
  near(history[0].releasedEnergy, state.elasticEnergy, 1.e-29, 1.e-10,
       "chatter prevention still accounts for released elastic energy once");
  for (int i = 0; i < 2; ++i) {
    nearVector(bodies[i].velocity, g::scale(force[i], step / mass), 2.e-11, 1.e-7,
               "open chatter endpoint satisfies the external linear impulse");
    near(bodies[i].omega[2], step * torque[i][2] / initial[i].inertiaBody[2],
         2.e-5, 1.e-6, "open chatter endpoint satisfies the external angular impulse");
  }
  // A later resolved compression must still reactivate this genuinely released
  // pair. This fixture has no current pair attraction, so its newly born rolling
  // cap is zero: the historical adhesive spring must never reappear.
  const std::vector<g::Vec3> compression{{load, 0., 0.}, {-load, 0., 0.}}, zero(2);
  const auto recontact = g::advanceParticles(bodies, compression, zero, step, step,
                                            s, &cache, &history);
  accepted(recontact);
  require(recontact.contactActivations > 0 && history[0].active &&
              history[0].normalLoad > 1000. * s.forceAbsoluteTolerance,
          "a later resolved compressive reaction reactivates the released pair");
  near(history[0].rollingCap, 0., 0., 0., "recontact uses its new zero-adhesion birth cap");
  nearVector(history[0].elasticRoll, {}, 0., 0., "resolved recontact cannot revive old rolling strain");
  nearVector(history[0].rollingTorque, {}, 0., 0., "resolved recontact cannot revive old rolling torque");
  near(history[0].releasedEnergy, 0., 0., 0., "resolved recontact cannot release old energy twice");
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

void threeBodyContactChain(bool staticReaction) {
  auto s = settings();
  if (staticReaction) {
    // The unchanged static-reaction assertion allows 1e-14 N.  Over successive
    // backward-Euler steps, position errors bounded by the gap resolution can
    // contribute up to 4*m*gapTolerance/dt^2 to the reaction (a second position
    // difference).  A 1e-16 m gap tolerance bounds that contribution by 8e-15 N.
    // The original gap tolerance is covered separately by discrete impulse.
    s.contactGapTolerance = 1.e-16;
  }
  const double distance = 2. * radius + s.rough.gap;
  std::vector<g::Body> bodies{sphere({8.e-6 - distance, 10.e-6, 10.e-6}),
                            sphere({8.e-6, 10.e-6, 10.e-6}),
                            sphere({8.e-6 + distance, 10.e-6, 10.e-6})};
  const auto initial = bodies;
  const std::vector<g::Vec3> force{{load, 0., 0.}, {}, {-load, 0., 0.}}, zero(3);
  g::PersistentContactState history;
  for (int k = 0; k < 3; ++k) {
    const auto previous = bodies;
    const auto previousHistory = history;
    const auto d = g::advanceParticles(bodies, force, zero, step, k * step, s,
                                     nullptr, &history);
    accepted(d);
    require(d.contacts == 2 && history.size() == 3,
            "three-body chain must support exactly two normal contacts");
    require(history[0].active && !history[1].active && history[2].active,
            "contact identity must follow ordered particle pairs");
    if (staticReaction) {
      near(history[0].normalLoad, load, 5.e-15, 5.e-6, "left chain reaction");
      near(history[2].normalLoad, load, 5.e-15, 5.e-6, "right chain reaction");
    }
    auto netForce = force;
    for (int i = 0; i < 2; ++i) {
      const int p = 2 * i;
      const auto& contact = history[p];
      require(contact.normalLoad > 0., "loaded chain contacts remain compressive");
      const auto pairForce = g::add(g::scale(contact.normal, -contact.normalLoad),
                                    contact.tangentForce);
      netForce[i] = g::add(netForce[i], pairForce);
      netForce[i+1] = g::sub(netForce[i+1], pairForce);
      nearVector(contact.elasticSlip, {}, 1.e-20, 0.,
                 "normal chain loading creates no tangential spring");
      nearVector(contact.elasticRoll, {}, 1.e-14, 0.,
                 "normal chain loading creates no rolling spring");
      if (k > 0)
        near(contact.rollingCap, previousHistory[p].rollingCap, 0., 0.,
             "persistent chain contact retains its birth rolling cap");
    }
    for (int i = 0; i < 3; ++i) {
      const auto inertialForce = g::scale(g::sub(bodies[i].velocity,
          previous[i].velocity), bodies[i].mass / step);
      double forceReference = 0.;
      for (int axis = 0; axis < 3; ++axis)
        forceReference = std::max(forceReference,
            std::max(std::abs(netForce[i][axis]), std::abs(inertialForce[axis])));
      nearVector(inertialForce, netForce[i], s.forceAbsoluteTolerance +
          s.relativeTolerance * forceReference, 0., "chain discrete impulse balance");
    }
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
void nonlinearWorkAccounting(bool failLinearSolve) {
  auto s = settings();
  s.maxNewtonIterations = failLinearSolve ? 60 : 2;
  TemporaryDiagnostics files;
  s.solverDiagnosticsPrefix = files.prefix();
  std::vector<std::pair<std::string, std::string>> overrides{
      {"-gr_npc_snes_max_it", "999"}, {"-gr_snes_max_it", "999"},
      {"-gr_snes_norm_schedule", "none"}, {"-gr_npc_snes_norm_schedule", "none"}};
  if (failLinearSolve) {
    overrides.emplace_back("-gr_npc_ksp_max_it", "1");
    overrides.emplace_back("-gr_npc_ksp_rtol", "1e-30");
    overrides.emplace_back("-gr_npc_ksp_atol", "1e-50");
  }
  ScopedPetscOptions options(overrides);
  auto bodies = pairAtGap(s.rough.gap);
  bodies[0].velocity = {1.e-3, 3.e-3, 0.};
  bodies[1].velocity = {-1.e-3, -3.e-3, 0.};
  bodies[0].omega = {10., 20., 30.};
  bodies[1].omega = {-20., 30., 10.};
  g::PersistentContactState history{preloadedState(s)};
  const auto oldHistory = history;
  std::vector<g::GapCache> cache(1);
  cache[0].normal = {1., 0., 0.};
  cache[0].valid = true;
  const auto oldCache = cache;
  const std::vector<g::Vec3> zero(2);
  std::vector<g::Body> output;
  g::ParticleStepDiagnostics d;
  std::string error;
  // Capture the real failed solve's structured diagnostic trace rather than
  // guessing work from the outer NGMRES iteration count or parsing prose.
  g::particle_detail::ParticleDiagnosticScope diagnostics(s);
  const bool success = g::particle_detail::dispatchImplicitStep(
      bodies, zero, zero, step, 0., s, cache, history, output, d, error,
      0., step, 1, 0);
  require(!success, "forced solver-budget or Krylov failure must not be accepted");
  require(!diagnostics.attempts.empty() && !diagnostics.attempts.back().trace.empty(),
          "failed nonlinear attempt must retain an actual work-count trace");
  const auto& trace = diagnostics.attempts.back().trace;
  int attemptedNewton = 0;
  bool failedKrylovRecorded = false;
  for (const auto& row : trace) {
    attemptedNewton = std::max(attemptedNewton, row.iteration);
    failedKrylovRecorded = failedKrylovRecorded || row.kspReason < 0;
    require(row.iteration <= s.maxNewtonIterations,
            "PETSc option overrides cannot bypass the physical solver's Newton budget");
  }
  if (failLinearSolve) {
    require(attemptedNewton == 1 && failedKrylovRecorded,
            "an unsuccessful Krylov solve still consumes and records one Newton attempt");
  } else {
    require(attemptedNewton == s.maxNewtonIterations,
            "nested nonlinear options cannot turn two allowed Newton solves into 999");
  }
  require(history.size() == oldHistory.size() && identicalContact(history[0], oldHistory[0]),
          "an exhausted nested solve cannot commit trial contact history");
  require(cache.size() == oldCache.size() && cache[0].valid == oldCache[0].valid &&
              cache[0].normal == oldCache[0].normal,
          "an exhausted nested solve cannot commit a trial geometric cache");
  require(output.empty(), "an unsuccessful nested solve cannot publish particle output");
}

void failedContactReleaseDoesNotCommit() {
  auto s = settings();
  // The normal-state update and its new angular balance exceed this budget.
  // Rejecting that outer step must undo every attempted contact release.
  s.maxNewtonIterations = 2;
  TemporaryDiagnostics files;
  s.solverDiagnosticsPrefix = files.prefix();
  auto bodies = pairAtGap(s.rough.gap);
  const auto initial = bodies;
  const auto state = zeroNormalRollingState(s);
  g::PersistentContactState history{state};
  std::vector<g::GapCache> cache(1);
  cache[0].valid = true;
  cache[0].normal = {1., 0., 0.};
  const auto oldCache = cache;
  const double separatingForce = mass * s.contactGapTolerance / (step * step);
  const std::vector<g::Vec3> force{{-separatingForce, 0., 0.},
                                  {separatingForce, 0., 0.}};
  const std::vector<g::Vec3> torque{g::scale(state.rollingTorque, -1.),
                                   state.rollingTorque};
  bool threw = false;
  try {
    g::advanceParticles(bodies, force, torque, step, 0., s, &cache, &history);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  require(threw, "limited-budget contact release must report its unfinished angular solve");
  require(identicalBody(bodies[0], initial[0]) && identicalBody(bodies[1], initial[1]),
          "rejected contact release preserves the accepted bodies exactly");
  require(history.size() == 1 && identicalContact(history[0], state),
          "rejected contact release cannot erase the accepted springs or birth cap");
  require(cache.size() == 1 && cache[0].normal == oldCache[0].normal &&
              cache[0].valid == oldCache[0].valid,
          "rejected contact release preserves the accepted geometric cache");
  const auto replay = g::particle_detail::readParticleReplay(
      files.prefix() + "_failure.dat");
  require(replay.contacts.size() == 1 && identicalContact(replay.contacts[0], state),
          "failed release replay stores the accepted history, not a trial release");
}

void recoveredRetryWritesNoDiagnostics() {
  auto s = settings();
  TemporaryDiagnostics files;
  s.solverDiagnosticsPrefix = files.prefix();
  // Keep the original four-Newton physical fixture. Acceleration can solve it
  // directly, so a second independent call injects a three-Newton budget, one
  // below the full interval's measured cost, to exercise actual retry/no-I/O.
  s.rough.enabled = false;
  s.maxNewtonIterations = 4;
  s.maxSubsteps = 32;
  auto rotor = sphere({8.e-6, 10.e-6, 10.e-6});
  rotor.axes = {radius, .7 * radius, .4 * radius};
  rotor.inertiaBody = {.13 * mass * radius * radius,
                      .232 * mass * radius * radius,
                      .298 * mass * radius * radius};
  rotor.omega = {1.5e5, 3.e5, -2.25e5};
  const auto initialMomentum = g::particle_detail::worldMomentum(rotor);
  const std::vector<g::Vec3> zero(1);
  // For zero applied torque, each accepted substep has
  // |delta L| <= subdt * torqueAbsoluteTolerance / (1-relativeTolerance).
  // The total bound is independent of the number of accepted substeps.
  const double angularAllowance = step * s.torqueAbsoluteTolerance /
      (1. - s.relativeTolerance) + 32. * std::numeric_limits<double>::epsilon() *
      g::norm(initialMomentum);
  for (const int budget : {4, 3}) {
    s.maxNewtonIterations = budget;
    std::vector<g::Body> bodies{rotor};
    const auto d = g::advanceParticles(bodies, zero, zero, step, 0., s);
    require(d.substeps >= 1 && d.substeps <= s.maxSubsteps,
            "rotor must converge within its unchanged particle subdivision limit");
    if (budget == 3)
      require(d.substeps > 1, "injected budget must actually trigger a recovered retry");
    require(d.maxForceResidualRatio <= 1. && d.maxTorqueResidualRatio <= 1. &&
            d.contactGapViolation <= s.contactGapTolerance,
            "rotor must satisfy all unchanged physical acceptance criteria");
    nearVector(g::particle_detail::worldMomentum(bodies[0]), initialMomentum,
               angularAllowance, 0., "free rotor preserves angular momentum through retries");
    nearVector(bodies[0].position, rotor.position, 0., 0.,
               "rotational retries cannot move a force-free particle centre");
    require(std::filesystem::is_empty(files.directory),
            "recoverable failed attempts must not write trace, outcome, or replay files");
    std::cout << "ROTATION VALIDATION: budget=" << budget << "; substeps=" << d.substeps
              << "; accepted_Newton=" << d.newtonIterations
              << "; accepted_Krylov=" << d.krylovIterations << '\n';
  }
}

void preloadedCollisionWithinProductionBudget() {
  auto s = settings();
  s.maxNewtonIterations = 60;
  s.maxSubsteps = 32;
  auto bodies = pairAtGap(s.rough.gap);
  bodies[0].velocity = {1.e-3, 3.e-3, 0.};
  bodies[1].velocity = {-1.e-3, -3.e-3, 0.};
  bodies[0].omega = {10., 20., 30.};
  bodies[1].omega = {-20., 30., 10.};
  g::PersistentContactState history{preloadedState(s)};
  const auto initial = bodies;
  const auto initialContact = history[0];
  const double initialEnergy = kineticEnergy(bodies) + history[0].elasticEnergy;
  const std::vector<g::Vec3> zero(2);
  const auto d = g::advanceParticles(bodies, zero, zero, step, 0., s, nullptr, &history);
  require(d.maxForceResidualRatio <= 1. && d.maxTorqueResidualRatio <= 1. &&
          d.contactGapViolation <= s.contactGapTolerance,
          "preloaded collision satisfies unchanged physical acceptance criteria");
  require(history.size() == 1 && history[0].active && history[0].normalLoad >= 0.,
          "preloaded collision retains its compressive rough contact");
  near(history[0].rollingCap, initialContact.rollingCap, 1.e-27, 1.e-12,
       "preloaded collision preserves the adhesive birth rolling cap");
  require(g::norm(history[0].tangentForce) <=
              s.rough.friction * history[0].normalLoad + s.forceAbsoluteTolerance &&
          g::norm(history[0].rollingTorque) <=
              history[0].rollingCap + s.torqueAbsoluteTolerance,
          "preloaded collision retains Coulomb and rolling limits");
  nearVector(momentum(bodies), momentum(initial), 1.e-22, 0.,
             "preloaded collision conserves total linear momentum");
  require(kineticEnergy(bodies) + history[0].elasticEnergy <= initialEnergy * (1. + 1.e-6),
          "passive preloaded collision cannot create kinetic plus elastic energy");
}

void savedParticleReplay(const std::string& input, bool requireFrictionCorrection) {
  // The real failure state is supplied explicitly and remains outside the test
  // source tree. Run exactly its saved substep, retaining every physical input,
  // convergence tolerance, and iteration budget without further subdivision.
  auto replay = g::particle_detail::readParticleReplay(input);
  const auto startingContacts = replay.contacts;
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
  require(success, "saved particle substep must converge: " + error);
  require(std::isfinite(d.maxForceResidualRatio) && d.maxForceResidualRatio <= 1. &&
          std::isfinite(d.maxTorqueResidualRatio) && d.maxTorqueResidualRatio <= 1. &&
          std::isfinite(d.contactGapViolation) &&
          d.contactGapViolation <= replay.settings.contactGapTolerance,
          "saved substep must satisfy its original force, torque, and gap tolerances");
  require(d.frictionBranchAttempts >= d.frictionBranchCorrections,
          "accepted friction corrections cannot exceed attempted corrections");
  if (requireFrictionCorrection)
    require(d.frictionBranchCorrections > 0,
            "saved fixture must exercise an accepted friction-branch correction");
  require(d.newtonIterations <= replay.settings.maxNewtonIterations,
          "friction-branch correction must remain within the original iteration budget");
  require(output.size() == replay.bodies.size(), "replay preserves particle count");
  for (const auto& body : output)
    require(g::finite(body.position) && g::finite(body.velocity) && g::finite(body.omega),
            "accepted replay state has finite positions and velocities");
  for (std::size_t p = 0; p < replay.contacts.size(); ++p) {
    const auto& contact = replay.contacts[p];
    const auto& old = startingContacts[p];
    if (!contact.active) {
      near(contact.normalLoad, 0., 0., 0., "open replayed contact has no stored reaction");
      nearVector(contact.elasticSlip, {}, 0., 0., "open replayed contact has no slip history");
      nearVector(contact.elasticRoll, {}, 0., 0., "open replayed contact has no rolling history");
      nearVector(contact.rollingTorque, {}, 0., 0., "open replayed contact has no rolling torque");
      near(contact.rollingCap, 0., 0., 0., "open replayed contact has no stale birth cap");
      near(contact.releasedEnergy, old.active ? old.elasticEnergy : 0.,
           1.e-29, 1.e-12, "replayed release accounts for the accepted history once");
      continue;
    }
    if (old.active)
      near(contact.rollingCap, old.rollingCap, 1.e-28, 1.e-12,
           "tentative normal-mask changes cannot reset a retained contact birth cap");
    near(contact.releasedEnergy, 0., 0., 0.,
         "retained replayed contact cannot commit a rejected release");
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
            << "; contact_state_updates=" << d.contactStateUpdates
            << "; contact_activations=" << d.contactActivations
            << "; contact_releases=" << d.contactReleases
            << "; force_ratio=" << d.maxForceResidualRatio
            << "; torque_ratio=" << d.maxTorqueResidualRatio
            << "; gap_violation_m=" << d.contactGapViolation << '\n';
}
#endif

void productionOblateContact(bool localAdhesion = false) {
  auto s = settings();
  s.pair.hamaker = 9.9e-20;
  if (localAdhesion) {
    s.pair.sigma = 4.197e-10;
    s.pair.roughnessGap = s.rough.gap;
    s.pair.localGap = 3.e-10;
    s.pair.localGapFraction = 1.;
    s.pair.localSwitchExcessGap = 2.e-9;
    s.pair.localCutoffExcessGap = 10.e-9;
    s.rough.tangentialStiffness = 80.;
  }
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
  if (localAdhesion)
    near(birthAttraction, 936.214833e-9, 2.e-12, 1.e-5,
         "local graphite adhesion reaches the approved FF force scale");
  const std::vector<g::Vec3> zero(2);
  std::vector<g::GapCache> cache;
  g::PersistentContactState history;
  double birthRollingCap = 0.;
  double momentumAllowance = 0.;
  for (int k = 0; k < 3; ++k) {
    const auto before = bodies;
    const auto d = g::advanceParticles(bodies, zero, zero, actualDt, k*actualDt, s,
                                     &cache, &history);
    if (!localAdhesion)
      require(d.substeps == 1, "production pair should solve one actual LB interval directly");
    else
      require(d.substeps > 0 && d.substeps <= s.maxSubsteps,
              "strong local adhesion stays within the production substep budget");
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
  bool requireFrictionCorrection = false;
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
    } else if (argument == "--require-friction-correction") {
      requireFrictionCorrection = true;
    }
    else runtimeArguments.push_back(argv[i]);
  }
  // Preserve the earlier case name as a physical replay alias. A different
  // globalization may solve the same state without invoking a branch predictor;
  // only the explicit --require-friction-correction flag requires that path.
  if (selected == "saved_friction_branch_replay") {
    selected = "saved_particle_replay";
  }
  if (requireFrictionCorrection && replayFixture.empty()) {
    std::cerr << "--require-friction-correction requires --replay-fixture FILE\n";
    return 1;
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
      {"zero_normal_load_retains_rolling_until_release", zeroNormalLoadRetainsRollingUntilRelease},
      {"recontact_creates_fresh_adhesive_history", recontactCreatesFreshAdhesiveHistory},
      {"released_adhesive_contact_stays_open_inside_gap_tolerance", releasedAdhesiveContactStaysOpenInsideGapTolerance},
      {"sliding_rolling_history_persistence", slidingRollingHistoryPersistence},
      {"integrated_elastic_rolling", [] { integratedRollingResponse(false); }},
      {"integrated_yielded_rolling", [] { integratedRollingResponse(true); }},
      {"rigid_rotation_transports_history", rigidRotationTransportsHistory},
      {"sliding_and_rolling_yield", slidingAndRollingYield},
      {"three_body_contact_chain", [] { threeBodyContactChain(true); }},
      {"three_body_contact_chain_discrete_impulse", [] { threeBodyContactChain(false); }},
      {"decreasing_normal_load_projects_stored_slip", decreasingNormalLoadProjectsStoredSlip},
      {"coupled_friction_load_reversal", coupledFrictionLoadReversal},
      {"failed_step_does_not_commit", failedStepDoesNotCommit},
#ifdef SLURRY_USE_PETSC
      {"npc_iteration_override_respects_shared_budget", [] { nonlinearWorkAccounting(false); }},
      {"failed_krylov_attempt_is_counted", [] { nonlinearWorkAccounting(true); }},
      {"failed_contact_release_does_not_commit", failedContactReleaseDoesNotCommit},
      {"recovered_retry_writes_no_diagnostics", recoveredRetryWritesNoDiagnostics},
      {"preloaded_collision_within_production_budget", preloadedCollisionWithinProductionBudget},
#endif
      {"production_oblate_contact", [] { productionOblateContact(); }},
      {"production_oblate_local_adhesion", [] { productionOblateContact(true); }}};
#ifdef SLURRY_USE_PETSC
  if (!replayFixture.empty())
    tests.emplace_back("saved_particle_replay", [replayFixture, requireFrictionCorrection] {
      savedParticleReplay(replayFixture, requireFrictionCorrection);
    });
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
