// Replay one recorded failed particle substep without running OpenLB or any
// preceding simulation time. Physical inputs and tolerances come from the dump.
#include "particleSubsteps.h"
#ifndef SLURRY_USE_PETSC
#error "Build the replay tool with SLURRY_USE_PETSC and link PETSc"
#endif
#include <petscsys.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace g = slurry::gr_re2::graphite;

int main(int argc, char** argv) {
  std::string input, prefix, backend = "petsc";
  int maximumNewton = 0;
  std::vector<char*> petscArguments{argv[0]};
  try {
    for (int i = 1; i < argc; ++i) {
      const std::string argument = argv[i];
      if (argument == "--help" || argument == "-h") {
        std::cout << "Usage: replay_particle_step INPUT_failure.dat [--max-newton N]\n"
                     "       [--backend petsc|legacy] [--diagnostics-prefix PATH] [-gr_... PETSc options]\n"
                     "Default: replay exactly the saved substep with its original physical inputs\n"
                     "and tolerances. No additional subdivision is performed.\n";
        return 0;
      }
      if (argument == "--max-newton" || argument == "--backend" || argument == "--diagnostics-prefix") {
        if (++i >= argc) throw std::invalid_argument("Missing value for " + argument);
        if (argument == "--max-newton") {
          std::size_t used = 0;
          maximumNewton = std::stoi(argv[i], &used);
          if (used != std::string(argv[i]).size() || maximumNewton < 1)
            throw std::invalid_argument("--max-newton must be a positive integer");
        } else if (argument == "--backend") {
          backend = argv[i];
          if (backend != "petsc" && backend != "legacy")
            throw std::invalid_argument("--backend must be petsc or legacy");
        } else {
          prefix = argv[i];
        }
      } else if (input.empty() && !argument.empty() && argument[0] != '-') {
        input = argument;
      } else {
        petscArguments.push_back(argv[i]);
      }
    }
    if (input.empty()) throw std::invalid_argument("A *_failure.dat input file is required");
  } catch (const std::exception& error) {
    std::cerr << "REPLAY ARGUMENT ERROR: " << error.what() << '\n';
    return 1;
  }
  int petscArgc = static_cast<int>(petscArguments.size());
  petscArguments.push_back(nullptr);
  char** petscArgv = petscArguments.data();
  if (PetscInitialize(&petscArgc, &petscArgv, nullptr, nullptr)) return 1;
  int result = 1;
  try {
    auto replay = g::particle_detail::readParticleReplay(input);
    if (replay.bodies.empty() || !(replay.dt > 0.) || !std::isfinite(replay.dt)
        || !std::isfinite(replay.time))
      throw std::invalid_argument("Replay requires particles and a finite positive substep duration");
    const int savedMaximum = replay.settings.maxNewtonIterations;
    replay.settings.solverBackend = backend;
    replay.settings.solverDiagnosticsPrefix = prefix.empty() ? input + "_replay" : prefix;
    if (maximumNewton) replay.settings.maxNewtonIterations = maximumNewton;
    std::cout << std::setprecision(17)
              << "Replay: " << input << '\n'
              << "Backend=" << backend << "; bodies=" << replay.bodies.size()
              << "; t=" << replay.time << " s; saved subdt=" << replay.dt << " s\n"
              << "Original outer dt=" << replay.outerDt
              << " s; subdivision_count=" << replay.count
              << "; substep_index=" << replay.substep << '\n';
    if (maximumNewton && maximumNewton != savedMaximum)
      std::cout << "Iteration budget override: " << savedMaximum << " -> " << maximumNewton << '\n';
    std::vector<g::Body> output;
    g::ParticleStepDiagnostics diagnostics;
    std::string error;
    const bool success = g::particle_detail::dispatchImplicitStep(
        replay.bodies, replay.force, replay.torque, replay.dt, replay.time,
        replay.settings, replay.cache, replay.contacts, output, diagnostics, error,
        replay.outerTime, replay.outerDt, replay.count, replay.substep);
    if (success) {
      std::cout << "REPLAY ACCEPTED: Newton=" << diagnostics.newtonIterations
                << "; Krylov=" << diagnostics.krylovIterations
                << "; force_ratio=" << diagnostics.maxForceResidualRatio
                << "; torque_ratio=" << diagnostics.maxTorqueResidualRatio
                << "; friction_branch_attempts=" << diagnostics.frictionBranchAttempts
                << "; friction_branch_corrections=" << diagnostics.frictionBranchCorrections
                << "; gap_violation_m=" << diagnostics.contactGapViolation << '\n';
      result = 0;
    } else {
      std::cout << "REPLAY FAILED (one saved substep): " << error << '\n'
                << "Diagnostics prefix: " << replay.settings.solverDiagnosticsPrefix << '\n';
      result = 2;
    }
  } catch (const std::exception& error) {
    std::cerr << "REPLAY ERROR: " << error.what() << '\n';
    result = 1;
  }
  const auto finalError = PetscFinalize();
  return finalError ? 1 : result;
}
