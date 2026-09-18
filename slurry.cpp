// SPDX-License-Identifier: GPL-2.0-or-later
// Dispatch only. Each case retains its original initialization and time loop.
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#ifdef SLURRY_USE_PETSC
#include <petscversion.h>
#endif
#include "slurry_registry.h"

int main(int argc, char** argv) {
  if (argc == 2 && std::string(argv[1]) == "--build-info") {
    std::cout << "{\"mpi_enabled\":";
#ifdef PARALLEL_MODE_MPI
    std::cout << "true";
#else
    std::cout << "false";
#endif
    std::cout << ",\"rough_contact\":true,\"revision\":\"slurry-petsc-contact-1\"";
#ifdef SLURRY_USE_PETSC
    std::cout << ",\"particle_solver\":\"petsc\",\"petsc_version\":\""
              << PETSC_VERSION_MAJOR << '.' << PETSC_VERSION_MINOR << '.' << PETSC_VERSION_SUBMINOR << '"';
#else
    std::cout << ",\"particle_solver\":\"legacy\"";
#endif
    std::cout << ",\"engines\":[";
    bool first = true;
    for (const auto& e : slurry_registry) {
      if (!first) std::cout << ',';
      std::cout << '"' << e.name << '"'; first = false;
    }
    std::cout << "]}\n";
    return 0;
  }
  std::string selected = std::getenv("SLURRY_ENGINE") ? std::getenv("SLURRY_ENGINE") : "";
  std::vector<char*> args{argv[0]};
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--engine") {
      if (i + 1 == argc) { std::cerr << "--engine requires a value\n"; return 2; }
      selected = argv[++i];
    } else args.push_back(argv[i]);
  }
  if (selected.empty()) {
    std::cout << "slurry --engine NAME [original solver arguments]\nAvailable:";
    for (const auto& e : slurry_registry) std::cout << ' ' << e.name;
    std::cout << "\nUse the common run_slurry_cpu.sbatch for prepared runs.\n";
    return argc == 1 || (argc == 2 && std::string(argv[1]) == "--help") ? 0 : 2;
  }
  for (const auto& e : slurry_registry) if (selected == e.name) {
    const int count = static_cast<int>(args.size());
    args.push_back(nullptr);
    return e.entry(count, args.data());
  }
  std::cerr << "Unknown engine: " << selected << ". Gr+CMC is not implemented yet.\n";
  return 2;
}
