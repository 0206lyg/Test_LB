/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef SLURRY_GR_RE2_GRAPHITE_QUICK_CONFIG_H
#define SLURRY_GR_RE2_GRAPHITE_QUICK_CONFIG_H
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

// SLURRY SCOPE BEGIN
namespace slurry { namespace gr_re2 {

struct Config {
  double shear_rate=100,box_x=1e-5,box_y=1e-5,box_z=1e-5,dx=5e-8;
  double diameter=3.3e-6,thickness=4e-7,rho_particle=2400,rho_fluid=997,dynamic_viscosity=.000890;
  double nu_lattice=.5,target_mach=.1,time_step_s=0,epsilon_cells=.5;
  double hamaker=.99e-19,sigma_lj=3e-9,switch_gap=400e-9,cutoff_gap=500e-9;
  // Old configurations retain the uncorrected RE2 law unless explicitly enabled.
  double local_gap=.3e-9,local_gap_fraction=0.;
  double local_switch_excess_gap=2e-9,local_cutoff_excess_gap=10e-9;
  bool surface_adhesion=false;
  double adhesion_work=.0219,adhesion_range=6.7e-10;
  double curvature_switch_gap=5e-9,curvature_cutoff_gap=2e-8;
  double end_strain=10,particle_tolerance=1e-4,lubrication_cutoff_cells=1.;
  bool rough_contact_enabled=true;
  double roughness_gap=2e-9,sliding_friction=.5,tangential_stiffness=9.;
  double rolling_length=100e-9,rolling_yield_angle=.01;
  double particle_force_absolute_tolerance=1e-15,particle_torque_absolute_tolerance=1.65e-21;
  double contact_gap_tolerance=1e-12;
  std::uint64_t max_steps=0,sample_every=20,vtk_every=0,checkpoint_every=0,checkpoint_keep=2;
  double checkpoint_seconds=21000.; // 350 minutes of wall time; step-based saves disabled.
  unsigned particle_max_substeps=32,particle_max_iterations=20,particle_max_krylov_iterations=120;
  std::string particle_solver="petsc";
  bool solver_diagnostics=true;
  std::string output_dir="run",particles_csv,restart_dir;
};
inline std::string strip(const std::string& s) {
  const auto a=s.find_first_not_of(" \t\r\n");
  return a==std::string::npos?"":s.substr(a,s.find_last_not_of(" \t\r\n")-a+1);
}
inline Config parseConfig(int argc,char**argv) {
  std::string file; std::uint64_t overrideSteps=0;
  for(int i=1;i<argc;++i){
    const std::string a=argv[i];
    if(i+1>=argc) throw std::runtime_error("Expected --config FILE [--max-steps N]");
    if(a=="--config")file=argv[++i];
    else if(a=="--max-steps")overrideSteps=std::stoull(argv[++i]);
    else throw std::runtime_error("Unknown argument: "+a);
  }
  if(file.empty())throw std::runtime_error("--config is required");
  std::ifstream in(file);if(!in)throw std::runtime_error("Cannot read "+file);
  Config c; std::string line;bool explicitLocal=false,explicitSurface=false;
  while(std::getline(in,line)){
    line=strip(line.substr(0,line.find('#')));if(line.empty())continue;
    const auto eq=line.find('=');if(eq==std::string::npos)throw std::runtime_error("Expected key=value: "+line);
    const auto key=strip(line.substr(0,eq)),v=strip(line.substr(eq+1));
    if(key=="local_gap"||key=="local_gap_fraction"||key=="local_switch_excess_gap"||key=="local_cutoff_excess_gap")
      explicitLocal=true;
    if(key=="adhesion_work"||key=="adhesion_range"||key=="curvature_switch_gap"||key=="curvature_cutoff_gap")
      explicitSurface=true;
#define REAL(k) if(key==#k){c.k=std::stod(v);if(!std::isfinite(c.k))throw std::runtime_error("Nonfinite config: " #k);continue;}
    REAL(shear_rate) REAL(box_x) REAL(box_y) REAL(box_z) REAL(dx) REAL(diameter) REAL(thickness)
    REAL(rho_particle) REAL(rho_fluid) REAL(dynamic_viscosity) REAL(nu_lattice) REAL(target_mach)
    REAL(time_step_s) REAL(epsilon_cells) REAL(hamaker) REAL(sigma_lj) REAL(switch_gap) REAL(cutoff_gap)
    REAL(local_gap) REAL(local_gap_fraction) REAL(local_switch_excess_gap) REAL(local_cutoff_excess_gap)
    REAL(adhesion_work) REAL(adhesion_range) REAL(curvature_switch_gap) REAL(curvature_cutoff_gap)
    REAL(checkpoint_seconds)
    REAL(end_strain) REAL(particle_tolerance) REAL(lubrication_cutoff_cells)
    REAL(roughness_gap) REAL(sliding_friction) REAL(tangential_stiffness)
    REAL(rolling_length) REAL(rolling_yield_angle)
    REAL(particle_force_absolute_tolerance) REAL(particle_torque_absolute_tolerance) REAL(contact_gap_tolerance)
#undef REAL
    if(key=="particle_solver"){
      if(v!="petsc"&&v!="legacy")throw std::runtime_error("particle_solver must be petsc or legacy");
      c.particle_solver=v;continue;
    }
    if(key=="solver_diagnostics"){
      if(v!="0"&&v!="1")throw std::runtime_error("solver_diagnostics must be 0 or 1");
      c.solver_diagnostics=v=="1";continue;
    }
    if(key=="surface_adhesion"){
      if(v!="0"&&v!="1")throw std::runtime_error("surface_adhesion must be 0 or 1");
      c.surface_adhesion=v=="1";continue;
    }
    if(key=="rough_contact_enabled"){
      if(v!="0"&&v!="1")throw std::runtime_error("rough_contact_enabled must be 0 or 1");
      c.rough_contact_enabled=v=="1";continue;
    }
#define INTEGER(k) if(key==#k){if(v.empty()||v[0]=='-')throw std::runtime_error("Negative config: " #k);c.k=std::stoull(v);continue;}
    INTEGER(max_steps) INTEGER(sample_every) INTEGER(vtk_every) INTEGER(checkpoint_every) INTEGER(checkpoint_keep) INTEGER(particle_max_substeps) INTEGER(particle_max_iterations) INTEGER(particle_max_krylov_iterations)
#undef INTEGER
    if(key=="output_dir"){c.output_dir=v;continue;}
    if(key=="particles_csv"){c.particles_csv=v;continue;}
    if(key=="restart_dir"){c.restart_dir=v;continue;}
    throw std::runtime_error("Unknown configuration key: "+key);
  }
  if(overrideSteps)c.max_steps=overrideSteps;
  if(c.checkpoint_seconds<0.||c.checkpoint_keep<1)
    throw std::runtime_error("Require checkpoint_seconds >= 0 and checkpoint_keep >= 1");
  if(!(c.shear_rate>0&&c.dx>0&&c.box_x>0&&c.box_y>0&&c.box_z>0&&c.diameter>0&&c.thickness>0
    &&c.thickness<=c.diameter&&c.rho_particle>0&&c.rho_fluid>0&&c.dynamic_viscosity>0&&c.nu_lattice>0
    &&c.target_mach>0&&c.time_step_s>=0&&c.epsilon_cells>0&&c.hamaker>=0&&c.sigma_lj>0
    &&c.switch_gap>c.sigma_lj&&c.cutoff_gap>c.switch_gap&&c.end_strain>0&&c.sample_every>0
    &&c.particle_max_substeps>0&&c.particle_max_iterations>0&&c.particle_max_krylov_iterations>0&&c.particle_tolerance>0&&c.particle_tolerance<1
    &&c.lubrication_cutoff_cells>=0&&c.roughness_gap>0&&c.roughness_gap<c.cutoff_gap
    &&c.sliding_friction>=0&&c.tangential_stiffness>0&&c.rolling_length>=0&&c.rolling_yield_angle>0
    &&c.particle_force_absolute_tolerance>0&&c.particle_torque_absolute_tolerance>0
    &&c.contact_gap_tolerance>0&&c.contact_gap_tolerance<c.roughness_gap))
    throw std::runtime_error("Invalid geometric/material/solver configuration");
  if(!(c.local_gap>0&&c.local_gap_fraction>=0&&c.local_gap_fraction<=1
       &&c.local_switch_excess_gap>=0&&c.local_cutoff_excess_gap>c.local_switch_excess_gap))
    throw std::runtime_error("Require local_gap > 0, 0 <= local_gap_fraction <= 1, and 0 <= local_switch_excess_gap < local_cutoff_excess_gap");
  if(c.surface_adhesion){
    if(explicitLocal)
      throw std::runtime_error("Surface adhesion replaces local-gap adhesion; remove all local_* inputs");
    if(!c.rough_contact_enabled)
      throw std::runtime_error("Surface adhesion requires rough_contact_enabled=1");
    if(!(c.adhesion_work>0.&&c.adhesion_range>0.
       &&c.roughness_gap+c.adhesion_range<=c.curvature_switch_gap
       &&c.curvature_switch_gap<c.curvature_cutoff_gap&&c.curvature_cutoff_gap<=c.switch_gap))
      throw std::runtime_error("Require positive adhesion_work/adhesion_range and roughness_gap + adhesion_range <= curvature_switch_gap < curvature_cutoff_gap <= switch_gap");
    const double wbg=c.hamaker/(12.*std::acos(-1.)*c.roughness_gap*c.roughness_gap)
        *(1.-std::pow(c.sigma_lj/c.roughness_gap,6)/30.);
    if(!(std::isfinite(wbg)&&wbg>=0.&&c.adhesion_work>=wbg))
      throw std::runtime_error("Surface adhesion requires nonnegative background work and adhesion_work >= background work at roughness_gap");
  } else if(explicitSurface) {
    throw std::runtime_error("Surface adhesion parameters require surface_adhesion=1");
  }
  if(c.local_gap_fraction>0){
    if(!c.rough_contact_enabled)
      throw std::runtime_error("Local-gap adhesion requires rough_contact_enabled=1");
    if(c.local_gap>c.roughness_gap)
      throw std::runtime_error("Local-gap adhesion requires local_gap <= roughness_gap");
    if(c.roughness_gap+c.local_cutoff_excess_gap>c.switch_gap)
      throw std::runtime_error("Local-gap adhesion requires roughness_gap + local_cutoff_excess_gap <= switch_gap");
  }
  for(double length:{c.box_x,c.box_y,c.box_z})
    if(length/c.dx<8||std::abs(length/c.dx-std::round(length/c.dx))>1e-7)
      throw std::runtime_error("Box lengths must be integer multiples of dx and at least eight cells");
  if(c.diameter+c.cutoff_gap>=.5*std::min({c.box_x,c.box_y,c.box_z}))
    throw std::runtime_error("diameter + cutoff_gap must be less than half the shortest box length");
  c.output_dir=std::filesystem::absolute(c.output_dir).string();
  if(!c.restart_dir.empty())c.restart_dir=std::filesystem::absolute(c.restart_dir).string();
  return c;
}
struct Units {
  double rePhysical,alpha,rhoFluid,rhoParticle,nu,dt,omega,reNumeric,stNumeric,mach;
  explicit Units(const Config& c){
    dt=c.time_step_s>0?c.time_step_s:c.target_mach*c.dx/(std::sqrt(3.)*.5*c.shear_rate*c.box_y);
    rhoFluid=c.dynamic_viscosity*dt/(c.nu_lattice*c.dx*c.dx);
    alpha=rhoFluid/c.rho_fluid;rhoParticle=alpha*c.rho_particle;
    nu=c.dynamic_viscosity/rhoFluid;omega=1/(.5+3*c.nu_lattice);
    rePhysical=c.rho_fluid*c.shear_rate*c.diameter*c.diameter/c.dynamic_viscosity;
    reNumeric=alpha*rePhysical;stNumeric=reNumeric*c.rho_particle/c.rho_fluid;
    mach=std::sqrt(3.)*.5*c.shear_rate*c.box_y*dt/c.dx;
    // Quick screening: report Mach/Re/St; do not impose a creeping-flow gate
    // or silently change the user-selected fluid time step.
    if(!(std::isfinite(dt)&&dt>0&&std::isfinite(rhoFluid)&&rhoFluid>0))
      throw std::runtime_error("Nonfinite time/unit mapping");
  }
};

} } // SLURRY SCOPE END
#endif
