/* SPDX-License-Identifier: GPL-2.0-or-later
 * Pure graphite suspension: resolved HLBM + full RE2 + rough sliding/rolling
 * contacts + Lees--Edwards shear.
 * Quick magnitude screening with explicitly selected fluid dt/Mach.
 */
#include <olb.h>
#include "quickConfig.h"
#include "periodicCoupling.h"
#include "leesEdwards.h"
#include "particleSubsteps.h"
#include "bulkStress.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

// SLURRY SCOPE BEGIN
namespace slurry { namespace gr_re2 {

using namespace olb;
using namespace olb::names;
using namespace olb::particles;
using T=double;
using D=descriptors::D3Q19<descriptors::POROSITY,descriptors::VELOCITY_NUMERATOR,
                         descriptors::VELOCITY_DENOMINATOR>;
using P=descriptors::ResolvedParticle3D;
using PS=ParticleSystem<T,P>;
using Problem=Case<NavierStokes,Lattice<T,D>>;
using U64=std::uint64_t;
namespace fs=std::filesystem;
using Clock=std::chrono::steady_clock;
static double seconds(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}

// OpenLB initializes MPI first. PETSc therefore does not own MPI finalization.
struct ParticleSolverRuntime {
  bool ownsPetsc=false;
  explicit ParticleSolverRuntime(const Config& c) {
    if(c.particle_solver!="petsc")return;
#ifdef SLURRY_USE_PETSC
    PetscBool initialized=PETSC_FALSE;
    if(PetscInitialized(&initialized))throw std::runtime_error("Cannot query PETSc initialization");
    if(!initialized){
      if(PetscInitializeNoArguments())throw std::runtime_error("PETSc initialization failed");
      ownsPetsc=true;
    }
#else
    throw std::runtime_error("PETSc particle solver requested, but this executable was not built with PETSc. Run build_slurry_cpu.sbatch.");
#endif
  }
  ~ParticleSolverRuntime(){
#ifdef SLURRY_USE_PETSC
    if(ownsPetsc)PetscFinalize();
#endif
  }
};

std::vector<std::array<T,7>> readParticles(const Config& c){
  std::ifstream in(c.particles_csv);if(!in)throw std::runtime_error("Cannot read particles_csv: "+c.particles_csv);
  std::string line;std::getline(in,line);
  if(strip(line)!="id,x_m,y_m,z_m,angle_x_deg,angle_y_deg,angle_z_deg")
    throw std::runtime_error("Particle CSV header mismatch");
  std::vector<std::array<T,7>> rows;
  while(std::getline(in,line)){
    if(strip(line).empty())continue;
    std::replace(line.begin(),line.end(),',',' ');
    std::istringstream row(line);std::array<T,7> r{};
    for(T& x:r)if(!(row>>x)||!std::isfinite(x))throw std::runtime_error("Invalid particle row");
    if(r[0]!=T(rows.size()))throw std::runtime_error("Particle IDs must be contiguous from zero");
    rows.push_back(r);
  }
  if(rows.empty())throw std::runtime_error("This suspension run requires particles");
  return rows;
}

class InitialShear : public AnalyticalF3D<T,T> {
  T gamma,height,conversion;
public:
  InitialShear(const Config& c,const Units& u):AnalyticalF3D<T,T>(3),gamma(c.shear_rate),height(c.box_y),conversion(u.dt/c.dx){}
  bool operator()(T out[],const T in[])override{
    out[0]=gamma*(in[1]-height/2)*conversion;out[1]=out[2]=0;return true;
  }
};
void prepare(Problem& problem,const Config& c,const Units& u){
  auto& g=problem.getGeometry();
  for(int b=0;b<g.getLoadBalancer().size();++b){auto& bg=g.getBlockGeometry(b);
    bg.forSpatialLocations([&](LatticeR<3> loc){bg.set(loc,1);});
  }
  g.communicate();g.updateStatistics(false);
  auto& l=problem.getLattice(NavierStokes{});
  l.setUnitConverter(c.dx,u.dt,c.diameter,c.shear_rate*c.diameter,u.nu,u.rhoFluid);
  olb::dynamics::set<PorousParticleBGKdynamics>(l,g.getMaterialIndicator({1}));
  l.setParameter<descriptors::OMEGA>(u.omega);
  auto domain=g.getMaterialIndicator({1});
  fields::set<descriptors::POROSITY>(l,domain,1.);
  fields::set<descriptors::VELOCITY_NUMERATOR>(l,domain,Vector<T,3>(0.));
  fields::set<descriptors::VELOCITY_DENOMINATOR>(l,domain,0.);
  AnalyticalConst3D<T,T> rho(1);InitialShear initial(c,u);
  l.iniEquilibrium(domain,rho,initial);l.defineRhoU(domain,rho,initial);l.initialize();
}

std::vector<graphite::Body> addParticles(PS& ps,const Config& c,const Units& u){
  std::vector<graphite::Body> bodies;
  const graphite::Vec3 axes{c.diameter/2,c.diameter/2,c.thickness/2};
  const T mass=4*std::acos(-1.)/3*axes[0]*axes[1]*axes[2]*u.rhoParticle;
  for(const auto& r:readParticles(c)){
    creators::addResolvedEllipsoid3D(ps,Vector<T,3>{r[1],r[2],r[3]},
      Vector<T,3>{axes[0],axes[1],axes[2]},c.epsilon_cells*c.dx,u.rhoParticle,
      Vector<T,3>{r[4],r[5],r[6]},Vector<T,3>{c.shear_rate*(r[2]-c.box_y/2),0,0});
    auto p=ps.get(ps.size()-1);
    const Vector<T,3> spin{0,0,-c.shear_rate/2};
    p.setField<descriptors::MOBILITY,descriptors::ANG_VELOCITY>(spin);
    graphite::Body b;b.axes=axes;b.mass=mass;
    b.position={r[1],r[2],r[3]};b.velocity={c.shear_rate*(r[2]-c.box_y/2),0,0};b.omega={0,0,spin[2]};
    auto rotation=access::getRotationMatrix(p);for(int k=0;k<9;++k)b.rotation[k]=rotation[k];
    for(int k=0;k<3;++k)b.inertiaBody[k]=mass/5*(axes[(k+1)%3]*axes[(k+1)%3]+axes[(k+2)%3]*axes[(k+2)%3]);
    bodies.push_back(b);
  }
  return bodies;
}
graphite::Vec3 euler(const graphite::Mat3& r){
  const T y=std::asin(std::clamp(-r[6],-1.,1.));
  if(std::abs(std::cos(y))>1e-10)return {std::atan2(r[7],r[8]),y,std::atan2(r[3],r[0])};
  return {0.,y,std::atan2(-r[1],r[4])};
}
void syncParticles(PS& ps,const std::vector<graphite::Body>& bodies){
  for(std::size_t i=0;i<bodies.size();++i){auto p=ps.get(i);const auto& b=bodies[i];
    p.setField<descriptors::GENERAL,descriptors::POSITION>(Vector<T,3>{b.position[0],b.position[1],b.position[2]});
    p.setField<descriptors::MOBILITY,descriptors::VELOCITY>(Vector<T,3>{b.velocity[0],b.velocity[1],b.velocity[2]});
    p.setField<descriptors::MOBILITY,descriptors::ANG_VELOCITY>(Vector<T,3>{b.omega[0],b.omega[1],b.omega[2]});
    Vector<T,9> rotation;for(int k=0;k<9;++k)rotation[k]=b.rotation[k];
    p.setField<descriptors::SURFACE,descriptors::ROT_MATRIX>(rotation);
    const auto angle=euler(b.rotation);p.setField<descriptors::SURFACE,descriptors::ANGLE>(Vector<T,3>{angle[0],angle[1],angle[2]});
  }
}
void writePoses(std::ofstream& out,const std::vector<graphite::Body>& bodies,U64 step,T time){
  if(!singleton::mpi().isMainProcessor())return;
  for(std::size_t i=0;i<bodies.size();++i){const auto& b=bodies[i];const auto angle=euler(b.rotation);
    out<<step<<','<<time<<','<<i;
    for(const auto& vector:{b.position,angle,b.velocity,b.omega})for(T v:vector)out<<','<<v;
    out<<'\n';
  }out.flush();
}

void simulate(const Config& c){
  const Units u(c);OstreamManager log(std::cout,"graphiteCouette3d");
  if(singleton::mpi().isMainProcessor())fs::create_directories(c.output_dir);
  singleton::mpi().barrier();singleton::directories().setOutputDir(c.output_dir+"/vtk/");
  IndicatorCuboid3D<T> indicator(Vector<T,3>{c.box_x-c.dx,c.box_y-c.dx,c.box_z-c.dx},Vector<T,3>(0.));
  Mesh<T,3> mesh(indicator,c.dx,singleton::mpi().getSize());mesh.setOverlap(2);
  mesh.getCuboidDecomposition().setPeriodicity(Vector<bool,3>(true,true,true));
  Problem::ParametersD parameters;Problem problem(parameters,mesh);prepare(problem,c,u);
  auto& l=problem.getLattice(NavierStokes{});auto& g=problem.getGeometry();auto& converter=l.getUnitConverter();
  PS ps;auto bodies=addParticles(ps,c,u);const graphite::Vec3 box{c.box_x,c.box_y,c.box_z};
  graphiteCoupling::CachedCoupling<T,D,P> coupling(g,l,converter,
    Vector<T,3>{c.box_x,c.box_y,c.box_z},c.shear_rate,Vector<T,3>{c.diameter/2,c.diameter/2,c.thickness/2});
  graphiteCoupling::LeesEdwardsBoundary<T,D> le(l,c.dx,u.dt,c.shear_rate,box);
  graphite::ParticleStepSettings solver;
  solver.box=box;solver.shearRate=c.shear_rate;
  solver.pair={c.hamaker,c.sigma_lj,c.switch_gap,c.cutoff_gap};
  solver.pair.roughnessGap=c.roughness_gap;
  solver.pair.localGap=c.local_gap;
  solver.pair.localGapFraction=c.local_gap_fraction;
  solver.pair.localSwitchExcessGap=c.local_switch_excess_gap;
  solver.pair.localCutoffExcessGap=c.local_cutoff_excess_gap;
  solver.nearField.viscosity=c.dynamic_viscosity;solver.nearField.matchingGap=c.lubrication_cutoff_cells*c.dx;
  solver.nearField.enabled=c.lubrication_cutoff_cells>0;
  solver.maxSubsteps=c.particle_max_substeps;solver.maxNewtonIterations=c.particle_max_iterations;solver.relativeTolerance=c.particle_tolerance;
  solver.maxKrylovIterations=c.particle_max_krylov_iterations;
  solver.solverBackend=c.particle_solver;
  if(c.solver_diagnostics)solver.solverDiagnosticsPrefix=(fs::path(c.output_dir)/
      ("particle_solver_rank"+std::to_string(singleton::mpi().getRank()))).string();
  solver.forceAbsoluteTolerance=c.particle_force_absolute_tolerance;
  solver.torqueAbsoluteTolerance=c.particle_torque_absolute_tolerance;
  solver.contactGapTolerance=c.contact_gap_tolerance;
  solver.rough.enabled=c.rough_contact_enabled;
  solver.rough.gap=c.roughness_gap;
  solver.rough.friction=c.sliding_friction;
  solver.rough.tangentialStiffness=c.tangential_stiffness;
  solver.rough.rollingLength=c.rolling_length;
  solver.rough.rollingYieldAngle=c.rolling_yield_angle;
  std::vector<graphite::GapCache> pairCache;
  graphite::PersistentContactState contacts;
  auto pair=graphite::evaluateParticleState(bodies,0.,solver,&pairCache,&contacts);
  std::vector<graphite::Vec3> angularAcceleration(bodies.size()),force(bodies.size()),torque(bodies.size());
  const long double desired=std::ceil(static_cast<long double>(c.end_strain)/(c.shear_rate*u.dt));
  if(desired>std::numeric_limits<U64>::max())throw std::runtime_error("Requested strain exceeds step counter capacity");
  const U64 endStep=static_cast<U64>(desired),stopStep=c.max_steps?std::min(endStep,c.max_steps):endStep;
  U64 step=0;
  log<<std::setprecision(12)<<"dt_s="<<u.dt<<" imposed_Mach="<<u.mach<<" physical_Re="<<u.rePhysical
      <<" numerical_Re="<<u.reNumeric<<" numerical_St="<<u.stNumeric<<" inertia_scale="<<u.alpha
      <<" particles="<<bodies.size()<<" steps_to_target_strain="<<endStep<<std::endl;
  log<<"rough_contact="<<c.rough_contact_enabled<<" roughness_gap_nm="<<c.roughness_gap*1.e9
     <<" sliding_friction="<<c.sliding_friction<<" tangential_stiffness_N_m="<<c.tangential_stiffness
     <<" rolling_length_nm="<<c.rolling_length*1.e9<<" rolling_yield_angle_rad="<<c.rolling_yield_angle
     <<" end_strain="<<c.end_strain<<std::endl;
  log<<"hamaker_J="<<c.hamaker<<" sigma_lj_nm="<<c.sigma_lj*1.e9
     <<" local_gap_nm="<<c.local_gap*1.e9<<" local_gap_fraction="<<c.local_gap_fraction
     <<" local_switch_excess_gap_nm="<<c.local_switch_excess_gap*1.e9
     <<" local_cutoff_excess_gap_nm="<<c.local_cutoff_excess_gap*1.e9
     <<" switch_gap_nm="<<c.switch_gap*1.e9<<" cutoff_gap_nm="<<c.cutoff_gap*1.e9<<std::endl;
  log<<"particle_solver="<<c.particle_solver<<" particle_tolerance="<<c.particle_tolerance
     <<" force_absolute_tolerance_N="<<c.particle_force_absolute_tolerance
     <<" torque_absolute_tolerance_N_m="<<c.particle_torque_absolute_tolerance
     <<" contact_gap_tolerance_m="<<c.contact_gap_tolerance
     <<" max_substeps="<<c.particle_max_substeps<<" max_newton_iterations="<<c.particle_max_iterations
     <<" max_krylov_iterations="<<c.particle_max_krylov_iterations
     <<" solver_diagnostics="<<c.solver_diagnostics<<std::endl;
  std::ofstream history,poses;
  if(singleton::mpi().isMainProcessor()){
    std::ofstream meta(fs::path(c.output_dir)/"mapping.json");
    meta<<std::setprecision(17)<<"{\n\"dt_s\":"<<u.dt<<",\n\"imposed_mach\":"<<u.mach<<",\n\"inertia_scale\":"<<u.alpha
      <<",\n\"rho_fluid_numeric_kg_m3\":"<<u.rhoFluid<<",\n\"rho_particle_numeric_kg_m3\":"<<u.rhoParticle
      <<",\n\"Re_physical\":"<<u.rePhysical<<",\n\"Re_numeric\":"<<u.reNumeric<<",\n\"St_numeric\":"<<u.stNumeric
      <<",\n\"particle_count\":"<<bodies.size()<<",\n\"steps_to_target_strain\":"<<endStep
      <<",\n\"interaction\":{\"hamaker_J\":"<<c.hamaker<<",\"sigma_lj_m\":"<<c.sigma_lj
      <<",\"switch_gap_m\":"<<c.switch_gap<<",\"cutoff_gap_m\":"<<c.cutoff_gap
      <<",\"local_gap_m\":"<<c.local_gap<<",\"local_gap_fraction\":"<<c.local_gap_fraction
      <<",\"local_switch_excess_gap_m\":"<<c.local_switch_excess_gap
      <<",\"local_cutoff_excess_gap_m\":"<<c.local_cutoff_excess_gap<<"}"
      <<",\n\"rough_contact\":{\"enabled\":"<<(c.rough_contact_enabled?"true":"false")
      <<",\"roughness_gap_m\":"<<c.roughness_gap<<",\"sliding_friction\":"<<c.sliding_friction
      <<",\"tangential_stiffness_N_m\":"<<c.tangential_stiffness
      <<",\"rolling_length_m\":"<<c.rolling_length<<",\"rolling_yield_angle_rad\":"<<c.rolling_yield_angle<<"}\n}\n";
    history.open(fs::path(c.output_dir)/"history.csv");poses.open(fs::path(c.output_dir)/"particles.csv");
    if(!history||!poses)throw std::runtime_error("Cannot create output CSV");
    history<<"step,time_s,strain,eta_bulk_Pa_s,eta_relative,stress_total_Pa,stress_fluid_Pa,stress_surface_Pa,stress_pair_attractive_Pa,stress_pair_repulsive_Pa,stress_lubrication_Pa,stress_acceleration_Pa,stress_fluid_reynolds_Pa,stress_particle_reynolds_Pa,stress_noninertial_Pa,stress_inertial_Pa,max_mach,particle_mach_bound,density_drift,porosity_volume_fraction,analytic_volume_fraction,min_gap_m,max_pair_force_N,potential_energy_J,active_pairs,particle_substeps,newton_iterations,krylov_iterations,residual_evaluations,wall_seconds,steps_per_second,fluid_seconds,map_seconds,coupling_seconds,particle_seconds,output_seconds,stress_contact_normal_Pa,stress_contact_tangential_Pa,contact_count,sliding_contact_count,rolling_contact_count,contact_dissipation_W,contact_elastic_energy_J,force_residual_ratio,torque_residual_ratio,contact_gap_violation_m,max_fluid_mach,fluid_density_drift\n";
    poses<<"step,time_s,id,x_m,y_m,z_m,angle_x_rad,angle_y_rad,angle_z_rad,vx_m_s,vy_m_s,vz_m_s,omega_x_s_inv,omega_y_s_inv,omega_z_s_inv\n";
    history<<std::setprecision(17);poses<<std::setprecision(17);
  }
  SuperVTMwriter3D<T> writer("graphite");SuperLatticePhysVelocity3D<T,D> velocity(l,converter);
  SuperLatticePhysPressure3D<T,D> pressure(l,converter);SuperLatticePhysExternalPorosity3D<T,D> porosity(l,converter);
  writer.addFunctor(velocity);writer.addFunctor(pressure);writer.addFunctor(porosity);if(c.vtk_every)writer.createMasterFile();
  const auto start=Clock::now();T fluidSeconds=0,mapSeconds=0,couplingSeconds=0,particleSeconds=0,outputSeconds=0;
  // Full communication must precede mapping: mapping fills LE image auxiliary
  // fields in all halos; ordinary periodic communication would overwrite them.
  auto mapAndCouple=[&](){
    l.setProcessingContext(ProcessingContext::Evaluation);auto begin=Clock::now();l.communicate();
    coupling.mapParticles(ps,step*u.dt);mapSeconds+=seconds(begin);
    begin=Clock::now();le.refreshForceGhosts(step*u.dt);coupling.coupleFluidToParticles(ps,step*u.dt);couplingSeconds+=seconds(begin);
  };
  auto sample=[&](){const auto begin=Clock::now();
    // Normal contact reactions and friction history come from the accepted
    // particle solve. They cannot be reconstructed from the final gap alone.
    const auto instantaneous=graphite::evaluateParticleState(bodies,step*u.dt,solver,&pairCache,&contacts);
    const auto m=graphite::measureBulk(l,g,converter,box,c.shear_rate,u.omega,u.rhoFluid,c.dx,u.dt,bodies,angularAcceleration,
      coupling.stressletSum(),instantaneous.attractiveMoment,instantaneous.repulsiveMoment,instantaneous.lubricationMoment,
      instantaneous.contactNormalMoment,instantaneous.contactTangentialMoment);
    const T elapsed=seconds(start);
    if(singleton::mpi().isMainProcessor()){
      history<<step<<','<<step*u.dt<<','<<step*u.dt*c.shear_rate<<','<<m.eta<<','<<m.etaRelative<<','<<m.stressTotal
        <<','<<m.stressFluid<<','<<m.stressSurface<<','<<m.stressPairAttractive<<','<<m.stressPairRepulsive<<','<<m.stressLubrication
        <<','<<m.stressAcceleration<<','<<m.stressFluidReynolds<<','<<m.stressParticleReynolds<<','<<m.stressNonInertial<<','<<m.stressInertial
        <<','<<m.mach<<','<<m.particleMachUpperBound<<','<<m.densityDrift<<','<<m.phi<<','<<m.phiAnalytic
        <<','<<pair.minGap<<','<<pair.maxForce<<','<<pair.energyAtEnd<<','<<pair.activePairs
        <<','<<pair.substeps<<','<<pair.newtonIterations<<','<<pair.krylovIterations<<','<<pair.residualEvaluations
        <<','<<elapsed<<','<<step/std::max(elapsed,1e-12)<<','<<fluidSeconds<<','<<mapSeconds<<','<<couplingSeconds<<','<<particleSeconds<<','<<outputSeconds
        <<','<<m.stressContactNormal<<','<<m.stressContactTangential
        <<','<<instantaneous.contacts<<','<<instantaneous.slidingContacts<<','<<instantaneous.rollingContacts
        <<','<<pair.contactDissipation<<','<<instantaneous.elasticContactEnergy
        <<','<<pair.maxForceResidualRatio<<','<<pair.maxTorqueResidualRatio<<','<<pair.contactGapViolation
        <<','<<m.fluidMach<<','<<m.fluidDensityDrift<<'\n';history.flush();
    }
    writePoses(poses,bodies,step,step*u.dt);
    log<<"step="<<step<<" strain="<<step*u.dt*c.shear_rate<<" stress_bulk_Pa="<<m.stressTotal
      <<" eta_bulk_Pa_s="<<m.eta<<" Ma="<<m.mach<<" fluid_Ma="<<m.fluidMach<<" min_gap_nm="<<pair.minGap*1e9
      <<" contacts="<<instantaneous.contacts<<" residual_ratio="<<std::max(pair.maxForceResidualRatio,pair.maxTorqueResidualRatio)
      <<" steps/s="<<step/std::max(elapsed,1e-12)<<std::endl;
    outputSeconds+=seconds(begin);
  };
  mapAndCouple();sample();
  while(step<stopStep){
    const auto& hydro=coupling.particleHydrodynamics();
    for(std::size_t i=0;i<bodies.size();++i){force[i]=hydro[i].force;torque[i]=hydro[i].torque;}
    std::vector<graphite::Vec3> oldOmega;oldOmega.reserve(bodies.size());for(const auto& b:bodies)oldOmega.push_back(b.omega);
    auto begin=Clock::now();pair=graphite::advanceParticles(bodies,force,torque,u.dt,step*u.dt,solver,&pairCache,&contacts);
    for(std::size_t i=0;i<bodies.size();++i)angularAcceleration[i]=graphite::scale(graphite::sub(bodies[i].omega,oldOmega[i]),1/u.dt);
    syncParticles(ps,bodies);particleSeconds+=seconds(begin);
    // Explicit resolved coupling uses the beginning-of-step particle mask.
    // Pair/lubrication integration is implicit within the unchanged LB step.
    begin=Clock::now();l.setProcessingContext(ProcessingContext::Simulation);
    l.collide();le.apply(step*u.dt);l.AndStream();fluidSeconds+=seconds(begin);
    ++step;mapAndCouple();
    if(step%c.sample_every==0||step==stopStep)sample();
    if(c.vtk_every&&step%c.vtk_every==0){
      const U64 frame=step/c.vtk_every;if(frame>std::numeric_limits<int>::max())throw std::runtime_error("Too many VTK frames");
      writer.write(static_cast<int>(frame));
    }
  }
  if(singleton::mpi().isMainProcessor()){
    std::ofstream status(fs::path(c.output_dir)/"status.json");
    status<<std::setprecision(17)<<"{\"status\":\""<<(step>=endStep?"COMPLETED":"MAX_STEPS")
      <<"\",\"step\":"<<step<<",\"time_s\":"<<step*u.dt<<",\"strain\":"<<step*u.dt*c.shear_rate
      <<",\"wall_seconds\":"<<seconds(start)<<"}\n";
  }
}
int runCase(int argc,char** argv){
  if(argc==2&&std::string(argv[1])=="--build-info"){
#ifdef PARALLEL_MODE_MPI
    std::cout<<"{\"mpi_enabled\":true,\"rough_contact\":true,\"local_gap_adhesion\":true,\"revision\":\"local-gap-adhesion-1\"}\n";
#else
    std::cout<<"{\"mpi_enabled\":false,\"rough_contact\":true,\"local_gap_adhesion\":true,\"revision\":\"local-gap-adhesion-1\"}\n";
#endif
    return 0;
  }
  initialize(&argc,&argv);
  try{
    if(argc==2&&std::string(argv[1])=="--help"){
      if(singleton::mpi().isMainProcessor())std::cout<<"graphiteCouette3d --config run.cfg [--max-steps N]\n";
      return 0;
    }
    const auto config=parseConfig(argc,argv);
    ParticleSolverRuntime particleRuntime(config);
    simulate(config);return 0;
  }catch(const std::exception& e){
    std::cerr<<"graphiteCouette3d: "<<e.what()<<std::endl;
#ifdef PARALLEL_MODE_MPI
    MPI_Abort(MPI_COMM_WORLD,2);
#endif
    return 2;
  }
}

} } // SLURRY SCOPE END
