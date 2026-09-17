// SPDX-License-Identifier: GPL-2.0-or-later
// CMC Cross fluid: periodic-x, node-to-node planar Couette rheometer.
// OpenLB 1.9r0. Geometry/setter idioms follow the OpenLB cavity2d/poiseuille2d
// examples (OpenLB authors, GPL-2.0-or-later). No OpenLB core changes required.
#include "crossModel.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// SLURRY SCOPE BEGIN
namespace slurry { namespace cmc {

using namespace olb;
using namespace olb::names;
using T = double;
using D = descriptors::D2Q9<descriptors::OMEGA>;
using MyCase = Case<NavierStokes,Lattice<T,D>>;

struct Options {
  T eta0, etaInf, lambda, m;
  T concentration=16, gamma=100, gap=1e-4, rho=1000, dt=0;
  int resolution=64, width=16, maxSteps=400000, checkEvery=500, vtkEvery=0;
  T convergenceTol=1e-5, errorTol=0.01;
  bool vtk=true;
};

T number(const std::string& s) {
  std::size_t n=0;
  const T v=std::stod(s,&n);
  if (n!=s.size() || !std::isfinite(v)) throw std::invalid_argument("Invalid number: "+s);
  return v;
}

Options parseOptions(int argc,char** argv) {
  Options o{};
  std::map<std::string,std::string> a;
  for (int i=1;i<argc;++i) {
    const std::string key=argv[i];
    if (key=="--no-vtk") { o.vtk=false; continue; }
    if (key=="--help") {
      if (singleton::mpi().isMainProcessor()) std::cout
        << "Cross Couette 2D (OpenLB 1.9r0)\n"
        << "Required: --eta0 Pa.s --eta-inf Pa.s --lambda s --m exponent\n"
        << "Optional: --concentration g/L --shear-rate 1/s --resolution 64\n"
        << " --gap 1e-4 --density 1000 --dt seconds --width-cells 16\n"
        << " --max-steps 400000 --check-every 500 --vtk-every 0 --no-vtk\n"
        << " --convergence-tol 1e-5 --error-tol 0.01\n";
      std::exit(0);
    }
    if (i+1==argc || key.rfind("--",0)!=0 || a.count(key))
      throw std::invalid_argument("Invalid/duplicate option: "+key);
    a[key]=argv[++i];
  }
  auto read=[&](std::string key,T fallback,bool required=false) {
    auto it=a.find(key);
    if (it==a.end()) {
      if (required) throw std::invalid_argument("Missing "+key);
      return fallback;
    }
    const T v=number(it->second); a.erase(it); return v;
  };
  auto integer=[&](std::string key,int fallback) {
    T v=read(key,fallback);
    if (v<0 || v>1e8 || v!=std::floor(v)) throw std::invalid_argument(key+" must be a nonnegative integer");
    return static_cast<int>(v);
  };
  o.eta0=read("--eta0",0,true); o.etaInf=read("--eta-inf",0,true);
  o.lambda=read("--lambda",0,true); o.m=read("--m",0,true);
  o.concentration=read("--concentration",o.concentration);
  o.gamma=read("--shear-rate",o.gamma); o.gap=read("--gap",o.gap);
  o.rho=read("--density",o.rho); o.dt=read("--dt",o.dt);
  o.resolution=integer("--resolution",o.resolution);
  o.width=integer("--width-cells",o.width);
  o.maxSteps=integer("--max-steps",o.maxSteps);
  o.checkEvery=integer("--check-every",o.checkEvery);
  o.vtkEvery=integer("--vtk-every",o.vtkEvery);
  o.convergenceTol=read("--convergence-tol",o.convergenceTol);
  o.errorTol=read("--error-tol",o.errorTol);
  if (!a.empty()) throw std::invalid_argument("Unknown option: "+a.begin()->first);
  if (!(o.eta0>0 && o.etaInf>=0 && o.etaInf<=o.eta0 && o.lambda>0 && o.m>0
      && o.gamma>0 && o.gap>0 && o.rho>0 && o.dt>=0 && o.concentration>=0
      && o.convergenceTol>0 && o.errorTol>0)) throw std::invalid_argument("Invalid physical parameters");
  if (o.resolution<16 || o.width<8 || o.checkEvery<1 || o.maxSteps<o.checkEvery*5)
    throw std::invalid_argument("Use resolution>=16, width>=8 and >=5 convergence checks");
  if (singleton::mpi().getSize()>o.resolution/4)
    throw std::invalid_argument("Too many MPI ranks for this small grid; use <= resolution/4");
  return o;
}

T crossViscosity(const Options& o,T g) {
  return o.etaInf+(o.eta0-o.etaInf)/(1+std::pow(o.lambda*g,o.m));
}

struct Diagnostics {
  std::vector<T> ux, uy, stress, omega, rho, localGamma, count;
  T stressMean=0, eta=0, gammaSlope=0, profileError=0, densityDrift=0;
  T maxMach=0, tauMin=1e30, tauMax=0, stressNonuniformity=0;
};

// Gather only core nodes, never overlap nodes. All MPI ranks obtain same rows.
Diagnostics measure(MyCase& problem,const Options& o,T dx,T dt) {
  auto& lattice=problem.getLattice(NavierStokes{});
  auto& geometry=problem.getGeometry();
  lattice.setProcessingContext(ProcessingContext::Evaluation);
  const int rows=o.resolution+1, columns=7;
  std::vector<T> packed(rows*columns,0), reduced(packed.size(),0);
  T localMaxMach=0, localDensityDrift=0, localTauMin=1e30, localTauMax=0;
  int invalid=0;
  const T pressureFactor=o.rho*(dx/dt)*(dx/dt);
  for (int b=0;b<lattice.getLoadBalancer().size();++b) {
    auto& bg=geometry.getBlockGeometry(b);
    auto& bl=lattice.getBlock(b);
    bg.forCoreSpatialLocations([&](LatticeR<2> loc) {
      const int material=bg.getMaterial(loc);
      if (material<1 || material>3) return;
      const auto position=bg.getPhysR(loc);
      const int j=std::lround(position[1]/dx);
      if (j<0 || j>=rows) { invalid=1; return; }
      auto cell=bl.get(loc);
      T rho, u[2], pi[3];
      cell.computeAllMomenta(rho,u,pi);
      const T omega=cell.template getField<descriptors::OMEGA>();
      if (!(rho>0 && omega>0 && omega<2 && std::isfinite(rho)
         && std::isfinite(u[0]) && std::isfinite(u[1])
         && std::isfinite(pi[0]) && std::isfinite(pi[1]) && std::isfinite(pi[2]))) {
        invalid=1; return;
      }
      const T piNorm=std::sqrt(pi[0]*pi[0]+2*pi[1]*pi[1]+pi[2]*pi[2]);
      const T gamma=3*omega*piNorm/(std::sqrt(T(2))*rho*dt);
      // Pre-collision non-equilibrium momentum flux -> viscous Cauchy stress.
      // This does NOT evaluate eta(Cross)*gamma as the measured stress.
      const T sigma=-(1-omega/2)*pi[1]*pressureFactor;
      if (!std::isfinite(gamma) || !std::isfinite(sigma)) { invalid=1; return; }
      T* row=&packed[j*columns];
      row[0]+=1; row[1]+=u[0]*dx/dt; row[2]+=u[1]*dx/dt;
      row[3]+=sigma; row[4]+=omega; row[5]+=rho; row[6]+=gamma;
      localDensityDrift=std::max(localDensityDrift,std::abs(rho-1));
      localMaxMach=std::max(localMaxMach,std::sqrt(T(3)*(u[0]*u[0]+u[1]*u[1])));
      localTauMin=std::min(localTauMin,1/omega); localTauMax=std::max(localTauMax,1/omega);
    });
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().allreduce(packed.data(),reduced.data(),static_cast<int>(packed.size()),MPI_SUM);
  singleton::mpi().reduceAndBcast(invalid,MPI_MAX);
  singleton::mpi().reduceAndBcast(localMaxMach,MPI_MAX);
  singleton::mpi().reduceAndBcast(localDensityDrift,MPI_MAX);
  singleton::mpi().reduceAndBcast(localTauMin,MPI_MIN);
  singleton::mpi().reduceAndBcast(localTauMax,MPI_MAX);
#else
  reduced=packed;
#endif
  if (invalid) throw std::runtime_error("Nonfinite populations, invalid density/omega, or invalid grid row");
  Diagnostics d;
  d.ux.resize(rows); d.uy.resize(rows); d.stress.resize(rows); d.omega.resize(rows);
  d.rho.resize(rows); d.localGamma.resize(rows); d.count.resize(rows);
  d.maxMach=localMaxMach; d.densityDrift=localDensityDrift;
  d.tauMin=localTauMin; d.tauMax=localTauMax;
  const T U=o.gamma*o.gap;
  for (int j=0;j<rows;++j) {
    const T* r=&reduced[j*columns];
    if (r[0]!=o.width) throw std::runtime_error("Geometry node count mismatch in row "+std::to_string(j));
    d.count[j]=r[0]; d.ux[j]=r[1]/r[0]; d.uy[j]=r[2]/r[0];
    d.stress[j]=r[3]/r[0]; d.omega[j]=r[4]/r[0]; d.rho[j]=r[5]/r[0];
    d.localGamma[j]=r[6]/r[0];
    d.profileError=std::max(d.profileError,std::abs(d.ux[j]-U*T(j)/o.resolution)/U);
  }
  const int lo=o.resolution/4, hi=3*o.resolution/4;
  T sumY=0, sumU=0, sumYY=0, sumYU=0;
  for (int j=lo;j<=hi;++j) {
    const T y=j*dx;
    d.stressMean+=d.stress[j]; sumY+=y; sumU+=d.ux[j]; sumYY+=y*y; sumYU+=y*d.ux[j];
  }
  const T n=hi-lo+1;
  d.stressMean/=n;
  d.gammaSlope=(n*sumYU-sumY*sumU)/(n*sumYY-sumY*sumY);
  d.eta=d.stressMean/o.gamma;
  for (int j=lo;j<=hi;++j)
    d.stressNonuniformity=std::max(d.stressNonuniformity,
      std::abs(d.stress[j]-d.stressMean)/std::max(std::abs(d.stressMean),T(1e-100)));
  return d;
}

void saveProfile(const Diagnostics& d,const Options& o,T dx,T dt) {
  if (!singleton::mpi().isMainProcessor()) return;
  std::ofstream f("profile.csv");
  f << std::setprecision(17)
    << "y_m,ux_m_s,uy_m_s,shear_stress_Pa,omega,rho_lattice,shear_rate_from_neq_s_inv,eta_from_omega_Pa_s\n";
  for (int j=0;j<=o.resolution;++j) {
    const T viscosity=o.rho*(1/d.omega[j]-.5)/3*dx*dx/dt;
    f<<j*dx<<','<<d.ux[j]<<','<<d.uy[j]<<','<<d.stress[j]<<','<<d.omega[j]
      <<','<<d.rho[j]<<','<<d.localGamma[j]<<','<<viscosity<<'\n';
  }
  if (!f) throw std::runtime_error("Could not write profile.csv");
}

void writeVTK(MyCase& problem,SuperVTMwriter2D<T>& writer,int step) {
  auto& lattice=problem.getLattice(NavierStokes{});
  lattice.setProcessingContext(ProcessingContext::Evaluation);
  writer.write(step);
}

int run(const Options& o) {
  const T dx=o.gap/o.resolution;
  const T etaTarget=crossViscosity(o,o.gamma);
  const T nuRef=std::sqrt(o.eta0*crossViscosity(o,1000))/o.rho;
  const T dt=o.dt>0 ? o.dt : .1*dx*dx/nuRef;
  const T U=o.gamma*o.gap, latticeU=U*dt/dx;
  if (!(dt>0 && std::isfinite(dt) && std::sqrt(T(3))*latticeU<.1))
    throw std::invalid_argument("Invalid dt or imposed lattice Mach >= 0.1; refine grid / adjust dt");
  OstreamManager log(std::cout,"CrossCouette");
  log << std::setprecision(12) << "c=" << o.concentration << " g/L, shear=" << o.gamma
      << " 1/s, etaTarget=" << etaTarget << " Pa.s, dx=" << dx << " m, dt=" << dt << " s" << std::endl;

  IndicatorCuboid2D<T> box({(o.width-1)*dx,o.gap},{0,0});
  Mesh<T,2> mesh(box,dx,singleton::mpi().getSize());
  mesh.setOverlap(2);
  mesh.getCuboidDecomposition().setPeriodicity({true,false});
  MyCase::ParametersD params;
  MyCase problem(params,mesh);
  auto& geometry=problem.getGeometry();
  // Explicit y-only materials: periodic x padding has the same physical fluid.
  for (int b=0;b<geometry.getLoadBalancer().size();++b) {
    auto& block=geometry.getBlockGeometry(b);
    block.forSpatialLocations([&](LatticeR<2> loc) {
      const T y=block.getPhysR(loc)[1];
      int material=0;
      if (std::abs(y)<.25*dx) material=2;
      else if (std::abs(y-o.gap)<.25*dx) material=3;
      else if (y>0 && y<o.gap) material=1;
      block.set(loc,material);
    });
  }
  geometry.communicate(); geometry.updateStatistics(false);

  auto& lattice=problem.getLattice(NavierStokes{});
  lattice.setUnitConverter(dx,dt,o.gap,U,nuRef,o.rho);
  auto& converter=lattice.getUnitConverter();
  converter.print();
  lattice.defineDynamics<cmc::CrossBGKdynamics>(geometry,1);
  boundary::set<boundary::LocalVelocity,cmc::CrossBGKdynamics>(lattice,geometry,2);
  boundary::set<boundary::LocalVelocity,cmc::CrossBGKdynamics>(lattice,geometry,3);
  lattice.setParameter<cmc::NU_ZERO>(o.eta0/o.rho*dt/(dx*dx));
  lattice.setParameter<cmc::NU_INF>(o.etaInf/o.rho*dt/(dx*dx));
  lattice.setParameter<cmc::LAMBDA>(o.lambda/dt);
  lattice.setParameter<cmc::M>(o.m);
  const T initialOmega=1/(.5+3*o.eta0/o.rho*dt/(dx*dx));
  lattice.setParameter<descriptors::OMEGA>(initialOmega);
  AnalyticalConst2D<T,T> rho(1), rest(0,0), omega(initialOmega);
  auto domain=geometry.getMaterialIndicator({1,2,3});
  lattice.iniEquilibrium(domain,rho,rest);
  lattice.defineRhoU(domain,rho,rest);
  lattice.defineField<descriptors::OMEGA>(domain,omega);
  lattice.initialize();

  singleton::directories().setOutputDir("./tmp/");
  SuperVTMwriter2D<T> writer("crossCouette2d");
  SuperLatticePhysVelocity2D<T,D> velocity(lattice,converter);
  SuperLatticePhysPressure2D<T,D> pressure(lattice,converter);
  SuperLatticeField2D<T,D,descriptors::OMEGA> omegaField(lattice);
  SuperLatticePhysViscosity2D<T,D> nuField(lattice,converter);
  nuField.getName()="kinematicViscosity_m2_s";
  writer.addFunctor(velocity); writer.addFunctor(pressure);
  writer.addFunctor(omegaField); writer.addFunctor(nuField);
  if (o.vtk) { writer.createMasterFile(); writeVTK(problem,writer,0); }

  std::ofstream history;
  if (singleton::mpi().isMainProcessor()) {
    history.open("history.csv");
    history << "step,time_s,eta_measured_Pa_s,shear_rate_measured_s_inv,velocity_change,stress_change,profile_relative_error,viscosity_relative_error,density_relative_drift,max_mach,tau_min,tau_max\n";
    history << std::setprecision(17);
    if (!history) throw std::runtime_error("Could not open history.csv");
  }
  const int rampSteps=std::max(1000,o.resolution*o.resolution);
  const int minSteps=2*rampSteps;
  Diagnostics previous, current;
  bool converged=false;
  int stable=0, finalStep=0;
  T velocityChange=1, stressChange=1;
  for (int step=1;step<=o.maxSteps;++step) {
    if (step<=rampSteps) {
      const T s=T(step)/rampSteps;
      AnalyticalConst2D<T,T> top(latticeU*(3*s*s-2*s*s*s),0);
      lattice.defineU(geometry,3,top);
    }
    lattice.setProcessingContext(ProcessingContext::Simulation);
    lattice.collideAndStream();
    if (step%o.checkEvery==0 || step==o.maxSteps) {
      current=measure(problem,o,dx,dt);
      if (!previous.ux.empty()) {
        velocityChange=0;
        for (std::size_t j=0;j<current.ux.size();++j)
          velocityChange=std::max(velocityChange,std::abs(current.ux[j]-previous.ux[j])/U);
        stressChange=std::abs(current.stressMean-previous.stressMean)/(etaTarget*o.gamma);
      }
      const T error=std::abs(current.eta/etaTarget-1);
      if (current.maxMach>=.1 || current.densityDrift>.01)
        throw std::runtime_error("Mach or density-drift gate exceeded");
      if (step>=minSteps && velocityChange<o.convergenceTol && stressChange<o.convergenceTol)
        ++stable;
      else stable=0;
      if (singleton::mpi().isMainProcessor()) {
        history<<step<<','<<step*dt<<','<<current.eta<<','<<current.gammaSlope
          <<','<<velocityChange<<','<<stressChange<<','<<current.profileError<<','<<error
          <<','<<current.densityDrift<<','<<current.maxMach<<','<<current.tauMin<<','<<current.tauMax<<'\n';
        history.flush();
        if (!history) throw std::runtime_error("Could not write history.csv");
      }
      if (step%(o.checkEvery*10)==0 || stable>=5 || step==o.maxSteps)
        log<<"step="<<step<<" eta="<<current.eta<<" relError="<<error
           <<" dU="<<velocityChange<<" dStress="<<stressChange<<" stableChecks="<<stable<<std::endl;
      previous=current;
      finalStep=step;
      if (stable>=5) { converged=true; break; }
    }
    if (o.vtk && o.vtkEvery>0 && step%o.vtkEvery==0) writeVTK(problem,writer,step);
  }
  const T error=std::abs(current.eta/etaTarget-1);
  const T gammaError=std::abs(current.gammaSlope/o.gamma-1);
  const bool pass=converged && error<=o.errorTol && gammaError<=o.errorTol
    && current.profileError<=o.errorTol && current.stressNonuniformity<=o.errorTol;
  saveProfile(current,o,dx,dt);
  if (o.vtk) writeVTK(problem,writer,finalStep);
  if (singleton::mpi().isMainProcessor()) {
    std::ofstream f("result.json");
    f<<std::setprecision(17)<<"{\n"
      <<"  \"status\": \""<<(pass?"PASS":"FAIL")<<"\",\n"
      <<"  \"converged\": "<<(converged?"true":"false")<<",\n"
      <<"  \"stress_method\": \"bulk_non_equilibrium_momentum_flux\",\n"
      <<"  \"concentration_g_L\": "<<o.concentration<<",\n"
      <<"  \"shear_rate_s_inv\": "<<o.gamma<<",\n"
      <<"  \"shear_rate_measured_s_inv\": "<<current.gammaSlope<<",\n"
      <<"  \"eta_target_Pa_s\": "<<etaTarget<<",\n"
      <<"  \"eta_measured_Pa_s\": "<<current.eta<<",\n"
      <<"  \"shear_stress_Pa\": "<<current.stressMean<<",\n"
      <<"  \"relative_error\": "<<error<<",\n"
      <<"  \"profile_relative_error\": "<<current.profileError<<",\n"
      <<"  \"stress_nonuniformity\": "<<current.stressNonuniformity<<",\n"
      <<"  \"density_relative_drift\": "<<current.densityDrift<<",\n"
      <<"  \"max_mach\": "<<current.maxMach<<",\n"
      <<"  \"tau_min\": "<<current.tauMin<<",\n"
      <<"  \"tau_max\": "<<current.tauMax<<",\n"
      <<"  \"steps\": "<<finalStep<<",\n"
      <<"  \"dt_s\": "<<dt<<",\n"
      <<"  \"dx_m\": "<<dx<<",\n"
      <<"  \"physical_time_s\": "<<finalStep*dt<<",\n"
      <<"  \"resolution\": "<<o.resolution<<",\n"
      <<"  \"width_cells\": "<<o.width<<",\n"
      <<"  \"gap_m\": "<<o.gap<<",\n"
      <<"  \"density_kg_m3\": "<<o.rho<<",\n"
      <<"  \"convergence_tolerance\": "<<o.convergenceTol<<",\n"
      <<"  \"acceptance_relative_error\": "<<o.errorTol<<"\n}\n";
    if (!f) throw std::runtime_error("Could not write result.json");
  }
  log<<(pass?"PASS":"FAIL")<<": c="<<o.concentration<<", shear="<<o.gamma
      <<", relative viscosity error="<<error<<std::endl;
  return pass?0:2;
}

int runCase(int argc,char** argv) {
  initialize(&argc,&argv);
  try { return run(parseOptions(argc,argv)); }
  catch (const std::exception& e) {
    std::cerr<<"ERROR [rank "<<singleton::mpi().getRank()<<"]: "<<e.what()<<'\n';
#ifdef PARALLEL_MODE_MPI
    MPI_Abort(MPI_COMM_WORLD,1);
#endif
    return 1;
  }
}

} } // SLURRY SCOPE END
