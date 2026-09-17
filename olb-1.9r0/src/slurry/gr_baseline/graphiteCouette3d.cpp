/* SPDX-License-Identifier: GPL-2.0-or-later
 * Resolved oblate graphite suspension, OpenLB 1.9r0.
 * Fluid: HLBM PorousParticleBGK. Rigid motion: Newton--Euler/Verlet.
 * Native overlap-volume elastic normal contact, zero tangential friction.
 * No attraction, gravity, Brownian forcing, CMC or lubrication correction.
 */
#include <olb.h>
#include "periodicCoupling.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// SLURRY SCOPE BEGIN
namespace slurry { namespace gr_baseline {

using namespace olb;
using namespace olb::names;
using namespace olb::particles;
using namespace olb::particles::contact;
using T=double;
using D=descriptors::D3Q19<descriptors::POROSITY,descriptors::VELOCITY_NUMERATOR,
  descriptors::VELOCITY_DENOMINATOR,descriptors::CONTACT_DETECTION>;
using P=descriptors::ResolvedParticleWithContact3D;
using PS=ParticleSystem<T,P>;
using PC=ParticleContactArbitraryFromOverlapVolume<T,3,true>;
using WC=WallContactArbitraryFromOverlapVolume<T,3,true>;
using CC=ContactContainer<T,PC,WC>;
using Problem=Case<NavierStokes,Lattice<T,D>>;
using U64=std::uint64_t;
namespace fs=std::filesystem;
static volatile std::sig_atomic_t stopRequested=0;
void requestStop(int) {stopRequested=1;}

struct Config {
  T shear_rate=10,box_x=1e-5,box_y=1e-5,box_z=1e-5,dx=5e-8;
  T diameter=3.3e-6,thickness=4e-7,rho_particle=2400,rho_fluid=997,dynamic_viscosity=.000890;
  T nu_lattice=.1,target_re=.001,epsilon_cells=.5,contact_enlargement_cells=.2;
  T contact_young=1000,contact_poisson=.3,contact_restitution=1,end_strain=1;
  U64 max_steps=0,sample_every=10000,vtk_every=0,checkpoint_every=1000000;
  unsigned contact_resolution=8;
  std::string output_dir="run",particles_csv,restart_dir;
};
std::string strip(std::string s) {
  const auto begin=s.find_first_not_of(" \t\r\n");
  if(begin==std::string::npos)return {};
  return s.substr(begin,s.find_last_not_of(" \t\r\n")-begin+1);
}
Config parse(int argc,char**argv) {
  std::string file,restart;U64 maxOverride=0;
  for(int i=1;i<argc;++i){std::string arg=argv[i];
    if(arg=="--help"){if(singleton::mpi().isMainProcessor())std::cout<<"graphiteCouette3d --config run.cfg [--max-steps N] [--restart checkpoint_dir]\n";std::exit(0);}
    if(i+1>=argc)throw std::runtime_error("Missing value for "+arg);
    if(arg=="--config")file=argv[++i];else if(arg=="--restart")restart=argv[++i];
    else if(arg=="--max-steps")maxOverride=std::stoull(argv[++i]);
    else throw std::runtime_error("Unknown argument "+arg);
  }
  if(file.empty())throw std::runtime_error("--config is required");
  std::ifstream in(file);if(!in)throw std::runtime_error("Cannot read "+file);
  Config c;std::string line;
  while(std::getline(in,line)){
    line=strip(line.substr(0,line.find('#')));if(line.empty())continue;
    auto eq=line.find('=');if(eq==std::string::npos)throw std::runtime_error("Expected key=value: "+line);
    auto key=strip(line.substr(0,eq)),v=strip(line.substr(eq+1));
#define REAL(k) if(key==#k){c.k=std::stod(v);continue;}
    REAL(shear_rate) REAL(box_x) REAL(box_y) REAL(box_z) REAL(dx) REAL(diameter) REAL(thickness)
    REAL(rho_particle) REAL(rho_fluid) REAL(dynamic_viscosity) REAL(nu_lattice) REAL(target_re)
    REAL(epsilon_cells) REAL(contact_enlargement_cells) REAL(contact_young) REAL(contact_poisson)
    REAL(contact_restitution) REAL(end_strain)
#undef REAL
#define INTEGER(k) if(key==#k){c.k=std::stoull(v);continue;}
    INTEGER(max_steps) INTEGER(sample_every) INTEGER(vtk_every) INTEGER(checkpoint_every) INTEGER(contact_resolution)
#undef INTEGER
#define STRING(k) if(key==#k){c.k=v;continue;}
    STRING(output_dir) STRING(particles_csv) STRING(restart_dir)
#undef STRING
    throw std::runtime_error("Unknown configuration key "+key);
  }
  if(maxOverride)c.max_steps=maxOverride;
  if(!restart.empty())c.restart_dir=restart;
  if(!(c.shear_rate>0&&c.dx>0&&c.box_x>0&&c.box_y>0&&c.box_z>0&&c.diameter>0&&c.thickness>0
    &&c.thickness<=c.diameter&&c.rho_particle>0&&c.rho_fluid>0&&c.dynamic_viscosity>0&&c.nu_lattice>0
    &&c.target_re>=0&&c.end_strain>0&&c.epsilon_cells>0&&c.contact_enlargement_cells>=0
    &&c.contact_young>0&&c.contact_poisson>-.99&&c.contact_poisson<.5&&c.sample_every>0
    &&c.contact_resolution>=4&&c.contact_restitution==1))
    throw std::runtime_error("Invalid config; this frictionless elastic baseline requires restitution=1");
  for(T length:{c.box_x,c.box_y,c.box_z})if(std::abs(length/c.dx-std::round(length/c.dx))>1e-8||length/c.dx<8)
    throw std::runtime_error("Box lengths must be integer multiples of dx, at least eight cells");
  c.output_dir=fs::absolute(c.output_dir).string();
  if(!c.restart_dir.empty())c.restart_dir=fs::absolute(c.restart_dir).string();
  return c;
}
struct Units {
  T rePhysical,alpha,rhoFluid,rhoParticle,nu,dt,omega,reNumeric,stNumeric,mach;
  Units(const Config& c){
    rePhysical=c.rho_fluid*c.shear_rate*c.diameter*c.diameter/c.dynamic_viscosity;
    alpha=c.target_re>0?c.target_re/rePhysical:1;
    rhoFluid=alpha*c.rho_fluid;rhoParticle=alpha*c.rho_particle;
    nu=c.dynamic_viscosity/rhoFluid;dt=c.nu_lattice*c.dx*c.dx/nu;
    omega=1/(.5+3*c.nu_lattice);reNumeric=alpha*rePhysical;
    stNumeric=reNumeric*c.rho_particle/c.rho_fluid;
    mach=std::sqrt(3.)*.5*c.shear_rate*c.box_y*dt/c.dx;
    if(!(std::isfinite(dt)&&dt>0&&mach<.1&&stNumeric<.1))
      throw std::runtime_error("Invalid mapping: dt nonfinite, imposed Mach>=0.1 or numerical St>=0.1");
  }
};
constexpr auto periodicity=[](){return Vector<bool,3>(true,false,true);};

std::vector<std::array<T,7>> readParticles(const Config&c){
  std::vector<std::array<T,7>> rows;if(c.particles_csv.empty())return rows;
  std::ifstream in(c.particles_csv);if(!in)throw std::runtime_error("Cannot open particles_csv");
  std::string line;std::getline(in,line);
  if(strip(line)!="id,x_m,y_m,z_m,angle_x_deg,angle_y_deg,angle_z_deg")
    throw std::runtime_error("Particle CSV header mismatch");
  while(std::getline(in,line)){if(strip(line).empty())continue;std::replace(line.begin(),line.end(),',',' ');
    std::istringstream row(line);std::array<T,7>a{};for(T&x:a)if(!(row>>x)||!std::isfinite(x))throw std::runtime_error("Invalid particle row");
    if(a[0]!=T(rows.size()))throw std::runtime_error("Particle IDs must be contiguous from zero");
    rows.push_back(a);
  }return rows;
}

void prepareGeometry(Problem&p,const Config&c){
  auto&g=p.getGeometry();
  for(int b=0;b<g.getLoadBalancer().size();++b){auto&bg=g.getBlockGeometry(b);
    bg.forSpatialLocations([&](LatticeR<3> loc){T y=bg.getPhysR(loc)[1];int m=0;
      if(std::abs(y)<.25*c.dx)m=2;else if(std::abs(y-c.box_y)<.25*c.dx)m=3;
      else if(y>0&&y<c.box_y)m=1;
      bg.set(loc,m);
    });
  }g.communicate();g.updateStatistics(false);
}
class InitialShear : public AnalyticalF3D<T,T> {
  T gamma,gap,conversion;
public:
  InitialShear(const Config&c,const Units&u):AnalyticalF3D<T,T>(3),gamma(c.shear_rate),gap(c.box_y),conversion(u.dt/c.dx){}
  bool operator()(T out[],const T in[])override{out[0]=gamma*(in[1]-gap/2)*conversion;out[1]=0;out[2]=0;return true;}
};
void prepareLattice(Problem&p,const Config&c,const Units&u){
  auto&l=p.getLattice(NavierStokes{});auto&g=p.getGeometry();
  l.setUnitConverter(c.dx,u.dt,c.diameter,c.shear_rate*c.diameter,u.nu,u.rhoFluid);
  olb::dynamics::set<PorousParticleBGKdynamics>(l,g.getMaterialIndicator({1}));
  boundary::set<boundary::LocalVelocity,BGKdynamics>(l,g,2);
  boundary::set<boundary::LocalVelocity,BGKdynamics>(l,g,3);
  l.setParameter<descriptors::OMEGA>(u.omega);
  auto domain=g.getMaterialIndicator({0,1,2,3});
  fields::set<descriptors::POROSITY>(l,domain,1.);
  fields::set<descriptors::VELOCITY_NUMERATOR>(l,domain,Vector<T,3>(0.));
  fields::set<descriptors::VELOCITY_DENOMINATOR>(l,domain,0.);
  AnalyticalConst3D<T,T> rho(1);InitialShear initial(c,u);
  l.iniEquilibrium(g.getMaterialIndicator({1,2,3}),rho,initial);
  l.defineRhoU(g.getMaterialIndicator({1,2,3}),rho,initial);
  // Initialize linear Couette velocity; non-equilibrium stress relaxes from equilibrium populations.
  const T uw=.5*c.shear_rate*c.box_y*u.dt/c.dx;
  AnalyticalConst3D<T,T> bottom(-uw,0,0),top(uw,0,0);
  l.defineU(g,2,bottom);l.defineU(g,3,top);
  auto&comm=l.getCommunicator(stage::PostPostProcess());
  comm.requestFields<descriptors::POROSITY,descriptors::VELOCITY_NUMERATOR,descriptors::VELOCITY_DENOMINATOR>();
  comm.requestOverlap(l.getOverlap());comm.exchangeRequests();
  l.initialize();
}
void addParticles(PS&ps,const Config&c,const Units&u){
  ps.defineDynamics<olb::particles::dynamics::VerletParticleDynamics<T,P>>();
  for(auto&r:readParticles(c)){
    creators::addResolvedEllipsoid3D(ps,Vector<T,3>{r[1],r[2],r[3]},
      Vector<T,3>{c.diameter/2,c.diameter/2,c.thickness/2},c.epsilon_cells*c.dx,u.rhoParticle,
      Vector<T,3>{r[4],r[5],r[6]},Vector<T,3>{c.shear_rate*(r[2]-c.box_y/2),0,0});
    auto particle=ps.get(ps.size()-1);
    particle.setField<descriptors::MOBILITY,descriptors::ANG_VELOCITY>(Vector<T,3>{0,0,-c.shear_rate/2});
    particle.setField<descriptors::MECHPROPERTIES,descriptors::MATERIAL>(0);
    particle.setField<descriptors::NUMERICPROPERTIES,descriptors::ENLARGEMENT_FOR_CONTACT>(c.contact_enlargement_cells*c.dx);
  }
}
std::vector<SolidBoundary<T,3>> makeWalls(const Config&c){
  std::vector<SolidBoundary<T,3>> walls;
  // Lateral extents extend beyond every periodic particle image; only y faces act as walls.
  walls.emplace_back(std::make_unique<IndicatorCuboid3D<T>>(
    Vector<T,3>{5*c.box_x,2*c.diameter,5*c.box_z},Vector<T,3>{-2*c.box_x,-2*c.diameter,-2*c.box_z}),2,0,0.);
  walls.emplace_back(std::make_unique<IndicatorCuboid3D<T>>(
    Vector<T,3>{5*c.box_x,2*c.diameter,5*c.box_z},Vector<T,3>{-2*c.box_x,c.box_y,-2*c.box_z}),3,0,0.);
  return walls;
}
void wrapParticles(PS&ps,const Config&c){
  for(std::size_t i=0;i<ps.size();++i){auto p=ps.get(i);auto x=access::getPosition(p);
    for(int k:{0,2}){T L=k==0?c.box_x:c.box_z;x[k]-=std::floor(x[k]/L)*L;}
    if(!std::isfinite(x[0])||!std::isfinite(x[1])||!std::isfinite(x[2])||x[1]<0||x[1]>c.box_y)
      throw std::runtime_error("Nonfinite particle state or particle centre outside walls");
    p.setField<descriptors::GENERAL,descriptors::POSITION>(x);
  }
}
struct Measure {T bottom=0,top=0,fluid=0,eta=0,mach=0,drift=0,phi=0,mismatch=0;};
Measure measure(Problem&p,const Config&c,const Units&u){
  auto&l=p.getLattice(NavierStokes{});auto&g=p.getGeometry();
  l.setProcessingContext(ProcessingContext::Evaluation);
  std::array<T,8>a{},sum{};T mach=0,drift=0;int bad=0;
  const T pressure=u.rhoFluid*(c.dx/u.dt)*(c.dx/u.dt);
  for(int b=0;b<l.getLoadBalancer().size();++b){auto&bg=g.getBlockGeometry(b);auto&bl=l.getBlock(b);
    bg.forCoreSpatialLocations([&](LatticeR<3> loc){int m=bg.getMaterial(loc);if(m<1||m>3)return;
      auto cell=bl.get(loc);T rho,vel[3],pi[6];cell.computeAllMomenta(rho,vel,pi);
      T porosity=cell.template getField<descriptors::POROSITY>();
      T stress=-(1-u.omega/2)*pi[1]*pressure;
      if(!(rho>0&&std::isfinite(stress)&&std::isfinite(vel[0])&&std::isfinite(vel[1])&&std::isfinite(vel[2]))){bad=1;return;}
      mach=std::max(mach,std::sqrt(3*(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2])));drift=std::max(drift,std::abs(rho-1));
      if(m==2){a[0]+=stress;a[1]+=1;}if(m==3){a[2]+=stress;a[3]+=1;}
      if(m==1){if(porosity>1-1e-8){a[4]+=stress;a[5]+=1;}a[6]+=1-porosity;a[7]+=1;}
    });
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().allreduce(a.data(),sum.data(),8,MPI_SUM);
  singleton::mpi().reduceAndBcast(mach,MPI_MAX);singleton::mpi().reduceAndBcast(drift,MPI_MAX);singleton::mpi().reduceAndBcast(bad,MPI_MAX);
#else
  sum=a;
#endif
  if(bad||sum[1]==0||sum[3]==0)throw std::runtime_error("Nonfinite lattice state or missing wall");
  Measure m;m.bottom=sum[0]/sum[1];m.top=sum[2]/sum[3];m.fluid=sum[5]>0?sum[4]/sum[5]:0;
  m.eta=(m.top+m.bottom)/(2*c.shear_rate);m.mach=mach;m.drift=drift;m.phi=sum[6]*c.dx*c.dx*c.dx/(c.box_x*c.box_y*c.box_z);
  m.mismatch=std::abs(m.top-m.bottom)/std::max(std::abs((m.top+m.bottom)/2),c.dynamic_viscosity*c.shear_rate);
  return m;
}
void writePoses(std::ofstream&out,PS&ps,U64 step,T time){
  if(!singleton::mpi().isMainProcessor())return;
  for(std::size_t i=0;i<ps.size();++i){auto p=ps.get(i);out<<step<<','<<time<<','<<i;
    for(auto v:{access::getPosition(p),access::getAngle(p),access::getVelocity(p),p.getField<descriptors::MOBILITY,descriptors::ANG_VELOCITY>()})for(unsigned k=0;k<3;++k)out<<','<<v[k];
    out<<'\n';}out.flush();
}
template<class V>void put(std::ostream&out,const V&v){out.write(reinterpret_cast<const char*>(&v),sizeof(v));}
template<class V>void get(std::istream&in,V&v){in.read(reinterpret_cast<char*>(&v),sizeof(v));if(!in)throw std::runtime_error("Truncated checkpoint");}
std::string signature(const Config&c,const Units&u,std::size_t count){
  std::ostringstream s;s<<std::setprecision(17)<<"GRHLBM2 "<<singleton::mpi().getSize()<<' '<<count;
  for(T v:{c.shear_rate,c.box_x,c.box_y,c.box_z,c.dx,c.diameter,c.thickness,c.rho_particle,c.rho_fluid,c.dynamic_viscosity,c.nu_lattice,c.target_re,c.epsilon_cells,c.contact_enlargement_cells,c.contact_young,c.contact_poisson,c.contact_restitution,u.dt})s<<' '<<v;
  s<<' '<<c.contact_resolution;return s.str();
}
template<class Contact>void saveContacts(std::ostream&out,std::vector<Contact>&contacts){
  U64 n=contacts.size();put(out,n);std::vector<std::uint8_t>b(Contact::getSerialSize());
  for(auto&contact:contacts){contact.serialize(b.data());out.write(reinterpret_cast<char*>(b.data()),b.size());}
}
template<class Contact>void loadContacts(std::istream&in,std::vector<Contact>&contacts){
  U64 n;get(in,n);if(n>10000000)throw std::runtime_error("Invalid contact count");contacts.resize(n);
  std::vector<std::uint8_t>b(Contact::getSerialSize());for(auto&contact:contacts){in.read(reinterpret_cast<char*>(b.data()),b.size());if(!in)throw std::runtime_error("Truncated contacts");contact.deserialize(b.data());}
}
std::string checkpoint(Problem&p,PS&ps,CC&contacts,const Config&c,const Units&u,U64 step,const std::array<T,3>&contactDiagnostics){
  auto&l=p.getLattice(NavierStokes{});l.setProcessingContext(ProcessingContext::Evaluation);
  fs::path final=fs::path(c.output_dir)/("checkpoint_"+std::to_string(step));
  fs::path staging=final.string()+".partial";
  if(singleton::mpi().isMainProcessor())fs::create_directories(staging);
  singleton::mpi().barrier();
  if(!l.template save<false>((staging/"lattice").string()))throw std::runtime_error("Lattice checkpoint save failed");
  if(singleton::mpi().isMainProcessor()){
    std::ofstream metadata(staging/"metadata.txt");metadata<<signature(c,u,ps.size())<<'\n'<<step<<'\n'<<std::setprecision(17)<<contactDiagnostics[0]<<' '<<contactDiagnostics[1]<<' '<<contactDiagnostics[2]<<'\n';metadata.close();
    std::ofstream out(staging/"particles.bin",std::ios::binary);
    for(std::size_t i=0;i<ps.size();++i){auto particle=ps.get(i);U64 n=particle.getSerialSize();put(out,n);std::vector<std::uint8_t>b(n);particle.serialize(b.data());out.write(reinterpret_cast<char*>(b.data()),n);}
    saveContacts(out,contacts.particleContacts);saveContacts(out,contacts.wallContacts);out.close();
    if(!out)throw std::runtime_error("Particle checkpoint save failed");
    if(fs::exists(final))fs::remove_all(final);
    fs::rename(staging,final);
    std::ofstream latest(fs::path(c.output_dir)/"latest_checkpoint.txt.tmp");latest<<final.string()<<'\n';latest.close();
    fs::rename(fs::path(c.output_dir)/"latest_checkpoint.txt.tmp",fs::path(c.output_dir)/"latest_checkpoint.txt");
  }singleton::mpi().barrier();return final.string();
}
U64 restore(Problem&p,PS&ps,CC&contacts,const Config&c,const Units&u,std::array<T,3>&contactDiagnostics){
  std::ifstream metadata(fs::path(c.restart_dir)/"metadata.txt");std::string saved;std::getline(metadata,saved);U64 step;metadata>>step>>contactDiagnostics[0]>>contactDiagnostics[1]>>contactDiagnostics[2];
  if(!metadata||saved!=signature(c,u,ps.size()))throw std::runtime_error("Checkpoint incompatible with config, MPI ranks or particle count");
  std::ifstream in(fs::path(c.restart_dir)/"particles.bin",std::ios::binary);if(!in)throw std::runtime_error("Missing particle checkpoint");
  for(std::size_t i=0;i<ps.size();++i){auto particle=ps.get(i);auto*surface=particle.getField<descriptors::SURFACE,descriptors::SINDICATOR>();
    U64 n;get(in,n);if(n!=particle.getSerialSize())throw std::runtime_error("Particle ABI mismatch");
    std::vector<std::uint8_t>b(n);in.read(reinterpret_cast<char*>(b.data()),n);if(!in)throw std::runtime_error("Truncated particle checkpoint");
    particle.deserialize(b.data());
    // Serialized pointers are not valid across processes. Restore the newly created shape pointer immediately.
    particle.setField<descriptors::SURFACE,descriptors::SINDICATOR>(surface);
  }loadContacts(in,contacts.particleContacts);loadContacts(in,contacts.wallContacts);
  auto&l=p.getLattice(NavierStokes{});
  if(!l.template load<false>((fs::path(c.restart_dir)/"lattice").string()))throw std::runtime_error("Lattice checkpoint load failed");
  return step;
}
void simulate(const Config&c){
  Units u(c);OstreamManager log(std::cout,"graphiteCouette3d");
  if(singleton::mpi().isMainProcessor())fs::create_directories(c.output_dir);
  singleton::mpi().barrier();
  singleton::directories().setOutputDir(c.output_dir+"/vtk/");
  IndicatorCuboid3D<T> box(Vector<T,3>{c.box_x-c.dx,c.box_y,c.box_z-c.dx},Vector<T,3>(0.));
  Mesh<T,3> mesh(box,c.dx,singleton::mpi().getSize());mesh.setOverlap(2);mesh.getCuboidDecomposition().setPeriodicity(periodicity());
  Problem::ParametersD parameters;Problem problem(parameters,mesh);prepareGeometry(problem,c);prepareLattice(problem,c,u);
  auto&l=problem.getLattice(NavierStokes{});auto&g=problem.getGeometry();auto&converter=l.getUnitConverter();
  PS ps;addParticles(ps,c,u);CC contacts;auto walls=makeWalls(c);
  ContactProperties<T,1> cp;cp.set(0,0,evalEffectiveYoungModulus(c.contact_young,c.contact_young,c.contact_poisson,c.contact_poisson),1.,0.,0.);
  auto mapParticles=[&](){
    // Elastic frictionless contact has no tangential/damping history. Rebuild to avoid stale periodic contact boxes.
    contacts.clearContacts();
    graphiteCoupling::mapParticles(ps,contacts,g,l,converter,walls,periodicity);
  };
  std::array<T,3>contactDiagnostics{};
  U64 step=0;if(c.restart_dir.empty())mapParticles();else step=restore(problem,ps,contacts,c,u,contactDiagnostics);
  const long double desired=std::ceil(static_cast<long double>(c.end_strain)/(c.shear_rate*u.dt));
  if(desired>std::numeric_limits<U64>::max())throw std::runtime_error("end_strain requires more than uint64 steps");
  U64 endStep=static_cast<U64>(desired),startStep=step;
  U64 stopStep=c.max_steps?std::min(endStep,c.max_steps):endStep;
  if(stopStep<step)throw std::runtime_error("max_steps must be at least checkpoint step");
  log<<std::setprecision(12)<<"physical_Re="<<u.rePhysical<<" numerical_Re="<<u.reNumeric<<" numerical_St="<<u.stNumeric
    <<" inertia_scale="<<u.alpha<<" rho_fluid_numeric="<<u.rhoFluid<<" rho_particle_numeric="<<u.rhoParticle<<" dt_s="<<u.dt
    <<" particles="<<ps.size()<<" target_step="<<endStep<<std::endl;
  if(singleton::mpi().isMainProcessor()){
    std::ofstream meta(fs::path(c.output_dir)/"mapping.json");meta<<std::setprecision(17)<<"{\n\"dt_s\":"<<u.dt<<",\n\"inertia_scale\":"<<u.alpha
      <<",\n\"rho_fluid_numeric_kg_m3\":"<<u.rhoFluid<<",\n\"rho_particle_numeric_kg_m3\":"<<u.rhoParticle
      <<",\n\"Re_physical\":"<<u.rePhysical<<",\n\"Re_numeric\":"<<u.reNumeric<<",\n\"St_numeric\":"<<u.stNumeric<<",\n\"particle_count\":"<<ps.size()<<"\n}\n";
  }
  std::ofstream history,poses;
  if(singleton::mpi().isMainProcessor()){
    bool append=!c.restart_dir.empty()&&fs::exists(fs::path(c.output_dir)/"history.csv");
    history.open(fs::path(c.output_dir)/"history.csv",append?std::ios::app:std::ios::out);
    if(!append)history<<"step,time_s,strain,eta_wall_Pa_s,eta_relative,stress_bottom_Pa,stress_top_Pa,stress_fluid_mean_Pa,max_mach,density_drift,porosity_volume_fraction,wall_stress_mismatch,wall_seconds,steps_per_second,wall_contact_force_bottom_x_N,wall_contact_force_top_x_N,max_contact_indentation_m\n";
    poses.open(fs::path(c.output_dir)/"particles.csv",append?std::ios::app:std::ios::out);
    if(!append)poses<<"step,time_s,id,x_m,y_m,z_m,angle_x_rad,angle_y_rad,angle_z_rad,vx_m_s,vy_m_s,vz_m_s,omega_x_s_inv,omega_y_s_inv,omega_z_s_inv\n";
    history<<std::setprecision(17);poses<<std::setprecision(17);
  }
  SuperVTMwriter3D<T> writer("graphite");SuperLatticePhysVelocity3D<T,D>velocity(l,converter);
  SuperLatticePhysPressure3D<T,D>pressure(l,converter);SuperLatticePhysExternalPorosity3D<T,D>porosity(l,converter);
  writer.addFunctor(velocity);writer.addFunctor(pressure);writer.addFunctor(porosity);if(c.vtk_every)writer.createMasterFile();
  T&wallForceBottom=contactDiagnostics[0];T&wallForceTop=contactDiagnostics[1];T&maxIndentation=contactDiagnostics[2];
  const auto start=std::chrono::steady_clock::now();
  auto sample=[&](){auto m=measure(problem,c,u);
    const T area=c.box_x*c.box_z;
    m.bottom+=wallForceBottom/area;m.top-=wallForceTop/area;
    m.eta=(m.top+m.bottom)/(2*c.shear_rate);
    m.mismatch=std::abs(m.top-m.bottom)/std::max(std::abs((m.top+m.bottom)/2),c.dynamic_viscosity*c.shear_rate);
    T elapsed=std::chrono::duration<T>(std::chrono::steady_clock::now()-start).count();
    if(singleton::mpi().isMainProcessor()){
      history<<step<<','<<step*u.dt<<','<<step*u.dt*c.shear_rate<<','<<m.eta<<','<<m.eta/c.dynamic_viscosity<<','<<m.bottom<<','<<m.top<<','<<m.fluid
        <<','<<m.mach<<','<<m.drift<<','<<m.phi<<','<<m.mismatch<<','<<elapsed<<','<<(step-startStep)/std::max(elapsed,1e-12)<<','<<wallForceBottom<<','<<wallForceTop<<','<<maxIndentation<<'\n';history.flush();}
    writePoses(poses,ps,step,step*u.dt);
    log<<"step="<<step<<" strain="<<step*u.dt*c.shear_rate<<" eta_wall="<<m.eta<<" Pa.s wall_mismatch="<<m.mismatch<<" max_Mach="<<m.mach<<std::endl;
    if(m.mach>=.1||m.drift>.05)throw std::runtime_error("Mach >=0.1 or density drift >5%; reduce time step and investigate");
  };
  sample();std::string latest;
  while(step<stopStep){
    l.setProcessingContext(ProcessingContext::Evaluation);
    graphiteCoupling::coupleFluidToParticles(ps,g,l,converter,periodicity());
    wallForceBottom=wallForceTop=maxIndentation=0;
    auto contactDiagnostic=[&](const PhysR<T,3>&position,const Vector<T,3>&normal,const Vector<T,3>&normalForce,
        const Vector<T,3>&tangentialForce,const std::array<std::size_t,2>&ids,T overlapVolume,T indentation,bool isWall){
      (void)position;(void)normal;(void)overlapVolume;
      maxIndentation=std::max(maxIndentation,indentation);
      if(isWall){const T fx=-(normalForce[0]+tangentialForce[0]);if(ids[1]==0)wallForceBottom+=fx;else wallForceTop+=fx;}
    };
    processContacts<T,P,PC,WC,ContactProperties<T,1>>(ps,walls,contacts,cp,g,c.contact_resolution,T(4/(3*std::sqrt(M_PI))),periodicity,contactDiagnostic);
    ps.process(u.dt);wrapParticles(ps,c);
    l.setProcessingContext(ProcessingContext::Simulation);l.collideAndStream();
    l.setProcessingContext(ProcessingContext::Evaluation);mapParticles();++step;
    if(step%c.sample_every==0||step==stopStep)sample();
    if(c.vtk_every&&step%c.vtk_every==0){const U64 frame=step/c.vtk_every;if(frame>std::numeric_limits<int>::max())throw std::runtime_error("Too many VTK frames");writer.write(static_cast<int>(frame));}
    int stop=stopRequested;
    if((step%c.sample_every==0||step==stopStep)&&singleton::mpi().isMainProcessor()&&fs::exists(fs::path(c.output_dir)/"STOP_REQUEST"))stop=1;
#ifdef PARALLEL_MODE_MPI
    // Synchronize signals at sample cadence to avoid an extra collective every time step.
    if(step%c.sample_every==0||step==stopStep)singleton::mpi().reduceAndBcast(stop,MPI_MAX);else stop=0;
#endif
    if((c.checkpoint_every&&step%c.checkpoint_every==0)||stop){latest=checkpoint(problem,ps,contacts,c,u,step,contactDiagnostics);if(stop)break;}
  }
  if(latest.empty()||latest!=fs::path(c.output_dir)/("checkpoint_"+std::to_string(step)))latest=checkpoint(problem,ps,contacts,c,u,step,contactDiagnostics);
  if(singleton::mpi().isMainProcessor()){
    std::ofstream status(fs::path(c.output_dir)/"status.json");status<<std::setprecision(17)<<"{\"status\":\""<<(step>=endStep?"COMPLETED":"CHECKPOINTED")
      <<"\",\"step\":"<<step<<",\"time_s\":"<<step*u.dt<<",\"strain\":"<<step*u.dt*c.shear_rate<<",\"checkpoint_dir\":\""<<latest<<"\"}\n";
  }
}
int runCase(int argc,char**argv){initialize(&argc,&argv);std::signal(SIGTERM,requestStop);std::signal(SIGINT,requestStop);std::signal(SIGUSR1,requestStop);
  try{simulate(parse(argc,argv));return 0;}catch(const std::exception&e){std::cerr<<"graphiteCouette3d: "<<e.what()<<std::endl;
#ifdef PARALLEL_MODE_MPI
    MPI_Abort(MPI_COMM_WORLD,2);
#endif
    return 2;}
}

} } // SLURRY SCOPE END
