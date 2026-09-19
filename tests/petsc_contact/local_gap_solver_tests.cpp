// Local adhesion must survive diagnostic replay and respect the domain of both
// Newton backends.  This test can also run without PETSc or OpenLB.
#include "particleSubsteps.h"
#include <chrono>
#include <filesystem>
#include <iostream>

namespace g=slurry::gr_re2::graphite;
namespace d=g::particle_detail;

namespace {
void require(bool condition,const char* message) {
  if(!condition)throw std::runtime_error(message);
}
g::ParticleStepSettings settings() {
  g::ParticleStepSettings s;
  s.pair.sigma=4.197e-10;
  s.pair.roughnessGap=s.rough.gap=2.e-9;
  s.pair.localGap=3.e-10;s.pair.localGapFraction=1.;
  s.pair.localSwitchExcessGap=2.e-9;s.pair.localCutoffExcessGap=10.e-9;
  s.rough.enabled=true;s.rough.tangentialStiffness=80.;
  s.shearRate=0.;s.box={100.e-6,100.e-6,100.e-6};
  return s;
}
std::vector<g::Body> bodies(const g::ParticleStepSettings& s) {
  std::vector<g::Body> b(2);
  b[1].position[2]=b[0].axes[2]+b[1].axes[2]+s.rough.gap;
  return b;
}
struct Scratch {
  std::filesystem::path path;
  Scratch() {
    path=std::filesystem::temp_directory_path()/
        ("gr_local_gap_replay_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directory(path);
  }
  ~Scratch(){std::error_code error;std::filesystem::remove_all(path,error);}
};

void trialDomainAndBirth() {
  auto s=settings();auto b=bodies(s);
  require(d::admissibleTrialGap(s.rough.gap,s),"Contact gap must remain admissible");
  require(d::admissibleTrialGap(1.75e-9,s),"Valid positive local trial gap rejected");
  require(!d::admissibleTrialGap(1.65e-9,s),"Negative local trial gap accepted");
  require(!d::admissibleTrialGap(g::minimumPairGap(s.pair),s),"Singular local trial gap accepted");
  auto legacy=s;legacy.pair.localGapFraction=0.;
  require(d::admissibleTrialGap(1.65e-9,legacy),"Disabled local correction changed legacy domain");
  require(!d::admissibleTrialGap(1.5e-9,legacy),"Legacy rough-contact domain changed");

  std::vector<g::Vec3> force(2),torque(2);std::vector<g::GapCache> cache(1);
  g::PersistentContactState contacts(1);std::vector<int> slots{0};
  const double L=b[0].axes[0],dt=1.e-8;
  d::Residual residual{b,force,torque,s,cache,contacts,slots,dt,0.,L,1.,1};
  d::Vector q(13,0.);d::Evaluation e;std::string error;
  const auto pair=g::evaluatePair(b[0],b[1],s.pair);
  const double adhesion=std::max(0.,g::dot(pair.forceI,pair.normal));
  q[12]=adhesion;
  require(residual(q,e,error),"Local contact residual failed at the physical contact gap");
  require(e.contacts[0].active,"New local adhesive contact was not activated");
  require(std::abs(e.contacts[0].rollingCap-s.rough.rollingLength*adhesion)
      <1.e-10*s.rough.rollingLength*adhesion,"Rolling birth did not use corrected adhesion");
  require(!contacts[0].active,"Residual mutated committed contact history");

  q[8]=(1.65e-9-s.rough.gap)/L;
  require(!residual(q,e,error),"Invalid local trial gap escaped legacy domain rejection");
  require(!error.empty(),"Rejected local trial lacks a domain diagnostic");
#ifdef SLURRY_USE_PETSC
  d::PetscParticleContext context(residual);
  require(!context.evaluate(q,e),"Invalid local trial gap escaped PETSc domain rejection");
  require(context.domainErrors==1,"PETSc did not record the local domain error");
#endif

  auto mismatch=s;mismatch.rough.gap=3.e-9;bool rejected=false;
  try {g::validateParticlePairSettings(mismatch);}catch(const std::exception&){rejected=true;}
  require(rejected,"Mismatched contact and local adhesion geometry accepted");
}

void replayVersions() {
  Scratch scratch;d::ParticleReplayInput input;
  input.settings=settings();input.dt=input.outerDt=1.e-8;input.bodies=bodies(input.settings);
  input.force.resize(2);input.torque.resize(2);input.contacts.resize(1);input.cache.resize(1);
  input.contacts[0].active=true;input.contacts[0].normalLoad=9.36e-7;
  input.contacts[0].rollingCap=9.36e-14;
  const auto v2=scratch.path/"v2.dat";d::writeParticleReplay(input,v2.string());
  const auto read=d::readParticleReplay(v2.string());
  const auto& a=input.settings.pair;const auto& b=read.settings.pair;
  require(a.roughnessGap==b.roughnessGap&&a.localGap==b.localGap
      &&a.localGapFraction==b.localGapFraction&&a.localSwitchExcessGap==b.localSwitchExcessGap
      &&a.localCutoffExcessGap==b.localCutoffExcessGap,"Replay lost local adhesion parameters");
  require(read.settings.rough.tangentialStiffness==80.,"Replay lost configured tangential stiffness");
  require(read.contacts[0].rollingCap==input.contacts[0].rollingCap,"Replay lost frozen rolling cap");
  const auto expected=g::evaluatePair(input.bodies[0],input.bodies[1],a);
  const auto actual=g::evaluatePair(read.bodies[0],read.bodies[1],b);
  require(expected.forceI==actual.forceI&&expected.energy==actual.energy,"Replay changed local pair interaction");

  // The v1 layout is the same payload without the added local-parameter line.
  std::ifstream source(v2);const auto v1=scratch.path/"v1.dat";std::ofstream target(v1);
  std::string line;int index=0;
  while(std::getline(source,line)) {
    if(index==0)target<<"GR_PARTICLE_REPLAY 1\n";
    else if(index!=5)target<<line<<'\n';
    ++index;
  }
  target.close();const auto old=d::readParticleReplay(v1.string());
  require(old.settings.pair.localGapFraction==0.,"Legacy replay unexpectedly enabled local adhesion");
  require(old.contacts[0].normalLoad==input.contacts[0].normalLoad,"Legacy replay contact data shifted");
}
}

int main() {
  try {trialDomainAndBirth();replayVersions();std::cout<<"Local gap solver/replay tests passed\n";return 0;}
  catch(const std::exception& error){std::cerr<<"Local gap solver/replay test failed: "<<error.what()<<'\n';return 1;}
}
