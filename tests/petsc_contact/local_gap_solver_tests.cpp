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
  auto free=s;free.pair.localGapFraction=0.;free.pair.freeCmcRepulsionPressure=52000.;
  g::validateParticlePairSettings(free);
  free.rough.enabled=false;rejected=false;
  try {g::validateParticlePairSettings(free);}catch(const std::exception&){rejected=true;}
  require(rejected,"Free-CMC repulsion accepted without rough contact");
  free.rough.enabled=true;free.rough.gap=3.e-9;rejected=false;
  try {g::validateParticlePairSettings(free);}catch(const std::exception&){rejected=true;}
  require(rejected,"Free-CMC repulsion accepted a mismatched roughness gap");
}

void replayVersions() {
  Scratch scratch;d::ParticleReplayInput input;
  input.settings=settings();input.dt=input.outerDt=1.e-8;input.bodies=bodies(input.settings);
  input.force.resize(2);input.torque.resize(2);input.contacts.resize(1);input.cache.resize(1);
  input.contacts[0].active=true;input.contacts[0].normalLoad=9.36e-7;
  input.contacts[0].rollingCap=9.36e-14;
  const auto v7=scratch.path/"v7.dat";d::writeParticleReplay(input,v7.string());
  const auto read=d::readParticleReplay(v7.string());
  const auto& a=input.settings.pair;const auto& b=read.settings.pair;
  require(a.roughnessGap==b.roughnessGap&&a.localGap==b.localGap
      &&a.localGapFraction==b.localGapFraction&&a.localSwitchExcessGap==b.localSwitchExcessGap
      &&a.localCutoffExcessGap==b.localCutoffExcessGap,"Replay lost local adhesion parameters");
  require(read.settings.rough.tangentialStiffness==80.,"Replay lost configured tangential stiffness");
  require(read.contacts[0].rollingCap==input.contacts[0].rollingCap,"Replay lost frozen rolling cap");
  require(read.settings.passMax==input.settings.passMax,"Replay lost maximum-iteration acceptance policy");
  const auto expected=g::evaluatePair(input.bodies[0],input.bodies[1],a);
  const auto actual=g::evaluatePair(read.bodies[0],read.bodies[1],b);
  require(expected.forceI==actual.forceI&&expected.energy==actual.energy,"Replay changed local pair interaction");

  require(!b.surfaceAdhesion,"Legacy replay unexpectedly enabled surface adhesion");
  // v1-v6 omit net CMC; v1-v5 omit inner CMC; v1-v3 omit outer CMC;
  // v1-v2 omit surface parameters;
  // v1 additionally omits the legacy local line.
  for(int version:{1,2,3,4,5,6}) {
    std::ifstream source(v7);const auto oldPath=scratch.path/("v"+std::to_string(version)+".dat");
    std::ofstream target(oldPath);std::string line;int index=0;
    while(std::getline(source,line)) {
      if(index==0)target<<"GR_PARTICLE_REPLAY "<<version<<'\n';
      else if(index!=10&&(version>=6||index!=9)&&(version>=5||index!=8)&&(version>=4||index!=7)
          &&(version>=3||index!=6)&&(version!=1||index!=5))target<<line<<'\n';
      ++index;
    }
    target.close();const auto old=d::readParticleReplay(oldPath.string());
    require(!old.settings.pair.surfaceAdhesion,"Legacy replay unexpectedly enabled surface adhesion");
    require(old.settings.pair.freeCmcRepulsionPressure==0.,"Legacy replay unexpectedly enabled free-CMC repulsion");
    require(old.settings.pair.freeCmcInnerRepulsionWork==0.,"Legacy replay unexpectedly enabled inner CMC repulsion");
    require(old.settings.pair.cmcNetBlend==0.,"Legacy replay unexpectedly enabled net CMC potential");
    require(old.settings.passMax==(version>=5?input.settings.passMax:0),
        "Legacy failure replay must retain its maximum-iteration policy");
    require(old.settings.pair.localGapFraction==(version==1?0.:a.localGapFraction),
            "Legacy replay changed local adhesion selection");
    require(old.contacts[0].normalLoad==input.contacts[0].normalLoad,"Legacy replay contact data shifted");
  }

  input.settings.pair.localGapFraction=0.;
  input.settings.pair.surfaceAdhesion=true;
  input.settings.pair.adhesionWork=.018;
  input.settings.pair.adhesionRange=5.e-10;
  input.settings.pair.curvatureSwitchGap=4.e-9;
  input.settings.pair.curvatureCutoffGap=18.e-9;
  input.settings.pair.freeCmcRepulsionPressure=52535.69484753685;
  input.settings.pair.freeCmcRepulsionLength=5.e-9;
  input.settings.pair.freeCmcInnerRepulsionWork=.007;
  input.settings.pair.freeCmcInnerRepulsionRange=4.e-10;
  input.settings.pair.freeCmcInnerRepulsionPower=2.1;
  const auto surface=scratch.path/"surface.dat";d::writeParticleReplay(input,surface.string());
  const auto restored=d::readParticleReplay(surface.string());const auto& p=restored.settings.pair;
  require(p.surfaceAdhesion&&p.adhesionWork==.018&&p.adhesionRange==5.e-10
      &&p.curvatureSwitchGap==4.e-9&&p.curvatureCutoffGap==18.e-9,
      "Replay lost surface adhesion parameters");
  require(p.freeCmcRepulsionPressure==input.settings.pair.freeCmcRepulsionPressure
      &&p.freeCmcRepulsionLength==input.settings.pair.freeCmcRepulsionLength,
      "Replay lost free-CMC repulsion parameters");
  require(p.freeCmcInnerRepulsionWork==input.settings.pair.freeCmcInnerRepulsionWork
      &&p.freeCmcInnerRepulsionRange==input.settings.pair.freeCmcInnerRepulsionRange
      &&p.freeCmcInnerRepulsionPower==input.settings.pair.freeCmcInnerRepulsionPower,
      "Replay lost inner CMC repulsion parameters");
  const auto before=g::evaluatePair(input.bodies[0],input.bodies[1],input.settings.pair);
  const auto after=g::evaluatePair(restored.bodies[0],restored.bodies[1],p);
  require(before.forceI==after.forceI&&before.energy==after.energy,
          "Replay changed surface adhesion interaction");
  // A genuine v6 inner-CMC file must retain its nonzero terms and use net blend 0.
  const auto oldInner=scratch.path/"inner_v6.dat";
  {
    std::ifstream source(surface);std::ofstream target(oldInner);std::string line;int index=0;
    while(std::getline(source,line)) {
      if(index==0)target<<"GR_PARTICLE_REPLAY 6\n";
      else if(index!=10)target<<line<<'\n';
      ++index;
    }
  }
  const auto oldInnerRead=d::readParticleReplay(oldInner.string());
  const auto oldInnerPair=g::evaluatePair(oldInnerRead.bodies[0],oldInnerRead.bodies[1],oldInnerRead.settings.pair);
  require(oldInnerRead.settings.pair.cmcNetBlend==0.
      &&oldInnerRead.settings.pair.freeCmcInnerRepulsionWork==p.freeCmcInnerRepulsionWork
      &&oldInnerPair.forceI==before.forceI&&oldInnerPair.energy==before.energy,
      "Version 6 inner-CMC replay changed its original interaction");

  input.settings.pair.cmcNetBlend=.63;
  input.settings.pair.cmcNetContactForce=1.7e-10;
  input.settings.pair.cmcNetBarrierForce=2.3e-11;
  input.settings.pair.cmcNetAttractionRange=1.1e-9;
  input.settings.pair.cmcNetRepulsionRange=6.3e-9;
  input.settings.pair.cmcNetReferenceLength=14.e-6;
  const auto net=scratch.path/"net.dat";d::writeParticleReplay(input,net.string());
  const auto netRead=d::readParticleReplay(net.string());
  const auto& netParameters=netRead.settings.pair;
  require(netParameters.cmcNetBlend==input.settings.pair.cmcNetBlend
      &&netParameters.cmcNetContactForce==input.settings.pair.cmcNetContactForce
      &&netParameters.cmcNetBarrierForce==input.settings.pair.cmcNetBarrierForce
      &&netParameters.cmcNetAttractionRange==input.settings.pair.cmcNetAttractionRange
      &&netParameters.cmcNetRepulsionRange==input.settings.pair.cmcNetRepulsionRange
      &&netParameters.cmcNetReferenceLength==input.settings.pair.cmcNetReferenceLength,
      "Replay lost net-CMC potential parameters");
  const auto netBefore=g::evaluatePair(input.bodies[0],input.bodies[1],input.settings.pair);
  const auto netAfter=g::evaluatePair(netRead.bodies[0],netRead.bodies[1],netParameters);
  require(netBefore.forceI==netAfter.forceI&&netBefore.energy==netAfter.energy,
      "Replay changed the blended net-CMC interaction");
  input.settings.pair.cmcNetBlend=0.;
  input.settings.pair.freeCmcRepulsionPressure=0.;
  input.settings.pair.freeCmcInnerRepulsionWork=0.;
  input.settings.pair.contactGap=input.settings.rough.gap=3.e-9;
  input.settings.pair.cohesionRetention=.2;input.settings.passMax=1;
  const auto coated=scratch.path/"coated.dat";d::writeParticleReplay(input,coated.string());
  const auto coatedRead=d::readParticleReplay(coated.string());
  require(coatedRead.settings.pair.contactGap==3.e-9
      &&coatedRead.settings.pair.cohesionRetention==.2&&coatedRead.settings.passMax==1,
      "Replay lost coated-contact law or maximum-iteration policy");
}

void netContactBirth() {
  auto s=settings();s.pair.localGapFraction=0.;s.pair.surfaceAdhesion=true;s.pair.cmcNetBlend=1.;
  auto b=bodies(s);
  std::vector<g::Vec3> force(2),torque(2);std::vector<g::GapCache> cache(1);
  g::PersistentContactState contacts(1);std::vector<int> slots{0};
  const double L=b[0].axes[0];
  d::Residual residual{b,force,torque,s,cache,contacts,slots,1.e-8,0.,L,1.,1};
  const auto pair=g::evaluatePair(b[0],b[1],s.pair);
  const double adhesion=std::max(0.,g::dot(pair.forceI,pair.normal));
  d::Vector q(13,0.);q[12]=adhesion;d::Evaluation e;std::string error;
  require(residual(q,e,error),"Net-CMC contact residual failed at h0");
  require(adhesion>0.&&e.contacts[0].active,"Net-CMC attraction did not permit true contact");
  require(std::abs(e.contacts[0].rollingCap-s.rough.rollingLength*adhesion)
      <=1.e-10*s.rough.rollingLength*adhesion,"Net-CMC rolling birth used a legacy force");
  require(!contacts[0].active,"Net-CMC residual mutated committed contact history");
  s.rough.gap=3.e-9;bool rejected=false;
  try {g::validateParticlePairSettings(s);}catch(const std::exception&){rejected=true;}
  require(rejected,"Net-CMC potential accepted a different physical contact plane");
}
}

int main() {
  try {trialDomainAndBirth();replayVersions();netContactBirth();std::cout<<"Local gap solver/replay tests passed\n";return 0;}
  catch(const std::exception& error){std::cerr<<"Local gap solver/replay test failed: "<<error.what()<<'\n';return 1;}
}
