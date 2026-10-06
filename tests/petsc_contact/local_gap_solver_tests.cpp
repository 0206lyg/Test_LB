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
  const auto v8=scratch.path/"v8.dat";d::writeParticleReplay(input,v8.string());
  const auto read=d::readParticleReplay(v8.string());
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
  // v1-v7 omit coordination; v1-v6 omit net CMC; v1-v5 omit inner CMC;
  // v1-v3 omit outer CMC;
  // v1-v2 omit surface parameters;
  // v1 additionally omits the legacy local line.
  for(int version:{1,2,3,4,5,6,7}) {
    std::ifstream source(v8);const auto oldPath=scratch.path/("v"+std::to_string(version)+".dat");
    std::ofstream target(oldPath);std::string line;int index=0;
    while(std::getline(source,line)) {
      if(index==0)target<<"GR_PARTICLE_REPLAY "<<version<<'\n';
      else if(index!=11&&(version>=7||index!=10)&&(version>=6||index!=9)&&(version>=5||index!=8)&&(version>=4||index!=7)
          &&(version>=3||index!=6)&&(version!=1||index!=5))target<<line<<'\n';
      ++index;
    }
    target.close();const auto old=d::readParticleReplay(oldPath.string());
    require(!old.settings.pair.surfaceAdhesion,"Legacy replay unexpectedly enabled surface adhesion");
    require(old.settings.pair.freeCmcRepulsionPressure==0.,"Legacy replay unexpectedly enabled free-CMC repulsion");
    require(old.settings.pair.freeCmcInnerRepulsionWork==0.,"Legacy replay unexpectedly enabled inner CMC repulsion");
    require(old.settings.pair.cmcNetBlend==0.,"Legacy replay unexpectedly enabled net CMC potential");
    require(!old.settings.pair.cmcCoordinationEnabled,"Legacy replay unexpectedly enabled coordination");
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
      else if(index!=10&&index!=11)target<<line<<'\n';
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
  const auto oldNet=scratch.path/"net_v7.dat";
  {
    std::ifstream source(net);std::ofstream target(oldNet);std::string line;int index=0;
    while(std::getline(source,line)) {
      if(index==0)target<<"GR_PARTICLE_REPLAY 7\n";
      else if(index!=11)target<<line<<'\n';
      ++index;
    }
  }
  const auto oldNetRead=d::readParticleReplay(oldNet.string());
  const auto oldNetPair=g::evaluatePair(oldNetRead.bodies[0],oldNetRead.bodies[1],oldNetRead.settings.pair);
  require(!oldNetRead.settings.pair.cmcCoordinationEnabled
      &&!oldNetRead.settings.rough.currentAdhesionRolling
      &&oldNetPair.forceI==netBefore.forceI&&oldNetPair.energy==netBefore.energy,
      "Version 7 replay changed the original net interaction or frozen rolling law");
  input.settings.pair.cmcCoordinationEnabled=true;
  input.settings.pair.cmcCoordinationStart=1.25;
  input.settings.pair.cmcCoordinationEnd=2.75;
  input.settings.pair.cmcCoordinationFloor=.13;
  const auto coordinated=scratch.path/"coordinated.dat";
  d::writeParticleReplay(input,coordinated.string());
  const auto coordinatedRead=d::readParticleReplay(coordinated.string());
  const auto& cp=coordinatedRead.settings.pair;
  require(cp.cmcCoordinationEnabled&&cp.cmcCoordinationStart==1.25
      &&cp.cmcCoordinationEnd==2.75&&cp.cmcCoordinationFloor==.13
      &&coordinatedRead.settings.rough.currentAdhesionRolling,
      "Replay lost coordination settings or the derived current-strength rolling law");
  input.settings.pair.cmcCoordinationEnabled=false;
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

void coordinatedResidualAndSnapshot() {
  auto s=settings();s.pair.localGapFraction=0.;s.pair.surfaceAdhesion=true;
  s.pair.cmcNetBlend=1.;s.pair.cmcCoordinationEnabled=true;s.nearField.enabled=false;
  s.box={10.e-6,10.e-6,10.e-6};s.shearRate=10.;
  constexpr double advancedTime=.099999999,dt=1.e-8,pi=3.14159265358979323846;
  constexpr double radius=.5e-6;
  const auto star=[&](double gap) {
    std::vector<g::Body> b(4);
    for(auto& body:b)body.axes={radius,radius,radius};
    b[0].position={5.e-6,9.7e-6,5.e-6};
    for(int i=1;i<4;++i) {
      const double angle=2.*pi*(i-1)/3.;
      b[i].position=g::add(b[0].position,{(2.*radius+gap)*std::cos(angle),
          (2.*radius+gap)*std::sin(angle),0.});
    }
    return b;
  };
  auto unwrapped=star(s.rough.gap+.35e-9),b=unwrapped;
  for(auto& body:b){d::wrap(body,advancedTime,s);body.velocity={};}
  std::vector<g::Vec3> zero(4);std::vector<g::GapCache> cache(6);
  g::PersistentContactState history(6);std::vector<int> slots(6,-1);
  d::Residual residual{b,zero,zero,s,cache,history,slots,dt,advancedTime-dt,radius,1.,0};
  d::Vector q(24,0.);d::Evaluation base;std::string error;
  require(residual(q,base,error),"Coordinated many-body residual failed near LE phase wrap");
  const auto physicalForce=[&](const d::Evaluation& e,std::size_t i) {
    const double scale=b[i].mass*radius/(dt*dt);
    return g::Vec3{-scale*e.residual[6*i],-scale*e.residual[6*i+1],-scale*e.residual[6*i+2]};
  };
  g::Vec3 totalForce{};
  for(std::size_t i=0;i<b.size();++i)totalForce=g::add(totalForce,physicalForce(base,i));
  require(g::norm(totalForce)<1.e-23,"Coordination residual lost equal-and-opposite force balance");
  const double epsilon=2.e-13;
  for(std::size_t i=0;i<b.size();++i)for(int axis=0;axis<2;++axis) {
    auto plus=q,minus=q;plus[6*i+axis]+=epsilon/radius;minus[6*i+axis]-=epsilon/radius;
    d::Evaluation ep,em;
    require(residual(plus,ep,error)&&residual(minus,em,error),"Coordination energy derivative left trial domain");
    const double force=physicalForce(base,i)[axis];
    const double derivative=(ep.diagnostic.energyAtEnd-em.diagnostic.energyAtEnd)/(2.*epsilon);
    require(std::abs(force+derivative)<1.e-18+2.e-5*std::abs(force),
        "Many-body residual force is not the gradient of its complete energy");
  }
  auto movedQ=q;movedQ[6]+=epsilon/radius;d::Evaluation moved;
  require(residual(movedQ,moved,error),"Neighbour perturbation residual failed");
  require(g::norm(g::sub(physicalForce(moved,2),physicalForce(base,2)))>1.e-18,
      "Residual omitted the force response through a third particle's coordination");
  const auto snapshot=g::evaluateParticleState(b,advancedTime,s);
  auto nonSheared=s;nonSheared.shearRate=0.;
  const auto reference=g::evaluateParticleState(unwrapped,0.,nonSheared);
  const auto close=[](double a,double c){return std::abs(a-c)<1.e-28+1.e-8*std::max(std::abs(a),std::abs(c));};
  require(close(base.diagnostic.energyAtEnd,snapshot.energyAtEnd)
      &&close(snapshot.energyAtEnd,reference.energyAtEnd),"Snapshot coordination energy disagrees with residual or LE images");
  for(int k=0;k<9;++k)
    require(close(base.diagnostic.pairMoment[k],snapshot.pairMoment[k])
        &&close(snapshot.pairMoment[k],reference.pairMoment[k]),
        "Snapshot coordination virial uses a different force or periodic branch");

  // Existing contacts must follow the corrected current strength as the
  // neighbourhood changes, including a previously frozen stronger rolling cap.
  b=star(s.rough.gap);s.shearRate=0.;
  cache.assign(6,{});history.assign(6,{});slots={0,1,2,-1,-1,-1};
  const auto corrected=d::coordinatedPairs(b,0.,s,cache);
  d::Vector contactQ(27,0.);
  for(std::size_t p=0;p<3;++p) {
    const auto isolated=g::evaluatePair(b[0],b[p+1],s.pair);
    auto& old=history[p];old.active=true;old.normal=isolated.normal;
    old.rollingCap=s.rough.rollingLength*std::max(0.,g::dot(isolated.forceI,isolated.normal));
    old.rollingStiffness=old.rollingCap/s.rough.rollingYieldAngle;
    old.elasticRoll={0.,0.,.003};
    old.elasticEnergy=.5*old.rollingStiffness*g::dot(old.elasticRoll,old.elasticRoll);
    const auto& current=corrected[p].result;
    contactQ[24+p]=std::max(0.,g::dot(current.forceI,current.normal));
  }
  d::Residual contactResidual{b,zero,zero,s,cache,history,slots,dt,0.,radius,1.,3};
  d::Evaluation contactState;
  require(contactResidual(contactQ,contactState,error),"Coordinated existing-contact residual failed");
  for(std::size_t p=0;p<3;++p) {
    const double expected=s.rough.rollingLength*contactQ[24+p];
    require(contactState.contacts[p].active&&expected<history[p].rollingCap
        &&std::abs(contactState.contacts[p].rollingCap-expected)<1.e-28+1.e-9*expected,
        "Existing rolling cap did not follow the current coordinated force");
  }
}
}

int main() {
  try {trialDomainAndBirth();replayVersions();netContactBirth();coordinatedResidualAndSnapshot();
    std::cout<<"Local gap solver/replay tests passed\n";return 0;}
  catch(const std::exception& error){std::cerr<<"Local gap solver/replay test failed: "<<error.what()<<'\n';return 1;}
}
