// Standalone checkpoint migration checks, without an OpenLB/MPI installation.
/* c++ -std=c++17 -O2 -Wall -Wextra -pedantic -I olb-1.9r0/src/slurry/gr_re2 \
     tests/checkpoint/state_version_tests.cpp -o /tmp/state_version_tests
   /tmp/state_version_tests */
#include <cstddef>
namespace olb {
class Serializer {
public:
  template<class T> explicit Serializer(T&) {}
  void computeSize() {}
  std::size_t getSize() const { return 0; }
  bool* getNextBlock(std::size_t&,bool) { return nullptr; }
};
namespace singleton {
struct SerialMpi {
  int getRank() const { return 0; }
  int getSize() const { return 1; }
  void barrier() const {}
};
inline SerialMpi& mpi() { static SerialMpi instance;return instance; }
}
}
#include "coupledCheckpoint.h"
#include <chrono>
#include <functional>
#include <iostream>

namespace c=slurry::gr_re2;
namespace cp=c::checkpoint_detail;
namespace g=c::graphite;

namespace {
void require(bool condition,const char* message) {
  if(!condition)throw std::runtime_error(message);
}
void rejects(const std::function<void()>& operation,const char* message) {
  try {operation();}catch(const std::runtime_error&){return;}
  throw std::runtime_error(message);
}
struct Scratch {
  cp::fs::path path=cp::fs::temp_directory_path()/
      ("slurry-checkpoint-state-"+std::to_string(
          std::chrono::high_resolution_clock::now().time_since_epoch().count()));
  Scratch(){cp::fs::create_directory(path);}
  ~Scratch(){std::error_code error;cp::fs::remove_all(path,error);}
};
// Independent description of the old on-disk order, frozen before the new
// diagnostic field.  Arrays keep this fixture independent of migration's names.
struct FrozenV1Diagnostic {
  double scalars[4]{};
  double moments[6][9]{};
  double moreScalars[5]{};
  int counts[12]{};
  std::size_t activePairs=0;
};
void legacyFile(const cp::fs::path& path,cp::U64 diagnosticSize=sizeof(FrozenV1Diagnostic),
                cp::U64 version=1) {
  cp::Output out(path);
  out.put(cp::U64(0x4752535441544531ULL));out.put(version);out.put(cp::byteOrder);
  for(cp::U64 size:{sizeof(g::Body),sizeof(g::GapCache),sizeof(g::RoughContactState),diagnosticSize})
    out.put(size);
  const std::string signature="old-pure-gr-signature";
  out.put(cp::U64(signature.size()));out.bytes(signature.data(),signature.size());
  out.put(cp::U64(42));out.put(cp::U64(1));
  FrozenV1Diagnostic d;
  for(int i=0;i<4;++i)d.scalars[i]=i+.125;
  for(int i=0;i<6;++i)for(int j=0;j<9;++j)d.moments[i][j]=10*i+j+.25;
  for(int i=0;i<5;++i)d.moreScalars[i]=i+.5;
  for(int i=0;i<12;++i)d.counts[i]=i+1;
  d.activePairs=13;out.put(d);
  out.put(std::array<double,6>{32.,33.,34.,35.,36.,37.});
  g::Body body;body.position={1.,2.,3.};out.put(body);
  out.put(g::Vec3{.1,.2,.3});out.finish();
}
void migrateLegacyAndRoundtrip(const cp::fs::path& directory) {
  const auto path=directory/"legacy.bin";
  legacyFile(path);
  auto state=cp::loadState(path,1,"old-pure-gr-signature");
  const auto& d=state.diagnostic;
  require(state.step==42&&state.maxIterationPassesTotal==0,"Legacy cumulative counter migration");
  require(d.maxIterationPasses==0,"Legacy per-step counter migration");
  require(d.minGap==.125&&d.maxForce==1.125&&d.energyAtEnd==2.125
      &&d.lubricationDissipation==3.125,"Legacy leading diagnostic scalars");
  const std::array<g::Mat3,6> moments{d.pairMoment,d.attractiveMoment,d.repulsiveMoment,
      d.lubricationMoment,d.contactNormalMoment,d.contactTangentialMoment};
  for(int i=0;i<6;++i)for(int j=0;j<9;++j)
    require(moments[i][j]==10*i+j+.25,"Legacy stress moments");
  require(d.contactDissipation==.5&&d.elasticContactEnergy==1.5
      &&d.maxForceResidualRatio==2.5&&d.maxTorqueResidualRatio==3.5
      &&d.contactGapViolation==4.5,"Legacy trailing diagnostic scalars");
  const std::array<int,12> counts{d.contacts,d.slidingContacts,d.rollingContacts,
      d.substeps,d.newtonIterations,d.krylovIterations,d.residualEvaluations,
      d.frictionBranchAttempts,d.frictionBranchCorrections,d.contactStateUpdates,
      d.contactActivations,d.contactReleases};
  for(int i=0;i<12;++i)require(counts[i]==i+1,"Legacy solver/contact counts");
  require(d.activePairs==13,"Legacy active pair count");
  require(state.timing==std::array<double,6>{32.,33.,34.,35.,36.,37.},"Legacy timing offset");
  require(state.bodies[0].position==g::Vec3{1.,2.,3.},"Legacy body offset");
  require(state.angularAcceleration[0]==g::Vec3{.1,.2,.3},"Legacy acceleration offset");

  state.diagnostic.maxIterationPasses=3;
  state.maxIterationPassesTotal=(cp::U64(1)<<40)+7;
  const auto current=directory/"current.bin";
  cp::saveState(current,state,"old-pure-gr-signature");
  const auto restored=cp::loadState(current,1,"old-pure-gr-signature");
  require(restored.diagnostic.maxIterationPasses==3,"Current per-step counter roundtrip");
  require(restored.maxIterationPassesTotal==state.maxIterationPassesTotal,
      "Current cumulative counter roundtrip must preserve 64-bit counts");
  require(restored.step==state.step&&restored.timing==state.timing
      &&restored.bodies[0].position==state.bodies[0].position
      &&restored.angularAcceleration==state.angularAcceleration,"Current state offsets");
  rejects([&]{cp::loadState(current,1,"new-physical-settings");},"Changed physics accepted");

  legacyFile(directory/"bad-size.bin",sizeof(FrozenV1Diagnostic)+8);
  rejects([&]{cp::loadState(directory/"bad-size.bin",1,"old-pure-gr-signature");},
      "Version 1 with incorrect diagnostics ABI accepted");
  legacyFile(directory/"bad-version.bin",sizeof(FrozenV1Diagnostic),99);
  rejects([&]{cp::loadState(directory/"bad-version.bin",1,"old-pure-gr-signature");},
      "Unsupported state version accepted");
}
void signatures() {
  c::Config cfg;c::Units units(cfg);
  const auto pure=cp::signature(cfg,units,2,1);
  require(pure.find("cmc_")==std::string::npos,"Inactive CMC changed pure signature");
  require(!cfg.current_adhesion_rolling&&pure.find("current_adhesion_rolling")==std::string::npos,
      "Default rolling mode changed the pure-Gr signature");
  cfg.pass_max=0;
  require(cp::signature(cfg,units,2,1)==pure,"pass_max policy became immutable physics");
  cfg.cmc_contact_version=1;cfg.cmc_contact_gap=4.e-9;cfg.cmc_cohesion_retention=.2;
  const auto coated=cp::signature(cfg,units,2,1);
  require(coated!=pure&&coated.find("cmc_contact_gap=")!=std::string::npos
      &&coated.find("cmc_cohesion_retention=")!=std::string::npos
      &&coated.find("cmc_contact_version=1\n")!=std::string::npos,
      "Active CMC physics missing from signature");
  cfg.cmc_cohesion_retention=0.;
  require(cp::signature(cfg,units,2,1)!=coated,"Changed screening did not invalidate checkpoint");
  cfg.cmc_contact_version=0;cfg.cmc_contact_gap=0.;cfg.cmc_cohesion_retention=1.;
  cfg.free_cmc_inner_repulsion_range=5.e-10;
  cfg.free_cmc_inner_repulsion_power=2.3;
  require(cp::signature(cfg,units,2,1)==pure,"Dormant inner repulsion changed pure checkpoint signature");
  cfg.free_cmc_inner_repulsion_work=.007;
  const auto inner=cp::signature(cfg,units,2,1);
  require(inner!=pure&&inner.find("free_cmc_inner_repulsion_version=1\n")!=std::string::npos
      &&inner.find("free_cmc_inner_repulsion_work=")!=std::string::npos,
      "Active inner repulsion missing from checkpoint signature");
  cfg.free_cmc_inner_repulsion_work=.008;
  require(cp::signature(cfg,units,2,1)!=inner,"Changed inner work accepted by checkpoint signature");
  cfg.free_cmc_inner_repulsion_work=.007;cfg.free_cmc_inner_repulsion_range=6.e-10;
  require(cp::signature(cfg,units,2,1)!=inner,"Changed inner range accepted by checkpoint signature");
  cfg.free_cmc_inner_repulsion_range=5.e-10;cfg.free_cmc_inner_repulsion_power=2.4;
  require(cp::signature(cfg,units,2,1)!=inner,"Changed inner power accepted by checkpoint signature");
  cfg.free_cmc_inner_repulsion_work=0.;
  cfg.cmc_net_contact_force=1.7e-10;cfg.cmc_net_barrier_force=2.3e-11;
  cfg.cmc_net_attraction_range=1.1e-9;cfg.cmc_net_repulsion_range=6.3e-9;
  cfg.cmc_net_reference_length=14.e-6;
  require(cp::signature(cfg,units,2,1)==pure,"Dormant net CMC changed the legacy checkpoint signature");
  cfg.cmc_net_blend=.63;
  const auto net=cp::signature(cfg,units,2,1);
  require(net!=pure&&net.find("cmc_net_potential_version=1\n")!=std::string::npos,
      "Active net CMC is missing from the checkpoint signature");
  for(auto field:{&c::Config::cmc_net_blend,&c::Config::cmc_net_contact_force,
      &c::Config::cmc_net_barrier_force,&c::Config::cmc_net_attraction_range,
      &c::Config::cmc_net_repulsion_range,&c::Config::cmc_net_reference_length}) {
    auto changed=cfg;changed.*field*=1.1;
    require(cp::signature(changed,units,2,1)!=net,"Changed net-CMC parameter accepted by checkpoint signature");
  }
  cfg.cmc_coordination_start=.7;cfg.cmc_coordination_end=2.3;cfg.cmc_coordination_floor=.2;
  require(cp::signature(cfg,units,2,1)==net,"Disabled coordination changed the previous net-CMC signature");
  cfg.cmc_coordination_enabled=true;
  const auto coordinated=cp::signature(cfg,units,2,1);
  require(coordinated!=net&&coordinated.find("cmc_coordination_version=1\n")!=std::string::npos
      &&coordinated.find("cmc_coordination_enabled=1\n")!=std::string::npos,
      "Active coordination is missing from the checkpoint signature");
  require(coordinated.find("current_adhesion_rolling")==std::string::npos,
      "Legacy implicit coordination rolling acquired a new fingerprint");
  for(auto field:{&c::Config::cmc_coordination_start,&c::Config::cmc_coordination_end,
      &c::Config::cmc_coordination_floor}) {
    auto changed=cfg;changed.*field*=1.1;
    require(cp::signature(changed,units,2,1)!=coordinated,
        "Changed coordination parameter accepted by checkpoint signature");
  }
  cfg.cmc_net_blend=0.;
  require(cp::signature(cfg,units,2,1)==pure,
      "Coordination with zero net blend changed the legacy checkpoint signature");
  cfg.current_adhesion_rolling=true;
  const auto rolling=cp::signature(cfg,units,2,1);
  require(rolling!=pure&&rolling.find("current_adhesion_rolling=1\n")!=std::string::npos
      &&rolling.find("current_adhesion_rolling_version=1\n")!=std::string::npos,
      "Independent current-adhesion rolling is missing from the signature");
  cfg.current_adhesion_rolling=false;
  require(cp::signature(cfg,units,2,1)==pure,"Disabled independent rolling changed the old signature");
  cfg.rough_contact_enabled=false;
  const auto contactDisabled=cp::signature(cfg,units,2,1);
  cfg.current_adhesion_rolling=true;
  require(cp::signature(cfg,units,2,1)==contactDisabled,
      "Independent rolling changed the signature with rough contact disabled");
}

void netConfiguration(const cp::fs::path& directory) {
  const std::string base="surface_adhesion=1\ncmc_net_blend=0.63\n"
      "cmc_net_contact_force=1.7e-10\ncmc_net_barrier_force=2.3e-11\n"
      "cmc_net_attraction_range=1.1e-9\ncmc_net_repulsion_range=6.3e-9\n"
      "cmc_net_reference_length=1.4e-5\n";
  const auto parse=[&](const std::string& text) {
    const auto path=directory/"net.cfg";{std::ofstream output(path);output<<text;}
    std::string program="config-test",option="--config",name=path.string();
    char* args[]={program.data(),option.data(),name.data()};
    return c::parseConfig(3,args);
  };
  const auto parsed=parse(base);
  require(parsed.cmc_net_blend==.63&&parsed.cmc_net_contact_force==1.7e-10
      &&parsed.cmc_net_barrier_force==2.3e-11&&parsed.cmc_net_attraction_range==1.1e-9
      &&parsed.cmc_net_repulsion_range==6.3e-9&&parsed.cmc_net_reference_length==14.e-6,
      "Net-CMC configuration fields were not parsed");
  for(const std::string invalid:{"cmc_net_blend=1.1\n","cmc_net_contact_force=0\n",
      "cmc_net_barrier_force=nan\n","cmc_net_attraction_range=0\n",
      "cmc_net_repulsion_range=4e-7\n","cmc_net_reference_length=0\n",
      "rough_contact_enabled=0\n","surface_adhesion=0\n","local_gap_fraction=0.1\n",
      "cmc_contact_gap=3e-9\n","cmc_cohesion_retention=0.5\n"})
    rejects([&]{parse(base+invalid);},"Invalid net-CMC configuration accepted");
  const auto inactive=parse("cmc_net_blend=0\n");
  require(!inactive.surface_adhesion,"Inactive net CMC unexpectedly changed the legacy model");
  require(!parsed.cmc_coordination_enabled&&parsed.cmc_coordination_start==1.
      &&parsed.cmc_coordination_end==2.&&parsed.cmc_coordination_floor==.1,
      "Coordination defaults changed the legacy configuration");
  const auto coordinated=parse(base+"cmc_coordination_enabled=1\ncmc_coordination_start=0.7\n"
      "cmc_coordination_end=2.3\ncmc_coordination_floor=0.2\n");
  require(coordinated.cmc_coordination_enabled&&coordinated.cmc_coordination_start==.7
      &&coordinated.cmc_coordination_end==2.3&&coordinated.cmc_coordination_floor==.2,
      "Coordination configuration fields were not parsed");
  for(const std::string invalid:{"cmc_coordination_enabled=2\n","cmc_coordination_enabled=true\n",
      "cmc_coordination_start=-0.1\n","cmc_coordination_start=2\n",
      "cmc_coordination_end=1\n","cmc_coordination_start=nan\n",
      "cmc_coordination_end=inf\n","cmc_coordination_floor=nan\n",
      "cmc_coordination_floor=-0.1\n","cmc_coordination_floor=1.1\n",
      "cmc_coordination_start=2.2250738585072014e-308\ncmc_coordination_end=2.225073858507202e-308\n"})
    rejects([&]{parse(base+invalid);},"Invalid coordination configuration accepted");
  require(parse(base+"cmc_coordination_floor=0\n").cmc_coordination_floor==0.
      &&parse(base+"cmc_coordination_floor=1\n").cmc_coordination_floor==1.,
      "Coordination floor endpoints must be valid");
  const auto dormant=parse("cmc_net_blend=0\ncmc_coordination_enabled=1\n");
  require(dormant.cmc_coordination_enabled&&!dormant.surface_adhesion,
      "Zero-blend coordination must remain dormant");
  require(!parsed.current_adhesion_rolling,"Independent rolling is enabled by default");
  require(parse(base+"current_adhesion_rolling=1\n").current_adhesion_rolling
      &&!parse(base+"current_adhesion_rolling=0\n").current_adhesion_rolling,
      "Independent current-adhesion rolling was not parsed");
  for(const std::string invalid:{"2","-1","true","false","1.0","nan",""})
    rejects([&]{parse(base+"current_adhesion_rolling="+invalid+"\n");},
        "Invalid independent rolling boolean accepted");
  const auto dormantRolling=parse("rough_contact_enabled=0\ncurrent_adhesion_rolling=1\n");
  require(dormantRolling.current_adhesion_rolling&&!dormantRolling.rough_contact_enabled,
      "Independent rolling must be allowed to remain dormant without rough contact");
}

void coordinationRestart(const cp::fs::path& directory) {
  c::Config cfg;cfg.surface_adhesion=true;cfg.cmc_net_blend=1.;
  c::Units units(cfg);
  const auto previous=cp::signature(cfg,units,1,1);
  cp::State state;state.bodies.resize(1);state.angularAcceleration.resize(1);state.step=17;
  const auto path=directory/"net-before-coordination.bin";
  cp::saveState(path,state,previous);
  cfg.cmc_coordination_enabled=true;
  const auto updated=cp::signature(cfg,units,1,1);
  rejects([&]{cp::loadState(path,1,updated);},
      "Previous net-potential state was accepted by the coordination model");
  const auto activePath=directory/"coordination.bin";
  cp::saveState(activePath,state,updated);
  require(cp::loadState(activePath,1,updated).step==17,"Coordination checkpoint did not roundtrip");
  rejects([&]{cp::loadState(activePath,1,previous);},
      "Coordination state was accepted with coordination disabled");
  cfg.cmc_coordination_enabled=false;cfg.cmc_coordination_floor=.8;
  require(cp::loadState(path,1,cp::signature(cfg,units,1,1)).step==17,
      "Disabled coordination rejected an unchanged net-potential checkpoint");
}

void currentRollingRestart(const cp::fs::path& directory) {
  c::Config cfg;cfg.surface_adhesion=true;cfg.cmc_net_blend=1.;
  c::Units units(cfg);
  const auto previous=cp::signature(cfg,units,1,1);
  cp::State state;state.bodies.resize(1);state.angularAcceleration.resize(1);state.step=23;
  const auto oldPath=directory/"net-before-current-rolling.bin";
  cp::saveState(oldPath,state,previous);
  cfg.current_adhesion_rolling=true;
  const auto updated=cp::signature(cfg,units,1,1);
  rejects([&]{cp::loadState(oldPath,1,updated);},
      "Birth-adhesion rolling checkpoint accepted by current-adhesion rolling mode");
  const auto activePath=directory/"net-current-rolling.bin";
  cp::saveState(activePath,state,updated);
  require(cp::loadState(activePath,1,updated).step==23,"Current-adhesion rolling checkpoint failed roundtrip");
  rejects([&]{cp::loadState(activePath,1,previous);},
      "Current-adhesion rolling checkpoint accepted by birth-adhesion rolling mode");
  cfg.current_adhesion_rolling=false;
  require(cp::loadState(oldPath,1,cp::signature(cfg,units,1,1)).step==23,
      "Disabled independent rolling rejected the previous net-potential checkpoint");
}
}
int main() {
  try {Scratch scratch;migrateLegacyAndRoundtrip(scratch.path);signatures();netConfiguration(scratch.path);
    coordinationRestart(scratch.path);currentRollingRestart(scratch.path);
    std::cout<<"PASS: checkpoint v1 migration, v2 counters, ABI/version guards, coordination/rolling config and physics signatures\n";
    return 0;
  }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
