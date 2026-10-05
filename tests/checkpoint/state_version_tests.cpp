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
}
}
int main() {
  try {Scratch scratch;migrateLegacyAndRoundtrip(scratch.path);signatures();
    std::cout<<"PASS: checkpoint v1 migration, v2 counters, ABI/version guards and physics signatures\n";
    return 0;
  }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
