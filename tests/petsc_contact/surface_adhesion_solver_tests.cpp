// Production PETSc regressions for the matched surface-adhesion potential.
// These isolate the particle/contact solve from the outer LB integration.
#include "particleSubsteps.h"
#include <petscsys.h>
#include <iostream>
#include <string>

namespace g=slurry::gr_re2::graphite;
namespace d=g::particle_detail;

namespace {
void require(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}
void near(double actual,double expected,double absolute,double relative,const std::string& message) {
  if(!std::isfinite(actual)||std::abs(actual-expected)>absolute+relative*std::abs(expected)) {
    std::ostringstream out;out<<std::setprecision(17)<<message<<": "<<actual<<" != "<<expected;
    throw std::runtime_error(out.str());
  }
}
g::ParticleStepSettings settings() {
  g::ParticleStepSettings s;
  s.solverBackend="petsc";s.box={100.e-6,100.e-6,100.e-6};s.shearRate=0.;
  s.pair.surfaceAdhesion=true;s.pair.localGapFraction=0.;
  s.pair.sigma=.4197e-9;s.pair.roughnessGap=s.rough.gap=2.e-9;
  s.pair.adhesionWork=.0219;s.pair.adhesionRange=.67e-9;
  s.pair.curvatureSwitchGap=5.e-9;s.pair.curvatureCutoffGap=20.e-9;
  s.nearField.enabled=false;s.rough.enabled=true;s.rough.friction=1.;
  s.rough.tangentialStiffness=80.;
  s.maxSubsteps=1;s.maxNewtonIterations=80;s.maxKrylovIterations=200;s.maxLineSearch=24;
  s.relativeTolerance=1.e-7;s.forceAbsoluteTolerance=1.e-14;
  s.torqueAbsoluteTolerance=1.e-21;s.contactGapTolerance=1.e-13;
  return s;
}
void setMass(g::Body& body) {
  body.mass=2.e-11;
  for(int k=0;k<3;++k)
    body.inertiaBody[k]=body.mass*(body.axes[(k+1)%3]*body.axes[(k+1)%3]
        +body.axes[(k+2)%3]*body.axes[(k+2)%3])/5.;
}
std::vector<g::Body> facePair(double gap) {
  std::vector<g::Body> bodies(2);
  bodies[0].position={20.e-6,20.e-6,20.e-6};
  bodies[1].position=g::add(bodies[0].position,{0.,0.,bodies[0].axes[2]+bodies[1].axes[2]+gap});
  for(auto& body:bodies)setMass(body);
  return bodies;
}
void accepted(const g::ParticleStepDiagnostics& info,const g::ParticleStepSettings& s) {
  require(info.substeps>=1&&info.substeps<=s.maxSubsteps,"Invalid accepted subdivision count");
  require(info.maxForceResidualRatio<=1.&&info.maxTorqueResidualRatio<=1.,
          "Accepted state violates physical force/torque tolerances");
  require(info.contactGapViolation<=s.contactGapTolerance,"Accepted state penetrates the contact constraint");
}
bool sameState(const g::RoughContactState& a,const g::RoughContactState& b) {
  return a.active==b.active&&a.normal==b.normal&&a.normalLoad==b.normalLoad
      &&a.elasticSlip==b.elasticSlip&&a.elasticRoll==b.elasticRoll
      &&a.rollingCap==b.rollingCap&&a.rollingStiffness==b.rollingStiffness
      &&a.elasticEnergy==b.elasticEnergy&&a.plasticSlipWork==b.plasticSlipWork
      &&a.plasticRollWork==b.plasticRollWork&&a.releasedEnergy==b.releasedEnergy;
}

void closedFaceBalanceAndTrialHistory() {
  auto s=settings();s.maxSubsteps=256;auto bodies=facePair(s.rough.gap);
  const std::vector<g::Vec3> zero(2);g::PersistentContactState history;
  std::vector<g::GapCache> cache;
  constexpr double dt=1.e-7;int maxSubdivision=0;
  for(int k=0;k<3;++k) {
    const auto info=g::advanceParticles(bodies,zero,zero,dt,k*dt,s,&cache,&history);
    accepted(info,s);maxSubdivision=std::max(maxSubdivision,info.substeps);
    const auto pair=g::evaluatePair(bodies[0],bodies[1],s.pair);
    require(history[0].active,"Adhesive FF contact was not activated");
    const double adhesion=g::dot(pair.forceI,pair.normal);
    require(adhesion>9.e-7&&adhesion<1.e-6,"Fixture no longer exercises the physical FF adhesive force");
    near(history[0].normalLoad,adhesion,1.e-12,1.e-3,
         "Constraint reaction must balance adhesion once");
    near(pair.gap,s.rough.gap,s.contactGapTolerance,0.,"Closed FF gap");
    near(g::norm(bodies[0].velocity),0.,s.contactGapTolerance*info.substeps/dt,0.,
         "Static contact drift exceeds the accepted gap-resolution speed");
  }

  // Probe a saturated Coulomb trial at a feasible, laterally moved FF contact.
  // An unaccepted trial must not commit its slip, plastic work, or birth data.
  const auto committed=history[0];const double length=bodies[0].axes[0];
  std::vector<int> slots{0};
  d::Residual residual{bodies,zero,zero,s,cache,history,slots,dt,0.,length,1.,1};
  d::Vector q(13,0.);q[6]=50.e-9/length;
  g::Body moved=bodies[1];moved.position[0]+=50.e-9;
  for(int iter=0;iter<4;++iter) {
    const auto gap=g::closestEllipsoidGap(bodies[0],moved);
    moved.position=g::add(moved.position,g::scale(gap.normal,s.rough.gap-gap.gap));
  }
  for(int k=0;k<3;++k)q[6+k]=(moved.position[k]-bodies[1].position[k])/length;
  q[12]=committed.normalLoad;
  d::Evaluation evaluation;std::string error;
  require(residual(q,evaluation,error),"Feasible sliding trial failed: "+error);
  require(evaluation.contacts[0].sliding,"Trial must exercise the Coulomb yield cap");
  near(g::norm(evaluation.contacts[0].tangentForce),s.rough.friction*q[12],1.e-18,1.e-12,
       "Coulomb cap must use mu*N, without adding adhesion twice");
  require(evaluation.contacts[0].plasticSlipWork>0.,"Saturated trial must produce candidate plastic slip work");
  require(sameState(history[0],committed),"Residual committed irreversible contact history");

  // The old shifted-gap model was singular below 1.7 nm. The new polynomial
  // continuation permits Newton trials there while retaining h>=h0 at acceptance.
  require(d::admissibleTrialGap(1.65e-9,s),"Obsolete shifted-gap domain remains active");
  q.assign(13,0.);q[8]=(1.65e-9-s.rough.gap)/length;q[12]=committed.normalLoad;
  d::PetscParticleContext context(residual);
  require(context.evaluate(q,evaluation),"PETSc rejected a nonsingular surface-adhesion trial");
  require(context.domainErrors==0,"Obsolete local-gap singularity registered a domain error");
  require(evaluation.diagnostic.contactGapViolation>0.,"Penetrating trial must still report constraint violation");
  require(sameState(history[0],committed),"Penetrating trial mutated committed contact history");
  std::cout<<"PASS: closed FF balance, single-count Coulomb load, transactional trial history/domain; max subdivisions="
           <<maxSubdivision<<'\n';
}

void contactCreationAndRelease() {
  auto s=settings();s.maxSubsteps=256;auto bodies=facePair(s.rough.gap+.15e-9);
  const std::vector<g::Vec3> zero(2);g::PersistentContactState history;
  std::vector<g::GapCache> cache;constexpr double dt=1.e-7;
  accepted(g::advanceParticles(bodies,zero,zero,dt,0.,s,&cache,&history),s);
  require(history[0].active,"Cohesive approach did not create a contact");
  near(g::closestEllipsoidGap(bodies[0],bodies[1]).gap,s.rough.gap,
       s.contactGapTolerance,0.,"Created contact gap");

  // A valid preloaded elastic state lets the release test check one-time energy
  // accounting as well as the normal complementarity transition.
  auto& old=history[0];old.elasticSlip={1.e-11,0.,0.};old.elasticRoll={0.,2.e-4,0.};
  old.tangentForce=g::scale(old.elasticSlip,-s.rough.tangentialStiffness);
  old.rollingTorque=g::scale(old.elasticRoll,-old.rollingStiffness);
  old.elasticEnergy=.5*s.rough.tangentialStiffness*g::dot(old.elasticSlip,old.elasticSlip)
      +.5*old.rollingStiffness*g::dot(old.elasticRoll,old.elasticRoll);
  const double stored=old.elasticEnergy;
  const auto pair=g::evaluatePair(bodies[0],bodies[1],s.pair);
  const double pull=2.*g::dot(pair.forceI,pair.normal);
  const std::vector<g::Vec3> force{g::scale(pair.normal,-pull),g::scale(pair.normal,pull)};
  const auto released=g::advanceParticles(bodies,force,zero,dt,dt,s,&cache,&history);
  accepted(released,s);
  require(!history[0].active,"Tensile separation did not release the contact");
  near(history[0].normalLoad,0.,1.e-20,0.,"Released normal reaction");
  near(g::norm(history[0].elasticSlip),0.,1.e-24,0.,"Release clears slip history");
  near(g::norm(history[0].elasticRoll),0.,1.e-20,0.,"Release clears roll history");
  near(released.contactDissipation*dt,stored,1.e-28,1.e-10,
       "Release accounts stored contact energy once across all accepted substeps");
  require(g::closestEllipsoidGap(bodies[0],bodies[1]).gap>s.rough.gap+s.pair.adhesionRange,
          "Released fixture must pass the cohesive cutoff");
  const auto open=g::advanceParticles(bodies,force,zero,dt,2.*dt,s,&cache,&history);
  accepted(open,s);
  near(open.contactDissipation*dt,0.,1.e-28,0.,"Stored contact energy was released twice");
  std::cout<<"PASS: cohesive contact creation and release across cutoff, one-time release work\n";
}

g::Vec3 support(const g::Body& body,const g::Vec3& n) {
  g::Vec3 a2{};for(int k=0;k<3;++k)a2[k]=body.axes[k]*body.axes[k];
  const auto q=g::rotatedDiagonal(body.rotation,a2);
  return g::scale(g::mul(q,n),1./std::sqrt(g::dot(n,g::mul(q,n))));
}
void offsetTiltedCutoffCrossings() {
  const auto s=settings();auto bodies=facePair(s.rough.gap);
  bodies[0].rotation=g::rotationIncrement({.005,.05,.03});
  bodies[1].rotation=g::rotationIncrement({.015,-.06,.11});
  const auto n=g::normalized(g::Vec3{.07,.03,1.});
  const double gap=s.rough.gap+.9*s.pair.adhesionRange;
  bodies[1].position=g::add(bodies[0].position,
      g::add(g::add(support(bodies[0],n),support(bodies[1],n)),g::scale(n,gap)));
  require(std::abs(bodies[1].position[0]-bodies[0].position[0])>3.e-7,
          "Tilted fixture must also have a substantial lateral offset");
  const std::vector<g::Vec3> zero(2);g::PersistentContactState history;
  std::vector<g::GapCache> cache;constexpr double dt=1.e-8;
  bodies[0].velocity=g::scale(n,-.005);bodies[1].velocity=g::scale(n,.005);
  accepted(g::advanceParticles(bodies,zero,zero,dt,0.,s,&cache,&history),s);
  auto pair=g::evaluatePair(bodies[0],bodies[1],s.pair);
  require(pair.gap>s.rough.gap+s.pair.adhesionRange,"Opening tilted pair did not cross cohesive cutoff");
  require(!history[0].active,"Open cohesive zone incorrectly generated contact history");
  auto background=s.pair;background.adhesionWork=g::surfaceBackgroundWork(background);
  const auto outside=g::evaluatePair(bodies[0],bodies[1],background);
  near(g::norm(g::sub(pair.forceI,outside.forceI)),0.,1.e-20,0.,"Cohesive force persists beyond range");

  bodies[0].velocity=g::scale(n,.005);bodies[1].velocity=g::scale(n,-.005);
  accepted(g::advanceParticles(bodies,zero,zero,dt,dt,s,&cache,&history),s);
  pair=g::evaluatePair(bodies[0],bodies[1],s.pair);
  require(pair.gap<s.rough.gap+s.pair.adhesionRange&&pair.gap>s.rough.gap,
          "Closing tilted pair did not return inside the open cohesive zone");
  require(!history[0].active,"Noncontact cohesive attraction activated contact history");
  const auto inside=g::evaluatePair(bodies[0],bodies[1],background);
  require(g::norm(g::sub(pair.forceI,inside.forceI))>1.e-8,
          "Closing pair lacks a resolved cohesive force");
  std::cout<<"PASS: offset/tilted moving pair crosses cohesive cutoff in both directions\n";
}
}

int main(int argc,char** argv) {
  if(PetscInitialize(&argc,&argv,nullptr,nullptr))return 1;
  int result=0;
  try {
    closedFaceBalanceAndTrialHistory();contactCreationAndRelease();offsetTiltedCutoffCrossings();
    std::cout<<"Surface adhesion PETSc solver tests passed\n";
  } catch(const std::exception& error) {
    std::cerr<<"Surface adhesion PETSc solver test failed: "<<error.what()<<'\n';result=1;
  }
  return PetscFinalize()?1:result;
}
