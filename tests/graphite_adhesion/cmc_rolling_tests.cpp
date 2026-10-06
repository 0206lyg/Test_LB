/* c++ -std=c++17 -O2 -Wall -Wextra -pedantic -I olb-1.9r0/src/slurry/gr_re2 \
     tests/graphite_adhesion/cmc_rolling_tests.cpp -o /tmp/cmc_rolling_tests */
#include "roughContact.h"
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace g=slurry::gr_re2::graphite;
namespace {
constexpr double adhesion=2.e-10;
const g::Vec3 normal{0.,0.,1.};
void require(bool value,const std::string& message){if(!value)throw std::runtime_error(message);}
void near(double actual,double expected,const std::string& message,double relative=2.e-12){
  if(!std::isfinite(actual)||std::abs(actual-expected)>1.e-33+relative*std::abs(expected)){
    std::ostringstream s;s.precision(17);s<<message<<": actual="<<actual<<", expected="<<expected;
    throw std::runtime_error(s.str());
  }
}
g::RoughContactSettings settings(bool current=true){
  g::RoughContactSettings s;s.enabled=true;s.currentAdhesionRolling=current;
  s.friction=0.;s.tangentialStiffness=0.;return s;
}
g::RoughContactState loaded(const g::Vec3& angle){
  const auto s=settings();g::RoughContactState old;old.active=true;old.normal=normal;
  old.rollingCap=s.rollingLength*adhesion;old.rollingStiffness=old.rollingCap/s.rollingYieldAngle;
  old.elasticRoll=angle;old.elasticEnergy=.5*old.rollingStiffness*g::dot(angle,angle);
  old.rollingTorque=g::scale(angle,-old.rollingStiffness);return old;
}
g::RoughContactResult evaluate(const g::RoughContactState& old,double force,
                              const g::RoughContactSettings& s=settings(),
                              const g::Vec3& increment={},bool engaged=true){
  g::Body a,b;b.position={0.,0.,.4e-6+s.gap};
  // Equal/opposite angular velocities give zero common frame rotation, so
  // increment is exactly the physical relative rolling angle in this test.
  a.omega=g::scale(increment,.5);b.omega=g::scale(increment,-.5);
  return g::roughContact(a,b,normal,{0.,0.,.2e-6},{0.,0.,-.2e-6},0.,force,1.,old,s,engaged);
}
bool same(const g::RoughContactResult& a,const g::RoughContactResult& b){
  const auto& x=a.candidateState;const auto& y=b.candidateState;
  return a.forceI==b.forceI&&a.torqueI==b.torqueI&&a.torqueJ==b.torqueJ
      &&a.rollingTorque==b.rollingTorque&&a.trialRoll==b.trialRoll
      &&a.plasticRollWork==b.plasticRollWork&&a.releasedEnergy==b.releasedEnergy
      &&a.elasticEnergy==b.elasticEnergy&&a.dTorqueDRollVelocity==b.dTorqueDRollVelocity
      &&a.rolling==b.rolling&&x.elasticRoll==y.elasticRoll
      &&x.rollingCap==y.rollingCap&&x.rollingStiffness==y.rollingStiffness;
}
void frozenCompatibility(){
  const auto old=loaded({.003,.004,0.});const auto s=settings(false);
  const auto reference=evaluate(old,adhesion,s,{.001,-.002,0.});
  for(double force:{0.,adhesion*.01,adhesion*1.e4,-adhesion})
    require(same(reference,evaluate(old,force,s,{.001,-.002,0.})),
            "disabled current law must ignore the current adhesion on an existing contact bitwise");
  const auto unchanged=evaluate(old,adhesion,settings(),{.001,-.002,0.});
  require(same(reference,unchanged),"unchanged current stiffness must use legacy arithmetic bitwise");
  const g::RoughContactState fresh;
  require(same(evaluate(fresh,adhesion,s,{.003,.004,0.}),
               evaluate(fresh,adhesion,settings(),{.003,.004,0.})),
          "current and birth-frozen laws must agree for a fresh contact");
}
void energyNeutralCycles(){
  auto old=loaded({.0006,.0008,0.});const double energy=old.elasticEnergy;
  for(double factor:{4.,.25,.5,1.,8.,.25,1.}){
    const auto r=evaluate(old,adhesion*factor);
    near(r.elasticEnergy,energy,"strength changes below capacity must preserve stored energy");
    near(r.releasedEnergy,0.,"no release on energy-admissible zero-motion cycles");
    near(r.plasticRollWork,0.,"no plastic work without rolling motion");
    near(r.candidateState.rollingCap,settings().rollingLength*adhesion*factor,"current cap");
    require(g::norm(r.rollingTorque)<=r.candidateState.rollingCap*(1.+1.e-14),"torque within current cap");
    old=r.candidateState;
  }
  const auto empty=loaded({});
  for(double factor:{0.,1.,20.,.01}){
    const auto r=evaluate(empty,adhesion*factor);
    require(r.elasticEnergy==0.&&r.releasedEnergy==0.&&r.plasticRollWork==0.
            &&r.rollingTorque==g::Vec3{},"zero-history changes must not create energy or dissipation");
  }
}
void weakeningAndRecovery(){
  const auto old=loaded({.008,0.,0.});const auto weak=evaluate(old,.25*adhesion);
  const double capacity=.5*weak.candidateState.rollingCap*settings().rollingYieldAngle;
  near(weak.elasticEnergy,capacity,"weakening respects new storage capacity");
  near(weak.releasedEnergy,old.elasticEnergy-capacity,"lost rolling energy is released exactly");
  near(weak.plasticRollWork,0.,"material weakening is not rolling plastic work");
  near(weak.candidateState.releasedEnergy,weak.releasedEnergy,"release committed in candidate state");
  near(g::norm(weak.rollingTorque),weak.candidateState.rollingCap,"weakened torque cap");
  const auto repeat=evaluate(weak.candidateState,.25*adhesion);
  near(repeat.releasedEnergy,0.,"unchanged zero-motion step does not release energy twice");
  near(repeat.plasticRollWork,0.,"unchanged zero-motion step does not dissipate");
  const auto strong=evaluate(weak.candidateState,4.*adhesion);
  near(strong.elasticEnergy,weak.elasticEnergy,"strength recovery cannot recreate released energy");
  near(strong.releasedEnergy,0.,"strength recovery is energy neutral");
  const auto zero=evaluate(strong.candidateState,0.);
  near(zero.releasedEnergy,strong.elasticEnergy,"zero cap releases the remaining energy");
  require(zero.candidateState.elasticRoll==g::Vec3{}&&zero.rollingTorque==g::Vec3{}
          &&zero.elasticEnergy==0.,"zero cap clears rolling history");
  const auto recovered=evaluate(zero.candidateState,adhesion);
  require(recovered.candidateState.elasticRoll==g::Vec3{}&&recovered.rollingTorque==g::Vec3{}
          &&recovered.elasticEnergy==0.&&recovered.releasedEnergy==0.,"recovery cannot resurrect zeroed history");
  const auto open=evaluate(weak.candidateState,adhesion,settings(),{},false);
  near(open.releasedEnergy,weak.elasticEnergy,"opening releases only energy still present");
}
void materialAndMechanicalEnergy(){
  for(double force:{.25*adhesion,4.*adhesion}){
    const auto old=loaded({.008,0.,0.});
    const auto material=evaluate(old,force);
    for(const g::Vec3 increment:{g::Vec3{.004,0.,0.},g::Vec3{-.004,.003,0.}}){
      const auto r=evaluate(old,force,settings(),increment);
      const double work=-g::dot(r.rollingTorque,increment);
      const auto elasticChange=g::sub(r.candidateState.elasticRoll,material.candidateState.elasticRoll);
      const double implicitLoss=.5*r.candidateState.rollingStiffness*g::dot(elasticChange,elasticChange);
      // Exact backward-Euler spring/return-map identity, with the material
      // release accounted separately from physical rolling plastic work.
      near(work+old.elasticEnergy,r.elasticEnergy+r.plasticRollWork+r.releasedEnergy+implicitLoss,
           "mechanical work plus old energy equals new energy and all losses");
      near(r.releasedEnergy,material.releasedEnergy,"material release precedes physical rolling");
      require(r.plasticRollWork>=0.&&r.releasedEnergy>=0.,"all losses nonnegative");
      require(g::norm(r.rollingTorque)<=r.candidateState.rollingCap*(1.+1.e-14),"final return map respects cap");
    }
  }
}
void transportedEnergy(){
  const auto old=loaded({.003,.004,0.});const auto s=settings();
  const g::Vec3 omega{.2,-.1,.3};const double dt=.5;
  const auto rotation=g::rotationIncrement(g::scale(omega,dt));
  g::Body a,b;b.position=g::mul(rotation,{0.,0.,.4e-6+s.gap});
  a.omega=b.omega=omega;a.velocity=g::cross(omega,a.position);b.velocity=g::cross(omega,b.position);
  const auto r=g::roughContact(a,b,g::mul(rotation,normal),g::mul(rotation,{0.,0.,.2e-6}),
      g::mul(rotation,{0.,0.,-.2e-6}),0.,4.*adhesion,dt,old,s);
  const auto expected=g::scale(g::mul(rotation,old.elasticRoll),.5);
  for(int k=0;k<3;++k)near(r.candidateState.elasticRoll[k],expected[k],"energy coordinate rotates objectively");
  near(r.elasticEnergy,old.elasticEnergy,"rigid transport and strengthening preserve energy");
  near(r.releasedEnergy,0.,"rigid transport has no material release");
  near(r.plasticRollWork,0.,"rigid transport has no plastic rolling");
}
}
int main(){
  const std::vector<std::pair<std::string,std::function<void()>>> tests={
    {"frozen compatibility",frozenCompatibility},{"energy neutral cycles",energyNeutralCycles},
    {"weakening and recovery",weakeningAndRecovery},{"material and mechanical energy",materialAndMechanicalEnergy},
    {"transported energy",transportedEnergy}
  };
  try{for(const auto& test:tests){test.second();std::cout<<"PASS "<<test.first<<'\n';}}
  catch(const std::exception& error){std::cerr<<"FAIL "<<error.what()<<'\n';return 1;}
}
