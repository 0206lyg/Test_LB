// Standalone checks of the conservative, finite-support effective CMC model.
/* c++ -std=c++17 -O2 -Wall -Wextra -pedantic \
     -I olb-1.9r0/src/slurry/gr_re2 \
     tests/graphite_adhesion/cmc_net_potential_tests.cpp -o /tmp/cmc_net_potential_tests
   /tmp/cmc_net_potential_tests */
#include "reSquaredPotential.h"

#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace g=slurry::gr_re2::graphite;
namespace {
constexpr double pi=3.1415926535897932384626433832795;
using Bodies=std::pair<g::Body,g::Body>;
void require(bool condition,const std::string& message){
  if(!condition)throw std::runtime_error(message);
}
void near(double actual,double expected,double absolute,double relative,const std::string& message){
  if(!std::isfinite(actual)||!std::isfinite(expected)
     ||std::abs(actual-expected)>absolute+relative*std::abs(expected)){
    std::ostringstream s;s<<std::setprecision(17)<<message<<": actual="<<actual<<", expected="<<expected;
    throw std::runtime_error(s.str());
  }
}
void nearVector(const g::Vec3& a,const g::Vec3& b,double absolute,double relative,const std::string& message){
  for(int k=0;k<3;++k)near(a[k],b[k],absolute,relative,message+" component "+std::to_string(k));
}
void rejects(const std::function<void()>& f,const std::string& message){
  try{f();}catch(const std::domain_error&){return;}
  throw std::runtime_error(message);
}
g::PairParameters parameters(double blend=1.){
  g::PairParameters p;p.sigma=.4197e-9;p.surfaceAdhesion=true;
  p.adhesionWork=.007666863221834307;
  // Nonzero legacy terms ensure full replacement really removes their tails.
  p.freeCmcRepulsionPressure=52143.97983709085;
  p.freeCmcInnerRepulsionWork=.0071170133084142865;
  p.cmcNetBlend=blend;return p;
}
g::Vec3 support(const g::Body& b,const g::Vec3& n){
  g::Vec3 a2{};for(int k=0;k<3;++k)a2[k]=b.axes[k]*b.axes[k];
  const auto qn=g::mul(g::rotatedDiagonal(b.rotation,a2),n);
  return g::scale(qn,1./std::sqrt(g::dot(n,qn)));
}
Bodies atGap(double h,int geometry=0){
  Bodies b;g::Vec3 n{0.,0.,1.};
  if(geometry==1)b.second.rotation=g::rotationIncrement({0.,pi/2.,0.});
  if(geometry==2)n={1.,0.,0.};
  if(geometry==3){
    b.first.axes={1.65e-6,1.33e-6,.20e-6};b.second.axes={1.45e-6,1.75e-6,.23e-6};
    b.first.rotation=g::rotationIncrement({.19,-.28,.08});
    b.second.rotation=g::rotationIncrement({-.11,.21,-.15});
    n=g::normalized({.18,-.12,1.});
  }
  if(geometry==4)n=g::normalized({.025,.013,1.});
  b.first.position={2.e-6,-3.e-6,1.e-6};
  b.second.position=g::add(b.first.position,
      g::add(g::add(support(b.first,n),support(b.second,n)),g::scale(n,h)));
  return b;
}
g::PairResult evaluate(const Bodies& b,const g::PairParameters& p){return g::evaluatePair(b.first,b.second,p);}
double outward(const g::PairResult& r){return -g::dot(r.forceI,r.normal);}
Bodies perturb(Bodies b,int k,double delta){
  if(k<3)b.second.position[k]+=delta;
  else{g::Vec3 angle{};angle[(k-3)%3]=delta;auto& body=k<6?b.first:b.second;
    body.rotation=g::rotateLaboratory(body.rotation,angle);}
  return b;
}
double derivative(const g::PairResult& r,int k){
  return k<3?r.forceI[k]:-(k<6?r.torqueI[k-3]:r.torqueJ[k-6]);
}
bool same(const g::PairResult& a,const g::PairResult& b){
  return a.ua==b.ua&&a.ur==b.ur&&a.energy==b.energy&&a.forceI==b.forceI
      &&a.torqueI==b.torqueI&&a.torqueJ==b.torqueJ
      &&a.forceAttractiveI==b.forceAttractiveI&&a.forceRepulsiveI==b.forceRepulsiveI
      &&a.torqueAttractiveI==b.torqueAttractiveI&&a.torqueAttractiveJ==b.torqueAttractiveJ
      &&a.torqueRepulsiveI==b.torqueRepulsiveI&&a.torqueRepulsiveJ==b.torqueRepulsiveJ;
}

void alignedPhysicalTargets(){
  const auto p=parameters();const double a=1.65e-6,c=.20e-6;
  const double lambda[]={a*a/c,2./std::sqrt((c/(a*a)+a/(c*c))*(c/(a*a)+1./a)),c};
  for(int geometry=0;geometry<3;++geometry){
    const double ratio=lambda[geometry]/p.cmcNetReferenceLength;
    const auto contact=evaluate(atGap(2.e-9,geometry),p);
    const auto saddle=evaluate(atGap(3.e-9,geometry),p);
    const auto peak=evaluate(atGap(6.e-9,geometry),p);
    near(outward(contact),-150.e-12*ratio,1.e-23,2.e-8,"specified contact pull-off force");
    near(contact.energy,14.e-21*ratio,1.e-31,2.e-8,"metastable contact energy relative to separated pair");
    near(saddle.energy,64.e-21*ratio,1.e-31,2.e-8,"specified approach work");
    near(saddle.energy-contact.energy,50.e-21*ratio,1.e-31,2.e-8,"specified contact escape work");
    near(outward(peak),20.e-12*ratio,1.e-23,2.e-8,"specified approach-force peak");
    near(peak.energy,32.e-21*ratio,1.e-31,2.e-8,"symmetric repulsive lobe half-work");
    nearVector(contact.torqueI,{},1.e-25,0.,"aligned first torque");
    nearVector(contact.torqueJ,{},1.e-25,0.,"aligned second torque");
    for(int k=0;k<100;++k)
      require(outward(evaluate(atGap((2.+.0099*k)*1.e-9,geometry),p))<0.,
              "actual contact is reached without an extra inner barrier");
    for(int k=1;k<100;++k)
      require(outward(evaluate(atGap((3.+.06*k)*1.e-9,geometry),p))>0.,
              "outer lobe has no secondary attractive basin");
  }
}

void supportAndContinuation(){
  const auto p=parameters();
  for(int geometry:{0,1,2,3,4})for(double h:{9.001e-9,20.e-9,450.e-9,501.e-9}){
    const auto r=evaluate(atGap(h,geometry),p);
    require(r.energy==0.&&r.ua==0.&&r.ur==0.&&r.forceI==g::Vec3{}
            &&r.torqueI==g::Vec3{}&&r.torqueJ==g::Vec3{},
            "full replacement has exactly zero interaction beyond its own cutoff");
  }
  const auto trial=evaluate(atGap(1.5e-9),p);
  near(outward(trial),-337.5e-12,1.e-23,2.e-8,"unclamped positive-gap Newton continuation");
  near(trial.energy,-104.75e-21,1.e-31,2.e-8,"continuation remains an energy gradient");
  require(g::effectiveContactGap(p)==2.e-9&&g::minimumPairGap(p)==0.,"contact plane and positive-gap domain unchanged");
  rejects([&]{evaluate(atGap(-.1e-9),p);},"overlap must remain rejected in full replacement");
  auto changed=p;changed.hamaker*=1.2;changed.adhesionWork*=1.1;
  changed.freeCmcRepulsionPressure*=7.;changed.freeCmcInnerRepulsionWork*=3.;
  for(double h:{1.8e-9,2.e-9,2.4e-9,4.e-9,8.e-9,20.e-9})
    require(same(evaluate(atGap(h,3),p),evaluate(atGap(h,3),changed)),
            "full net branch is independent of replaced legacy amplitudes");
}

void conservativeGeometry(){
  const auto p=parameters();
  for(int geometry:{3,4})for(double h:{1.8e-9,2.2e-9,2.7e-9,3.e-9,3.4e-9,6.e-9,8.5e-9,9.e-9}){
    const auto b=atGap(h,geometry);const auto r=evaluate(b,p);
    for(int k=0;k<9;++k){
      // The C2 joins have zero radial force: shorten the stencil there so
      // the cubic truncation error does not mask the small curvature force.
      const double joinScale=(h==3.e-9||h==9.e-9)?.2:1.;
      const double epsilon=joinScale*(k<3?1.e-13:1.e-7);
      const double numerical=(evaluate(perturb(b,k,epsilon),p).energy
          -evaluate(perturb(b,k,-epsilon),p).energy)/(2.*epsilon);
      near(derivative(r,k),numerical,k<3?2.e-20:2.e-26,8.e-5,"full geometry energy gradient "+std::to_string(k));
    }
    const auto moment=g::sub(g::add(r.torqueI,r.torqueJ),
        g::cross(g::sub(b.second.position,b.first.position),r.forceI));
    nearVector(moment,{},2.e-26+2.e-10*g::norm(r.torqueI),0.,"internal angular momentum balance");
    const auto swapped=evaluate({b.second,b.first},p);
    near(swapped.energy,r.energy,1.e-30,2.e-8,"particle-exchange energy");
    nearVector(swapped.forceI,g::scale(r.forceI,-1.),1.e-21,2.e-7,"particle-exchange force");
    nearVector(swapped.torqueI,r.torqueJ,1.e-27,2.e-7,"particle-exchange torque");
    const auto rotation=g::rotationIncrement({.71,-.36,.43});auto moved=b;
    for(auto* body:{&moved.first,&moved.second}){
      body->position=g::add(g::mul(rotation,body->position),{4.e-6,-2.e-6,3.e-6});
      body->rotation=g::multiply(rotation,body->rotation);
    }
    const auto rotated=evaluate(moved,p);
    near(rotated.energy,r.energy,1.e-30,2.e-8,"rigid-motion invariant energy");
    nearVector(rotated.forceI,g::mul(rotation,r.forceI),1.e-21,2.e-7,"rigid-motion covariant force");
    nearVector(rotated.torqueI,g::mul(rotation,r.torqueI),1.e-27,2.e-7,"rigid-motion covariant torque");
  }
}

void blending(){
  auto old=parameters(0.),dormant=old,net=parameters();
  dormant.cmcNetContactForce*=4.;dormant.cmcNetBarrierForce*=3.;
  dormant.cmcNetAttractionRange*=2.;dormant.cmcNetRepulsionRange*=2.;
  dormant.cmcNetReferenceLength*=.7;
  for(double h:{1.8e-9,2.3e-9,4.e-9,8.e-9,20.e-9,450.e-9}){
    const auto b=atGap(h,3);const auto a=evaluate(b,old),z=evaluate(b,net);
    require(same(a,evaluate(b,dormant)),"zero blend ignores all new scales bitwise");
    for(double weight:{.1,.5,.9}){
      auto p=old;p.cmcNetBlend=weight;const auto r=evaluate(b,p);
      near(r.ua,(1.-weight)*a.ua+weight*z.ua,1.e-32,2.e-14,"attractive energy blend");
      near(r.ur,(1.-weight)*a.ur+weight*z.ur,1.e-32,2.e-14,"repulsive energy blend");
      nearVector(r.forceI,g::add(g::scale(a.forceI,1.-weight),g::scale(z.forceI,weight)),1.e-23,2.e-12,"force blend");
      nearVector(r.torqueI,g::add(g::scale(a.torqueI,1.-weight),g::scale(z.torqueI,weight)),1.e-29,2.e-12,"first torque blend");
      nearVector(r.torqueJ,g::add(g::scale(a.torqueJ,1.-weight),g::scale(z.torqueJ,weight)),1.e-29,2.e-12,"second torque blend");
    }
  }
  for(bool local:{false,true}){
    old.surfaceAdhesion=false;old.freeCmcRepulsionPressure=0.;old.freeCmcInnerRepulsionWork=0.;
    old.localGapFraction=local?.1:0.;dormant=old;dormant.cmcNetContactForce*=3.;
    require(same(evaluate(atGap(2.3e-9,3),old),evaluate(atGap(2.3e-9,3),dormant)),
            "zero blend preserves pure RE2 and legacy local-gap models");
  }
}

void smoothJoins(){
  const auto p=parameters();
  for(double join:{3.e-9,9.e-9}){
    double previousForce=0.,previousStiffness=0.;
    for(double epsilon:{2.e-12,1.e-12,.5e-12}){
      const auto left=evaluate(atGap(join-epsilon),p),right=evaluate(atGap(join+epsilon),p);
      const auto center=evaluate(atGap(join),p);
      near(left.energy,right.energy,2.e-27,0.,"C0 join energy");
      const double force=std::max(std::abs(outward(left)),std::abs(outward(right)));
      const double stiffness=std::max(std::abs(outward(left)-outward(center)),
                                     std::abs(outward(right)-outward(center)))/epsilon;
      if(previousForce>0.){
        require(force<.27*previousForce,"force approaches its zero join quadratically");
        require(stiffness<.54*previousStiffness,"force derivative approaches zero continuously");
      }
      previousForce=force;previousStiffness=stiffness;
    }
  }
}

void validation(){
  const auto base=parameters();
  const std::vector<std::function<void(g::PairParameters&)>> invalid={
    [](auto& p){p.cmcNetBlend=-.1;},[](auto& p){p.cmcNetBlend=1.1;},
    [](auto& p){p.cmcNetBlend=std::numeric_limits<double>::quiet_NaN();},
    [](auto& p){p.cmcNetContactForce=-1.;},[](auto& p){p.cmcNetContactForce=0.;},
    [](auto& p){p.cmcNetBarrierForce=0.;},[](auto& p){p.cmcNetBarrierForce=std::numeric_limits<double>::infinity();},
    [](auto& p){p.cmcNetAttractionRange=0.;},[](auto& p){p.cmcNetRepulsionRange=-1.;},
    [](auto& p){p.cmcNetReferenceLength=0.;},[](auto& p){p.cmcNetRepulsionRange=p.switchGap;},
    [](auto& p){p.surfaceAdhesion=false;},[](auto& p){p.localGapFraction=.1;},
    [](auto& p){p.cohesionRetention=.5;},[](auto& p){p.contactGap=2.1e-9;},
    [](auto& p){p.cmcNetAttractionRange=std::numeric_limits<double>::denorm_min();}
  };
  for(const auto& change:invalid){auto p=base;change(p);rejects([&]{g::validatePairParameters(p);},"invalid net parameters accepted");}
  auto p=base;p.contactGap=p.roughnessGap;g::validatePairParameters(p);
}
}
int main(){
  const std::vector<std::pair<std::string,std::function<void()>>> tests={
    {"aligned physical targets",alignedPhysicalTargets},{"support and continuation",supportAndContinuation},
    {"conservative geometry",conservativeGeometry},{"blending",blending},
    {"smooth C2 joins",smoothJoins},{"validation",validation}
  };
  try{for(const auto& test:tests){test.second();std::cout<<"PASS "<<test.first<<'\n';}}
  catch(const std::exception& error){std::cerr<<"FAIL "<<error.what()<<'\n';return 1;}
}
