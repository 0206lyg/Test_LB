/* c++ -std=c++17 -O2 -Wall -Wextra -pedantic -I olb-1.9r0/src/slurry/gr_re2 \
     tests/graphite_adhesion/cmc_coordination_tests.cpp -o /tmp/cmc_coordination_tests */
#include "cmcCoordination.h"
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace g=slurry::gr_re2::graphite;
namespace {
void require(bool condition,const std::string& message){if(!condition)throw std::runtime_error(message);}
void near(double actual,double expected,double absolute,double relative,const std::string& message){
  if(!std::isfinite(actual)||!std::isfinite(expected)||std::abs(actual-expected)>absolute+relative*std::abs(expected)){
    std::ostringstream s;s<<std::setprecision(17)<<message<<": actual="<<actual<<", expected="<<expected;
    throw std::runtime_error(s.str());
  }
}
void nearVector(const g::Vec3& a,const g::Vec3& b,double absolute,double relative,const std::string& message){
  for(int k=0;k<3;++k)near(a[k],b[k],absolute,relative,message+" component "+std::to_string(k));
}
g::PairParameters parameters(double blend=1.){
  g::PairParameters p;p.surfaceAdhesion=true;p.sigma=.4197e-9;
  p.adhesionWork=.007666863221834307;p.freeCmcRepulsionPressure=52143.97983709085;
  p.freeCmcInnerRepulsionWork=.0071170133084142865;p.cmcNetBlend=blend;
  p.cmcNetContactForce=450.e-12;p.cmcCoordinationEnabled=true;return p;
}
g::Vec3 support(const g::Body& b,const g::Vec3& n){
  g::Vec3 a2{};for(int k=0;k<3;++k)a2[k]=b.axes[k]*b.axes[k];
  const auto qn=g::mul(g::rotatedDiagonal(b.rotation,a2),n);
  return g::scale(qn,1./std::sqrt(g::dot(n,qn)));
}
std::vector<g::Body> star(const std::vector<double>& gaps,bool oblique=false){
  std::vector<g::Body> b(gaps.size()+1);
  b[0].position={1.e-6,-2.e-6,3.e-6};
  const g::Vec3 directions[]={{0.,0.,1.},{0.,0.,-1.},{1.,0.,0.}};
  // Three separated edge neighbours avoid neighbour-neighbour overlap even
  // for thin, independently tilted platelets.
  const g::Vec3 obliqueDirections[]={{1.,0.,0.},{-.5,std::sqrt(.75),0.},{-.5,-std::sqrt(.75),0.}};
  if(oblique){b[0].axes={1.65e-6,1.4e-6,.2e-6};b[0].rotation=g::rotationIncrement({.09,-.07,.04});}
  const auto common=g::rotationIncrement({.13,.08,-.11});
  for(std::size_t k=1;k<b.size();++k){
    if(oblique){
      b[k].axes={1.55e-6+.02e-6*k,1.4e-6+.03e-6*k,.21e-6};
      b[k].rotation=g::rotationIncrement({.02*static_cast<double>(k),-.025*static_cast<double>(k),.015});
    }
    const auto n=oblique?g::mul(common,obliqueDirections[k-1]):directions[k-1];
    b[k].position=g::add(b[0].position,g::add(g::add(support(b[0],n),support(b[k],n)),g::scale(n,gaps[k-1])));
  }
  return b;
}
struct State {
  std::vector<g::CoordinationPair> pairs;
  std::vector<g::Vec3> force,torque;
  double energy=0.;
};
State evaluate(const std::vector<g::Body>& b,const g::PairParameters& p,bool apply=true,bool reverse=false){
  State s;s.force.resize(b.size());s.torque.resize(b.size());
  for(std::size_t i=0;i<b.size();++i)for(std::size_t j=i+1;j<b.size();++j)
    s.pairs.push_back({i,j,g::evaluatePair(b[i],b[j],p)});
  if(reverse)std::reverse(s.pairs.begin(),s.pairs.end());
  if(apply)g::applyCmcCoordination(b.size(),s.pairs,p);
  for(const auto& pair:s.pairs){
    const auto& r=pair.result;s.energy+=r.energy;
    s.force[pair.i]=g::add(s.force[pair.i],r.forceI);s.force[pair.j]=g::sub(s.force[pair.j],r.forceI);
    s.torque[pair.i]=g::add(s.torque[pair.i],r.torqueI);s.torque[pair.j]=g::add(s.torque[pair.j],r.torqueJ);
  }
  return s;
}
bool samePhysics(const g::PairResult& a,const g::PairResult& b){
  return a.ua==b.ua&&a.ur==b.ur&&a.energy==b.energy&&a.forceI==b.forceI&&a.torqueI==b.torqueI&&a.torqueJ==b.torqueJ
      &&a.forceAttractiveI==b.forceAttractiveI&&a.forceRepulsiveI==b.forceRepulsiveI
      &&a.torqueAttractiveI==b.torqueAttractiveI&&a.torqueAttractiveJ==b.torqueAttractiveJ
      &&a.torqueRepulsiveI==b.torqueRepulsiveI&&a.torqueRepulsiveJ==b.torqueRepulsiveJ;
}

void isolatedAndDisabled(){
  auto p=parameters();
  for(double h:{1.8e-9,2.e-9,2.4e-9,3.e-9,6.e-9,20.e-9}){
    const auto b=star({h},true);const auto raw=evaluate(b,p,false),corrected=evaluate(b,p);
    require(samePhysics(raw.pairs[0].result,corrected.pairs[0].result),"isolated pair is unchanged bitwise");
  }
  const auto b=star({2.1e-9,2.2e-9,2.3e-9},true);
  p.cmcCoordinationEnabled=false;
  auto raw=evaluate(b,p,false),corrected=evaluate(b,p);
  for(std::size_t k=0;k<raw.pairs.size();++k){
    require(samePhysics(raw.pairs[k].result,corrected.pairs[k].result),"disabled correction is bitwise inert");
    require(!corrected.pairs[k].result.cmcCoordinationApplied,"disabled correction does not mark metadata");
  }
  p.cmcCoordinationEnabled=true;p.cmcNetBlend=0.;
  raw=evaluate(b,p,false);corrected=evaluate(b,p);
  for(std::size_t k=0;k<raw.pairs.size();++k)
    require(samePhysics(raw.pairs[k].result,corrected.pairs[k].result),"zero net blend preserves all legacy interactions");
  p=parameters();p.cmcCoordinationFloor=1.;
  raw=evaluate(b,p,false);corrected=evaluate(b,p);
  for(std::size_t k=0;k<raw.pairs.size();++k)
    require(samePhysics(raw.pairs[k].result,corrected.pairs[k].result),"unit floor is exactly inert");
}

void coordinationTargets(){
  const auto p=parameters();
  // At full occupancy two neighbours retain the isolated strength; adding a
  // third weakens all incident bonds symmetrically, not a history-selected one.
  const auto two=star({1.999e-9,1.999e-9});
  const auto twoRaw=evaluate(two,p,false),twoCorrected=evaluate(two,p);
  for(std::size_t k=0;k<twoRaw.pairs.size();++k)
    require(samePhysics(twoRaw.pairs[k].result,twoCorrected.pairs[k].result),"two neighbours retain full attraction");
  const auto three=star({1.999e-9,1.999e-9,1.999e-9});
  const auto raw=evaluate(three,p,false),corrected=evaluate(three,p);
  for(std::size_t k=0;k<raw.pairs.size();++k){
    const auto& before=raw.pairs[k];const auto& after=corrected.pairs[k].result;
    if(before.i==0){
      near(after.ua,.1*before.result.ua,1.e-31,2.e-14,"degree-three central attenuation");
      nearVector(after.forceAttractiveI,g::scale(before.result.forceAttractiveI,.1),1.e-23,2.e-12,"plateau attraction scales by floor");
    }
    require(after.ur==before.result.ur&&after.forceRepulsiveI==before.result.forceRepulsiveI,
            "coordination leaves repulsion unchanged");
  }
  // A complete four-body graph at full occupancy exercises both endpoint
  // factors without depending on a realizable close-packed oblate geometry.
  std::vector<g::CoordinationPair> graph;
  for(std::size_t i=0;i<4;++i)for(std::size_t j=i+1;j<4;++j){
    g::PairResult r;r.cmcCoordinationOccupancy=1.;r.ua=r.energy=r.cmcNetAttractiveEnergy=-1.e-20;
    graph.push_back({i,j,r});
  }
  g::applyCmcCoordination(4,graph,p);
  for(const auto& pair:graph)near(pair.result.energy,-1.e-22,1.e-34,2.e-14,"both endpoint factors multiply");
}

void gradientCase(std::vector<g::Body> b,g::PairParameters p){
  const auto s=evaluate(b,p);const auto raw=evaluate(b,p,false);
  require(std::abs(s.energy-raw.energy)>1.e-24,"fixture exercises nontrivial coordination");
  bool occupancyDerivative=false;
  for(const auto& pair:s.pairs)occupancyDerivative|=pair.result.cmcCoordinationGapDerivative!=0.;
  require(occupancyDerivative,"fixture exercises occupancy gradients");
  for(std::size_t i=0;i<b.size();++i)for(int k=0;k<6;++k){
    const double epsilon=k<3?2.e-14:2.e-8;auto plus=b,minus=b;
    if(k<3){plus[i].position[k]+=epsilon;minus[i].position[k]-=epsilon;}
    else{g::Vec3 angle{};angle[k-3]=epsilon;plus[i].rotation=g::rotateLaboratory(plus[i].rotation,angle);
      angle[k-3]=-epsilon;minus[i].rotation=g::rotateLaboratory(minus[i].rotation,angle);}
    const double numerical=-(evaluate(plus,p).energy-evaluate(minus,p).energy)/(2.*epsilon);
    near(k<3?s.force[i][k]:s.torque[i][k-3],numerical,k<3?3.e-18:3.e-24,2.e-4,
         "many-body energy derivative body "+std::to_string(i)+" coordinate "+std::to_string(k));
  }
  g::Vec3 force{},moment{};
  for(std::size_t i=0;i<b.size();++i){force=g::add(force,s.force[i]);moment=g::add(moment,g::add(s.torque[i],g::cross(b[i].position,s.force[i])));}
  nearVector(force,{},1.e-24,0.,"total force balance");nearVector(moment,{},2.e-24,0.,"total angular momentum balance");
  // The corrected antisymmetric pair decomposition must reproduce affine
  // work; this is the same displacement/force contraction used in the virial.
  double predicted=0.;
  for(const auto& pair:s.pairs){
    const auto dr=g::sub(b[pair.j].position,b[pair.i].position);
    predicted+=pair.result.forceI[0]*dr[1];
  }
  const double epsilon=2.e-8;auto plus=b,minus=b;
  for(std::size_t i=0;i<b.size();++i){plus[i].position[0]+=epsilon*b[i].position[1];minus[i].position[0]-=epsilon*b[i].position[1];}
  const double numerical=(evaluate(plus,p).energy-evaluate(minus,p).energy)/(2.*epsilon);
  near(predicted,numerical,3.e-24,3.e-4,"pair virial equals affine energy derivative");
}
void conservativeManyBody(){
  gradientCase(star({2.13e-9,2.24e-9,2.39e-9},true),parameters());
  gradientCase(star({2.13e-9,2.24e-9,2.39e-9},true),parameters(.4));
  auto p=parameters();p.cmcCoordinationStart=.2;p.cmcCoordinationEnd=1.;
  gradientCase(star({2.18e-9,2.33e-9},true),p);
}

void permutationAndObjectivity(){
  const auto p=parameters();const auto b=star({2.13e-9,2.24e-9,2.39e-9},true);
  const auto original=evaluate(b,p),reverse=evaluate(b,p,true,true);
  near(reverse.energy,original.energy,1.e-31,2.e-14,"edge-order invariant energy");
  for(std::size_t i=0;i<b.size();++i){
    nearVector(reverse.force[i],original.force[i],1.e-23,2.e-12,"edge-order invariant force");
    nearVector(reverse.torque[i],original.torque[i],1.e-29,2.e-12,"edge-order invariant torque");
  }
  const std::size_t order[]={2,0,3,1};std::vector<g::Body> reordered;
  for(auto i:order)reordered.push_back(b[i]);
  const auto permuted=evaluate(reordered,p);
  near(permuted.energy,original.energy,2.e-28,2.e-8,"body-label invariant energy");
  for(std::size_t i=0;i<b.size();++i){
    nearVector(permuted.force[i],original.force[order[i]],2.e-18,2.e-7,"body-label invariant force");
    nearVector(permuted.torque[i],original.torque[order[i]],2.e-24,2.e-7,"body-label invariant torque");
  }
  const auto rotation=g::rotationIncrement({.27,-.16,.21});auto moved=b;
  for(auto& body:moved){body.position=g::add(g::mul(rotation,body.position),{4.e-6,-2.e-6,3.e-6});body.rotation=g::multiply(rotation,body.rotation);}
  const auto transformed=evaluate(moved,p);near(transformed.energy,original.energy,2.e-28,2.e-8,"rigid-motion energy");
  for(std::size_t i=0;i<b.size();++i){
    nearVector(transformed.force[i],g::mul(rotation,original.force[i]),2.e-18,2.e-7,"rigid-motion force");
    nearVector(transformed.torque[i],g::mul(rotation,original.torque[i]),2.e-24,2.e-7,"rigid-motion torque");
  }
}

void smoothMembershipAndCutoff(){
  const auto p=parameters();
  for(int k=1;k<=100;++k){
    const auto s=evaluate(star({p.roughnessGap+k*2.e-23}),p);
    const double q=s.pairs[0].result.cmcCoordinationOccupancy;
    require(q>=0.&&q<=1.,"roundoff near contact cannot exceed unit occupancy");
  }
  for(int exponent=-24;exponent<=-18;++exponent)for(double factor:{1.,2.3,7.1}){
    const auto s=evaluate(star({p.roughnessGap+factor*std::pow(10.,exponent)}),p);
    const double q=s.pairs[0].result.cmcCoordinationOccupancy;
    require(q>=0.&&q<=1.,"sub-picometre openings retain bounded occupancy");
  }
  for(int exponent=-24;exponent<=-3;++exponent)for(double factor:{1.,2.3,7.1}){
    const auto g0=g::coordination_detail::attenuation(p.cmcCoordinationStart+factor*std::pow(10.,exponent),p);
    require(g0.value>=p.cmcCoordinationFloor&&g0.value<=1.&&g0.derivative<=0.,
            "attenuation near transition start stays bounded and nonincreasing");
  }
  for(double join:{2.e-9,3.e-9}){
    double lastValue=0.,lastSlope=0.;
    for(double epsilon:{2.e-12,1.e-12,.5e-12}){
      const auto left=evaluate(star({join-epsilon}),p,false).pairs[0].result;
      const auto right=evaluate(star({join+epsilon}),p,false).pairs[0].result;
      const double value=std::abs(left.cmcCoordinationOccupancy-right.cmcCoordinationOccupancy);
      const double slope=std::max(std::abs(left.cmcCoordinationGapDerivative),std::abs(right.cmcCoordinationGapDerivative));
      if(lastValue>0.)require(value<.14*lastValue,"occupancy reaches endpoint cubically");
      if(lastSlope>0.)require(slope<.27*lastSlope,"occupancy derivative reaches endpoint quadratically");
      lastValue=value;lastSlope=slope;
    }
  }
  const auto b=star({2.1e-9,2.2e-9,3.001e-9});const auto withFar=evaluate(b,p);
  const auto without=evaluate(std::vector<g::Body>(b.begin(),b.begin()+3),p);
  for(std::size_t i=0;i<3;++i){
    // The third particle's outer repulsion remains physical; only the net
    // attractive branch and its coordination correction vanish at 3 nm.
    g::Vec3 attractive{};
    for(const auto& pair:withFar.pairs){if(pair.i==i)attractive=g::add(attractive,pair.result.forceAttractiveI);if(pair.j==i)attractive=g::sub(attractive,pair.result.forceAttractiveI);}
    g::Vec3 expected{};
    for(const auto& pair:without.pairs){if(pair.i==i)expected=g::add(expected,pair.result.forceAttractiveI);if(pair.j==i)expected=g::sub(expected,pair.result.forceAttractiveI);}
    nearVector(attractive,expected,1.e-24,2.e-13,"outside attraction range does not enter coordination");
  }
}

void invalidAndInactive(){
  auto p=parameters();
  const std::vector<std::function<void(g::PairParameters&)>> invalid={
    [](auto& q){q.cmcCoordinationStart=-1.;},[](auto& q){q.cmcCoordinationEnd=q.cmcCoordinationStart;},
    [](auto& q){q.cmcCoordinationFloor=-.1;},[](auto& q){q.cmcCoordinationFloor=1.1;},
    [](auto& q){q.cmcCoordinationEnd=std::numeric_limits<double>::infinity();}
  };
  for(const auto& change:invalid){auto bad=p;change(bad);bool caught=false;try{g::validatePairParameters(bad);}catch(const std::domain_error&){caught=true;}require(caught,"invalid coordination parameters rejected");}
  auto state=evaluate(star({2.2e-9}),p);bool caught=false;
  try{g::applyCmcCoordination(2,state.pairs,p);}catch(const std::logic_error&){caught=true;}
  require(caught,"double application is rejected");
  std::vector<g::CoordinationPair> emptyFar{{0,1,{}}};g::applyCmcCoordination(2,emptyFar,p);
  require(emptyFar[0].result.energy==0.&&emptyFar[0].result.forceI==g::Vec3{},"inactive far-pair entries are accepted");
  p.cmcCoordinationFloor=0.;g::validatePairParameters(p);
}
}
int main(){
  const std::vector<std::pair<std::string,std::function<void()>>> tests={
    {"isolated and disabled compatibility",isolatedAndDisabled},{"coordination strength targets",coordinationTargets},
    {"many-body energy and virial derivatives",conservativeManyBody},{"permutation and objectivity",permutationAndObjectivity},
    {"smooth membership and cutoff",smoothMembershipAndCutoff},{"validation and inactive pairs",invalidAndInactive}
  };
  try{for(const auto& test:tests){test.second();std::cout<<"PASS "<<test.first<<'\n';}}
  catch(const std::exception& error){std::cerr<<"FAIL "<<error.what()<<'\n';return 1;}
}
