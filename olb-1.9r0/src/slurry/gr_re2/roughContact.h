#ifndef SLURRY_GR_RE2_GRAPHITE_ROUGH_CONTACT_H
#define SLURRY_GR_RE2_GRAPHITE_ROUGH_CONTACT_H

#include "particleMath.h"

// SLURRY SCOPE BEGIN
namespace slurry { namespace gr_re2 {

namespace graphite {

// Effective asperity-contact law.  The particle solver supplies the normal
// reaction; no normal penalty spring or overlap-dependent force is added here.
struct RoughContactSettings {
  bool enabled=false;
  double gap=2.e-9;
  double friction=.5;
  double tangentialStiffness=9.;              // N/m
  double rollingLength=100.e-9;               // m; not the asperity height
  double rollingYieldAngle=.01;              // rad
};

// One state per persistent, ordered particle-ID pair.  LE image changes must
// not change that identity.  Only an accepted particle substep commits a state.
struct RoughContactState {
  bool active=false,sliding=false,rolling=false;
  Vec3 normal{},elasticSlip{},elasticRoll{};
  double rollingCap=0.,rollingStiffness=0.;
  double normalLoad=0.,elasticEnergy=0.;
  double plasticSlipWork=0.,plasticRollWork=0.,releasedEnergy=0.,stepDuration=0.;
  Vec3 tangentForce{},rollingTorque{},leverI{},leverJ{};
};

struct RoughContactResult {
  Vec3 forceI{},normalForce{},tangentForce{},torqueI{},torqueJ{},rollingTorque{};
  Vec3 leverI{},leverJ{};
  Vec3 trialSlip{},trialRoll{};
  RoughContactState candidateState{};
  bool active=false,sliding=false,rolling=false;
  double plasticSlipWork=0.,plasticRollWork=0.,releasedEnergy=0.,elasticEnergy=0.;
  // Derivatives at fixed contact geometry and transported history.  They are
  // exact on the current stick/slip branch, suitable for contact preconditioning
  // and a semismooth generalized tangent.  Rotation of the history and changing
  // geometry must be handled by the outer geometric Newton linearization.
  Mat3 dForceDSlipVelocity{},dTorqueDRollVelocity{};
  Vec3 normalForceDerivative{},tangentForceLoadDerivative{};
};

namespace rough_detail {
inline Vec3 tangent(const Vec3& v,const Vec3& n) {
  return sub(v,scale(n,dot(n,v)));
}
inline Mat3 tangentProjector(const Vec3& n) {
  Mat3 p=identity3;for(int i=0;i<3;++i)for(int j=0;j<3;++j)p[3*i+j]-=n[i]*n[j];return p;
}
inline Mat3 alignNormals(const Vec3& oldNormal,const Vec3& newNormal) {
  const Vec3 a=normalized(oldNormal),b=normalized(newNormal);
  const Vec3 axis=cross(a,b);
  const double sine=norm(axis),cosine=std::max(-1.,std::min(1.,dot(a,b)));
  if(sine>1.e-14)return rotationIncrement(scale(axis,std::atan2(sine,cosine)/sine));
  if(cosine>0.)return identity3;
  // Antiparallel normals have no unique minimum rotation; choose a deterministic
  // tangent axis.  Normal evolution through ordinary accepted steps is smooth.
  constexpr double pi=3.1415926535897932384626433832795;
  const Vec3 axisPi=normalized(cross(a,std::abs(a[0])<.8?Vec3{1.,0.,0.}:Vec3{0.,1.,0.}));
  return rotationIncrement(scale(axisPi,pi));
}
inline Mat3 historyTransport(const RoughContactState& old,const Vec3& normal,
                             const Vec3& averageOmega,double dt) {
  // First carry the frame with common particle rotation, including spin about
  // the contact normal.  Then minimally correct any geometric normal change.
  // For a common finite rigid rotation this reproduces that rotation exactly.
  const Mat3 common=rotationIncrement(scale(averageOmega,dt));
  const Mat3 correction=alignNormals(mul(common,old.normal),normal);
  return multiply(correction,common);
}
struct ReturnMap {
  Vec3 load{},elastic{};
  Mat3 velocityDerivative{};
  Vec3 capDerivative{};
  double plasticWork=0.,energy=0.;
  bool yielded=false;
};
inline ReturnMap returnMap(const Vec3& trial,double stiffness,double cap,
                           const Vec3& normal,double dt) {
  ReturnMap out;
  if(!(stiffness>0.))return out;
  const Vec3 trialLoad=scale(trial,-stiffness);
  const double magnitude=norm(trialLoad);
  if(!(cap>0.)) {out.yielded=magnitude>0.;return out;}
  const Mat3 p=tangentProjector(normal);
  if(magnitude<=cap) {
    out.load=trialLoad;out.elastic=trial;
    for(int k=0;k<9;++k)out.velocityDerivative[k]=-stiffness*dt*p[k];
  } else {
    out.yielded=true;
    const Vec3 direction=scale(trialLoad,1./magnitude);
    out.load=scale(direction,cap);out.elastic=scale(out.load,-1./stiffness);
    out.capDerivative=direction;
    const double factor=-stiffness*dt*cap/magnitude;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
      out.velocityDerivative[3*i+j]=factor*(p[3*i+j]-direction[i]*direction[j]);
    out.plasticWork=cap*norm(sub(trial,out.elastic));
  }
  out.energy=.5*stiffness*dot(out.elastic,out.elastic);
  return out;
}
} // namespace rough_detail

inline RoughContactResult roughContact(
    const Body& bi,const Body& bj,const Vec3& normal,
    const Vec3& surfaceLeverI,const Vec3& surfaceLeverJ,double normalLoad,
    double adhesiveBirthForce,double dt,const RoughContactState& old,
    const RoughContactSettings& settings,bool engaged=true) {
  RoughContactResult out;
  if(!settings.enabled || !engaged) {
    // Release of stored contact elasticity is explicitly accounted as a model
    // loss.  It is not added again as pair-force mechanical work.
    if(old.active)out.releasedEnergy=old.elasticEnergy;
    out.candidateState.releasedEnergy=out.releasedEnergy;
    out.candidateState.stepDuration=dt;
    return out;
  }
  if(!(settings.gap>0.) || !(settings.friction>=0.) ||
     !(settings.tangentialStiffness>=0.) || !(settings.rollingLength>=0.) ||
     !(settings.rollingYieldAngle>0.) || !(dt>=0.) || !std::isfinite(dt) ||
     !std::isfinite(normalLoad) || !std::isfinite(adhesiveBirthForce))
    throw std::invalid_argument("Invalid rough-contact settings or loads");
  const Vec3 n=normalized(normal);
  const Vec3 point=scale(add(add(bi.position,surfaceLeverI),
                             add(bj.position,surfaceLeverJ)),.5);
  out.leverI=sub(point,bi.position);out.leverJ=sub(point,bj.position);
  const Vec3 vi=add(bi.velocity,cross(bi.omega,out.leverI));
  const Vec3 vj=add(bj.velocity,cross(bj.omega,out.leverJ));
  const Vec3 slip=rough_detail::tangent(sub(vi,vj),n);
  const Vec3 roll=rough_detail::tangent(sub(bi.omega,bj.omega),n);

  RoughContactState state;
  if(old.active) {
    state=old;
    const Mat3 transport=rough_detail::historyTransport(old,n,scale(add(bi.omega,bj.omega),.5),dt);
    state.elasticSlip=rough_detail::tangent(mul(transport,old.elasticSlip),n);
    state.elasticRoll=rough_detail::tangent(mul(transport,old.elasticRoll),n);
  } else {
    state.rollingCap=settings.rollingLength*std::max(0.,adhesiveBirthForce);
    state.rollingStiffness=state.rollingCap/settings.rollingYieldAngle;
  }
  out.trialSlip=add(state.elasticSlip,scale(slip,dt));
  out.trialRoll=add(state.elasticRoll,scale(roll,dt));
  const auto sliding=rough_detail::returnMap(out.trialSlip,
      settings.tangentialStiffness,settings.friction*std::max(0.,normalLoad),n,dt);
  const auto rolling=rough_detail::returnMap(out.trialRoll,
      state.rollingStiffness,state.rollingCap,n,dt);

  // n points from I to J.  Positive N therefore repels I in direction -n.
  // Negative trial reactions are algebraic Newton values only: the active-set
  // solver must release or reject them before committing an accepted state.
  out.normalForce=scale(n,-normalLoad);
  out.tangentForce=sliding.load;out.rollingTorque=rolling.load;
  out.forceI=add(out.normalForce,out.tangentForce);
  out.torqueI=add(cross(out.leverI,out.forceI),out.rollingTorque);
  out.torqueJ=scale(add(cross(out.leverJ,out.forceI),out.rollingTorque),-1.);
  out.sliding=sliding.yielded;out.rolling=rolling.yielded;out.active=true;
  out.plasticSlipWork=sliding.plasticWork;out.plasticRollWork=rolling.plasticWork;
  out.elasticEnergy=sliding.energy+rolling.energy;
  out.dForceDSlipVelocity=sliding.velocityDerivative;
  out.dTorqueDRollVelocity=rolling.velocityDerivative;
  out.normalForceDerivative=scale(n,-1.);
  out.tangentForceLoadDerivative=scale(sliding.capDerivative,normalLoad>0.?settings.friction:0.);
  state.active=true;state.normal=n;state.elasticSlip=sliding.elastic;state.elasticRoll=rolling.elastic;
  state.normalLoad=normalLoad;state.tangentForce=out.tangentForce;state.rollingTorque=out.rollingTorque;
  state.leverI=out.leverI;state.leverJ=out.leverJ;state.elasticEnergy=out.elasticEnergy;
  state.sliding=out.sliding;state.rolling=out.rolling;
  state.plasticSlipWork=out.plasticSlipWork;state.plasticRollWork=out.plasticRollWork;
  state.releasedEnergy=0.;state.stepDuration=dt;
  out.candidateState=state;
  return out;
}

struct RoughContactIncrement {Vec3 forceI{},torqueI{},torqueJ{};};
// Apply the frozen-geometry generalized contact tangent to arbitrary velocity
// and normal-load perturbations.  Equal/opposite force and common-point torque
// conventions are identical to those used for the nonlinear force itself.
inline RoughContactIncrement roughContactTangentAction(
    const RoughContactResult& base,const Vec3& dVelocityI,const Vec3& dOmegaI,
    const Vec3& dVelocityJ,const Vec3& dOmegaJ,double dNormalLoad) {
  RoughContactIncrement out;if(!base.active)return out;
  const Vec3 dSlip=sub(add(dVelocityI,cross(dOmegaI,base.leverI)),
                       add(dVelocityJ,cross(dOmegaJ,base.leverJ)));
  out.forceI=add(mul(base.dForceDSlipVelocity,dSlip),
      scale(add(base.normalForceDerivative,base.tangentForceLoadDerivative),dNormalLoad));
  const Vec3 dRolling=mul(base.dTorqueDRollVelocity,sub(dOmegaI,dOmegaJ));
  out.torqueI=add(cross(base.leverI,out.forceI),dRolling);
  out.torqueJ=scale(add(cross(base.leverJ,out.forceI),dRolling),-1.);
  return out;
}

} // namespace graphite

} } // SLURRY SCOPE END
#endif
