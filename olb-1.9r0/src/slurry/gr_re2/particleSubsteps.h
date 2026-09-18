#ifndef SLURRY_GR_RE2_GRAPHITE_PARTICLE_SUBSTEPS_H
#define SLURRY_GR_RE2_GRAPHITE_PARTICLE_SUBSTEPS_H

#include "reSquaredPotential.h"
#include "nearFieldResistance.h"
#include "roughContact.h"
#include <functional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <fstream>
#include <iomanip>
#include <deque>
#ifdef SLURRY_USE_PETSC
#include <petscsnes.h>
#endif

// SLURRY SCOPE BEGIN
namespace slurry { namespace gr_re2 {

namespace graphite {

struct ParticleStepSettings {
  PairParameters pair;
  NearFieldSettings nearField;
  RoughContactSettings rough;
  Vec3 box{10.e-6,10.e-6,10.e-6};
  double shearRate=100.;
  int maxSubsteps=32;
  int maxNewtonIterations=12;
  int maxKrylovIterations=120;
  std::string solverBackend="petsc";
  std::string solverDiagnosticsPrefix;
  int maxLineSearch=14;
  double relativeTolerance=1.e-7;
  double forceAbsoluteTolerance=1.e-15;
  double torqueAbsoluteTolerance=1.65e-21;
  double contactGapTolerance=1.e-12;
  double finiteDifferenceStep=1.e-7;
};

struct ParticleStepDiagnostics {
  double minGap=std::numeric_limits<double>::infinity();
  double maxForce=0.,energyAtEnd=0.,lubricationDissipation=0.;
  // Symmetric Cauchy stress moments, units N*m. Divide by total box volume.
  // Pair moments are -sym(F_i outer (x_i-x_j_image)).
  Mat3 pairMoment{},attractiveMoment{},repulsiveMoment{},lubricationMoment{};
  Mat3 contactNormalMoment{},contactTangentialMoment{};
  double contactDissipation=0.,elasticContactEnergy=0.;
  double maxForceResidualRatio=0.,maxTorqueResidualRatio=0.,contactGapViolation=0.;
  int contacts=0,slidingContacts=0,rollingContacts=0;
  int substeps=0,newtonIterations=0,krylovIterations=0,residualEvaluations=0;
  int frictionBranchAttempts=0,frictionBranchCorrections=0;
  int contactStateUpdates=0,contactActivations=0,contactReleases=0;
  std::size_t activePairs=0;
};

using PersistentContactState=std::vector<RoughContactState>;

namespace particle_detail {
using Vector=std::vector<double>;
inline double inner(const Vector& a,const Vector& b) {
  double s=0.;for(std::size_t k=0;k<a.size();++k)s+=a[k]*b[k];return s;
}
inline double length(const Vector& v) {return std::sqrt(inner(v,v));}
inline double maxAbs(const Vector& v) {
  double s=0.;for(double x:v)s=std::max(s,std::abs(x));return s;
}
inline void addMoment(Mat3& m,const Vec3& force,const Vec3& rij,double weight=1.) {
  for(int a=0;a<3;++a)for(int b=0;b<3;++b)
    m[3*a+b]-=.5*weight*(force[a]*rij[b]+force[b]*rij[a]);
}
inline void addMatrix(Mat3& dst,const Mat3& src,double weight) {
  for(int k=0;k<9;++k)dst[k]+=weight*src[k];
}
inline double radius(const Body& b) {return *std::max_element(b.axes.begin(),b.axes.end());}
inline Vec3 worldMomentum(const Body& b) {
  return mul(rotatedDiagonal(b.rotation,b.inertiaBody),b.omega);
}

// Find the nearest image in the sliding lattice.  Rounding y alone is not a
// Euclidean nearest-image rule in a sheared box, so inspect adjacent y images.
inline Body closestImage(const Body& i,const Body& j,double time,
                         const ParticleStepSettings& s) {
  const double shift=std::remainder(s.shearRate*s.box[1]*time,s.box[0]);
  const int ny0=static_cast<int>(std::llround((i.position[1]-j.position[1])/s.box[1]));
  double best=std::numeric_limits<double>::infinity();Body chosen=j;
  for(int ny=ny0-1;ny<=ny0+1;++ny) {
    Body image=j;
    image.position[1]+=ny*s.box[1];
    image.position[0]+=ny*shift;
    image.velocity[0]+=ny*s.shearRate*s.box[1];
    image.position[0]+=std::nearbyint((i.position[0]-image.position[0])/s.box[0])*s.box[0];
    image.position[2]+=std::nearbyint((i.position[2]-image.position[2])/s.box[2])*s.box[2];
    const double d2=dot(sub(image.position,i.position),sub(image.position,i.position));
    if(d2<best){best=d2;chosen=image;}
  }
  return chosen;
}

inline void wrap(Body& b,double time,const ParticleStepSettings& s) {
  const double ny=std::floor(b.position[1]/s.box[1]);
  b.position[1]-=ny*s.box[1];
  b.position[0]-=ny*std::remainder(s.shearRate*s.box[1]*time,s.box[0]);
  b.velocity[0]-=ny*s.shearRate*s.box[1];
  b.position[0]-=std::floor(b.position[0]/s.box[0])*s.box[0];
  b.position[2]-=std::floor(b.position[2]/s.box[2])*s.box[2];
}

struct PairLinearization {
  std::size_t i=0,j=0,index=0;
  int slot=-1;
  double gap=0.;
  Vec3 normal{},leverI{},leverJ{};
  Mat3 resistance{},rollingResistance{};
  Vec3 forcePerNormalLoad{},torqueIPerNormalLoad{},torqueJPerNormalLoad{};
  RoughContactResult contact;
};

struct Evaluation {
  Vector residual;
  std::vector<Body> bodies;
  std::vector<GapCache> cache;
  PersistentContactState contacts;
  std::vector<double> gaps,normalLoads,forceReference,torqueReference;
  std::vector<PairLinearization> pairs;
  ParticleStepDiagnostics diagnostic;
};

struct Residual {
  const std::vector<Body>& old;
  const std::vector<Vec3>& externalForce;
  const std::vector<Vec3>& externalTorque;
  const ParticleStepSettings& settings;
  const std::vector<GapCache>& startingCache;
  const PersistentContactState& startingContacts;
  const std::vector<int>& activeSlot;
  double dt,time,lengthScale,reactionScale;
  int activeCount;
  mutable int evaluations=0;
  bool complementarity=false;
  const std::vector<unsigned char>* engagementOverride=nullptr;

  bool operator()(const Vector& q,Evaluation& e,std::string& error) const {
    ++evaluations;
    try {
      const std::size_t n=old.size(),np=n*(n-1)/2;
      e=Evaluation{};e.bodies=old;e.cache=startingCache;
      e.contacts.resize(np);e.gaps.assign(np,std::numeric_limits<double>::infinity());
      e.normalLoads.assign(np,0.);e.forceReference.resize(n);e.torqueReference.resize(n);
      e.residual.assign(6*n+activeCount,0.);
      std::vector<Vec3> force=externalForce,torque=externalTorque;
      for(std::size_t i=0;i<n;++i) {
        const Vec3 displacement{q[6*i]*lengthScale,q[6*i+1]*lengthScale,q[6*i+2]*lengthScale};
        const Vec3 angle{q[6*i+3],q[6*i+4],q[6*i+5]};
        e.bodies[i].position=add(old[i].position,displacement);
        e.bodies[i].velocity=scale(displacement,1./dt);
        e.bodies[i].omega=scale(angle,1./dt);
        e.bodies[i].rotation=rotateLaboratory(old[i].rotation,angle);
        if(!finite(e.bodies[i].position)||!finite(e.bodies[i].velocity)||!finite(e.bodies[i].omega))
          throw std::domain_error("Non-finite particle trial state");
      }
      std::size_t index=0;
      const double range=std::max(settings.pair.cutoffGap,
          std::max(settings.nearField.matchingGap,settings.rough.enabled?settings.rough.gap:0.));
      for(std::size_t i=0;i<n;++i)for(std::size_t j=i+1;j<n;++j,++index) {
        const Body image=closestImage(e.bodies[i],e.bodies[j],time+dt,settings);
        const Vec3 rij=sub(e.bodies[i].position,image.position);
        const double bound=radius(e.bodies[i])+radius(image)+range;
        if(dot(rij,rij)>bound*bound && activeSlot[index]<0 && !startingContacts[index].active)continue;
        PairResult p=evaluatePair(e.bodies[i],image,settings.pair,&e.cache[index]);
        if(!std::isfinite(p.gap)) {
          const auto gap=closestEllipsoidGap(e.bodies[i],image,&e.cache[index]);
          p.gap=gap.gap;p.normal=gap.normal;p.leverI=gap.leverI;p.leverJ=gap.leverJ;
        }
        if(!(p.gap>0.) || !std::isfinite(p.gap))
          throw std::domain_error("Particle trial state has a non-positive ellipsoid gap");
        if(settings.rough.enabled && p.gap<.8*settings.rough.gap)
          throw std::domain_error("Particle trial exceeds the rough-contact Newton gap domain");
        e.gaps[index]=p.gap;
        e.diagnostic.minGap=std::min(e.diagnostic.minGap,p.gap);
        force[i]=add(force[i],p.forceI);force[j]=sub(force[j],p.forceI);
        torque[i]=add(torque[i],p.torqueI);torque[j]=add(torque[j],p.torqueJ);
        e.diagnostic.maxForce=std::max(e.diagnostic.maxForce,norm(p.forceI));
        e.diagnostic.energyAtEnd+=p.energy;
        e.diagnostic.activePairs+=p.active?1:0;
        addMoment(e.diagnostic.pairMoment,p.forceI,rij);
        addMoment(e.diagnostic.attractiveMoment,p.forceAttractiveI,rij);
        addMoment(e.diagnostic.repulsiveMoment,p.forceRepulsiveI,rij);
        const auto lub=nearFieldResistance(e.bodies[i],image,p.gap,p.normal,p.leverI,p.leverJ,settings.nearField);
        force[i]=add(force[i],lub.forceI);force[j]=sub(force[j],lub.forceI);
        torque[i]=add(torque[i],lub.torqueI);torque[j]=add(torque[j],lub.torqueJ);
        addMoment(e.diagnostic.lubricationMoment,lub.forceI,rij);
        e.diagnostic.lubricationDissipation+=lub.dissipation;
        PairLinearization linear;
        linear.i=i;linear.j=j;linear.index=index;linear.slot=activeSlot[index];linear.gap=p.gap;
        linear.normal=p.normal;
        const Vec3 point=scale(add(add(e.bodies[i].position,p.leverI),add(image.position,p.leverJ)),.5);
        linear.leverI=sub(point,e.bodies[i].position);linear.leverJ=sub(point,image.position);
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)
          linear.resistance[3*a+b]=lub.normalResistance*p.normal[a]*p.normal[b]
            +lub.tangentialResistance*((a==b?1.:0.)-p.normal[a]*p.normal[b]);
        if(settings.rough.enabled) {
          e.diagnostic.contactGapViolation=std::max(e.diagnostic.contactGapViolation,
              std::max(0.,settings.rough.gap-p.gap));
          const bool candidate=linear.slot>=0;
          const double normalLoad=candidate?q[6*n+linear.slot]*reactionScale:0.;
          e.normalLoads[index]=normalLoad;
          const bool active=candidate && (!complementarity || (engagementOverride
              ? (*engagementOverride)[index]!=0
              : p.gap<=settings.rough.gap+settings.contactGapTolerance
                && (normalLoad>0. || startingContacts[index].active)));
          double adhesiveBirthForce=0.;
          if(active && !startingContacts[index].active) {
            Body birthImage=image;
            birthImage.position=add(birthImage.position,scale(p.normal,settings.rough.gap-p.gap));
            const auto birth=evaluatePair(e.bodies[i],birthImage,settings.pair);
            adhesiveBirthForce=std::max(0.,dot(birth.forceI,birth.normal));
          }
          auto contact=roughContact(e.bodies[i],image,p.normal,p.leverI,p.leverJ,
              normalLoad,adhesiveBirthForce,dt,startingContacts[index],settings.rough,active);
          // A candidate may be open: only its algebraic normal multiplier enters
          // momentum. It must not acquire tangential or rolling history.
          if(candidate && !active) {
            contact.normalForce=scale(p.normal,-normalLoad);contact.forceI=contact.normalForce;
            contact.leverI=linear.leverI;contact.leverJ=linear.leverJ;
            contact.torqueI=cross(linear.leverI,contact.forceI);
            contact.torqueJ=scale(cross(linear.leverJ,contact.forceI),-1.);
          }
          e.contacts[index]=contact.candidateState;
          linear.contact=contact;
          if(candidate) {
            force[i]=add(force[i],contact.forceI);force[j]=sub(force[j],contact.forceI);
            torque[i]=add(torque[i],contact.torqueI);torque[j]=add(torque[j],contact.torqueJ);
            addMoment(e.diagnostic.contactNormalMoment,contact.normalForce,rij);
            addMoment(e.diagnostic.contactTangentialMoment,contact.tangentForce,rij);
            e.diagnostic.maxForce=std::max(e.diagnostic.maxForce,norm(contact.forceI));
            e.diagnostic.contacts+=active?1:0;
            e.diagnostic.slidingContacts+=contact.sliding?1:0;
            e.diagnostic.rollingContacts+=contact.rolling?1:0;
            e.residual[6*n+linear.slot]=(p.gap-settings.rough.gap)/lengthScale;
          }
          e.diagnostic.contactDissipation+=(contact.plasticSlipWork+contact.plasticRollWork+contact.releasedEnergy)/dt;
          e.diagnostic.elasticContactEnergy+=contact.elasticEnergy;
          // Filled below with the analytic constitutive tangents supplied by the
          // return map.  The Newton product also retains geometry derivatives.
          if(candidate) {
            for(int k=0;k<9;++k) {
              linear.resistance[k]-=contact.dForceDSlipVelocity[k];
              linear.rollingResistance[k]-=contact.dTorqueDRollVelocity[k];
            }
            linear.forcePerNormalLoad=add(scale(p.normal,-1.),contact.tangentForceLoadDerivative);
            linear.torqueIPerNormalLoad=cross(linear.leverI,linear.forcePerNormalLoad);
            linear.torqueJPerNormalLoad=scale(cross(linear.leverJ,linear.forcePerNormalLoad),-1.);
          }
        }
        if(lub.active || linear.slot>=0)e.pairs.push_back(linear);
      }
      for(std::size_t i=0;i<n;++i) {
        const double forceScale=old[i].mass*lengthScale/(dt*dt);
        const double inertiaMax=*std::max_element(old[i].inertiaBody.begin(),old[i].inertiaBody.end());
        const double torqueScale=inertiaMax/(dt*dt);
        const Vec3 inertialForce=scale(sub(e.bodies[i].velocity,old[i].velocity),old[i].mass/dt);
        const Vec3 inertialTorque=scale(sub(worldMomentum(e.bodies[i]),worldMomentum(old[i])),1./dt);
        const Vec3 rf=sub(inertialForce,force[i]),rt=sub(inertialTorque,torque[i]);
        double fError=0.,tError=0.,fCurrent=0.,tCurrent=0.;
        for(int c=0;c<3;++c) {
          e.residual[6*i+c]=rf[c]/forceScale;e.residual[6*i+3+c]=rt[c]/torqueScale;
          fError=std::max(fError,std::abs(rf[c]));tError=std::max(tError,std::abs(rt[c]));
          fCurrent=std::max(fCurrent,std::max(std::abs(inertialForce[c]),std::abs(force[i][c])));
          tCurrent=std::max(tCurrent,std::max(std::abs(inertialTorque[c]),std::abs(torque[i][c])));
        }
        e.forceReference[i]=fCurrent;e.torqueReference[i]=tCurrent;
        e.diagnostic.maxForceResidualRatio=std::max(e.diagnostic.maxForceResidualRatio,
            fError/(settings.forceAbsoluteTolerance+settings.relativeTolerance*fCurrent));
        e.diagnostic.maxTorqueResidualRatio=std::max(e.diagnostic.maxTorqueResidualRatio,
            tError/(settings.torqueAbsoluteTolerance+settings.relativeTolerance*tCurrent));
      }
      for(double r:e.residual)if(!std::isfinite(r))throw std::domain_error("Non-finite particle residual");
      return true;
    } catch(const std::exception& ex) {error=ex.what();return false;}
  }
};

// Unrestarted GMRES for a matrix-free Newton correction; the dimension is at
// 6N plus the active normal reactions; no dense global matrix is assembled.
inline bool gmres(const std::function<bool(const Vector&,Vector&)>& apply,
                  const Vector& rhs,int maxIterations,double tolerance,
                  Vector& solution,int& iterations,double* trueResidual=nullptr,double* estimatedResidual=nullptr) {
  const std::size_t n=rhs.size();const int m=std::min<int>(maxIterations,n);
  const double beta=length(rhs);solution.assign(n,0.);
  if(beta<=tolerance)return true;
  std::vector<Vector> v(m+1,Vector(n,0.));
  std::vector<Vector> h(m+1,Vector(m,0.));
  Vector cs(m),sn(m),g(m+1,0.);g[0]=beta;
  for(std::size_t k=0;k<n;++k)v[0][k]=rhs[k]/beta;
  int used=0;
  for(int j=0;j<m;++j) {
    Vector w;
    if(!apply(v[j],w))return false;
    ++iterations;
    for(int k=0;k<=j;++k) {
      h[k][j]=inner(w,v[k]);for(std::size_t c=0;c<n;++c)w[c]-=h[k][j]*v[k][c];
    }
    // A second orthogonalization pass helps close-contact ill conditioning.
    for(int k=0;k<=j;++k) {
      const double correction=inner(w,v[k]);h[k][j]+=correction;
      for(std::size_t c=0;c<n;++c)w[c]-=correction*v[k][c];
    }
    h[j+1][j]=length(w);
    if(h[j+1][j]>1.e-30)for(std::size_t c=0;c<n;++c)v[j+1][c]=w[c]/h[j+1][j];
    for(int k=0;k<j;++k) {
      const double t=cs[k]*h[k][j]+sn[k]*h[k+1][j];
      h[k+1][j]=-sn[k]*h[k][j]+cs[k]*h[k+1][j];h[k][j]=t;
    }
    const double d=std::hypot(h[j][j],h[j+1][j]);
    if(!(d>0.))break;
    cs[j]=h[j][j]/d;sn[j]=h[j+1][j]/d;
    h[j][j]=d;h[j+1][j]=0.;
    g[j+1]=-sn[j]*g[j];g[j]*=cs[j];used=j+1;
    if(std::abs(g[j+1])<=.1*tolerance)break;
  }
  if(used==0)return false;
  Vector y(used,0.);
  for(int j=used-1;j>=0;--j) {
    double s=g[j];for(int k=j+1;k<used;++k)s-=h[j][k]*y[k];
    if(std::abs(h[j][j])<1.e-30)return false;
    y[j]=s/h[j][j];
  }
  for(int j=0;j<used;++j)for(std::size_t c=0;c<n;++c)solution[c]+=y[j]*v[j][c];
  if(!std::all_of(solution.begin(),solution.end(),[](double a){return std::isfinite(a);}))return false;
  Vector check;if(!apply(solution,check))return false;
  for(std::size_t k=0;k<n;++k)check[k]-=rhs[k];
  if(trueResidual)*trueResidual=length(check);
  if(estimatedResidual)*estimatedResidual=std::abs(g[used]);
  return length(check)<=tolerance;
}

using Block6=std::array<double,36>;
inline bool inverse6(Block6 a,Block6& inverse) {
  inverse.fill(0.);for(int k=0;k<6;++k)inverse[6*k+k]=1.;
  for(int k=0;k<6;++k) {
    int pivot=k;for(int r=k+1;r<6;++r)if(std::abs(a[6*r+k])>std::abs(a[6*pivot+k]))pivot=r;
    if(!(std::abs(a[6*pivot+k])>1.e-25))return false;
    if(pivot!=k)for(int c=0;c<6;++c){std::swap(a[6*k+c],a[6*pivot+c]);std::swap(inverse[6*k+c],inverse[6*pivot+c]);}
    const double scale=1./a[6*k+k];for(int c=0;c<6;++c){a[6*k+c]*=scale;inverse[6*k+c]*=scale;}
    for(int r=0;r<6;++r)if(r!=k){const double m=a[6*r+k];for(int c=0;c<6;++c){a[6*r+c]-=m*a[6*k+c];inverse[6*r+c]-=m*inverse[6*k+c];}}
  }
  return true;
}
using Six=std::array<double,6>;
inline Six multiply6(const Block6& a,const Six& x) {
  Six y{};for(int r=0;r<6;++r)for(int c=0;c<6;++c)y[r]+=a[6*r+c]*x[c];return y;
}
inline double dot6(const Six& a,const Six& b){double s=0.;for(int k=0;k<6;++k)s+=a[k]*b[k];return s;}
struct ContactBlock {std::size_t i,j;int slot;Six bi{},bj{},ci{},cj{};double schur=0.;};

// Invert local 6x6 inertia + near-field/contact resistance blocks, then a
// diagonal approximation to the normal-constraint Schur complement.  No dense
// 6N-by-6N finite-difference matrix is assembled.
struct BlockPreconditioner {
  std::vector<Block6> inverse;
  std::vector<ContactBlock> contact;
  std::size_t n=0;
  bool build(const Evaluation& e,const Residual& residual) {
    n=e.bodies.size();std::vector<Block6> blocks(n);inverse.resize(n);contact.clear();
    const double dt=residual.dt,L=residual.lengthScale;
    for(std::size_t i=0;i<n;++i) {
      for(int c=0;c<3;++c)blocks[i][6*c+c]=1.;
      const double im=*std::max_element(e.bodies[i].inertiaBody.begin(),e.bodies[i].inertiaBody.end());
      const auto inertia=rotatedDiagonal(e.bodies[i].rotation,e.bodies[i].inertiaBody);
      for(int r=0;r<3;++r)for(int c=0;c<3;++c)blocks[i][6*(r+3)+c+3]=inertia[3*r+c]/im;
    }
    for(const auto& p:e.pairs) {
      for(int side=0;side<2;++side) {
        const std::size_t i=side?p.j:p.i;const Vec3 lever=side?p.leverJ:p.leverI;
        const double fs=e.bodies[i].mass*L/(dt*dt);
        const double ts=*std::max_element(e.bodies[i].inertiaBody.begin(),e.bodies[i].inertiaBody.end())/(dt*dt);
        for(int col=0;col<6;++col) {
          Vec3 dv{},dw{};if(col<3)dv[col]=L/dt;else dw[col-3]=1./dt;
          const Vec3 force=mul(p.resistance,add(dv,cross(dw,lever)));
          const Vec3 torque=add(cross(lever,force),mul(p.rollingResistance,dw));
          for(int c=0;c<3;++c){blocks[i][6*c+col]+=force[c]/fs;blocks[i][6*(c+3)+col]+=torque[c]/ts;}
        }
      }
    }
    for(std::size_t i=0;i<n;++i)if(!inverse6(blocks[i],inverse[i]))return false;
    for(const auto& p:e.pairs)if(p.slot>=0) {
      ContactBlock c;c.i=p.i;c.j=p.j;c.slot=p.slot;
      const double fi=e.bodies[p.i].mass*L/(dt*dt),fj=e.bodies[p.j].mass*L/(dt*dt);
      const double ti=*std::max_element(e.bodies[p.i].inertiaBody.begin(),e.bodies[p.i].inertiaBody.end())/(dt*dt);
      const double tj=*std::max_element(e.bodies[p.j].inertiaBody.begin(),e.bodies[p.j].inertiaBody.end())/(dt*dt);
      const Vec3 gi=scale(cross(p.leverI,p.normal),-1./L),gj=scale(cross(p.leverJ,p.normal),1./L);
      for(int k=0;k<3;++k) {
        c.bi[k]=-p.forcePerNormalLoad[k]*residual.reactionScale/fi;
        c.bj[k]=p.forcePerNormalLoad[k]*residual.reactionScale/fj;
        c.bi[k+3]=-p.torqueIPerNormalLoad[k]*residual.reactionScale/ti;
        c.bj[k+3]=-p.torqueJPerNormalLoad[k]*residual.reactionScale/tj;
        c.ci[k]=-p.normal[k];c.cj[k]=p.normal[k];c.ci[k+3]=gi[k];c.cj[k+3]=gj[k];
      }
      c.schur=dot6(c.ci,multiply6(inverse[c.i],c.bi))+dot6(c.cj,multiply6(inverse[c.j],c.bj));
      if(!std::isfinite(c.schur) || std::abs(c.schur)<1.e-25)return false;
      contact.push_back(c);
    }
    return true;
  }
  Vector operator()(const Vector& r) const {
    Vector out(r.size(),0.);
    for(std::size_t i=0;i<n;++i){Six x{};for(int k=0;k<6;++k)x[k]=r[6*i+k];const auto y=multiply6(inverse[i],x);for(int k=0;k<6;++k)out[6*i+k]=y[k];}
    Vector corrected=r;
    for(const auto& c:contact) {
      Six zi{},zj{};for(int k=0;k<6;++k){zi[k]=out[6*c.i+k];zj[k]=out[6*c.j+k];}
      const double lambda=(dot6(c.ci,zi)+dot6(c.cj,zj)-r[6*n+c.slot])/c.schur;
      out[6*n+c.slot]=lambda;
      for(int k=0;k<6;++k){corrected[6*c.i+k]-=c.bi[k]*lambda;corrected[6*c.j+k]-=c.bj[k]*lambda;}
    }
    for(std::size_t i=0;i<n;++i){Six x{};for(int k=0;k<6;++k)x[k]=corrected[6*i+k];const auto y=multiply6(inverse[i],x);for(int k=0;k<6;++k)out[6*i+k]=y[k];}
    return out;
  }
};

inline bool physicallyConverged(const Evaluation& e,const ParticleStepSettings& s) {
  if(e.diagnostic.maxForceResidualRatio>1. || e.diagnostic.maxTorqueResidualRatio>1.)return false;
  if(s.rough.enabled) {
    if(e.diagnostic.contactGapViolation>s.contactGapTolerance)return false;
    for(std::size_t p=0;p<e.contacts.size();++p)if(e.contacts[p].active) {
      if(e.contacts[p].normalLoad<0. || std::abs(e.gaps[p]-s.rough.gap)>s.contactGapTolerance)return false;
      if(norm(e.contacts[p].tangentForce)>s.rough.friction*e.contacts[p].normalLoad+s.forceAbsoluteTolerance)return false;
      if(norm(e.contacts[p].rollingTorque)>e.contacts[p].rollingCap+s.torqueAbsoluteTolerance)return false;
    }
  }
  return true;
}

inline bool implicitStep(const std::vector<Body>& old,const std::vector<Vec3>& force,
                         const std::vector<Vec3>& torque,double dt,double time,
                         const ParticleStepSettings& settings,std::vector<GapCache>& cache,
                         PersistentContactState& contacts,std::vector<Body>& output,
                         ParticleStepDiagnostics& diagnostic,std::string& error) {
  double lengthScale=0.;for(const auto& b:old)lengthScale=std::max(lengthScale,radius(b));
  const std::size_t n=old.size(),np=n*(n-1)/2;
  const double reactionScale=old.front().mass*lengthScale/(dt*dt);
  std::vector<int> activeSlot(np,-1);std::vector<double> birthReaction(np,0.);int activeCount=0;
  Vector q(6*n,0.);
  for(std::size_t i=0;i<n;++i)for(int c=0;c<3;++c) {
    q[6*i+c]=old[i].velocity[c]*dt/lengthScale;
    q[6*i+3+c]=old[i].omega[c]*dt;
  }
  if(settings.rough.enabled) {
    std::size_t index=0;
    for(std::size_t i=0;i<n;++i)for(std::size_t j=i+1;j<n;++j,++index) {
      const Body image=closestImage(old[i],old[j],time,settings);
      const Vec3 rij=sub(old[i].position,image.position);
      const double bound=radius(old[i])+radius(image)+settings.pair.cutoffGap;
      if(dot(rij,rij)>bound*bound && !contacts[index].active)continue;
      const auto gap=closestEllipsoidGap(old[i],image,&cache[index]);
      if(gap.gap<settings.rough.gap-settings.contactGapTolerance) {
        error="Initial accepted particle state violates the rough contact gap";return false;
      }
      const Vec3 vi=add(old[i].velocity,cross(old[i].omega,gap.leverI));
      const Vec3 vj=add(image.velocity,cross(image.omega,gap.leverJ));
      const double predicted=gap.gap+dt*dot(sub(vj,vi),gap.normal);
      if((contacts[index].active && gap.gap<=settings.rough.gap+settings.contactGapTolerance)
          || predicted<=settings.rough.gap) {
        activeSlot[index]=activeCount++;
        // An initial Newton guess only: the reaction remains a solved unknown.
        // Starting a newly adhesive contact at N=0 creates a Coulomb corner
        // and can stall line search while the gap moves to its constraint.
        Body birthImage=image;birthImage.position=add(image.position,scale(gap.normal,settings.rough.gap-gap.gap));
        const auto birth=evaluatePair(old[i],birthImage,settings.pair);
        birthReaction[index]=std::max(0.,dot(birth.forceI,birth.normal));
      }
    }
  }
  q.resize(6*n+activeCount,0.);
  for(std::size_t p=0;p<np;++p)if(activeSlot[p]>=0)q[6*n+activeSlot[p]]=(contacts[p].active?std::max(0.,contacts[p].normalLoad):birthReaction[p])/reactionScale;
  Residual residual{old,force,torque,settings,cache,contacts,activeSlot,dt,time,lengthScale,reactionScale,activeCount};
  Evaluation base;bool predictorValid=false;
  for(int attempt=0;attempt<12;++attempt) {
    if(residual(q,base,error)){predictorValid=true;break;}
    for(std::size_t k=0;k<6*n;++k)q[k]*=.5;
  }
  if(!predictorValid){for(std::size_t k=0;k<6*n;++k)q[k]=0.;if(!residual(q,base,error))return false;}
  int newton=0,krylov=0;
  for(;newton<settings.maxNewtonIterations;++newton) {
    // A fixed active set is used throughout each Krylov solve.  Reactions and
    // gap violations update it only at a completed outer Newton iteration.
    if(settings.rough.enabled) {
      std::vector<int> revised(np,-1);int count=0;bool changed=false,clamped=false;
      for(std::size_t p=0;p<np;++p) {
        bool on=activeSlot[p]>=0;
        if(on && q[6*n+activeSlot[p]]*reactionScale < -settings.forceAbsoluteTolerance)on=false;
        else if(!on && base.gaps[p]<settings.rough.gap-1.e-15)on=true;
        if(on)revised[p]=count++;
        changed|=(on!=(activeSlot[p]>=0));
      }
      if(changed) {
        Vector revisedQ(6*n+count,0.);std::copy(q.begin(),q.begin()+6*n,revisedQ.begin());
        for(std::size_t p=0;p<np;++p)if(revised[p]>=0) {
          if(activeSlot[p]>=0)revisedQ[6*n+revised[p]]=std::max(0.,q[6*n+activeSlot[p]]);
          else {
            const auto found=std::find_if(base.pairs.begin(),base.pairs.end(),
                [&](const PairLinearization& a){return a.index==p;});
            if(found!=base.pairs.end()) {
              Body image=closestImage(base.bodies[found->i],base.bodies[found->j],time+dt,settings);
              image.position=add(image.position,scale(found->normal,settings.rough.gap-found->gap));
              const auto birth=evaluatePair(base.bodies[found->i],image,settings.pair);
              revisedQ[6*n+revised[p]]=std::max(0.,dot(birth.forceI,birth.normal))/reactionScale;
            }
          }
        }
        activeSlot.swap(revised);activeCount=count;residual.activeCount=count;q.swap(revisedQ);
        if(!residual(q,base,error))return false;
      }
      for(std::size_t p=0;p<np;++p)if(activeSlot[p]>=0 && q[6*n+activeSlot[p]]<0.) {
        q[6*n+activeSlot[p]]=0.;clamped=true;
      }
      if(clamped && !residual(q,base,error))return false;
    }
    if(physicallyConverged(base,settings)) {
      output=std::move(base.bodies);cache=std::move(base.cache);contacts=std::move(base.contacts);diagnostic=base.diagnostic;
      diagnostic.newtonIterations=newton;diagnostic.krylovIterations=krylov;
      diagnostic.residualEvaluations=residual.evaluations;
      for(auto& b:output)wrap(b,time+dt,settings);
      return true;
    }
    BlockPreconditioner preconditioner;
    if(!preconditioner.build(base,residual)) {error="Particle block/contact preconditioner is singular";return false;}
    auto jacobian=[&](const Vector& direction,Vector& product)->bool {
      double magnitude=0.;for(std::size_t k=0;k<6*n;++k)magnitude+=direction[k]*direction[k];magnitude=std::sqrt(magnitude);
      product.assign(q.size(),0.);
      if(magnitude>0.) {
        double bodyMagnitude=0.;for(std::size_t k=0;k<6*n;++k)bodyMagnitude+=q[k]*q[k];
        double epsilon=settings.finiteDifferenceStep*(1.+std::sqrt(bodyMagnitude))/magnitude;
        Vector trial=q;Evaluation displaced;std::string derivativeError;bool valid=false;
        for(int attempt=0;attempt<8;++attempt) {
          for(std::size_t k=0;k<6*n;++k)trial[k]=q[k]+epsilon*direction[k];
          if(residual(trial,displaced,derivativeError)){valid=true;break;}
          epsilon*=-.5;
        }
        if(!valid){error=derivativeError;return false;}
        for(std::size_t k=0;k<q.size();++k)product[k]=(displaced.residual[k]-base.residual[k])/epsilon;
        // Linearize the ENTIRE return map on one semismooth branch. The
        // trial-vector difference retains changing contact geometry and common
        // spin/history transport. Replacing just a frozen-geometry part leaves
        // a second, inconsistent material branch inside the geometric Jv.
        auto returnDerivative=[](const Vec3& trial,const Vec3& change,double stiffness,
                                 double cap,double dStiffness,double dCap)->Vec3 {
          const double r=norm(trial);
          if(cap<=0.)return r>0.?scale(trial,-dCap/r):Vec3{};
          if(stiffness*r<=cap)return add(scale(change,-stiffness),scale(trial,-dStiffness));
          const Vec3 direction=scale(trial,1./r);
          return add(scale(sub(change,scale(direction,dot(direction,change))),-cap/r),
                     scale(direction,-dCap));
        };
        for(const auto& p:base.pairs)if(p.slot>=0) {
          const auto found=std::find_if(displaced.pairs.begin(),displaced.pairs.end(),
              [&](const PairLinearization& a){return a.index==p.index;});
          if(found==displaced.pairs.end()){error="Active contact missing from geometric derivative";return false;}
          const auto& contact=p.contact;const auto& next=found->contact;
          const Vec3 dSlip=scale(sub(next.trialSlip,contact.trialSlip),1./epsilon);
          const Vec3 dRoll=scale(sub(next.trialRoll,contact.trialRoll),1./epsilon);
          const Vec3 dNormal=scale(sub(found->normal,p.normal),1./epsilon);
          const Vec3 dLeverI=scale(sub(found->leverI,p.leverI),1./epsilon);
          const Vec3 dLeverJ=scale(sub(found->leverJ,p.leverJ),1./epsilon);
          const auto& state=contact.candidateState;const auto& stateNext=next.candidateState;
          const Vec3 dFt=returnDerivative(contact.trialSlip,dSlip,settings.rough.tangentialStiffness,
              settings.rough.friction*std::max(0.,state.normalLoad),0.,0.);
          const Vec3 dMr=returnDerivative(contact.trialRoll,dRoll,state.rollingStiffness,state.rollingCap,
              (stateNext.rollingStiffness-state.rollingStiffness)/epsilon,
              (stateNext.rollingCap-state.rollingCap)/epsilon);
          const Vec3 dF=add(scale(dNormal,-state.normalLoad),dFt);
          const Vec3 dTi=add(add(cross(dLeverI,contact.forceI),cross(p.leverI,dF)),dMr);
          const Vec3 dTj=scale(add(add(cross(dLeverJ,contact.forceI),cross(p.leverJ,dF)),dMr),-1.);
          const Vec3 correctF=sub(dF,scale(sub(next.forceI,contact.forceI),1./epsilon));
          const Vec3 correctTi=sub(dTi,scale(sub(next.torqueI,contact.torqueI),1./epsilon));
          const Vec3 correctTj=sub(dTj,scale(sub(next.torqueJ,contact.torqueJ),1./epsilon));
          const double fi=old[p.i].mass*lengthScale/(dt*dt),fj=old[p.j].mass*lengthScale/(dt*dt);
          const double ti=*std::max_element(old[p.i].inertiaBody.begin(),old[p.i].inertiaBody.end())/(dt*dt);
          const double tj=*std::max_element(old[p.j].inertiaBody.begin(),old[p.j].inertiaBody.end())/(dt*dt);
          for(int c=0;c<3;++c) {
            product[6*p.i+c]-=correctF[c]/fi;product[6*p.j+c]+=correctF[c]/fj;
            product[6*p.i+c+3]-=correctTi[c]/ti;product[6*p.j+c+3]-=correctTj[c]/tj;
          }
        }
      }
      // Exact normal-load and friction-cap columns.  In particular, the hard
      // constraint is never approximated by a normal penalty stiffness.
      for(const auto& c:preconditioner.contact) {
        const double d=direction[6*n+c.slot];
        for(int k=0;k<6;++k){product[6*c.i+k]+=c.bi[k]*d;product[6*c.j+k]+=c.bj[k]*d;}
      }
      return true;
    };
    auto apply=[&](const Vector& direction,Vector& product)->bool {
      Vector physical;if(!jacobian(direction,physical))return false;product=preconditioner(physical);return true;
    };
    Vector rhs=base.residual;for(double& x:rhs)x=-x;rhs=preconditioner(rhs);Vector correction;
    const double linearTolerance=std::max(1.e-16,.05*length(rhs));
    double linearResidual=std::numeric_limits<double>::infinity(),estimatedResidual=linearResidual;
    if(!gmres(apply,rhs,settings.maxKrylovIterations,linearTolerance,correction,krylov,&linearResidual,&estimatedResidual)) {
      std::ostringstream message;message<<"Particle Newton GMRES failed its true linear residual criterion: actual="
        <<linearResidual<<", Arnoldi="<<estimatedResidual<<", target="<<linearTolerance<<", rhs="<<length(rhs);
      error=message.str();return false;
    }
    const double before=length(base.residual);bool accepted=false;
    Vector trial=q;Evaluation next;double fraction=1.;
    for(int line=0;line<settings.maxLineSearch;++line,fraction*=.5) {
      for(std::size_t k=0;k<q.size();++k)trial[k]=q[k]+fraction*correction[k];
      if(residual(trial,next,error) && (length(next.residual)<before*(1.-1.e-4*fraction)
          || physicallyConverged(next,settings))) {
        q=trial;base=std::move(next);accepted=true;break;
      }
    }
    if(!accepted){error="Particle Newton line search failed to reduce the admissible-gap residual";return false;}
  }
  if(physicallyConverged(base,settings)) {
    output=std::move(base.bodies);cache=std::move(base.cache);contacts=std::move(base.contacts);diagnostic=base.diagnostic;
    diagnostic.newtonIterations=newton;diagnostic.krylovIterations=krylov;
    diagnostic.residualEvaluations=residual.evaluations;
    for(auto& b:output)wrap(b,time+dt,settings);
    return true;
  }
  std::ostringstream message;message<<"Particle implicit solve exceeded "<<settings.maxNewtonIterations
      <<" Newton iterations; force residual ratio="<<base.diagnostic.maxForceResidualRatio
      <<", torque residual ratio="<<base.diagnostic.maxTorqueResidualRatio
      <<", gap violation="<<base.diagnostic.contactGapViolation<<" m";
  error=message.str();return false;
}

#include "particlePetscSolver.h"

} // namespace particle_detail

inline ParticleStepDiagnostics evaluateParticleState(
    const std::vector<Body>& bodies,double time,const ParticleStepSettings& settings,
    std::vector<GapCache>* persistentCache=nullptr,
    const PersistentContactState* persistentContacts=nullptr) {
  ParticleStepDiagnostics diagnostic;
  const std::size_t pairs=bodies.size()*(bodies.size()-1)/2;
  std::vector<GapCache> localCache;
  if(!persistentCache){localCache.resize(pairs);persistentCache=&localCache;}
  if(persistentCache->size()!=pairs)persistentCache->resize(pairs);
  if(persistentContacts && !persistentContacts->empty() && persistentContacts->size()!=pairs)
    throw std::invalid_argument("Contact history size does not match ordered particle pairs");
  std::size_t index=0;
  for(std::size_t i=0;i<bodies.size();++i)for(std::size_t j=i+1;j<bodies.size();++j,++index) {
    if(settings.rough.enabled && persistentContacts && !persistentContacts->empty()) {
      const auto& contact=(*persistentContacts)[index];
      if(contact.stepDuration>0.)diagnostic.contactDissipation+=(contact.plasticSlipWork+contact.plasticRollWork+contact.releasedEnergy)/contact.stepDuration;
    }
    const Body image=particle_detail::closestImage(bodies[i],bodies[j],time,settings);
    const Vec3 rij=sub(bodies[i].position,image.position);
    const double range=std::max(settings.pair.cutoffGap,settings.nearField.matchingGap);
    const double bound=particle_detail::radius(bodies[i])+particle_detail::radius(image)+range;
    if(dot(rij,rij)>bound*bound)continue;
    auto pair=evaluatePair(bodies[i],image,settings.pair,&(*persistentCache)[index]);
    if(!std::isfinite(pair.gap)) {
      const auto gap=closestEllipsoidGap(bodies[i],image,&(*persistentCache)[index]);
      pair.gap=gap.gap;pair.normal=gap.normal;pair.leverI=gap.leverI;pair.leverJ=gap.leverJ;
    }
    if(!(pair.gap>0.))throw std::domain_error("Initial particle state has overlapping ellipsoids");
    diagnostic.minGap=std::min(diagnostic.minGap,pair.gap);
    diagnostic.maxForce=std::max(diagnostic.maxForce,norm(pair.forceI));
    diagnostic.energyAtEnd+=pair.energy;diagnostic.activePairs+=pair.active?1:0;
    particle_detail::addMoment(diagnostic.pairMoment,pair.forceI,rij);
    particle_detail::addMoment(diagnostic.attractiveMoment,pair.forceAttractiveI,rij);
    particle_detail::addMoment(diagnostic.repulsiveMoment,pair.forceRepulsiveI,rij);
    const auto lub=nearFieldResistance(bodies[i],image,pair.gap,pair.normal,pair.leverI,pair.leverJ,settings.nearField);
    particle_detail::addMoment(diagnostic.lubricationMoment,lub.forceI,rij);
    diagnostic.lubricationDissipation+=lub.dissipation;
    if(settings.rough.enabled && persistentContacts && !persistentContacts->empty()) {
      const auto& contact=(*persistentContacts)[index];
      diagnostic.contactGapViolation=std::max(diagnostic.contactGapViolation,std::max(0.,settings.rough.gap-pair.gap));
      if(contact.active) {
        particle_detail::addMoment(diagnostic.contactNormalMoment,scale(contact.normal,-contact.normalLoad),rij);
        particle_detail::addMoment(diagnostic.contactTangentialMoment,contact.tangentForce,rij);
        diagnostic.elasticContactEnergy+=contact.elasticEnergy;
        ++diagnostic.contacts;diagnostic.slidingContacts+=contact.sliding?1:0;diagnostic.rollingContacts+=contact.rolling?1:0;
        diagnostic.maxForce=std::max(diagnostic.maxForce,norm(add(scale(contact.normal,-contact.normalLoad),contact.tangentForce)));
      }
    }
  }
  return diagnostic;
}

// Advance one fixed outer LB interval.  Resolved hydrodynamic loads are frozen
// over this interval; only particle work is subdivided.  There is no global LB
// timestep alteration, checkpoint, restart, force clipping, or attraction ramp.
inline ParticleStepDiagnostics advanceParticles(
    std::vector<Body>& bodies,const std::vector<Vec3>& hydroForce,
    const std::vector<Vec3>& hydroTorque,double dt,double time,
    const ParticleStepSettings& settings,std::vector<GapCache>* persistentCache=nullptr,
    PersistentContactState* persistentContacts=nullptr) {
  if(bodies.empty())return {};
  if(hydroForce.size()!=bodies.size() || hydroTorque.size()!=bodies.size())
    throw std::invalid_argument("Particle hydrodynamic load count does not match particle count");
  if(!(dt>0.) || !std::isfinite(dt) || settings.maxSubsteps<1 || settings.maxNewtonIterations<1
      || settings.maxKrylovIterations<1 || !(settings.relativeTolerance>0.)
      || !(settings.forceAbsoluteTolerance>0.) || !(settings.torqueAbsoluteTolerance>0.)
      || !(settings.contactGapTolerance>0.))
    throw std::invalid_argument("Invalid implicit particle step settings");
  for(double l:settings.box)if(!(l>0.))throw std::invalid_argument("Invalid periodic box length");
  for(const Body& b:bodies) {
    if(!(b.mass>0.) || !std::isfinite(b.mass))throw std::invalid_argument("Particle mass must be positive");
    for(double a:b.axes)if(!(a>0.))throw std::invalid_argument("Particle semiaxes must be positive");
    for(double inertia:b.inertiaBody)if(!(inertia>0.))throw std::invalid_argument("Particle principal inertia must be positive");
  }
  if(settings.rough.enabled && !persistentContacts)
    throw std::invalid_argument("Rough contacts require persistent history passed to advanceParticles");
  const std::size_t pairs=bodies.size()*(bodies.size()-1)/2;
  if(persistentContacts && !persistentContacts->empty() && persistentContacts->size()!=pairs)
    throw std::invalid_argument("Contact history size does not match ordered particle pairs");
  PersistentContactState originalContacts=(persistentContacts && persistentContacts->size()==pairs)
      ?*persistentContacts:PersistentContactState(pairs);
  std::vector<GapCache> originalCache=(persistentCache && persistentCache->size()==pairs)
      ?*persistentCache:std::vector<GapCache>(pairs);
  std::string failure;
  particle_detail::ParticleDiagnosticScope diagnosticScope(settings);
  for(int count=1;count<=settings.maxSubsteps;count*=2) {
    auto working=bodies;auto cache=originalCache;auto contacts=originalContacts;ParticleStepDiagnostics total;
    total.substeps=count;const double subdt=dt/count;bool okay=true;
    for(int substep=0;substep<count;++substep) {
      std::vector<Body> next;ParticleStepDiagnostics d;
      if(!particle_detail::dispatchImplicitStep(working,hydroForce,hydroTorque,subdt,time+substep*subdt,
                                       settings,cache,contacts,next,d,failure,time,dt,count,substep)) {okay=false;break;}
      working=std::move(next);
      total.minGap=std::min(total.minGap,d.minGap);total.maxForce=std::max(total.maxForce,d.maxForce);
      total.energyAtEnd=d.energyAtEnd;total.activePairs=d.activePairs;
      total.elasticContactEnergy=d.elasticContactEnergy;total.contacts=d.contacts;
      total.slidingContacts=d.slidingContacts;total.rollingContacts=d.rollingContacts;
      total.maxForceResidualRatio=std::max(total.maxForceResidualRatio,d.maxForceResidualRatio);
      total.maxTorqueResidualRatio=std::max(total.maxTorqueResidualRatio,d.maxTorqueResidualRatio);
      total.contactGapViolation=std::max(total.contactGapViolation,d.contactGapViolation);
      total.newtonIterations+=d.newtonIterations;total.krylovIterations+=d.krylovIterations;
      total.frictionBranchAttempts+=d.frictionBranchAttempts;
      total.frictionBranchCorrections+=d.frictionBranchCorrections;
      total.contactStateUpdates+=d.contactStateUpdates;
      total.contactActivations+=d.contactActivations;
      total.contactReleases+=d.contactReleases;
      total.residualEvaluations+=d.residualEvaluations;
      const double weight=1./count;
      total.lubricationDissipation+=weight*d.lubricationDissipation;
      total.contactDissipation+=weight*d.contactDissipation;
      particle_detail::addMatrix(total.pairMoment,d.pairMoment,weight);
      particle_detail::addMatrix(total.attractiveMoment,d.attractiveMoment,weight);
      particle_detail::addMatrix(total.repulsiveMoment,d.repulsiveMoment,weight);
      particle_detail::addMatrix(total.lubricationMoment,d.lubricationMoment,weight);
      particle_detail::addMatrix(total.contactNormalMoment,d.contactNormalMoment,weight);
      particle_detail::addMatrix(total.contactTangentialMoment,d.contactTangentialMoment,weight);
    }
    if(okay) {
      // Recoverable retries are discarded with the bounded in-memory trace.
      // Routine successful timesteps produce no diagnostic filesystem traffic.
      bodies=std::move(working);if(persistentCache)*persistentCache=std::move(cache);
      if(persistentContacts)*persistentContacts=std::move(contacts);
      return total;
    }
    if(count>settings.maxSubsteps/2)break;
  }
  diagnosticScope.flushFailure();
  particle_detail::writeOuterOutcome(settings,time,dt,settings.maxSubsteps,false);
  std::ostringstream message;message<<"Particle solve failed within the configured "<<settings.maxSubsteps
      <<" particle-only subdivisions at t="<<time<<" s and fixed dt_LB="<<dt<<" s: "<<failure
      <<". The global LB timestep was not changed.";
  throw std::runtime_error(message.str());
}

} // namespace graphite

} } // SLURRY SCOPE END
#endif
