#ifndef SLURRY_GR_RE2_GRAPHITE_RE_SQUARED_POTENTIAL_H
#define SLURRY_GR_RE2_GRAPHITE_RE_SQUARED_POTENTIAL_H
#include "ellipsoidGap.h"
#include "contactCurvature.h"
#include <limits>

// SLURRY SCOPE BEGIN
namespace slurry { namespace gr_re2 {

namespace graphite {
struct PairParameters {
  double hamaker=.99e-19,sigma=3e-9,switchGap=400e-9,cutoffGap=500e-9;
  // The bare rough-contact reference remains at roughnessGap. Local adhesion replaces a
  // fraction of the near-contact RE2 energy with the same interaction at
  // d=localGap+(h-roughnessGap). The fraction is an effective contribution,
  // not a measured real-contact area fraction. Zero preserves legacy RE2.
  double roughnessGap=2e-9,localGap=3e-10,localGapFraction=0.;
  double localSwitchExcessGap=2e-9,localCutoffExcessGap=10e-9;
  // Version 1: local-curvature Derjaguin background, matched to far RE2,
  // plus a finite-range Dugdale surface traction. Legacy configurations and
  // replay fixtures keep the old evaluator unless explicitly enabled.
  bool surfaceAdhesion=false;
  double adhesionWork=.0219,adhesionRange=.67e-9;
  double curvatureSwitchGap=5e-9,curvatureCutoffGap=20e-9;
  // CMC contact version 1: actual carbon gap still determines the RE2
  // background, while finite-range adhesion opens from the coating contact
  // plane. Zero contactGap selects the bare reference for exact compatibility.
  // Retention screens both attractive energies, including their geometry
  // derivatives. It is an effective cohesion closure, not an EDL force law.
  double contactGap=0.,cohesionRetention=1.;
  // Legacy exponential free-CMC pressure (Pa) and decay length (m). Active
  // legacy pressure cannot be combined with coated-contact cohesion screening.
  double freeCmcRepulsionPressure=0.,freeCmcRepulsionLength=5e-9;
  // Compact free-CMC overlap repulsion about the unchanged rough-contact gap.
  // U=pi*lambda_D*W*delta/p * [1-(h-h0)/delta]_+^p. W is an effective
  // planar work (J/m^2), supplied separately from adsorbed-CMC adhesion.
  // Zero work preserves all previous pair models, including replay fixtures.
  double freeCmcInnerRepulsionWork=0.,freeCmcInnerRepulsionRange=.67e-9;
  double freeCmcInnerRepulsionPower=2.1;
  // Effective net CMC interaction: replace, rather than add to, the complete
  // legacy pair energy. Forces are reference aligned-face amplitudes (N);
  // local curvature/reference length supplies the geometry dependence.
  // The concentration-dependent blend is supplied by the caller and is held
  // constant while differentiating a pair. Zero preserves legacy arithmetic.
  double cmcNetBlend=0.;
  double cmcNetContactForce=1.5e-10,cmcNetBarrierForce=2.e-11;
  double cmcNetAttractionRange=1.e-9,cmcNetRepulsionRange=6.e-9;
  double cmcNetReferenceLength=13.6125e-6;
};
inline double effectiveContactGap(const PairParameters&p){
  return p.contactGap>0.?p.contactGap:p.roughnessGap;
}
// Planar work already present in the geometric-gap LJ background. This is a
// surface energy (J/m^2), not the complete curved-pair separation work (J).
inline double surfaceBackgroundWork(const PairParameters&p){
  constexpr double pi=3.1415926535897932384626433832795;
  return p.hamaker/(12.*pi*p.roughnessGap*p.roughnessGap)
      *(1.-std::pow(p.sigma/p.roughnessGap,6)/30.);
}
inline void validatePairParameters(const PairParameters&p){
  if(!std::isfinite(p.hamaker)||!(p.hamaker>=0.)||!std::isfinite(p.sigma)||!(p.sigma>0.)
     ||!std::isfinite(p.switchGap)||!(p.switchGap>=0.)||!std::isfinite(p.cutoffGap)||!(p.cutoffGap>p.switchGap))
    throw std::domain_error("Invalid RE2 parameters");
  if(!std::isfinite(p.localGapFraction)||p.localGapFraction<0.||p.localGapFraction>1.)
    throw std::domain_error("RE2 local gap fraction must be finite and between zero and one");
  if(!std::isfinite(p.contactGap)||p.contactGap<0.
     ||(p.contactGap>0.&&(!std::isfinite(p.roughnessGap)||!(p.roughnessGap>0.)
                         ||p.contactGap<p.roughnessGap)))
    throw std::domain_error("CMC contact gap must be zero (bare reference) or finite and at least h0>0");
  if(!std::isfinite(p.cohesionRetention)||p.cohesionRetention<0.||p.cohesionRetention>1.)
    throw std::domain_error("CMC cohesion retention must be finite and between zero and one");
  if(!std::isfinite(p.cmcNetBlend)||p.cmcNetBlend<0.||p.cmcNetBlend>1.
     ||!std::isfinite(p.cmcNetContactForce)||!(p.cmcNetContactForce>0.)
     ||!std::isfinite(p.cmcNetBarrierForce)||!(p.cmcNetBarrierForce>0.)
     ||!std::isfinite(p.cmcNetAttractionRange)||!(p.cmcNetAttractionRange>0.)
     ||!std::isfinite(p.cmcNetRepulsionRange)||!(p.cmcNetRepulsionRange>0.)
     ||!std::isfinite(p.cmcNetReferenceLength)||!(p.cmcNetReferenceLength>0.))
    throw std::domain_error("CMC net potential requires blend in [0,1] and finite forces/ranges/reference length>0");
  if(p.cmcNetBlend>0.){
    const double saddle=p.roughnessGap+p.cmcNetAttractionRange;
    const double cutoff=saddle+p.cmcNetRepulsionRange;
    const double attractionEnergy=p.cmcNetContactForce*(p.cmcNetAttractionRange/3.);
    const double barrierEnergy=(8./15.)*p.cmcNetBarrierForce*p.cmcNetRepulsionRange;
    const double outerCoefficient=16.*p.cmcNetBarrierForce*p.cmcNetRepulsionRange;
    if(!p.surfaceAdhesion||p.localGapFraction!=0.||p.cohesionRetention!=1.
       ||effectiveContactGap(p)!=p.roughnessGap
       ||!std::isfinite(p.roughnessGap)||!(p.roughnessGap>0.)
       ||!std::isfinite(saddle)||!(saddle>p.roughnessGap)
       ||!std::isfinite(cutoff)||!(cutoff>saddle)||!(cutoff<=p.switchGap)
       ||!std::isfinite(1./p.cmcNetAttractionRange)
       ||!std::isfinite(1./p.cmcNetRepulsionRange)
       ||!std::isfinite(1./p.cmcNetReferenceLength)
       ||!std::isfinite(attractionEnergy)||!std::isfinite(barrierEnergy)
       ||!(attractionEnergy>0.)||!(barrierEnergy>0.)
       ||!std::isfinite(outerCoefficient)||!(outerCoefficient>0.)
       ||!std::isfinite(attractionEnergy/p.cmcNetReferenceLength)
       ||!std::isfinite(barrierEnergy/p.cmcNetReferenceLength))
      throw std::domain_error("Active CMC net potential requires unshifted surface adhesion, no legacy local-gap/screening, cutoff<=far switch, and representable scales");
  }
  if((effectiveContactGap(p)!=p.roughnessGap||p.cohesionRetention!=1.)
     &&p.localGapFraction>0.)
    throw std::domain_error("Coated-contact cohesion cannot be combined with legacy local-gap adhesion");
  if(!std::isfinite(p.freeCmcRepulsionPressure)||p.freeCmcRepulsionPressure<0.
     ||!std::isfinite(p.freeCmcRepulsionLength)||!(p.freeCmcRepulsionLength>0.))
    throw std::domain_error("Free-CMC repulsion requires finite pressure>=0 and finite length>0");
  if(p.freeCmcRepulsionPressure>0.){
    if(effectiveContactGap(p)!=p.roughnessGap||p.cohesionRetention!=1.)
      throw std::domain_error("Legacy free-CMC exponential cannot be combined with coated-contact cohesion");
    const double work=p.freeCmcRepulsionPressure*p.freeCmcRepulsionLength;
    const double energyPerLength=work*p.freeCmcRepulsionLength;
    if(!std::isfinite(p.roughnessGap)||!(p.roughnessGap>0.)
       ||!std::isfinite(1./p.freeCmcRepulsionLength)
       ||!std::isfinite(p.roughnessGap/p.freeCmcRepulsionLength)
       ||!std::isfinite(work)||!(work>0.)
       ||!std::isfinite(energyPerLength)||!(energyPerLength>0.))
      throw std::domain_error("Active free-CMC repulsion requires h0>0 and representable pressure-length scales");
  }
  if(!std::isfinite(p.freeCmcInnerRepulsionWork)||p.freeCmcInnerRepulsionWork<0.
     ||!std::isfinite(p.freeCmcInnerRepulsionRange)||!(p.freeCmcInnerRepulsionRange>0.)
     ||!std::isfinite(p.freeCmcInnerRepulsionPower)||!(p.freeCmcInnerRepulsionPower>2.))
    throw std::domain_error("Compact free-CMC repulsion requires finite work>=0, range>0, and power>2");
  if(p.freeCmcInnerRepulsionWork>0.){
    if(effectiveContactGap(p)!=p.roughnessGap||p.cohesionRetention!=1.)
      throw std::domain_error("Compact free-CMC repulsion cannot be combined with coated-contact cohesion");
    const double energyPerLength=p.freeCmcInnerRepulsionWork*
        (p.freeCmcInnerRepulsionRange/p.freeCmcInnerRepulsionPower);
    const double pressure=(p.freeCmcInnerRepulsionPower-1.)*
        (p.freeCmcInnerRepulsionWork/p.freeCmcInnerRepulsionRange);
    if(!std::isfinite(p.roughnessGap)||!(p.roughnessGap>0.)
       ||!std::isfinite(1./p.freeCmcInnerRepulsionRange)
       ||!std::isfinite(p.roughnessGap/p.freeCmcInnerRepulsionRange)
       ||!(p.roughnessGap+p.freeCmcInnerRepulsionRange>p.roughnessGap)
       ||!(p.roughnessGap+p.freeCmcInnerRepulsionRange<=p.switchGap)
       ||!std::isfinite(energyPerLength)||!(energyPerLength>0.)
       ||!std::isfinite(pressure)||!(pressure>0.))
      throw std::domain_error("Active compact free-CMC repulsion requires h0>0, cutoff<=far switch, and representable work-range scales");
  }
  if(p.localGapFraction>0.
     &&(!std::isfinite(p.roughnessGap)||!std::isfinite(p.localGap)||!(p.localGap>0.)||!(p.localGap<=p.roughnessGap)
        ||!std::isfinite(p.localSwitchExcessGap)||!(p.localSwitchExcessGap>=0.)
        ||!std::isfinite(p.localCutoffExcessGap)||!(p.localCutoffExcessGap>p.localSwitchExcessGap)
        ||!(p.roughnessGap+p.localCutoffExcessGap<=p.switchGap)))
    throw std::domain_error("Invalid RE2 local adhesion gaps: require 0<D0<=h0, 0<=local switch<local cutoff, and h0+local cutoff<=far switch");
  if(p.surfaceAdhesion){
    if(p.localGapFraction!=0.)
      throw std::domain_error("Surface adhesion cannot be combined with legacy local-gap adhesion");
    if(!std::isfinite(p.roughnessGap)||!(p.roughnessGap>0.)
       ||!std::isfinite(p.adhesionWork)||!(p.adhesionWork>0.)
       ||!std::isfinite(p.adhesionRange)||!(p.adhesionRange>0.)
       ||!std::isfinite(p.curvatureSwitchGap)||!std::isfinite(p.curvatureCutoffGap)
       ||!(effectiveContactGap(p)+p.adhesionRange<=p.curvatureSwitchGap)
       ||!(p.curvatureSwitchGap<p.curvatureCutoffGap)||!(p.curvatureCutoffGap<=p.switchGap))
      throw std::domain_error("Invalid surface adhesion: require h0>0, range>0, Hc+range<=curvature switch<curvature cutoff<=far switch");
    const double background=surfaceBackgroundWork(p);
    if(!std::isfinite(background)||background<0.||p.adhesionWork<background)
      throw std::domain_error("Surface adhesion work must be at least the nonnegative planar background work; repulsive screening is a separate model");
  }
}
// Strict lower bound on the geometric gap. Only the legacy local-gap model
// has an extra shifted-gap singularity. Surface adhesion uses a polynomial
// continuation for trial h<Hc; the accepted rough constraint remains h>=Hc.
inline double minimumPairGap(const PairParameters&p){
  return !p.surfaceAdhesion&&p.localGapFraction>0.?std::max(0.,p.roughnessGap-p.localGap):0.;
}
struct PairResult {
  Vec3 forceI{},torqueI{},torqueJ{};
  Vec3 forceAttractiveI{},forceRepulsiveI{};
  Vec3 torqueAttractiveI{},torqueAttractiveJ{},torqueRepulsiveI{},torqueRepulsiveJ{};
  double ua=0.,ur=0.,energy=0.,gap=std::numeric_limits<double>::infinity();
  Vec3 normal{},leverI{},leverJ{};
  bool active=false;
};
namespace re2_detail {
// Nine forward derivatives: relative centre displacement (3), laboratory
// infinitesimal rotation of I (3), laboratory infinitesimal rotation of J (3).
// Gap derivatives use the envelope theorem, so no perturbed closest-gap solves
// are required in the production force/torque evaluation.
struct AD {
  double v=0.;std::array<double,9> d{};
  AD()=default;AD(double value):v(value){}
};
inline AD operator+(const AD& a,const AD& b){AD c(a.v+b.v);for(int k=0;k<9;++k)c.d[k]=a.d[k]+b.d[k];return c;}
inline AD operator-(const AD& a,const AD& b){AD c(a.v-b.v);for(int k=0;k<9;++k)c.d[k]=a.d[k]-b.d[k];return c;}
inline AD operator-(const AD& a){AD c(-a.v);for(int k=0;k<9;++k)c.d[k]=-a.d[k];return c;}
inline AD operator*(const AD& a,const AD& b){AD c(a.v*b.v);for(int k=0;k<9;++k)c.d[k]=a.d[k]*b.v+a.v*b.d[k];return c;}
inline AD operator/(const AD& a,const AD& b){AD c(a.v/b.v);for(int k=0;k<9;++k)c.d[k]=(a.d[k]-c.v*b.d[k])/b.v;return c;}
inline AD sqrt(const AD& a){AD c(std::sqrt(a.v));for(int k=0;k<9;++k)c.d[k]=a.d[k]/(2.*c.v);return c;}
inline AD exp(const AD& a){AD c(std::exp(a.v));for(int k=0;k<9;++k)c.d[k]=c.v*a.d[k];return c;}
inline AD power(AD a,int n){AD b(1.);for(;n;n>>=1,a=a*a)if(n&1)b=b*a;return b;}
// The compact branch calls this only for nonnegative opening and power>2.
// At its cutoff both the derivative and its first derivative vanish.
inline AD realPower(const AD& a,double exponent){
  AD c(std::pow(a.v,exponent));
  const double slope=exponent*std::pow(a.v,exponent-1.);
  for(int k=0;k<9;++k)c.d[k]=slope*a.d[k];
  return c;
}
using AVec=std::array<AD,3>;using AMat=std::array<AD,9>;
inline AD dot(const AVec&a,const AVec&b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline AVec mv(const AMat&a,const AVec&b){AVec c{};for(int i=0;i<3;++i)for(int k=0;k<3;++k)c[i]=c[i]+a[3*i+k]*b[k];return c;}
inline AD determinant(const AMat&m){return m[0]*(m[4]*m[8]-m[5]*m[7])-m[1]*(m[3]*m[8]-m[5]*m[6])+m[2]*(m[3]*m[7]-m[4]*m[6]);}
inline AMat inverse(const AMat&m){
  AMat c{m[4]*m[8]-m[5]*m[7],m[2]*m[7]-m[1]*m[8],m[1]*m[5]-m[2]*m[4],
    m[5]*m[6]-m[3]*m[8],m[0]*m[8]-m[2]*m[6],m[2]*m[3]-m[0]*m[5],
    m[3]*m[7]-m[4]*m[6],m[1]*m[6]-m[0]*m[7],m[0]*m[4]-m[1]*m[3]};
  const AD det=determinant(m);for(auto& x:c)x=x/det;return c;
}
inline AMat rotation(const Mat3&r,int offset){
  AMat m;for(int i=0;i<9;++i)m[i]=AD(r[i]);
  for(int axis=0;axis<3;++axis){Vec3 e{};e[axis]=1.;for(int col=0;col<3;++col){const Vec3 v{r[col],r[3+col],r[6+col]};const Vec3 dv=cross(e,v);for(int row=0;row<3;++row)m[3*row+col].d[offset+axis]=dv[row];}}
  return m;
}
inline AMat shape(const AMat&r,const Vec3&d){AMat q{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)q[3*i+j]=q[3*i+j]+r[3*i+k]*d[k]*r[3*j+k];return q;}
inline AD orientationLength(const Body&bi,const Body&bj,const AVec&rh,double sigma){
  const AMat ri=rotation(bi.rotation,3),rj=rotation(bj.rotation,6);
  const double di=bi.axes[0]*bi.axes[1]*bi.axes[2],dj=bj.axes[0]*bj.axes[1]*bj.axes[2];
  Vec3 ai2{},aj2{},aiInv{},ajInv{};for(int k=0;k<3;++k){ai2[k]=bi.axes[k]*bi.axes[k];aj2[k]=bj.axes[k]*bj.axes[k];aiInv[k]=1./ai2[k];ajInv[k]=1./aj2[k];}
  const AMat qi=shape(ri,ai2),qj=shape(rj,aj2),ii=shape(ri,aiInv),ij=shape(rj,ajInv);
  const AD pi=1./sqrt(dot(rh,mv(ii,rh))),pj=1./sqrt(dot(rh,mv(ij,rh)));
  AMat h{},b{};for(int k=0;k<9;++k){h[k]=qi[k]/pi+qj[k]/pj;b[k]=qi[k]*(sigma/di)+qj[k]*(sigma/dj);}
  const AD chi=2.*dot(rh,mv(inverse(b),rh));
  const AD eta=(di/(pi*pi)+dj/(pj*pj))/sqrt(determinant(h)/(pi+pj));
  return eta*chi*sigma;
}
inline AD branch(const AD&h,const AD&ell,const Body&i,const Body&j,const PairParameters&p,bool repulsive){
  const double m=repulsive?2025.:-36.,o=repulsive?45./56.:3.,shapeScale=repulsive?std::cbrt(60.):2.;
  AD product(1.);for(int k=0;k<3;++k){product=product*(i.axes[k]/(i.axes[k]+h/shapeScale));product=product*(j.axes[k]/(j.axes[k]+h/shapeScale));}
  AD u=(p.hamaker/m)*(1.+o*ell/h)*product;
  if(repulsive)u=u*power(p.sigma/h,6);
  return u;
}
// Integral of the planar LJ surface energy over a local quadratic gap. Both
// attraction and repulsion use the same actual contact curvature.
inline AD derjaguinBranch(const AD&h,const AD&ell,const PairParameters&p,bool repulsive){
  if(repulsive)return (p.hamaker/2520.)*(ell/h)*power(p.sigma/h,6);
  return (-p.hamaker/12.)*(ell/h);
}
inline AD smoothSwitch(const AD&t){return 1.-10.*power(t,3)+15.*power(t,4)-6.*power(t,5);}
inline void unpack(const AD&u,Vec3&forceI,Vec3&torqueI,Vec3&torqueJ){for(int k=0;k<3;++k){forceI[k]=u.d[k];torqueI[k]=-u.d[3+k];torqueJ[k]=-u.d[6+k];}}
}
inline PairResult evaluatePair(const Body&bi,const Body&bj,const PairParameters&p=PairParameters{},GapCache*cache=nullptr){
  validatePairParameters(p);
  PairResult result;
  const Vec3 dr=sub(bj.position,bi.position);const double distance=norm(dr);
  const double bound=*std::max_element(bi.axes.begin(),bi.axes.end())+*std::max_element(bj.axes.begin(),bj.axes.end());
  if(distance>bound+p.cutoffGap)return result;
  const GapResult gap=closestEllipsoidGap(bi,bj,cache);
  result.gap=gap.gap;result.normal=gap.normal;result.leverI=gap.leverI;result.leverJ=gap.leverJ;
  if(gap.gap>=p.cutoffGap)return result;
  if(!(gap.gap>0.))throw std::domain_error("RE2 evaluated at overlapping ellipsoids: trial particle step must remain disjoint");
  if(!(gap.gap>minimumPairGap(p)))throw std::domain_error("RE2 evaluated outside positive local-gap domain: trial particle step must respect the local adhesion gap");
  if(!(distance>0.))throw std::domain_error("Coincident RE2 centres");
  using namespace re2_detail;
  AD h(gap.gap);const Vec3 gi=scale(cross(gap.leverI,gap.normal),-1.),gj=cross(gap.leverJ,gap.normal);
  for(int k=0;k<3;++k){h.d[k]=gap.normal[k];h.d[3+k]=gi[k];h.d[6+k]=gj[k];}
  AD ua,ur;
  // At full replacement do not evaluate a large legacy attraction and then
  // subtract it: the screened net model owns the entire pair energy.
  if(p.cmcNetBlend!=1.){
    AVec r{};for(int k=0;k<3;++k){r[k]=AD(dr[k]);r[k].d[k]=1.;}
    const AD length=sqrt(re2_detail::dot(r,r));AVec rh{};for(int k=0;k<3;++k)rh[k]=r[k]/length;
    const AD ell=orientationLength(bi,bj,rh,p.sigma);
    ua=branch(h,ell,bi,bj,p,false);ur=branch(h,ell,bi,bj,p,true);
    const bool localBackground=p.surfaceAdhesion&&gap.gap<p.curvatureCutoffGap;
    const bool innerRepulsion=p.freeCmcInnerRepulsionWork>0.
        &&gap.gap-p.roughnessGap<p.freeCmcInnerRepulsionRange;
    AD localLength;
    if(localBackground||p.freeCmcRepulsionPressure>0.||innerRepulsion){
      const auto curvature=contactCurvature(bi,bj,gap);
      localLength=AD(curvature.length);localLength.d=curvature.lengthDerivative;
    }
    if(localBackground){
      const AD nearA=derjaguinBranch(h,localLength,p,false);
      const AD nearR=derjaguinBranch(h,localLength,p,true);
      if(gap.gap<=p.curvatureSwitchGap){ua=nearA;ur=nearR;}
      else{
        const AD sw=smoothSwitch((h-p.curvatureSwitchGap)/(p.curvatureCutoffGap-p.curvatureSwitchGap));
        ua=ua+sw*(nearA-ua);ur=ur+sw*(nearR-ur);
      }
      const AD s=h-effectiveContactGap(p);
      if(s.v<p.adhesionRange){
        constexpr double pi=3.1415926535897932384626433832795;
        const double excessWork=p.adhesionWork-surfaceBackgroundWork(p);
        const AD opening=1.-s/p.adhesionRange;
        // G=pi*lambda_D=2*pi/sqrt(det_t K). The quadratic pair energy is
        // the Derjaguin integral of phi_coh=-DeltaW*(1-s/range)_+.
        // Its energy AND force vanish at the cutoff. For Newton trial points
        // s<0 use the same analytic polynomial; accepted states obey h>=Hc.
        // DeltaW retains the bare h0 benchmark when Hc moves with adsorption.
        // Differentiating localLength includes moving-contact force and torque.
        ua=ua-(.5*pi*excessWork*p.adhesionRange)*localLength*opening*opening;
      }
    }
    if(p.localGapFraction>0.&&gap.gap-p.roughnessGap<p.localCutoffExcessGap){
      const AD s=h-p.roughnessGap,d=s+p.localGap;
      // Differentiate the entire energy, including the local switch. Replacing
      // both branches avoids double counting and keeps forces and torques
      // conservative through the blend as well as at changing orientations.
      AD weight(p.localGapFraction);
      if(s.v>p.localSwitchExcessGap)
        weight=weight*smoothSwitch((s-p.localSwitchExcessGap)/(p.localCutoffExcessGap-p.localSwitchExcessGap));
      ua=ua+weight*(branch(d,ell,bi,bj,p,false)-ua);
      ur=ur+weight*(branch(d,ell,bi,bj,p,true)-ur);
    }
    if(p.freeCmcRepulsionPressure>0.){
      constexpr double pi=3.1415926535897932384626433832795;
      const double decay=p.freeCmcRepulsionLength;
      const double energyPerLength=(p.freeCmcRepulsionPressure*decay)*decay;
      // U=pi*lambda_D*P*ell^2*exp[-(H-h0)/ell]. Differentiate both the gap
      // and actual contact curvature; this term extends beyond the background's
      // curvature cutoff. Keep the same analytic continuation for H<h0 trials.
      ur=ur+localLength*(pi*energyPerLength)*exp(-(h-p.roughnessGap)/decay);
    }
    if(innerRepulsion){
      constexpr double pi=3.1415926535897932384626433832795;
      const AD opening=1.-(h-p.roughnessGap)/p.freeCmcInnerRepulsionRange;
      const double energyPerLength=p.freeCmcInnerRepulsionWork*
          (p.freeCmcInnerRepulsionRange/p.freeCmcInnerRepulsionPower);
      // Differentiate the complete energy, including moving-contact curvature.
      // Positive-gap Newton trials below h0 use the same analytic continuation;
      // neither the physical contact constraint nor the adhesion plane moves.
      ur=ur+localLength*(pi*energyPerLength)*realPower(opening,p.freeCmcInnerRepulsionPower);
    }
    // Screen the complete attractive energy before its force/torque derivatives
    // are unpacked. Avoid extra arithmetic in the unchanged pure-Gr branch.
    if(p.cohesionRetention!=1.)ua=ua*p.cohesionRetention;
    if(gap.gap>p.switchGap){const AD t=(h-p.switchGap)/(p.cutoffGap-p.switchGap);const AD sw=smoothSwitch(t);ua=ua*sw;ur=ur*sw;}
  }
  if(p.cmcNetBlend>0.){
    AD netA,netR;
    const double saddle=p.roughnessGap+p.cmcNetAttractionRange;
    const double cutoff=saddle+p.cmcNetRepulsionRange;
    if(gap.gap<cutoff){
      const auto curvature=contactCurvature(bi,bj,gap);
      AD localLength(curvature.length);localLength.d=curvature.lengthDerivative;
      const AD geometry=localLength/p.cmcNetReferenceLength;
      if(gap.gap<saddle){
        // The same cubic continues analytically below h0 for Newton trials;
        // it neither shifts nor softens the accepted hard-contact constraint.
        const AD opening=1.-(h-p.roughnessGap)/p.cmcNetAttractionRange;
        netA=-(p.cmcNetContactForce*(p.cmcNetAttractionRange/3.))*power(opening,3);
      }
      if(gap.gap<=saddle){
        netR=AD((8./15.)*p.cmcNetBarrierForce*p.cmcNetRepulsionRange);
      }else{
        // Integral of 16*FB*y^2*(1-y)^2, evaluated in z=1-y to
        // avoid cancellation at the outer cutoff. U, dU and d2U join.
        const AD z=(cutoff-h)/p.cmcNetRepulsionRange;
        netR=(16.*p.cmcNetBarrierForce*p.cmcNetRepulsionRange)*
            (power(z,3)/3.-power(z,4)/2.+power(z,5)/5.);
      }
      netA=geometry*netA;netR=geometry*netR;
    }
    if(p.cmcNetBlend==1.){ua=netA;ur=netR;}
    else{
      ua=(1.-p.cmcNetBlend)*ua+p.cmcNetBlend*netA;
      ur=(1.-p.cmcNetBlend)*ur+p.cmcNetBlend*netR;
    }
  }
  unpack(ua,result.forceAttractiveI,result.torqueAttractiveI,result.torqueAttractiveJ);
  unpack(ur,result.forceRepulsiveI,result.torqueRepulsiveI,result.torqueRepulsiveJ);
  result.forceI=add(result.forceAttractiveI,result.forceRepulsiveI);
  result.torqueI=add(result.torqueAttractiveI,result.torqueRepulsiveI);
  result.torqueJ=add(result.torqueAttractiveJ,result.torqueRepulsiveJ);
  result.ua=ua.v;result.ur=ur.v;result.energy=ua.v+ur.v;result.active=true;
  if(!std::isfinite(result.energy)||!finite(result.forceI)||!finite(result.torqueI)||!finite(result.torqueJ))throw std::runtime_error("Nonfinite RE2 force or energy");
  return result;
}
} // namespace graphite

} } // SLURRY SCOPE END
#endif
