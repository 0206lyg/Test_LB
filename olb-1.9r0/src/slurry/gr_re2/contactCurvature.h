#ifndef SLURRY_GR_RE2_GRAPHITE_CONTACT_CURVATURE_H
#define SLURRY_GR_RE2_GRAPHITE_CONTACT_CURVATURE_H

#include "ellipsoidGap.h"

// SLURRY SCOPE BEGIN
namespace slurry { namespace gr_re2 {
namespace graphite {

struct ContactCurvatureResult {
  double rootDet=0.;
  double length=0.; // Derjaguin length: 2/sqrt(det_t(K_I+K_J)).
  // Relative centre displacement (3), laboratory rotation of I (3), and
  // laboratory rotation of J (3), in the same order as the pair-energy AD.
  std::array<double,9> lengthDerivative{};
};

namespace contact_curvature_detail {
inline Mat3 angularDerivative(const Mat3& m,int axis) {
  Vec3 angle{};angle[axis]=1.;
  const Mat3 skew{0.,-angle[2],angle[1],angle[2],0.,-angle[0],-angle[1],angle[0],0.};
  Mat3 result=multiply(skew,m);const Mat3 right=multiply(m,skew);
  for(int k=0;k<9;++k)result[k]-=right[k];
  return result;
}
inline Mat3 adjugate(const Mat3& m) {
  return {m[4]*m[8]-m[5]*m[7],m[2]*m[7]-m[1]*m[8],m[1]*m[5]-m[2]*m[4],
          m[5]*m[6]-m[3]*m[8],m[0]*m[8]-m[2]*m[6],m[2]*m[3]-m[0]*m[5],
          m[3]*m[7]-m[4]*m[6],m[1]*m[6]-m[0]*m[7],m[0]*m[4]-m[1]*m[3]};
}
inline Mat3 adjugateDerivative(const Mat3& m,const Mat3& dm) {
  const auto difference=[&](int a,int b,int c,int d) {
    return dm[a]*m[b]+m[a]*dm[b]-dm[c]*m[d]-m[c]*dm[d];
  };
  return {difference(4,8,5,7),difference(2,7,1,8),difference(1,5,2,4),
          difference(5,6,3,8),difference(0,8,2,6),difference(2,3,0,5),
          difference(3,7,4,6),difference(1,6,0,7),difference(0,4,1,3)};
}
}

// Curvatures are evaluated at the actual closest surface points, not along
// the centre line. For Q=R diag(a_k^2) R^T the tangent curvature is the tangent
// restriction of sqrt(n.Q.n) Q^{-1}, exactly as in nfCurvature().
//
// The closest normal also changes with translation and rotation. Differentiate
// r-p_I(n)-p_J(n)=h*n and project onto the tangent plane:
//   [H_I+H_J+h I]_t dn = [dr-partial(p_I)-partial(p_J)]_t,
// where H=Q/s-(Qn)(Qn)^T/s^3 and s=sqrt(n.Q.n). This is an analytic 2x2
// implicit solve; production evaluation never perturbs the closest-gap solver.
inline ContactCurvatureResult contactCurvature(
    const Body& bi,const Body& bj,const GapResult& gap) {
  if(!(gap.gap>=0.)||!std::isfinite(gap.gap))
    throw std::domain_error("Contact curvature requires disjoint ellipsoids");
  const Vec3 n=normalized(gap.normal);
  const Vec3 reference=std::abs(n[0])<.8?Vec3{1.,0.,0.}:Vec3{0.,1.,0.};
  const Vec3 t1=normalized(cross(n,reference)),t2=cross(n,t1);
  Vec3 ai2{},aj2{},aiInv{},ajInv{};
  for(int k=0;k<3;++k) {
    if(!(bi.axes[k]>0.)||!(bj.axes[k]>0.))
      throw std::domain_error("Contact curvature requires positive ellipsoid axes");
    ai2[k]=bi.axes[k]*bi.axes[k];aj2[k]=bj.axes[k]*bj.axes[k];
    aiInv[k]=1./ai2[k];ajInv[k]=1./aj2[k];
  }
  const Mat3 qi=rotatedDiagonal(bi.rotation,ai2),qj=rotatedDiagonal(bj.rotation,aj2);
  const Mat3 ii=rotatedDiagonal(bi.rotation,aiInv),ij=rotatedDiagonal(bj.rotation,ajInv);
  const Vec3 vi=mul(qi,n),vj=mul(qj,n);
  const double si=std::sqrt(dot(n,vi)),sj=std::sqrt(dot(n,vj));
  Mat3 curvature{},normalHessian{};
  for(int k=0;k<3;++k)for(int l=0;l<3;++l) {
    const int kl=3*k+l;
    curvature[kl]=si*ii[kl]+sj*ij[kl];
    normalHessian[kl]=qi[kl]/si-vi[k]*vi[l]/(si*si*si)
                       +qj[kl]/sj-vj[k]*vj[l]/(sj*sj*sj)
                       +(k==l?gap.gap:0.);
  }
  const double b11=dot(t1,mul(normalHessian,t1));
  const double b12=dot(t1,mul(normalHessian,t2));
  const double b22=dot(t2,mul(normalHessian,t2));
  const double bdet=b11*b22-b12*b12;
  if(!(b11>0.)||!(b22>0.)||!(bdet>0.)||!std::isfinite(bdet))
    throw std::domain_error("Degenerate closest-normal curvature derivative");

  // det_t K=n.adj(K).n is invariant under tangent-basis changes. Differentiating
  // this expression includes the change of the tangent plane itself; keeping
  // a numerical tangent basis fixed when differentiating K would omit it.
  const Mat3 adj=contact_curvature_detail::adjugate(curvature);
  const Vec3 adjN=mul(adj,n);
  const double det=dot(n,adjN);
  if(!(det>0.)||!std::isfinite(det))
    throw std::domain_error("Degenerate ellipsoid contact curvature");
  ContactCurvatureResult result;
  result.rootDet=std::sqrt(det);result.length=2./result.rootDet;
  for(int dof=0;dof<9;++dof) {
    Mat3 dqi{},dqj{},dii{},dij{};Vec3 dr{};
    if(dof<3)dr[dof]=1.;
    else if(dof<6) {
      dqi=contact_curvature_detail::angularDerivative(qi,dof-3);
      dii=contact_curvature_detail::angularDerivative(ii,dof-3);
    } else {
      dqj=contact_curvature_detail::angularDerivative(qj,dof-6);
      dij=contact_curvature_detail::angularDerivative(ij,dof-6);
    }
    const Vec3 dvi=mul(dqi,n),dvj=mul(dqj,n);
    const double partialSi=.5*dot(n,dvi)/si,partialSj=.5*dot(n,dvj)/sj;
    const Vec3 partialPi=sub(scale(dvi,1./si),scale(vi,partialSi/(si*si)));
    const Vec3 partialPj=sub(scale(dvj,1./sj),scale(vj,partialSj/(sj*sj)));
    const Vec3 rhs=sub(sub(dr,partialPi),partialPj);
    const double rhs1=dot(t1,rhs),rhs2=dot(t2,rhs);
    const Vec3 dn=add(scale(t1,(b22*rhs1-b12*rhs2)/bdet),
                      scale(t2,(b11*rhs2-b12*rhs1)/bdet));
    const double dsi=partialSi+dot(dn,vi)/si,dsj=partialSj+dot(dn,vj)/sj;
    Mat3 dc{};
    for(int k=0;k<9;++k)dc[k]=dsi*ii[k]+si*dii[k]+dsj*ij[k]+sj*dij[k];
    const Mat3 dadj=contact_curvature_detail::adjugateDerivative(curvature,dc);
    const double ddet=2.*dot(dn,adjN)+dot(n,mul(dadj,n));
    result.lengthDerivative[dof]=-.5*result.length*ddet/det;
    if(!std::isfinite(result.lengthDerivative[dof]))
      throw std::runtime_error("Nonfinite contact-curvature derivative");
  }
  return result;
}

} // namespace graphite
} } // SLURRY SCOPE END
#endif
