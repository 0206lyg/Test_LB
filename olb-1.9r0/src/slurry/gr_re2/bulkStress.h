#ifndef SLURRY_GR_RE2_GRAPHITE_BULK_STRESS_H
#define SLURRY_GR_RE2_GRAPHITE_BULK_STRESS_H

#include <olb.h>
#include "particleMath.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <vector>

// SLURRY SCOPE BEGIN
namespace slurry { namespace gr_re2 {

namespace graphite {

// Cauchy stress convention: positive xy stress for positive imposed shear.
// All returned stress components are Pa. Moments supplied below are N*m.
struct BulkMeasure {
  double stressFluid=0., stressSurface=0.;
  double stressPairAttractive=0., stressPairRepulsive=0., stressLubrication=0.;
  double stressContactNormal=0., stressContactTangential=0.;
  double stressAcceleration=0., stressFluidReynolds=0., stressParticleReynolds=0.;
  double stressNonInertial=0., stressInertial=0., stressTotal=0.;
  double eta=0., etaRelative=0., solventViscosity=0.;
  double mach=0., particleMachUpperBound=0., densityDrift=0.;
  // Pure exterior cells only (porosity >= 0.99), in addition to the legacy
  // all-cell extrema which include fictitious fluid inside resolved particles.
  double fluidMach=0., fluidDensityDrift=0.;
  double phi=0., phiAnalytic=0., fluidVolume=0., sampledVolume=0.;
};

namespace bulk_detail {
inline Mat3 crossMatrix(const Vec3& w) {
  return {0.,-w[2],w[1],w[2],0.,-w[0],-w[1],w[0],0.};
}
inline double symmetricXY(const Mat3& a) { return .5*(a[1]+a[3]); }

struct RigidMoments {
  double accelerationXY=0.; // positive integral sym(r tensor rho*a)
  double reynoldsXY=0.;     // positive integral rho*u'_x*u'_y
};

// Exact rigid-ellipsoid volume integrals, without an extra voxel traversal.
// Q = integral rho*r*r^T dV; A = cross(alpha)+cross(omega)^2.
// The reference velocity is the imposed affine shear at every spatial point,
// including points inside the body, hence D = cross(omega)-grad(U_affine).
inline RigidMoments rigidMoments(const Body& body, const Vec3& alpha,
                                double gamma, double height) {
  const Vec3 diagonal{body.mass*body.axes[0]*body.axes[0]/5.,
                      body.mass*body.axes[1]*body.axes[1]/5.,
                      body.mass*body.axes[2]*body.axes[2]/5.};
  const Mat3 q=rotatedDiagonal(body.rotation,diagonal);
  const Mat3 w=crossMatrix(body.omega);
  Mat3 a=multiply(w,w);
  const Mat3 adot=crossMatrix(alpha);
  for (int k=0;k<9;++k) a[k]+=adot[k];
  Mat3 d=w;
  d[1]-=gamma;
  const Mat3 aq=multiply(a,q);
  const Mat3 dqd=multiply(multiply(d,q),transpose(d));
  const double cx=body.velocity[0]-gamma*(body.position[1]-.5*height);
  return {symmetricXY(aq),body.mass*cx*body.velocity[1]+dqd[1]};
}
} // namespace bulk_detail

// Bulk decomposition follows the surface-traction / acceleration / Reynolds
// formulation of Haddadi & Morris (2014), Eqs. 16--19:
// https://arxiv.org/abs/1403.7784 . The mass/density used here are the numerical
// values actually integrated, so inertia of the accelerated quick model stays
// visible; no small-Re assumption or Ma rejection is applied by this routine.
//
// Native HLBM Wen surface-MEA is treated as fluid-on-particle traction, as in
// Trunk et al. (2021), DOI:10.3390/computation9020011. It is NOT the negative
// volume integral of the collision forcing. Therefore no additional fictitious
// interior-fluid momentum derivative is added to this surface force/moment.
// The diffuse surface introduces the spatial approximation of this HLBM model.
//
// IMPORTANT INPUT CONTRACT:
// * Call after mapping CURRENT particle geometry: collision resets POROSITY.
// * hydroStressletSum is already globally summed +sym(lever tensor dF_hydro),
//   packed {xx,xy,xz,yy,yz,zz}, with dF_hydro the force ON the particle.
// * All pair/contact moments are already globally summed NEGATIVE virials,
//   -sum_{i<j}(r_ij tensor F_ij), NOT divided by box volume. The center branch
//   r_ij=x_i-x_j must be the SAME Lees--Edwards image used to evaluate F_ij.
// * Bodies, angular accelerations, and all supplied moments are replicated
//   global values. Only local fluid integrals are MPI-reduced here.
// * Fluid stress is weighted by actual exterior-fluid fraction. Neither an
//   extra solvent eta*gamma nor fictitious interior-fluid stress is added.
template<class LATTICE, class GEOMETRY, class CONVERTER>
BulkMeasure measureBulk(
    LATTICE& lattice, GEOMETRY& geometry, const CONVERTER& converter,
    const Vec3& box, double gamma, double omegaLB, double rhoNum,
    double dx, double dt, const std::vector<Body>& bodies,
    const std::vector<Vec3>& angularAcceleration,
    const std::array<double,6>& hydroStressletSum,
    const Mat3& pairAttractiveMoment, const Mat3& pairRepulsiveMoment,
    const Mat3& lubricationMoment,
    const Mat3& contactNormalMoment=Mat3{},
    const Mat3& contactTangentialMoment=Mat3{}) {
  (void)converter; // Units below are explicit and match numerical density.
  if (angularAcceleration.size()!=bodies.size())
    throw std::invalid_argument("Bulk stress requires one angular acceleration per body");
  const double volume=box[0]*box[1]*box[2];
  if (!(volume>0. && dx>0. && dt>0. && rhoNum>0.) || gamma==0.)
    throw std::invalid_argument("Invalid units or zero shear rate for bulk viscosity");
  const double cellVolume=dx*dx*dx;
  const double velocityUnit=dx/dt;
  const double pressureUnit=rhoNum*velocityUnit*velocityUnit;
  const double viscousFactor=-(1.-omegaLB/2.)*pressureUnit;
  // stress integral, fluid Reynolds integral, solid volume, fluid volume,
  // sampled volume. Exclude halo nodes using forCoreSpatialLocations.
  std::array<double,5> local{},global{};
  double mach=0.,densityDrift=0.,fluidMach=0.,fluidDensityDrift=0.;
  int bad=0;
  lattice.setProcessingContext(olb::ProcessingContext::Evaluation);
  for (int b=0;b<lattice.getLoadBalancer().size();++b) {
    auto& bg=geometry.getBlockGeometry(b);
    auto& bl=lattice.getBlock(b);
    bg.forCoreSpatialLocations([&](olb::LatticeR<3> loc) {
      if (bg.getMaterial(loc)!=1) return;
      auto cell=bl.get(loc);
      using Descriptor=typename std::decay_t<decltype(cell)>::descriptor_t;
      double rho,raw[3],velocity[3],pi[6];
      // Use the raw velocity of the BGK equilibrium for the EDM viscous
      // non-equilibrium moment. Exposed PorousParticle moments use a shifted
      // velocity without the corresponding off-diagonal stress correction.
      olb::lbm<Descriptor>::computeRhoU(cell,rho,raw);
      olb::lbm<Descriptor>::computeStress(cell,rho,raw,pi);
      cell.computeU(velocity); // native half-forcing velocity for transport
      const double epsilon=cell.template getField<olb::descriptors::POROSITY>();
      const double stress=viscousFactor*pi[1];
      if (!std::isfinite(rho) || !std::isfinite(epsilon) ||
          !std::isfinite(stress) || !std::isfinite(velocity[0]) ||
          !std::isfinite(velocity[1]) || !std::isfinite(velocity[2])) {
        bad=1; return;
      }
      olb::Vector<double,3> position;
      bg.getPhysR(position,loc);
      const double peculiarX=velocityUnit*velocity[0]-gamma*(position[1]-.5*box[1]);
      const double peculiarY=velocityUnit*velocity[1];
      local[0]+=epsilon*stress*cellVolume;
      local[1]+=epsilon*rhoNum*rho*peculiarX*peculiarY*cellVolume;
      local[2]+=(1.-epsilon)*cellVolume;
      local[3]+=epsilon*cellVolume;
      local[4]+=cellVolume;
      const double cellMach=std::sqrt(3.*(velocity[0]*velocity[0]+
                       velocity[1]*velocity[1]+velocity[2]*velocity[2]));
      mach=std::max(mach,cellMach);
      densityDrift=std::max(densityDrift,std::abs(rho-1.));
      if(epsilon>=.99) {
        fluidMach=std::max(fluidMach,cellMach);
        fluidDensityDrift=std::max(fluidDensityDrift,std::abs(rho-1.));
      }
    });
  }
#ifdef PARALLEL_MODE_MPI
  olb::singleton::mpi().allreduce(local.data(),global.data(),5,MPI_SUM);
  olb::singleton::mpi().reduceAndBcast(mach,MPI_MAX);
  olb::singleton::mpi().reduceAndBcast(densityDrift,MPI_MAX);
  olb::singleton::mpi().reduceAndBcast(fluidMach,MPI_MAX);
  olb::singleton::mpi().reduceAndBcast(fluidDensityDrift,MPI_MAX);
  olb::singleton::mpi().reduceAndBcast(bad,MPI_MAX);
#else
  global=local;
#endif
  if (bad) throw std::runtime_error("Nonfinite lattice state in bulk stress measurement");
  BulkMeasure result;
  result.stressFluid=global[0]/volume;
  result.stressFluidReynolds=-global[1]/volume;
  result.stressSurface=hydroStressletSum[1]/volume;
  result.stressPairAttractive=bulk_detail::symmetricXY(pairAttractiveMoment)/volume;
  result.stressPairRepulsive=bulk_detail::symmetricXY(pairRepulsiveMoment)/volume;
  result.stressLubrication=bulk_detail::symmetricXY(lubricationMoment)/volume;
  result.stressContactNormal=bulk_detail::symmetricXY(contactNormalMoment)/volume;
  result.stressContactTangential=bulk_detail::symmetricXY(contactTangentialMoment)/volume;
  result.phi=global[2]/volume;
  result.fluidVolume=global[3];
  result.sampledVolume=global[4];
  result.mach=mach;
  result.densityDrift=densityDrift;
  result.fluidMach=fluidMach;
  result.fluidDensityDrift=fluidDensityDrift;
  for (std::size_t i=0;i<bodies.size();++i) {
    const auto& body=bodies[i];
    const auto moments=bulk_detail::rigidMoments(body,angularAcceleration[i],gamma,box[1]);
    result.stressAcceleration-=moments.accelerationXY/volume;
    result.stressParticleReynolds-=moments.reynoldsXY/volume;
    result.phiAnalytic+=(4.*std::acos(-1.)/3.)*body.axes[0]*body.axes[1]*body.axes[2]/volume;
    const double maxAxis=*std::max_element(body.axes.begin(),body.axes.end());
    // Conservative surface speed bound, separate from measured fluid Ma.
    result.particleMachUpperBound=std::max(result.particleMachUpperBound,
      std::sqrt(3.)*(norm(body.velocity)+maxAxis*norm(body.omega))/velocityUnit);
  }
  result.stressNonInertial=result.stressFluid+result.stressSurface+
    result.stressPairAttractive+result.stressPairRepulsive+result.stressLubrication+
    result.stressContactNormal+result.stressContactTangential;
  result.stressInertial=result.stressAcceleration+result.stressFluidReynolds+
    result.stressParticleReynolds;
  result.stressTotal=result.stressNonInertial+result.stressInertial;
  result.eta=result.stressTotal/gamma;
  const double nuLB=(1./omegaLB-.5)/3.;
  result.solventViscosity=rhoNum*nuLB*dx*dx/dt;
  result.etaRelative=result.eta/result.solventViscosity;
  if (!std::isfinite(result.stressTotal) || !std::isfinite(result.etaRelative) ||
      !std::isfinite(result.particleMachUpperBound))
    throw std::runtime_error("Nonfinite particle state or stress in bulk measurement");
  return result;
}

} // namespace graphite

} } // SLURRY SCOPE END
#endif
