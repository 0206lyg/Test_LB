#ifndef SLURRY_GR_RE2_GRAPHITE_CMC_COORDINATION_H
#define SLURRY_GR_RE2_GRAPHITE_CMC_COORDINATION_H

#include "reSquaredPotential.h"
#include <vector>

// SLURRY SCOPE BEGIN
namespace slurry { namespace gr_re2 { namespace graphite {

struct CoordinationPair {
  std::size_t i=0,j=0;
  PairResult result{};
};

namespace coordination_detail {
struct Attenuation {double value=1.,derivative=0.;};
inline Attenuation attenuation(double others,const PairParameters& p){
  if(others<=p.cmcCoordinationStart)return {};
  if(others>=p.cmcCoordinationEnd)return {p.cmcCoordinationFloor,0.};
  const double width=p.cmcCoordinationEnd-p.cmcCoordinationStart;
  const double t=(others-p.cmcCoordinationStart)/width,z=1.-t;
  const double sw=t<=.5?1.-t*t*t*(10.-15.*t+6.*t*t):z*z*z*(10.-15.*z+6.*z*z);
  return {p.cmcCoordinationFloor+(1.-p.cmcCoordinationFloor)*sw,
          -30.*(1.-p.cmcCoordinationFloor)*t*t*z*z/width};
}
struct EdgeDerivative {
  double factor=1.,atI=0.,atJ=0.;
};
}

// Call once on fresh pair evaluations for the same complete body state and
// periodic images. Each unordered body pair must appear exactly once; inactive
// default PairResults are allowed. The caller retains the same i->j image
// displacement for the corrected pair virial. This routine is O(bodies+pairs).
//
// E_A=sum_e A_e*g(z_i-q_e)*g(z_j-q_e), z_i=sum_{e incident i} q_e.
// A_e is ONLY the blend-weighted net attraction. All legacy attraction and all
// repulsion retain their original energy and derivatives. Reverse accumulation
// includes derivatives through every other bond's occupancy; freezing the
// coordination number would not be conservative.
inline void applyCmcCoordination(std::size_t bodyCount,
                                 std::vector<CoordinationPair>& pairs,
                                 const PairParameters& p){
  if(!cmcCoordinationActive(p))return;
  validatePairParameters(p);
  std::vector<double> coordination(bodyCount,0.),sensitivity(bodyCount,0.);
  std::vector<coordination_detail::EdgeDerivative> edges(pairs.size());
  for(const auto& pair:pairs){
    if(pair.i>=bodyCount||pair.j>=bodyCount||pair.i==pair.j)
      throw std::domain_error("Invalid body indices in CMC coordination assembly");
    const auto& r=pair.result;
    if(r.cmcCoordinationApplied)
      throw std::logic_error("CMC coordination requires fresh pair results, not repeated correction");
    if(!std::isfinite(r.cmcCoordinationOccupancy)||r.cmcCoordinationOccupancy<0.
       ||r.cmcCoordinationOccupancy>1.||!std::isfinite(r.cmcCoordinationGapDerivative)
       ||!std::isfinite(r.cmcNetAttractiveEnergy))
      throw std::domain_error("Nonfinite or invalid CMC coordination pair metadata");
    coordination[pair.i]+=r.cmcCoordinationOccupancy;
    coordination[pair.j]+=r.cmcCoordinationOccupancy;
  }
  for(std::size_t k=0;k<pairs.size();++k){
    const auto& pair=pairs[k];const auto& r=pair.result;
    const auto gi=coordination_detail::attenuation(coordination[pair.i]-r.cmcCoordinationOccupancy,p);
    const auto gj=coordination_detail::attenuation(coordination[pair.j]-r.cmcCoordinationOccupancy,p);
    auto& edge=edges[k];edge.factor=gi.value*gj.value;
    edge.atI=r.cmcNetAttractiveEnergy*gi.derivative*gj.value;
    edge.atJ=r.cmcNetAttractiveEnergy*gi.value*gj.derivative;
    sensitivity[pair.i]+=edge.atI;sensitivity[pair.j]+=edge.atJ;
  }
  for(std::size_t k=0;k<pairs.size();++k){
    auto& pair=pairs[k];auto& r=pair.result;const auto& edge=edges[k];
    // Own q_e is excluded at BOTH endpoints. Its only influence is on the
    // energies of the other bonds sharing those endpoints.
    const double occupancyEnergyDerivative=(sensitivity[pair.i]-edge.atI)
                                          +(sensitivity[pair.j]-edge.atJ);
    const double gapDerivative=occupancyEnergyDerivative*r.cmcCoordinationGapDerivative;
    const double change=edge.factor-1.;
    if(change!=0.||gapDerivative!=0.){
      const Vec3 force=add(scale(r.cmcNetAttractiveForceI,change),scale(r.normal,gapDerivative));
      const Vec3 torqueI=add(scale(r.cmcNetAttractiveTorqueI,change),
                            scale(cross(r.leverI,r.normal),gapDerivative));
      const Vec3 torqueJ=sub(scale(r.cmcNetAttractiveTorqueJ,change),
                            scale(cross(r.leverJ,r.normal),gapDerivative));
      r.forceAttractiveI=add(r.forceAttractiveI,force);r.forceI=add(r.forceI,force);
      r.torqueAttractiveI=add(r.torqueAttractiveI,torqueI);r.torqueI=add(r.torqueI,torqueI);
      r.torqueAttractiveJ=add(r.torqueAttractiveJ,torqueJ);r.torqueJ=add(r.torqueJ,torqueJ);
      const double energy=change*r.cmcNetAttractiveEnergy;
      r.ua+=energy;r.energy+=energy;
      if(!std::isfinite(r.energy)||!finite(r.forceI)||!finite(r.torqueI)||!finite(r.torqueJ))
        throw std::runtime_error("Nonfinite CMC coordination force or energy");
    }
    r.cmcCoordinationApplied=true;
  }
}

} } } // SLURRY SCOPE END
#endif
