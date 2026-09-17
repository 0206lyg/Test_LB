/* SPDX-License-Identifier: GPL-2.0-or-later
 * Application-local Cartesian periodic-image treatment for OpenLB 1.9r0.
 *
 * The native resolved ParticleSystem wrappers evaluate the original particle
 * and one combined ghost.  A particle crossing both x and z boundaries needs
 * the two singly shifted images as well.  This header leaves the OpenLB source
 * untouched and calls the native single-image kernels for every image.
 */
#ifndef SLURRY_GR_BASELINE_GRAPHITE_PERIODIC_COUPLING_H
#define SLURRY_GR_BASELINE_GRAPHITE_PERIODIC_COUPLING_H

#include <olb.h>
#include <limits>
#include <stdexcept>
#include <vector>

// SLURRY SCOPE BEGIN
namespace slurry { namespace gr_baseline {

namespace graphiteCoupling {

// Including the nearest image on each side also covers lattice overlap cells
// and contact-detection shells.  A non-intersecting image is rejected by the
// native block/particle intersection test before any lattice-cell traversal.
template <typename T, unsigned DIM>
std::vector<olb::PhysR<T,DIM>> periodicImages(
    const olb::PhysR<T,DIM>& position,
    const olb::PhysR<T,DIM>& minimum,
    const olb::PhysR<T,DIM>& maximum,
    const olb::Vector<bool,DIM>& periodic) {
  std::vector<olb::PhysR<T,DIM>> images{position};
  for (unsigned k=0; k<DIM; ++k) {
    if (!periodic[k]) continue;
    const T length=maximum[k]-minimum[k];
    if (!(length>0)) throw std::runtime_error("Invalid periodic domain length");
    const auto count=images.size();
    for (std::size_t i=0; i<count; ++i) {
      auto before=images[i],after=images[i];
      before[k]-=length; after[k]+=length;
      images.push_back(before); images.push_back(after);
    }
  }
  return images;
}

template <typename T,typename P>
struct RestorePosition {
  olb::particles::Particle<T,P>& particle;
  olb::PhysR<T,P::d> original;
  explicit RestorePosition(olb::particles::Particle<T,P>& p)
    : particle(p),original(olb::particles::access::getPosition(p)) {}
  ~RestorePosition() {
    particle.template setField<olb::descriptors::GENERAL,
                               olb::descriptors::POSITION>(original);
  }
};

template <typename T,typename D,typename P,typename PC,typename WC,typename F>
void mapParticles(
    olb::particles::ParticleSystem<T,P>& particles,
    olb::particles::contact::ContactContainer<T,PC,WC>& contacts,
    olb::SuperGeometry<T,D::d>& geometry,
    olb::SuperLattice<T,D>& lattice,
    const olb::UnitConverter<T,D>& converter,
    std::vector<olb::SolidBoundary<T,D::d>>& walls,
    F periodicity) {
  static_assert(D::d==3 && P::d==3,"This application requires 3D particles");
  const auto minimum=olb::particles::communication::getCuboidMin<T,D::d>(
      geometry.getCuboidDecomposition());
  const auto maximum=olb::particles::communication::getCuboidMax<T,D::d>(
      geometry.getCuboidDecomposition(),minimum);

  // Rebuild both core and overlap fields.  Plain-BGK wall cells do not execute
  // PorousParticle's contact-field reset and would otherwise retain stale IDs.
  olb::particles::resetSuperParticleField(geometry,lattice);
  contacts.cleanContacts();
  for (std::size_t i=0; i<particles.size(); ++i) {
    auto particle=particles.get(i);
    RestorePosition<T,P> guard(particle);
    const auto images=periodicImages<T,D::d>(
        guard.original,minimum,maximum,periodicity());
    for (const auto& position: images) {
      particle.template setField<olb::descriptors::GENERAL,
                                 olb::descriptors::POSITION>(position);
      for (int b=0; b<lattice.getLoadBalancer().size(); ++b) {
        // This overload is the SINGLE-image kernel.  Keep true periodicity in
        // contact-coordinate unification; do not call the outer ghost wrapper.
        olb::particles::setBlockParticleField(
            geometry.getBlockGeometry(b),lattice.getBlock(b),converter,
            particles,contacts,i,particle,walls,minimum,maximum,periodicity);
      }
    }
  }
  olb::particles::contact::communicateContacts<T,PC,WC>(contacts);
}

template <typename T,typename D,typename P>
void coupleFluidToParticles(
    olb::particles::ParticleSystem<T,P>& particles,
    const olb::SuperGeometry<T,D::d>& geometry,
    olb::SuperLattice<T,D>& lattice,
    const olb::UnitConverter<T,D>& converter,
    const olb::Vector<bool,D::d>& periodic) {
  static_assert(D::d==3 && P::d==3,"This application requires 3D particles");
  constexpr std::size_t stride=7; // Fx,Fy,Fz,Tx,Ty,Tz,voxel count
  if (particles.size()>std::size_t(std::numeric_limits<int>::max())/stride)
    throw std::runtime_error("Too many particles for one MPI force reduction");
  if (particles.size()==0) return;
  const auto minimum=olb::particles::communication::getCuboidMin<T,D::d>(
      geometry.getCuboidDecomposition());
  const auto maximum=olb::particles::communication::getCuboidMax<T,D::d>(
      geometry.getCuboidDecomposition(),minimum);
  std::vector<T> local(stride*particles.size(),T{}),global(local.size(),T{});
  for (int b=0; b<lattice.getLoadBalancer().size(); ++b) {
    olb::BlockLatticeMomentumExchangeForce<T,D,P> kernel(
        lattice.getBlock(b),geometry.getBlockGeometry(b),particles,converter);
    for (std::size_t i=0; i<particles.size(); ++i) {
      auto particle=particles.get(i);
      if (!olb::particles::access::isValid(particle)) continue;
      RestorePosition<T,P> guard(particle);
      const auto images=periodicImages<T,D::d>(
          guard.original,minimum,maximum,periodic);
      for (const auto& position: images) {
        particle.template setField<olb::descriptors::GENERAL,
                                   olb::descriptors::POSITION>(position);
        // evaluate() is public and performs exactly one image.  Its torque
        // lever is relative to that same image centre, preserving rotations.
        kernel.evaluate(local.data(),particle,static_cast<int>(i));
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  olb::singleton::mpi().allreduce(local.data(),global.data(),
      static_cast<int>(local.size()),MPI_SUM);
#else
  global=local;
#endif
  for (std::size_t i=0; i<particles.size(); ++i) {
    auto particle=particles.get(i);
    const std::size_t offset=stride*i;
    particle.template setField<olb::descriptors::FORCING,
                               olb::descriptors::FORCE>(
        olb::Vector<T,3>{global[offset],global[offset+1],global[offset+2]});
    particle.template setField<olb::descriptors::FORCING,
                               olb::descriptors::TORQUE>(
        olb::Vector<T,3>{global[offset+3],global[offset+4],global[offset+5]});
  }
}

} // namespace graphiteCoupling

} } // SLURRY SCOPE END
#endif
