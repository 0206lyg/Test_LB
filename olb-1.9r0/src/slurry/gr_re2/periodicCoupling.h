/* SPDX-License-Identifier: GPL-2.0-or-later
 * Cached application-local HLBM coupling, OpenLB 1.9r0.
 * Native smooth indicator and Wen momentum exchange are retained. Each image
 * evaluates solid fraction once in its rotated ellipsoid AABB, then mapping,
 * force, torque and traction moment reuse that data and its surface-link list.
 * LE images: (x+nx*Lx+ny*gamma*Ly*t,y+ny*Ly,z+nz*Lz), vx+ny*gamma*Ly.
 */
#ifndef SLURRY_GR_RE2_GRAPHITE_PERIODIC_COUPLING_H
#define SLURRY_GR_RE2_GRAPHITE_PERIODIC_COUPLING_H
#include <olb.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

// SLURRY SCOPE BEGIN
namespace slurry { namespace gr_re2 {

namespace graphiteCoupling {
template <typename T,typename D,typename P>
class CachedCoupling {
  static_assert(D::d==3 && P::d==3,"CachedCoupling requires 3D particles");
  static_assert(D::q<=32,"Surface-link mask requires q<=32");
  using V=olb::Vector<T,3>;
  using L=olb::LatticeR<3>;
  using PS=olb::particles::ParticleSystem<T,P>;
  using Particle=olb::particles::Particle<T,P>;
public:
  struct Hydrodynamics {
    std::array<T,3> force{},torque{};
    // Symmetric fluid-on-particle traction moment, N m: xx,xy,xz,yy,yz,zz.
    // Enters bulk stress with PLUS sign, together with EXTERIOR fluid stress.
    // This is surface MEA, not a volumetric forcing reaction. Do not add an
    // IBM displaced-fluid acceleration correction to this force definition.
    std::array<T,6> stresslet{};
    T sampledVolume=0;
  };
  struct Statistics {
    std::size_t sdfEvaluations=0,mappedCells=0;
    std::size_t boundaryCells=0,boundaryLinks=0,activeImages=0;
  };
private:
  struct SurfaceCell { L location; std::uint32_t mask=0; };
  struct ImageCache {
    int block=0;std::size_t particle=0;
    V center{},halfExtent{};T velocityJump=0;
    L lower{},upper{},size{};
    std::vector<T> fraction;
    std::vector<L> mapped;
    std::vector<SurfaceCell> surface;
    std::size_t index(int x,int y,int z) const {
      return (std::size_t(x-lower[0])*std::size_t(size[1])
             +std::size_t(y-lower[1]))*std::size_t(size[2])
             +std::size_t(z-lower[2]);
    }
    T phi(int x,int y,int z) const{return fraction[index(x,y,z)];}
  };
  olb::SuperGeometry<T,3>& _geometry;
  olb::SuperLattice<T,D>& _lattice;
  const olb::UnitConverter<T,D>& _converter;
  V _box,_axes;
  T _gamma,_mapTime=std::numeric_limits<T>::quiet_NaN();
  std::vector<ImageCache> _images;
  std::size_t _activeImages=0;
  std::vector<T> _local,_global,_localVolumes;
  std::vector<Hydrodynamics> _hydro;
  std::array<T,6> _stressletSum{};
  Statistics _statistics{};
  static constexpr std::size_t stride=13;
  static V cross(const V&a,const V&b) {
    return V{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
  }
  static void addMoment(T*s,const V&r,const V&f) {
    s[0]+=r[0]*f[0];s[1]+=T(.5)*(r[0]*f[1]+r[1]*f[0]);
    s[2]+=T(.5)*(r[0]*f[2]+r[2]*f[0]);s[3]+=r[1]*f[1];
    s[4]+=T(.5)*(r[1]*f[2]+r[2]*f[1]);s[5]+=r[2]*f[2];
  }
  V halfExtent(Particle&particle) const {
    using namespace olb::descriptors;
    const auto rotation=olb::particles::access::getRotationMatrix(particle);
    const T shell=T(.5)*particle.template getField<SURFACE,SINDICATOR>()->getEpsilon();
    // Native sdf::ellipsoid is a lower-bound distance approximation outside
    // the ellipsoid. Its diffuse support is therefore not necessarily inside
    // the exact Minkowski sum (ellipsoid + epsilon/2 sphere). For k0=|x/a|,
    // native d=(k0-1)/|n/a| and |n/a|<=1/min(a); d<=shell implies
    // k0<=1+shell/min(a). This scaled ellipsoid is a guaranteed support bound.
    const T supportScale=T(1)+shell/std::min({_axes[0],_axes[1],_axes[2]});
    V result;
    for(int i=0;i<3;++i) {
      T square=0;
      for(int j=0;j<3;++j)square+=rotation[3*i+j]*rotation[3*i+j]*_axes[j]*_axes[j];
      result[i]=std::sqrt(square)*supportScale;
    }
    return result;
  }
  void resetMappedFields() {
    // Sparse reset handles stale overlap fields and repeated mapping without
    // another whole-domain sweep. Core fields are also reset by HLBM collision.
    for(std::size_t k=0;k<_activeImages;++k) {
      auto&cache=_images[k];auto&block=_lattice.getBlock(cache.block);
      for(const auto&location:cache.mapped) {
        auto cell=block.get(location);
        olb::particles::resetAllParticleRelatedFields<D,decltype(cell),T>(cell);
      }
    }
  }
  void buildImage(Particle&particle,std::size_t id,int blockId,
                  const V&center,const V&extent,T velocityJump) {
    using namespace olb::descriptors;
    auto&bg=_geometry.getBlockGeometry(blockId);auto&bl=_lattice.getBlock(blockId);
    const V origin=bg.getOrigin();const auto blockSize=bg.getExtent();
    const int padding=bg.getPadding();const T dx=bg.getDeltaR();
    L mapLower,mapUpper,forceLower,forceUpper;
    const T nativeRadius=olb::particles::access::getRadius(particle);
    for(int k=0;k<3;++k) {
      mapLower[k]=std::max(-padding,int(std::floor((center[k]-extent[k]-origin[k])/dx)));
      mapUpper[k]=std::min(blockSize[k]-1+padding,int(std::ceil((center[k]+extent[k]-origin[k])/dx)));
      if(mapLower[k]>mapUpper[k])return;
      // Match the original force functor's circumscribed-cube traversal even
      // when the conservative SDF cache extends beyond that traversal.
      forceLower[k]=std::max(0,int(std::floor((center[k]-nativeRadius-origin[k])/dx)));
      forceUpper[k]=std::min(blockSize[k]-1,int(std::ceil((center[k]+nativeRadius-origin[k])/dx)));
    }
    if(_activeImages==_images.size())_images.emplace_back();
    auto&cache=_images[_activeImages++];cache.block=blockId;cache.particle=id;
    cache.center=center;cache.halfExtent=extent;cache.velocityJump=velocityJump;
    cache.mapped.clear();cache.surface.clear();
    // One SDF halo contains all D3Q19 neighbors. Only core centers generate
    // force links, so no lattice access beyond its allocated overlap occurs.
    for(int k=0;k<3;++k){cache.lower[k]=mapLower[k]-1;cache.upper[k]=mapUpper[k]+1;
      cache.size[k]=cache.upper[k]-cache.lower[k]+1;}
    const std::size_t count=std::size_t(cache.size[0])*cache.size[1]*cache.size[2];
    cache.fraction.resize(count);
    auto indicator=particle.template getField<SURFACE,SINDICATOR>();
    const auto inverseRotation=olb::util::invertRotationMatrix<T,3>(
        particle.template getField<SURFACE,ROT_MATRIX>());
    const T epsilon=indicator->getEpsilon();
    // Same native signed distance and smoothing as evalSolidVolumeFraction.
    // Rotation inversion is common to the image, not repeated for every cell.
    for(int x=cache.lower[0];x<=cache.upper[0];++x)
      for(int y=cache.lower[1];y<=cache.upper[1];++y)
        for(int z=cache.lower[2];z<=cache.upper[2];++z) {
          const V phys=origin+V{T(x)*dx,T(y)*dx,T(z)*dx};
          const V body=olb::util::executeRotation<T,3,true>(phys,inverseRotation,center);
          const T distance=indicator->signedDistance(body);
          T phi=0;olb::sdf::evalSolidVolumeFraction(&phi,distance,epsilon);
          cache.fraction[cache.index(x,y,z)]=phi;
        }
    _statistics.sdfEvaluations+=count;
    V velocity=olb::particles::access::getVelocity(particle);velocity[0]+=velocityJump;
    const V angular=olb::particles::access::getAngularVelocity(particle);
    const T velocityScale=_converter.getPhysDeltaT()/dx;
    for(int x=mapLower[0];x<=mapUpper[0];++x)
      for(int y=mapLower[1];y<=mapUpper[1];++y)
        for(int z=mapLower[2];z<=mapUpper[2];++z) {
          const T phi=cache.phi(x,y,z);if(olb::util::nearZero(phi))continue;
          const L location{x,y,z};if(bg.getMaterial(location)!=1)continue;
          const V lever=origin+V{T(x)*dx,T(y)*dx,T(z)*dx}-center;
          const V latticeVelocity=(velocity+cross(angular,lever))*velocityScale;
          auto cell=bl.get(location);
          cell.template getFieldPointer<VELOCITY_NUMERATOR>()+=phi*latticeVelocity;
          const T denominator=cell.template getField<VELOCITY_DENOMINATOR>();
          cell.template setField<VELOCITY_DENOMINATOR>(denominator+phi);
          const T porosity=cell.template getField<POROSITY>();
          cell.template setField<POROSITY>(porosity*(T(1)-phi));
          cache.mapped.push_back(location);++_statistics.mappedCells;
          bool core=true;
          for(int k=0;k<3;++k)core=core&&location[k]>=0&&location[k]<blockSize[k];
          if(!core)continue;
          _localVolumes[id]+=phi*dx*dx*dx;
          bool forceRegion=true;
          for(int k=0;k<3;++k)forceRegion=forceRegion&&location[k]>=forceLower[k]&&location[k]<=forceUpper[k];
          if(!forceRegion)continue;
          std::uint32_t mask=0;
          for(int pop=1;pop<D::q;++pop) {
            const auto c=olb::descriptors::c<D>(pop);
            if(phi==T(1)&&cache.phi(x+c[0],y+c[1],z+c[2])==T(1))continue;
            mask|=std::uint32_t(1)<<pop;++_statistics.boundaryLinks;
          }
          if(mask){cache.surface.push_back({location,mask});++_statistics.boundaryCells;}
        }
  }
public:
  CachedCoupling(olb::SuperGeometry<T,3>&geometry,olb::SuperLattice<T,D>&lattice,
                 const olb::UnitConverter<T,D>&converter,const V&box,T gamma,
                 const V&semiaxes)
    :_geometry(geometry),_lattice(lattice),_converter(converter),
     _box(box),_axes(semiaxes),_gamma(gamma) {
    for(int k=0;k<3;++k)if(!(_box[k]>0&&_axes[k]>0))
      throw std::runtime_error("Invalid CachedCoupling box or ellipsoid axes");
    for(int b=0;b<lattice.getLoadBalancer().size();++b)
      if(geometry.getBlockGeometry(b).getPadding()<1)
        throw std::runtime_error("Momentum exchange requires one overlap cell");
  }
  void mapParticles(PS&particles,T time) {
    resetMappedFields();_activeImages=0;_statistics={};
    _localVolumes.assign(particles.size(),T(0));
    T shift=std::fmod(_gamma*_box[1]*time,_box[0]);if(shift<0)shift+=_box[0];
    for(std::size_t id=0;id<particles.size();++id) {
      auto particle=particles.get(id);if(!olb::particles::access::isValid(particle))continue;
      const V position=olb::particles::access::getPosition(particle);
      const V extent=halfExtent(particle);
      for(int b=0;b<_lattice.getLoadBalancer().size();++b) {
        const auto&bg=_geometry.getBlockGeometry(b);const T dx=bg.getDeltaR();
        const int padding=bg.getPadding();const V lower=bg.getOrigin()-V(T(padding)*dx);
        const auto n=bg.getExtent();
        const V upper=bg.getOrigin()+V{T(n[0]-1+padding)*dx,T(n[1]-1+padding)*dx,T(n[2]-1+padding)*dx};
        const int y0=int(std::ceil((lower[1]-extent[1]-position[1])/_box[1]));
        const int y1=int(std::floor((upper[1]+extent[1]-position[1])/_box[1]));
        const int z0=int(std::ceil((lower[2]-extent[2]-position[2])/_box[2]));
        const int z1=int(std::floor((upper[2]+extent[2]-position[2])/_box[2]));
        for(int ny=y0;ny<=y1;++ny) {
          const T shiftedX=position[0]+T(ny)*shift;
          const int x0=int(std::ceil((lower[0]-extent[0]-shiftedX)/_box[0]));
          const int x1=int(std::floor((upper[0]+extent[0]-shiftedX)/_box[0]));
          for(int nx=x0;nx<=x1;++nx)for(int nz=z0;nz<=z1;++nz)
            buildImage(particle,id,b,V{shiftedX+T(nx)*_box[0],position[1]+T(ny)*_box[1],
                       position[2]+T(nz)*_box[2]},extent,T(ny)*_gamma*_box[1]);
        }
      }
    }
    _statistics.activeImages=_activeImages;_mapTime=time;
  }
  // Caller first communicates populations, then refreshes outgoing LE force
  // ghost populations. No ordinary-periodic communication occurs inside here.
  void coupleFluidToParticles(PS&particles,T time) {
    if(!std::isfinite(_mapTime)||time!=_mapTime)
      throw std::runtime_error("Map particle fields at the force-evaluation time first");
    if(particles.size()>std::size_t(std::numeric_limits<int>::max())/stride)
      throw std::runtime_error("Too many particles for one MPI coupling reduction");
    _local.assign(stride*particles.size(),T(0));_global.resize(_local.size());
    for(std::size_t id=0;id<particles.size();++id)_local[stride*id+12]=_localVolumes[id];
    const T forceScale=_converter.getConversionFactorForce();
    const T velocityScale=_converter.getPhysDeltaT()/_converter.getPhysDeltaX();
    for(std::size_t k=0;k<_activeImages;++k) {
      const auto&cache=_images[k];auto particle=particles.get(cache.particle);
      auto&bl=_lattice.getBlock(cache.block);const auto&bg=_geometry.getBlockGeometry(cache.block);
      V velocity=olb::particles::access::getVelocity(particle);velocity[0]+=cache.velocityJump;
      const V angular=olb::particles::access::getAngularVelocity(particle);
      T*result=_local.data()+stride*cache.particle;
      for(const auto&boundary:cache.surface) {
        const V lever=bg.getPhysR(boundary.location)-cache.center;
        const V pVelocity=(velocity+cross(angular,lever))*velocityScale;
        auto inner=bl.get(boundary.location);V force(T(0));
        std::uint32_t mask=boundary.mask;
        while(mask) {
          const int pop=__builtin_ctz(mask);mask&=mask-1;
          const auto c=olb::descriptors::c<D>(pop);
          const T f1=bl.get(boundary.location+c)[pop];
          const T f2=inner[olb::descriptors::opposite<D>(pop)];
          for(int d=0;d<3;++d)
            force[d]-=forceScale*(f1*(T(c[d])-pVelocity[d])+f2*(T(c[d])+pVelocity[d]));
        }
        const V torque=cross(lever,force);
        for(int d=0;d<3;++d){result[d]+=force[d];result[3+d]+=torque[d];}
        addMoment(result+6,lever,force);
      }
    }
#ifdef PARALLEL_MODE_MPI
    if(!_local.empty())olb::singleton::mpi().allreduce(_local.data(),_global.data(),
        static_cast<int>(_local.size()),MPI_SUM);
#else
    _global=_local;
#endif
    _hydro.resize(particles.size());_stressletSum.fill(T(0));
    for(std::size_t id=0;id<particles.size();++id) {
      auto particle=particles.get(id);auto&h=_hydro[id];const T*r=_global.data()+stride*id;
      for(int d=0;d<3;++d){h.force[d]=r[d];h.torque[d]=r[3+d];}
      for(int d=0;d<6;++d){h.stresslet[d]=r[6+d];_stressletSum[d]+=h.stresslet[d];}
      h.sampledVolume=r[12];
      particle.template setField<olb::descriptors::FORCING,olb::descriptors::FORCE>(V{h.force[0],h.force[1],h.force[2]});
      particle.template setField<olb::descriptors::FORCING,olb::descriptors::TORQUE>(V{h.torque[0],h.torque[1],h.torque[2]});
    }
  }
  const std::vector<Hydrodynamics>&particleHydrodynamics()const{return _hydro;}
  const std::array<T,6>&stressletSum()const{return _stressletSum;}
  const Statistics&statistics()const{return _statistics;}
};
} // namespace graphiteCoupling

} } // SLURRY SCOPE END
#endif
