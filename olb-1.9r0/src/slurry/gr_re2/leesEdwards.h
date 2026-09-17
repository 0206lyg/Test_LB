#ifndef SLURRY_GR_RE2_GRAPHITE_LEES_EDWARDS_H
#define SLURRY_GR_RE2_GRAPHITE_LEES_EDWARDS_H

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

// SLURRY SCOPE BEGIN
namespace slurry { namespace gr_re2 {

// App-local OpenLB 1.9 sliding periodic populations. The surrounding geometry
// remains periodic in x,y,z. This class replaces the ordinary y-wrap halos
// AFTER their normal communication. No OpenLB source modification is required.
//
// Image convention: r_image = r + (nx*Lx + ny*s(t), ny*Ly, nz*Lz),
//                   u_image = u + (ny*gamma*Ly, 0, 0).
// Only the five D3Q19 populations actually required on each y plane are sent.
// Nonequilibrium populations are retained; an equilibrium difference provides
// the Galilean velocity jump. Fractional x displacement uses conservative
// periodic linear interpolation. The full equilibrium difference implements
// the rho/momentum/second-moment boost of Adhikari et al. (2005), Eqs. 31--34,
// https://arxiv.org/abs/cond-mat/0503175, with higher kinetic modes unchanged.
namespace graphiteCoupling {

template <typename T, typename DESCRIPTOR>
class LeesEdwardsBoundary {
  static_assert(DESCRIPTOR::d == 3 && DESCRIPTOR::q == 19,
                "This boundary is implemented for D3Q19");

  struct PlaneCell {
    int block, x, y, z, gx, gz, side;
  };

  olb::SuperLattice<T,DESCRIPTOR>& lattice_;
  T dx_, dt_, shearRate_, ly_, lx_;
  std::array<int,3> n_;
  std::array<std::array<int,5>,2> populations_{};
  std::vector<PlaneCell> sources_, ghosts_;
  std::vector<T> plane_;
  bool participant_ = false;
#ifdef PARALLEL_MODE_MPI
  MPI_Comm planeCommunicator_ = MPI_COMM_NULL;
#endif

  static int wrap(int x, int n) {
    const int r = x % n;
    return r < 0 ? r+n : r;
  }

  std::size_t index(int side, int x, int z, int p) const {
    return (((std::size_t(side)*n_[0]+x)*n_[2]+z)*5+p);
  }

  void remap(T time, bool outgoing) {
    if (!participant_) return;
    std::fill(plane_.begin(), plane_.end(), T{});
    const T velocityJump = shearRate_*ly_*dt_/dx_;

    // A bottom donor (side=0) supplies the top image with +velocityJump;
    // an upper donor supplies the bottom image with -velocityJump.
    for (const auto& source : sources_) {
      auto cell = lattice_.getBlock(source.block).get(source.x,source.y,source.z);
      T rho = T{1};
      T momentum[3] = {T{},T{},T{}};
      for (int p=0; p<DESCRIPTOR::q; ++p) {
        const T f = cell[p];
        rho += f;
        for (int d=0; d<3; ++d) {
          momentum[d] += T(olb::descriptors::c<DESCRIPTOR>(p,d))*f;
        }
      }
      const T u[3] = {momentum[0]/rho,momentum[1]/rho,momentum[2]/rho};
      const T du = source.side == 0 ? velocityJump : -velocityJump;
      // populations_[0] has cy=-1, populations_[1] has cy=+1.
      const auto& pop = populations_[outgoing ? 1-source.side : source.side];
      for (int k=0; k<5; ++k) {
        const int p = pop[k];
        plane_[index(source.side,source.gx,source.gz,k)] =
          cell[p] + equilibriumBoost(p,rho,u,du);
      }
    }
#ifdef PARALLEL_MODE_MPI
    olb::singleton::mpi().allReduceVect(plane_, MPI_SUM, planeCommunicator_);
#endif

    const T shift = phase(time)/dx_;
    for (const auto& ghost : ghosts_) {
      // side=0 is lower ghost (takes upper donor at x+s); side=1 is upper
      // ghost (takes lower donor at x-s). Include the x,z corner halos.
      const T donorX = T(ghost.gx) + (ghost.side == 0 ? shift : -shift);
      const T floorX = std::floor(donorX);
      const int x0 = wrap(static_cast<int>(floorX),n_[0]);
      const int x1 = x0+1 == n_[0] ? 0 : x0+1;
      const T fraction = donorX-floorX;
      const int donorSide = 1-ghost.side;
      const auto& pop = populations_[outgoing ? ghost.side : donorSide];
      auto cell = lattice_.getBlock(ghost.block).get(ghost.x,ghost.y,ghost.z);
      for (int k=0; k<5; ++k) {
        cell[pop[k]] =
          (T{1}-fraction)*plane_[index(donorSide,x0,ghost.gz,k)]
          + fraction*plane_[index(donorSide,x1,ghost.gz,k)];
      }
    }
  }

public:
  // lengths are the periodic cell lengths, NOT the final node coordinates:
  // nodes span [origin, origin+length-dx] with length/dx integer.
  LeesEdwardsBoundary(olb::SuperLattice<T,DESCRIPTOR>& lattice,
                     T dx, T dt, T shearRate,
                     const std::array<T,3>& lengths,
                     const std::array<T,3>& origin = {T{},T{},T{}})
    : lattice_(lattice), dx_(dx), dt_(dt), shearRate_(shearRate),
      ly_(lengths[1]), lx_(lengths[0]) {
    if (!(dx>0 && dt>0)) throw std::invalid_argument("LE requires positive dx and dt");
    for (int d=0; d<3; ++d) {
      n_[d] = static_cast<int>(std::llround(lengths[d]/dx));
      if (n_[d]<2 || std::abs(T(n_[d])*dx-lengths[d])>dx*T(1e-6)) {
        throw std::invalid_argument("LE lengths must contain an integer number of grid intervals");
      }
    }
    for (int sign=0; sign<2; ++sign) {
      int count=0;
      for (int p=0; p<DESCRIPTOR::q; ++p) {
        if (olb::descriptors::c<DESCRIPTOR>(p,1)==(sign==0 ? -1 : 1)) {
          if (count>=5) throw std::logic_error("Unexpected D3Q19 crossing population count");
          populations_[sign][count++]=p;
        }
      }
      if (count!=5) throw std::logic_error("Unexpected D3Q19 crossing population count");
    }

    auto& load = lattice_.getLoadBalancer();
    auto& decomposition = lattice_.getCuboidDecomposition();
    for (int block=0; block<load.size(); ++block) {
      const auto& cuboid = decomposition.get(load.glob(block));
      const auto extent = cuboid.getExtent();
      const auto blockOrigin = cuboid.getOrigin();
      std::array<int,3> offset;
      for (int d=0; d<3; ++d) {
        offset[d]=static_cast<int>(std::llround((blockOrigin[d]-origin[d])/dx));
        if (offset[d]<0 || offset[d]+extent[d]>n_[d]) {
          throw std::invalid_argument("LE cuboids do not match the specified periodic cell");
        }
      }
      if (lattice_.getBlock(block).getPadding()<1) {
        throw std::invalid_argument("LE requires at least one halo layer");
      }
      for (int side=0; side<2; ++side) {
        const bool touches = side==0 ? offset[1]==0 : offset[1]+extent[1]==n_[1];
        if (!touches) continue;
        participant_=true;
        const int sourceY = side==0 ? 0 : extent[1]-1;
        const int ghostY = side==0 ? -1 : extent[1];
        for (int x=0; x<extent[0]; ++x) for (int z=0; z<extent[2]; ++z) {
          sources_.push_back({block,x,sourceY,z,offset[0]+x,offset[2]+z,side});
        }
        for (int x=-1; x<=extent[0]; ++x) for (int z=-1; z<=extent[2]; ++z) {
          ghosts_.push_back({block,x,ghostY,z,offset[0]+x,wrap(offset[2]+z,n_[2]),side});
        }
      }
    }
#ifdef PARALLEL_MODE_MPI
    MPI_Comm_split(MPI_COMM_WORLD,participant_ ? 0 : MPI_UNDEFINED,
                   olb::singleton::mpi().getRank(),&planeCommunicator_);
#endif
    if (participant_) {
      const std::size_t size=std::size_t(2)*n_[0]*n_[2]*5;
      if (size>std::size_t(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("LE plane exceeds MPI count range");
      }
      plane_.resize(size);
      // One-time ownership check catches duplicated endpoint planes and a
      // geometry constructed with L/dx+1 nodes, before modifying populations.
      std::vector<int> coverage(std::size_t(2)*n_[0]*n_[2],0);
      for (const auto& s : sources_) {
        ++coverage[(std::size_t(s.side)*n_[0]+s.gx)*n_[2]+s.gz];
      }
#ifdef PARALLEL_MODE_MPI
      olb::singleton::mpi().allReduceVect(coverage,MPI_SUM,planeCommunicator_);
#endif
      if (std::any_of(coverage.begin(),coverage.end(),[](int x){return x!=1;})) {
        throw std::invalid_argument("LE boundary planes must have exactly one owner per node");
      }
    }
  }

  LeesEdwardsBoundary(const LeesEdwardsBoundary&) = delete;
  LeesEdwardsBoundary& operator=(const LeesEdwardsBoundary&) = delete;

  ~LeesEdwardsBoundary() {
#ifdef PARALLEL_MODE_MPI
    int finalized=0;
    MPI_Finalized(&finalized);
    if (!finalized && planeCommunicator_!=MPI_COMM_NULL) MPI_Comm_free(&planeCommunicator_);
#endif
  }

  T phase(T time) const {
    T shift=std::fmod(shearRate_*ly_*time,lx_);
    return shift<T{} ? shift+lx_ : shift;
  }

  // Call between lattice.collide() and lattice.AndStream(). collide() already
  // performed ordinary propagation-halo communication; time is the SAME
  // time as the donor distributions/mask (t_n, before streaming). Do not
  // communicate again between this call and AndStream().
  void apply(T time) { remap(time,false); }

  // Call AFTER lattice.communicate() and BEFORE surface momentum exchange.
  // Its ghost population directions are the opposite set from apply():
  // OpenLB ME reads neighbor[iPop] where c_i points out of the core cell.
  void refreshForceGhosts(T time) { remap(time,true); }

  std::size_t communicatedValues() const { return plane_.size(); }

  static T equilibriumBoost(int p, T rho, const T u[3], T du) {
    const T csInv=olb::descriptors::invCs2<T,DESCRIPTOR>();
    T cu=T{};
    for (int d=0; d<3; ++d) cu += T(olb::descriptors::c<DESCRIPTOR>(p,d))*u[d];
    const T cdu=T(olb::descriptors::c<DESCRIPTOR>(p,0))*du;
    return olb::descriptors::t<T,DESCRIPTOR>(p)*rho*
      (csInv*cdu + T(0.5)*csInv*csInv*(T{2}*cu*cdu+cdu*cdu)
       -T(0.5)*csInv*(T{2}*u[0]*du+du*du));
  }
};

} // namespace graphiteCoupling

} } // SLURRY SCOPE END
#endif
