/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Florian Kaiser, Stephan Simonis
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::descriptors;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierCauchy,
  Lattice<double,
    descriptors::D2Q8<
      tag::MRT>>
>;

namespace olb::parameters {
  struct KAPPA  : public descriptors::FIELD_BASE<1> { };
}

// TO DO
class Ellipse2D {
  Vector<MyCase::value_t,2> center1;
  Vector<MyCase::value_t,2> center2;
  Vector<MyCase::value_t,2> center3;

  public:
    Ellipse2D(const MyCase::value_t charL) :
      center1(Vector<MyCase::value_t,2>(.75 * charL, .75 * charL)),
      center2(Vector<MyCase::value_t,2>(.95 * charL, .6 * charL)),
      center3(Vector<MyCase::value_t,2>(.5 * charL, .9 * charL)),
      ellipse1(IndicatorEllipse2D<MyCase::value_t>(center1, .693 * charL, .548 * charL, -20)),
      ellipse2(IndicatorEllipse2D<MyCase::value_t>(center2, .1 * charL, .134 * charL)),
      ellipse3(IndicatorEllipse2D<MyCase::value_t>(center3, .187 * charL, .1 * charL))
    {};

    IndicatorEllipse2D<MyCase::value_t> ellipse1;
    IndicatorEllipse2D<MyCase::value_t> ellipse2;
    IndicatorEllipse2D<MyCase::value_t> ellipse3;
};

template <typename T>
class ForceField2D : public AnalyticalF2D<T, T> {
  protected:
    MyCase& aCase;

  public:
    ForceField2D(MyCase& myCase)
      : AnalyticalF2D<T, T>(2),
        aCase(myCase) {};

    bool operator()(T output[], const T input[]) override
    {
      auto& lattice = aCase.getLattice(NavierCauchy{});
      T pi = std::numbers::pi_v<double>;
      using DESCRIPTOR = MyCase::descriptor_t_of<NavierCauchy>;

      const auto& converter = lattice.getUnitConverter();
      T dx = converter.getPhysDeltaX();
      T dt = converter.getPhysDeltaT();
      T kappa = converter.getDampingFactor();
      T lambda = converter.getLatticeLambda();
      T mu = converter.getLatticeShearModulus();
      T charU = converter.getCharPhysDisplacement();

      T x = input[0];
      T y = input[1];

      T omega_11 = 1. / (mu / descriptors::invCs2<T,DESCRIPTOR>() + 0.5);
      T omega_d  = 1. / (2 * mu / (1 - descriptors::invCs2<T,DESCRIPTOR>()) + 0.5);
      T omega_s  = 1. / (2 * (mu + lambda) / (1 + descriptors::invCs2<T,DESCRIPTOR>()) + 0.5);

      T tau_11 = 1. / omega_11 - 0.5;
      T tau_s = 1. / omega_d - 0.5;
      T tau_p = 1. / omega_s - 0.5;
      T tau_12 = 0.5;

      T d1 = -1. / 4. - 1. / 2. * descriptors::invCs2<T,DESCRIPTOR>() * tau_12 * tau_s + pow(tau_s, 2) / 2.
          -  1. / 2. * descriptors::invCs2<T,DESCRIPTOR>() * pow(tau_s, 2) + 1. / 2. * descriptors::invCs2<T,DESCRIPTOR>() * tau_12 * tau_p
          +  pow(tau_p, 2) / 2. + 1. / 2. * descriptors::invCs2<T,DESCRIPTOR>() * pow(tau_p, 2);
      T d2 = -descriptors::invCs2<T,DESCRIPTOR>() / 2. + descriptors::invCs2<T,DESCRIPTOR>() * pow(tau_11, 2) + descriptors::invCs2<T,DESCRIPTOR>() * tau_11 * tau_12
          + 1. / 2. * descriptors::invCs2<T,DESCRIPTOR>() * tau_12 * tau_s - pow(tau_s, 2) / 2.
          + 1. / 2. * descriptors::invCs2<T,DESCRIPTOR>() * pow(tau_s, 2) + 1. / 2. * descriptors::invCs2<T,DESCRIPTOR>() * tau_12 * tau_p
          + pow(tau_p, 2) / 2. + 1. / 2. * descriptors::invCs2<T,DESCRIPTOR>() * pow(tau_p, 2);
      T d3 = -descriptors::invCs2<T,DESCRIPTOR>() / 4.
          +  descriptors::invCs2<T,DESCRIPTOR>() * pow(tau_11, 2)
          +  descriptors::invCs2<T,DESCRIPTOR>() * tau_11 * tau_12;

      T mu_phys     =     mu * (dx * dx * kappa) / dt;
      T lambda_phys = lambda * (dx * dx * kappa) / dt;


      output[0] = bx(0, x, y, dx, mu_phys, lambda_phys, charU, pi) * dt / kappa
                + (d1 * d2bx_dx2(0, x, y, dx, mu_phys, lambda_phys, charU, pi)
                + d2 * d2by_dxdy(0, x, y, dx, mu_phys, lambda_phys, charU, pi)
                + d3 * d2bx_dy2(0, x, y, dx, mu_phys, lambda_phys, charU, pi)) * dx * dx * dt / kappa;

      output[1] = by(0, x, y, dx, mu_phys, lambda_phys, charU, pi) * dt / kappa
                + (d1 * d2by_dy2(0, x, y, dx, mu_phys, lambda_phys, charU, pi)
                + d2 * d2bx_dxdy(0, x, y, dx, mu_phys, lambda_phys, charU, pi)
                + d3 * d2by_dx2(0, x, y, dx, mu_phys, lambda_phys, charU, pi)) * dx * dx * dt / kappa;

      return true;
    };

    T bx(T t, T x, T y, T dx, T mu, T lambda, T charU, T pi)
    {
      return   0.0036 * mu * pi * pi * x * x * util::sin(2 * pi * x * y)
            + 0.0036 * pi * pi * y * y * (lambda + 2 * mu) * util::sin(2 * pi * x * y)
            + 0.0028 * pi * x * (lambda + mu) * util::sin(2 * pi * y);
    }

    T d2bx_dx2(T t, T x, T y, T dx, T mu, T lambda, T charU, T pi)
    {
      return   pi * pi * (-0.0144 * mu * pi * pi * x * x * y * y * util::sin(2 * pi * x * y)
            + 0.0288 * mu * pi * x * y * util::cos(2 * pi * x * y)
            + 0.0072 * mu * util::sin(2 * pi * x * y)
            - 0.0144 * pi * pi * y * y * y * y * (lambda + 2 * mu) * util::sin(2 * pi * x * y));
    }

    T d2by_dxdy(T t, T x, T y, T dx, T mu, T lambda, T charU, T pi)
    {
      return    pi * pi * (-0.0112 * pi * x * (lambda + 2 * mu) * util::sin(2 * pi * y)
            + (lambda + mu) * (-0.0144 * pi * pi * x * x * y * y * util::sin(2 * pi * x * y)
            + 0.0288 * pi * x * y * util::cos(2 * pi * x * y)
            + 0.0072 * util::sin(2 * pi * x * y)));
    }

    T d2bx_dy2(T t, T x, T y, T dx, T mu, T lambda, T charU, T pi)
    {
      return   pi * pi * (-0.0144 * mu * pi * pi * x * x * x * x * util::sin(2 * pi * x * y)
            - 0.0144 * pi * pi * x * x * y * y * (lambda + 2 * mu) * util::sin(2 * pi * x * y)
            + 0.0288 * pi * x * y * (lambda + 2 * mu) * util::cos(2 * pi * x * y)
            - 0.0112 * pi * x * (lambda + mu) * util::sin(2 * pi * y)
            + 0.0072 * (lambda + 2 * mu) * util::sin(2 * pi * x * y));
    }

    T by(T t, T x, T y, T dx, T mu, T lambda, T charU, T pi)
    {
      return - 0.0014 * mu * cos(2 * pi * y)
            + 0.0028 * pi * pi * (lambda + 2 * mu) * (x * x + 1) * cos(2 * pi * y)
            - pi * (lambda + mu) * (-0.0036 * pi * x * y * sin(2 * pi * x * y)
            + 0.0018 * cos(2 * pi * x * y));
    }

    T d2by_dy2(T t, T x, T y, T dx, T mu, T lambda, T charU, T pi)
    {
      return   pi*pi*(0.0056*mu*cos(2*pi*y)
            - 0.0112*pi*pi*(lambda + 2*mu)*(x*x + 1)*cos(2*pi*y)
            - pi*x*x*(lambda + mu)*(0.0144*pi*x*y*sin(2*pi*x*y)
            - 0.0216*cos(2*pi*x*y)));
    }

    T d2bx_dxdy(T t, T x, T y, T dx, T mu, T lambda, T charU, T pi)
    {
      return   pi*pi*(-0.0144*mu*pi*pi*x*x*x*y*sin(2*pi*x*y)
            + 0.0216*mu*pi*x*x*cos(2*pi*x*y)
            - 0.0144*pi*pi*x*y*y*y*(lambda + 2*mu)*sin(2*pi*x*y)
            + 0.0216*pi*y*y*(lambda + 2*mu)*cos(2*pi*x*y)
            + 0.0056*(lambda + mu)*cos(2*pi*y));
    }

    T d2by_dx2(T t, T x, T y, T dx, T mu, T lambda, T charU, T pi)
    {
      return   pi*pi*(-pi*y*y*(lambda + mu)*(0.0144*pi*x*y*sin(2*pi*x*y)
            - 0.0216*cos(2*pi*x*y))
            + 0.0056*(lambda + 2*mu)*cos(2*pi*y));
    }
};

template <typename T>
class ManufacturedSolutionU2D : public AnalyticalF2D<T, T> {
  protected:

  public:
  ManufacturedSolutionU2D(MyCase& myCase) : AnalyticalF2D<T, T>(2) {};

  bool operator()(T output[], const T input[]) override
  {
    T pi = std::numbers::pi_v<double>;

    T x = input[0];
    T y = input[1];

    output[0] = 9. / 10000. * util::sin(2. * pi * x * y);
    output[1] = 7. / 10000. * util::cos(2. * pi * y) * (x * x + 1);

    return true;
  };
};

template <typename T>
class ManufacturedSolutionStress2D : public AnalyticalF2D<T, T> {
  protected:
    MyCase& aCase;

  public:
    ManufacturedSolutionStress2D(MyCase& myCase) : AnalyticalF2D<T, T>(4),
    aCase(myCase) {};

    bool operator()(T output[], const T input[]) override
    {
      T pi = std::numbers::pi_v<double>;

      const auto& converter = aCase.getLattice(NavierCauchy{}).getUnitConverter();
      T dx = converter.getPhysDeltaX();
      T dt = converter.getPhysDeltaT();
      T kappa = converter.getDampingFactor();
      T lambda = converter.getLatticeLambda();
      T mu = converter.getLatticeShearModulus();

      T x = input[0];
      T y = input[1];

      T latticeFactor = dt / (kappa * dx);

      // xx
      output[0] = latticeFactor * ((2. * mu + lambda) * dux_dx(0, x, y, pi)
                                            + lambda  * duy_dy(0, x, y, pi));
      // xy
      output[1] = latticeFactor * (mu * (  dux_dy(0, x, y, pi)
                                        + duy_dx(0, x, y, pi)  ));

      // yx
      output[2] = latticeFactor * (mu * (  dux_dy(0, x, y, pi)
                + duy_dx(0, x, y, pi)  ));

      // yy
      output[3] = latticeFactor * ((2. * mu + lambda) * duy_dy(0, x, y, pi)
                                            + lambda  * dux_dx(0, x, y, pi));

      return true;
    };

    T dux_dx(T t, T x, T y, T pi) {
      return 0.0018 * pi * y * util::cos(2. * pi * x * y);
    }

    T dux_dy(T t, T x, T y, T pi) {
      return 0.0018 * pi * x * util::cos(2 * pi * x * y);
    }

    T duy_dx(T t, T x, T y, T pi) {
      return 0.0014 * x * util::cos(2 * pi * y);
    }

    T duy_dy(T t, T x, T y, T pi) {
      return -0.0014 * pi * (x * x + 1.) * util::sin(2 * pi * y);
    }
};

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX = extent[0] / parameters.get<parameters::RESOLUTION>();
  const Vector origin{physDeltaX / 2., physDeltaX / 2.};
  IndicatorCuboid2D<T> cuboid(extent, origin);

  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({false, false});
  return mesh;
}

/// @brief Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers will be used to assign physics to lattice nodes
void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  geometry.rename( 0, 4 );

  auto physCharLength = params.get<parameters::PHYS_CHAR_LENGTH>();
  Ellipse2D Ellipses(physCharLength);
  auto& ellipse1 = Ellipses.ellipse1;
  auto& ellipse2 = Ellipses.ellipse2;
  auto& ellipse3 = Ellipses.ellipse3;

  geometry.rename( 4, 1, ellipse1 );
  geometry.rename( 1, 2, ellipse2 );
  geometry.rename( 1, 3, ellipse3 );

  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}
//new
void prepareLattice(MyCase& myCase) {
  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& lattice = myCase.getLattice(NavierCauchy{});

  using DESCRIPTOR = MyCase::descriptor_t_of<NavierCauchy>;

  auto extent = params.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX = extent[0] / params.get<parameters::RESOLUTION>();
  const T physDeltaT = physDeltaX * physDeltaX;
  const T charDisplacement = params.get<parameters::PHYS_CHAR_DISPLACEMENT>();
  const T ELattice = params.get<parameters::YOUNGS_MODULUS>();
  const T physPoissonRatio = params.get<parameters::POISSON_RATIO>();
  const T kappa = params.get<parameters::KAPPA>();
  const T charLength = params.get<parameters::PHYS_CHAR_LENGTH>();

  lattice.setUnitConverter<LinElaUnitConverter<T, DESCRIPTOR>>(
    (T) physDeltaX,
    (T) physDeltaT,
    (T) charLength,
    (T) charDisplacement,
    (T) ELattice * (physDeltaX * physDeltaX * kappa) / physDeltaT,
    (T) physPoissonRatio,
    (T) kappa
  );

  const auto& converter = lattice.getUnitConverter();
  converter.print();

  auto bulkIndicator = geometry.getMaterialIndicator(1);

  lattice.defineDynamics<BoolakeeLinearElasticityBoundary>( bulkIndicator );

  const T physCharLength = params.get<parameters::PHYS_CHAR_LENGTH>();
  Ellipse2D Ellipses(physCharLength);
  auto& ellipse1 = Ellipses.ellipse1;
  auto& ellipse2 = Ellipses.ellipse2;
  auto& ellipse3 = Ellipses.ellipse3;
  setBoolakeeNeumannBoundary<T,DESCRIPTOR,BoolakeeNeumannPostProcessor<T,DESCRIPTOR>>(lattice, geometry.getMaterialIndicator( 4 ), bulkIndicator, ellipse1, true);
  setBoolakeeDirichletBoundary<T,DESCRIPTOR,BoolakeeDirichletPostProcessor<T,DESCRIPTOR>>(lattice, geometry.getMaterialIndicator( 2 ), bulkIndicator, ellipse2);
  setBoolakeeDirichletBoundary<T,DESCRIPTOR,BoolakeeDirichletPostProcessor<T,DESCRIPTOR>>(lattice, geometry.getMaterialIndicator( 3 ), bulkIndicator, ellipse3);

  {
    auto& communicator = lattice.getCommunicator(stage::PostCollide());
    communicator.template requestField<descriptors::PREVIOUS_CELL>();
    communicator.template requestField<descriptors::SOLID_DISTANCE_FIELD>();
    communicator.template requestField<descriptors::BOUNDARY_COORDS_X>();
    communicator.template requestField<descriptors::BOUNDARY_COORDS_Y>();
    communicator.template requestField<descriptors::BARED_MOMENT_VECTOR>();
    communicator.requestOverlap(1);
    communicator.exchangeRequests();
  }

  std::vector<T>          iniPop = {0., 0., 0., 0., 0., 0., 0., 0.};
  std::vector<T>          iniDisp = {0., 0.};
  std::vector<T>          iniStress = {0., 0., 0.};

  AnalyticalConst2D<T, T> initialPopulationF(iniPop);
  AnalyticalConst2D<T, T> initialDispF(iniDisp);
  AnalyticalConst2D<T, T> initialStressF(iniStress);

  // dx, dt, theta, mü, lambda, kappa, uChar, epsilon
  T magic[8] = {converter.getConversionFactorLength(),
                converter.getConversionFactorTime(),
                descriptors::invCs2<T, DESCRIPTOR>(),
                converter.getLatticeShearModulus(),
                converter.getLatticeLambda(),
                converter.getDampingFactor(),
                converter.getCharPhysVelocity(),
                converter.getEpsilon()};

  lattice.setParameter<descriptors::MAGIC_SOLID>(magic);

  const T omega_11  = 1. / (      converter.getLatticeShearModulus()                                 /       invCs2<T,DESCRIPTOR>()  + 0.5);
  const T omega_d   = 1. / (2. *  converter.getLatticeShearModulus()                                 / (1. - invCs2<T,DESCRIPTOR>()) + 0.5);
  const T omega_s   = 1. / (2. * (converter.getLatticeShearModulus() + converter.getLatticeLambda()) / (1. + invCs2<T,DESCRIPTOR>()) + 0.5);

  const T tau_12 = 0.5;
  const T tau_22 = 0.5;

  const T tau_21 = tau_12;

  const T omega_12  = 1. / (tau_12 + 0.5);
  const T omega_21  = 1. / (tau_21 + 0.5);
  const T omega_22  = 1. / (tau_22 + 0.5);

  lattice.setParameter<descriptors::OMEGA_SOLID>({omega_11, omega_s, omega_d, omega_12, omega_21, omega_22});

  ForceField2D<T> force(myCase);
  lattice.defineField<FORCE>(geometry, 1, force);
  lattice.defineField<POPULATION>(geometry, 1, initialPopulationF);

  lattice.defineField<DISP_SOLID>(geometry, 1, initialDispF);
  lattice.defineField<SIGMA_SOLID>(geometry, 1, initialStressF);

  lattice.defineField<PREVIOUS_CELL>(geometry, 2, initialPopulationF);
  lattice.defineField<PREVIOUS_CELL>(geometry, 3, initialPopulationF);
  lattice.defineField<PREVIOUS_CELL>(geometry, 4, initialPopulationF);

  T neumann_constants[3] = {
    (2. * (1. - descriptors::invCs2<T,DESCRIPTOR>()) * (converter.getLatticeBulkModulus() - converter.getLatticeShearModulus())) / (descriptors::invCs2<T,DESCRIPTOR>() * (1. - descriptors::invCs2<T,DESCRIPTOR>() - 4. * converter.getLatticeShearModulus())),
     2. * converter.getLatticeShearModulus() / (descriptors::invCs2<T,DESCRIPTOR>() - 2. * converter.getLatticeShearModulus()),
     4. * converter.getLatticeShearModulus() / (1. - descriptors::invCs2<T,DESCRIPTOR>() - 4. * converter.getLatticeShearModulus())
  };

  lattice.setParameter<descriptors::NEUMANN_SOLID_C>(neumann_constants);

  /// Make the lattice ready for simulation
  lattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

/// Set initial condition for primal variables (velocity and density)
/// @param myCase The Case instance which keeps the simulation data
/// @note Be careful: initial values have to be set using lattice units
void setInitialValues(MyCase& myCase)
{
  // Nothing to do here, because simulation is initialized with 0
}

/// Update boundary values at times (and external fields, if they exist)
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
/// @note Be careful: boundary values have to be set using lattice units
void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{
  // Nothing to do here, because simulation does not depend on time
}

/// Compute simulation results at times
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  OstreamManager clout(std::cout, "getResults");
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierCauchy>;
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierCauchy{});
  auto& converter = lattice.getUnitConverter();
  auto& geometry = myCase.getGeometry();

  const std::size_t iTlog = parameters.get<parameters::IT_LOG_PSEUDO_TIME>();
  const std::size_t iTvtk = parameters.get<parameters::IT_VTK_PSEUDO_TIME>();
  const std::size_t iTMax = converter.getLatticeTime(parameters.get<parameters::MAX_PHYS_T>());

  SuperVTMwriter2D<T> vtmWriter("ellipseCombo");

  ManufacturedSolutionU2D<T> dispSol(myCase);
  SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR>      dispSolLattice(dispSol, lattice);

  ManufacturedSolutionStress2D<T> stressSol(myCase);
  SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR>      stressSolLattice(stressSol, lattice);

  // Fields for error calc
  SuperLatticeField2D<T, DESCRIPTOR, olb::descriptors::DISP_SOLID>     disp(lattice);
  SuperLatticeField2D<T, DESCRIPTOR, olb::descriptors::SIGMA_SOLID>    stress(lattice);

  if (iT == 0) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, DESCRIPTOR>   cuboid(lattice);
    SuperLatticeRank2D<T, DESCRIPTOR>     rank(lattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();
  }

  // Writes the VTK files
  if (iT%iTvtk == 0 && iT >= 0) {
    vtmWriter.addFunctor(disp,             "numerical disp");
    vtmWriter.addFunctor(dispSolLattice,   "analytical disp");
    vtmWriter.addFunctor(stress,           "numerical stress");
    vtmWriter.addFunctor(stressSolLattice, "analytical stress");

    lattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }

  /// Print some (numerical and computational) statistics
  if (iT%iTlog == 0 || iT % iTMax == 0) {
    T   l2UResult[2]        = {T(), T()};
    T   lInfUResult[2]      = {T(), T()};
    T   l2StressResult[2]   = {T(), T()};
    T   lInfStressResult[2] = {T(), T()};
    int tmp[]               = {int()};

    auto indicatorF = geometry.getMaterialIndicator(1);
    SuperRelativeErrorL2Norm2D<T>   relUErrorL2Norm(lattice, disp, dispSol, indicatorF);
    SuperRelativeErrorLinfNorm2D<T> relUErrorLinfNorm(lattice, disp, dispSol, indicatorF);
    SuperRelativeErrorL2Norm2D<T>   relStressErrorL2Norm( lattice, stress, stressSol, indicatorF );
    SuperRelativeErrorLinfNorm2D<T> relStressErrorLinfNorm( lattice, stress, stressSol, indicatorF );

    relUErrorL2Norm(l2UResult, tmp);
    relUErrorLinfNorm(lInfUResult, tmp);
    relStressErrorL2Norm(l2StressResult, tmp);
    relStressErrorLinfNorm(lInfStressResult, tmp);

    clout << "N\t" << "L2 U Error\t" << "LInf U Error\t" << "L2 Stress Error\t" << "LInf Stress Error" << std::endl;
    clout << converter.getResolution() << "\t" << l2UResult[0] << "\t" << lInfUResult[0] << "\t" << l2StressResult[0] << "\t\t" << lInfStressResult[0] << std::endl;
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.print(iT);
  }
}

/// @brief Execute simulation: set initial values and run time loop
/// @param myCase The Case instance which keeps the simulation data
void simulate(MyCase& myCase) {
  OstreamManager clout(std::cout, "simulate");
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  const T physMaxT = parameters.get<parameters::MAX_PHYS_T>();

  const std::size_t iTmax = myCase.getLattice(NavierCauchy{}).getUnitConverter().getLatticeTime(physMaxT);
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    timer.update(iT);
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    myCase.getLattice(NavierCauchy{}).collideAndStream();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT);
  }

  timer.stop();
  timer.printSummary();
}

/// Setup and run a simulation
int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<KAPPA                    >(  1.          );
    myCaseParameters.set<RESOLUTION               >( 40           );
    myCaseParameters.set<YOUNGS_MODULUS           >(  0.11        );
    myCaseParameters.set<parameters::POISSON_RATIO>(  0.7         );
    myCaseParameters.set<PHYS_CHAR_DISPLACEMENT   >(  1.0         );
    myCaseParameters.set<PHYS_CHAR_LENGTH         >(  1.0         );
    myCaseParameters.set<DOMAIN_EXTENT            >( { 1.5, 1.5 } );
    myCaseParameters.set<MAX_PHYS_T               >( 30           );
    myCaseParameters.set<IT_LOG_PSEUDO_TIME       >(100           );
    myCaseParameters.set<IT_VTK_PSEUDO_TIME       >(100           );
  }
  myCaseParameters.fromCLI(argc, argv);

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);
}
