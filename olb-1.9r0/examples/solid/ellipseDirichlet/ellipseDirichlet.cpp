/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Stephan Simonis, Florian Kaiser
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

#include "olb.h"

using namespace olb;
using namespace olb::names;
using namespace olb::descriptors;

using MyCase = Case<
  NavierCauchy, Lattice<double, D2Q8<tag::MRT>>
>;

namespace olb::parameters {
  struct KAPPA : public descriptors::FIELD_BASE<1> { };
}

class ellipse2D {
  Vector<MyCase::value_t,2> center1;
  Vector<MyCase::value_t,2> center2;
  Vector<MyCase::value_t,2> center3;

  public:
    ellipse2D(const MyCase::value_t charL) :
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

template<typename T, typename DESCRIPTOR>
class ForceField2D : public AnalyticalF2D<T, T> {
  protected:
    MyCase& myCase;
    const T pi = std::numbers::pi_v<T>;

  public:
    ForceField2D(MyCase& myCase) : AnalyticalF2D<T, T>(2), myCase(myCase) {};

    bool operator()(T output[], const T input[]) override
    {
      constexpr T theta = invCs2<T, DESCRIPTOR>();
      auto& lattice = myCase.getLattice(NavierCauchy{});
      const auto& converter = lattice.getUnitConverter();
      T dx = converter.getPhysDeltaX();
      T dt = converter.getPhysDeltaT();
      T kappa = converter.getDampingFactor();
      T lambda = converter.getLatticeLambda();
      T mu = converter.getLatticeShearModulus();
      T charU = converter.getCharPhysDisplacement();

      T x = input[0];
      T y = input[1];

      T omega_11 = 1. / (mu / theta + 0.5);
      T omega_d  = 1. / (2 * mu / (1 - theta) + 0.5);
      T omega_s  = 1. / (2 * (mu + lambda) / (1 + theta) + 0.5);

      T tau_11 = 1. / omega_11 - 0.5;
      T tau_s = 1. / omega_d - 0.5;
      T tau_p = 1. / omega_s - 0.5;
      T tau_12 = 0.5;

      // These are correct
      T d1 = -1. / 4. - 1. / 2. * theta * tau_12 * tau_s + pow(tau_s, 2) / 2.
          -  1. / 2. * theta * pow(tau_s, 2) + 1. / 2. * theta * tau_12 * tau_p
          +  pow(tau_p, 2) / 2. + 1. / 2. * theta * pow(tau_p, 2);
      T d2 = -theta / 2. + theta * pow(tau_11, 2) + theta * tau_11 * tau_12
          + 1. / 2. * theta * tau_12 * tau_s - pow(tau_s, 2) / 2.
          + 1. / 2. * theta * pow(tau_s, 2) + 1. / 2. * theta * tau_12 * tau_p
          + pow(tau_p, 2) / 2. + 1. / 2. * theta * pow(tau_p, 2);
      T d3 = -theta / 4.
          +  theta * pow(tau_11, 2)
          +  theta * tau_11 * tau_12;

      T mu_phys     =     mu * (dx * dx * kappa) / dt;
      T lambda_phys = lambda * (dx * dx * kappa) / dt;


      output[0] = bx(0, x, y, dx, mu_phys, lambda_phys, charU) * dt / kappa
                + (d1 * d2bx_dx2(0, x, y, dx, mu_phys, lambda_phys, charU)
                + d2 * d2by_dxdy(0, x, y, dx, mu_phys, lambda_phys, charU)
                + d3 * d2bx_dy2(0, x, y, dx, mu_phys, lambda_phys, charU)) * dx * dx * dt / kappa;

      output[1] = by(0, x, y, dx, mu_phys, lambda_phys, charU) * dt / kappa
                + (d1 * d2by_dy2(0, x, y, dx, mu_phys, lambda_phys, charU)
                + d2 * d2bx_dxdy(0, x, y, dx, mu_phys, lambda_phys, charU)
                + d3 * d2by_dx2(0, x, y, dx, mu_phys, lambda_phys, charU)) * dx * dx * dt / kappa;

      return true;
    };

    T bx(T t, T x, T y, T dx, T mu, T lambda, T charU)
    {
      return   0.0036 * mu * pi * pi * x * x * util::sin(2 * pi * x * y)
            + 0.0036 * pi * pi * y * y * (lambda + 2 * mu) * util::sin(2 * pi * x * y)
            + 0.0028 * pi * x * (lambda + mu) * util::sin(2 * pi * y);
    }

    T d2bx_dx2(T t, T x, T y, T dx, T mu, T lambda, T charU)
    {
      return   pi * pi * (-0.0144 * mu * pi * pi * x * x * y * y * util::sin(2 * pi * x * y)
            + 0.0288 * mu * pi * x * y * util::cos(2 * pi * x * y)
            + 0.0072 * mu * util::sin(2 * pi * x * y)
            - 0.0144 * pi * pi * y * y * y * y * (lambda + 2 * mu) * util::sin(2 * pi * x * y));
    }

    T d2by_dxdy(T t, T x, T y, T dx, T mu, T lambda, T charU)
    {
      return    pi * pi * (-0.0112 * pi * x * (lambda + 2 * mu) * util::sin(2 * pi * y)
            + (lambda + mu) * (-0.0144 * pi * pi * x * x * y * y * util::sin(2 * pi * x * y)
            + 0.0288 * pi * x * y * util::cos(2 * pi * x * y)
            + 0.0072 * util::sin(2 * pi * x * y)));
    }

    T d2bx_dy2(T t, T x, T y, T dx, T mu, T lambda, T charU)
    {
      return   pi * pi * (-0.0144 * mu * pi * pi * x * x * x * x * util::sin(2 * pi * x * y)
            - 0.0144 * pi * pi * x * x * y * y * (lambda + 2 * mu) * util::sin(2 * pi * x * y)
            + 0.0288 * pi * x * y * (lambda + 2 * mu) * util::cos(2 * pi * x * y)
            - 0.0112 * pi * x * (lambda + mu) * util::sin(2 * pi * y)
            + 0.0072 * (lambda + 2 * mu) * util::sin(2 * pi * x * y));
    }

    T by(T t, T x, T y, T dx, T mu, T lambda, T charU)
    {
      return - 0.0014 * mu * cos(2 * pi * y)
            + 0.0028 * pi * pi * (lambda + 2 * mu) * (x * x + 1) * cos(2 * pi * y)
            - pi * (lambda + mu) * (-0.0036 * pi * x * y * sin(2 * pi * x * y)
            + 0.0018 * cos(2 * pi * x * y));
    }

    T d2by_dy2(T t, T x, T y, T dx, T mu, T lambda, T charU)
    {
      return   pi*pi*(0.0056*mu*cos(2*pi*y)
            - 0.0112*pi*pi*(lambda + 2*mu)*(x*x + 1)*cos(2*pi*y)
            - pi*x*x*(lambda + mu)*(0.0144*pi*x*y*sin(2*pi*x*y)
            - 0.0216*cos(2*pi*x*y)));
    }

    T d2bx_dxdy(T t, T x, T y, T dx, T mu, T lambda, T charU)
    {
      return   pi*pi*(-0.0144*mu*pi*pi*x*x*x*y*sin(2*pi*x*y)
            + 0.0216*mu*pi*x*x*cos(2*pi*x*y)
            - 0.0144*pi*pi*x*y*y*y*(lambda + 2*mu)*sin(2*pi*x*y)
            + 0.0216*pi*y*y*(lambda + 2*mu)*cos(2*pi*x*y)
            + 0.0056*(lambda + mu)*cos(2*pi*y));
    }

    T d2by_dx2(T t, T x, T y, T dx, T mu, T lambda, T charU)
    {
      return   pi*pi*(-pi*y*y*(lambda + mu)*(0.0144*pi*x*y*sin(2*pi*x*y)
            - 0.0216*cos(2*pi*x*y))
            + 0.0056*(lambda + 2*mu)*cos(2*pi*y));
    }
};

template<typename T>
class ManufacturedSolutionU2D : public AnalyticalF2D<T, T> {

  protected:
    MyCase& myCase;

  public:
    ManufacturedSolutionU2D(MyCase& _myCase) : AnalyticalF2D<T, T>(2), myCase(_myCase){};

    bool operator()(T output[], const T input[]) override
    {
      constexpr T pi = std::numbers::pi_v<T>;
      T x = input[0];
      T y = input[1];

      // lattice units
      output[0] = 9. / 10000. * util::sin(2. * pi * x * y);
      output[1] = 7. / 10000. * util::cos(2. * pi * y) * (x * x + 1);

      return true;
    };
};

template<typename T>
class ManufacturedSolutionStress2D : public AnalyticalF2D<T, T> {
  protected:
    MyCase& myCase;
    const T pi = std::numbers::pi_v<T>;

  public:
    ManufacturedSolutionStress2D(MyCase& _myCase) : AnalyticalF2D<T, T>(4), myCase(_myCase) {};

    bool operator()(T output[], const T input[]) override
    {
      auto& lattice = myCase.getLattice(NavierCauchy{});
      const auto& converter = lattice.getUnitConverter();
      const T dx = converter.getPhysDeltaX();
      const T dt = converter.getPhysDeltaT();
      const T kappa = converter.getDampingFactor();
      const T lambda = converter.getLatticeLambda();
      const T mu = converter.getLatticeShearModulus();

      T x = input[0];
      T y = input[1];

      T latticeFactor = dt / (kappa * dx);

      // xx
      output[0] = latticeFactor * ((2. * mu + lambda) * dux_dx(0, x, y)
                                            + lambda  * duy_dy(0, x, y));
      // xy
      output[1] = latticeFactor * (mu * (  dux_dy(0, x, y)
                                        + duy_dx(0, x, y)  ));

      // yx
      output[2] = latticeFactor * (mu * (  dux_dy(0, x, y)
                + duy_dx(0, x, y)  ));

      // yy
      output[3] = latticeFactor * ((2. * mu + lambda) * duy_dy(0, x, y)
                                            + lambda  * dux_dx(0, x, y));

      return true;
    };

    T dux_dx(T t, T x, T y) {
      return 0.0018 * pi * y * util::cos(2. * pi * x * y);
    }

    T dux_dy(T t, T x, T y) {
      return 0.0018 * pi * x * util::cos(2 * pi * x * y);
    }

    T duy_dx(T t, T x, T y) {
      return 0.0014 * x * util::cos(2 * pi * y);
    }

    T duy_dy(T t, T x, T y) {
      return -0.0014 * pi * (x * x + 1.) * util::sin(2 * pi * y);
    }
};

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  const int noOfCuboids = singleton::mpi().getSize();

  const T conversionFactorLength = params.get<parameters::DOMAIN_EXTENT>()[0] / params.get<parameters::RESOLUTION>();
  const T origin = conversionFactorLength / 3.0;
  const T charL = params.get<parameters::PHYS_CHAR_LENGTH>();

  IndicatorCuboid2D<T> cuboid({params.get<parameters::DOMAIN_EXTENT>()[0] * charL,
    params.get<parameters::DOMAIN_EXTENT>()[1] * charL}, {origin, origin});
  CuboidDecomposition<T, 2> cuboidDecomposition(cuboid, conversionFactorLength, noOfCuboids);

  const T physDeltaX = params.get<parameters::DOMAIN_EXTENT>()[0]/ params.get<parameters::RESOLUTION>();

  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({false,false});
  return mesh;
}

void prepareGeometry( MyCase& myCase, ellipse2D& ellipseCase)
{
  auto& geometry = myCase.getGeometry();
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  geometry.rename( 0, 2 );
  geometry.rename( 2, 1, ellipseCase.ellipse1);
  geometry.rename( 1, 3, ellipseCase.ellipse2);
  geometry.rename( 1, 4, ellipseCase.ellipse3);

  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase, ellipse2D& ellipseCase, const int bulkNum = 1) {
  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& NClattice = myCase.getLattice(NavierCauchy{});

  using NCDESCRIPTOR = MyCase::descriptor_t_of<NavierCauchy>;

  // TODO for now, to be combined with unit converter refactor
  const T plateLength = 1.5;
  const T physDeltaX = plateLength / params.get<parameters::RESOLUTION>();
  const T physDeltaT = physDeltaX * physDeltaX;
  const T charDisplacement = params.get<parameters::PHYS_CHAR_DISPLACEMENT>();
  const T ELattice = params.get<parameters::YOUNGS_MODULUS>();
  const T physPoissonRatio = params.get<parameters::POISSON_RATIO>();
  const T kappa = params.get<parameters::KAPPA>();
  const T charLength = params.get<parameters::PHYS_CHAR_LENGTH>();

  constexpr T theta = invCs2<T, NCDESCRIPTOR>();

  clout << "physDeltaX=" << physDeltaX << std::endl;
  clout << "physDeltaT=" << physDeltaT << std::endl;
  clout << "charDisplacement=" << charDisplacement << std::endl;
  clout << "charLength=" << charLength << std::endl;
  clout << "ELattice=" << ELattice << std::endl;
  clout << "physPoissonRatio=" << physPoissonRatio << std::endl;
  clout << "kappa=" << kappa << std::endl;

  NClattice.setUnitConverter<LinElaUnitConverter<T, NCDESCRIPTOR>>(
    (T) physDeltaX, // physDeltaX
    (T) physDeltaT, // physDeltaT
    (T) charLength,  // charPhysLength
    (T) charDisplacement, // charPhysDisplacement
    (T) ELattice * (physDeltaX * physDeltaX * kappa) / physDeltaT,  // physViscosity
    (T) physPoissonRatio, // physViscosity
    (T) kappa // physThermalConductivity
  );

  const auto& converter = NClattice.getUnitConverter();
  converter.print();

  auto bulkIndicator = geometry.getMaterialIndicator({ bulkNum });

  NClattice.defineDynamics<BoolakeeLinearElasticityBoundary>( bulkIndicator );

  setBoolakeeDirichletBoundary<T,NCDESCRIPTOR,BoolakeeDirichletPostProcessor<T,NCDESCRIPTOR>>(NClattice, geometry.getMaterialIndicator( 2 ), bulkIndicator, ellipseCase.ellipse1);
  setBoolakeeDirichletBoundary<T,NCDESCRIPTOR,BoolakeeDirichletPostProcessor<T,NCDESCRIPTOR>>(NClattice, geometry.getMaterialIndicator( 3 ), bulkIndicator, ellipseCase.ellipse2);
  setBoolakeeDirichletBoundary<T,NCDESCRIPTOR,BoolakeeDirichletPostProcessor<T,NCDESCRIPTOR>>(NClattice, geometry.getMaterialIndicator( 4 ), bulkIndicator, ellipseCase.ellipse3);

  {
    auto& communicator = NClattice.getCommunicator(stage::PostCollide());
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
                theta,
                converter.getLatticeShearModulus(),
                converter.getLatticeLambda(),
                converter.getDampingFactor(),
                converter.getCharPhysVelocity(),
                converter.getEpsilon()};

  const T omega_11  = 1. / (      converter.getLatticeShearModulus()                                 /       theta  + 0.5);
  const T omega_d   = 1. / (2. *  converter.getLatticeShearModulus()                                 / (1. - theta) + 0.5);
  const T omega_s   = 1. / (2. * (converter.getLatticeShearModulus() + converter.getLatticeLambda()) / (1. + theta) + 0.5);

  const T tau_12 = 0.5;
  const T tau_22 = 0.5;

  const T tau_21 = tau_12;

  const T omega_12  = 1. / (tau_12 + 0.5);
  const T omega_21  = 1. / (tau_21 + 0.5);
  const T omega_22  = 1. / (tau_22 + 0.5);

  std::vector<T> allOmegas   = {omega_11, omega_s, omega_d, omega_12, omega_21, omega_22};

  NClattice.setParameter<descriptors::MAGIC_SOLID>(magic);
  NClattice.setParameter<descriptors::OMEGA_SOLID>(allOmegas);

  ForceField2D<T, NCDESCRIPTOR> force(myCase);
  NClattice.defineField<FORCE>(geometry, bulkNum, force);
  NClattice.defineField<POPULATION>(geometry, bulkNum, initialPopulationF);

  NClattice.defineField<DISP_SOLID>(geometry, bulkNum, initialDispF);
  NClattice.defineField<SIGMA_SOLID>(geometry, bulkNum, initialStressF);

  NClattice.defineField<PREVIOUS_CELL>(geometry, 2, initialPopulationF);
  NClattice.defineField<PREVIOUS_CELL>(geometry, 3, initialPopulationF);
  NClattice.defineField<PREVIOUS_CELL>(geometry, 4, initialPopulationF);

  /// Make the lattice ready for simulation
  NClattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

std::vector<MyCase::value_t> getResults( MyCase& myCase,
                           int iT,
                           int maxIt)
{
  using T = MyCase::value_t;
  using NCDESCRIPTOR = MyCase::descriptor_t_of<NavierCauchy>;
  auto& NClattice = myCase.getLattice(NavierCauchy{});
  const auto& converter = NClattice.getUnitConverter();
  auto& geometry = myCase.getGeometry();

  OstreamManager clout(std::cout, "getResults");

  SuperVTMwriter2D<T> vtmWriter("ellipseDirichlet");

  CSV<T> csvWriter("TEMP_CSV");

  SuperGeometryF<T,2> materials(geometry);
  vtmWriter.addFunctor( materials );

  if (iT == 0) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, NCDESCRIPTOR>   cuboid(NClattice);
    SuperLatticeRank2D<T, NCDESCRIPTOR>     rank(NClattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();
  }

  SuperLatticePhysVelocity2D<T, NCDESCRIPTOR>          velocityPhys(NClattice, NClattice.getUnitConverter());
  SuperLatticeVelocity2D<T, NCDESCRIPTOR>              velocity(NClattice);
  SuperLatticeDensity2D<T, NCDESCRIPTOR>               density(NClattice);
  SuperLatticeFpop2D<T, NCDESCRIPTOR>                  population(NClattice);
  SuperLatticeField2D<T, NCDESCRIPTOR, FORCE>          forceField(NClattice);

  ManufacturedSolutionU2D<T> dispSol(myCase);
  SuperLatticeFfromAnalyticalF2D<T, NCDESCRIPTOR>      dispSolLattice(dispSol, NClattice);

  ManufacturedSolutionStress2D<T> stressSol(myCase);
  SuperLatticeFfromAnalyticalF2D<T, NCDESCRIPTOR>      stressSolLattice(stressSol, NClattice);

  // Fields for error calc
  SuperLatticeField2D<T, NCDESCRIPTOR, DISP_SOLID>     moments(NClattice);
  SuperLatticeField2D<T, NCDESCRIPTOR, SIGMA_SOLID>    stress(NClattice);
  auto indicatorF = geometry.getMaterialIndicator(1);

  vtmWriter.addFunctor(population,       "population");
  vtmWriter.addFunctor(moments,          "numerical disp");
  vtmWriter.addFunctor(forceField,       "force");
  vtmWriter.addFunctor(dispSolLattice,   "analytical disp");
  vtmWriter.addFunctor(stressSolLattice, "analytical stress");
  vtmWriter.addFunctor(stress,           "numerical stress");

  vtmWriter.write(iT);

  T   l2UResult[2]        = {T(), T()};
  T   lInfUResult[2]      = {T(), T()};
  T   l2StressResult[2]   = {T(), T()};
  T   lInfStressResult[2] = {T(), T()};
  int tmp[]               = {int()};

  SuperRelativeErrorL2Norm2D<T>   relUErrorL2Norm(NClattice, moments, dispSol, indicatorF);
  SuperRelativeErrorLinfNorm2D<T> relUErrorLinfNorm(NClattice, moments, dispSol, indicatorF);
  SuperRelativeErrorL2Norm2D<T>   relStressErrorL2Norm( NClattice, stress, stressSol, indicatorF );
  SuperRelativeErrorLinfNorm2D<T> relStressErrorLinfNorm( NClattice, stress, stressSol, indicatorF );

  relUErrorL2Norm(l2UResult, tmp);
  relUErrorLinfNorm(lInfUResult, tmp);
  relStressErrorL2Norm(l2StressResult, tmp);
  relStressErrorLinfNorm(lInfStressResult, tmp);

  csvWriter.writeDataFile(iT, l2UResult[0],  "l2UErr");
  csvWriter.writeDataFile(iT, lInfUResult[0], "lInfUErr");
  csvWriter.writeDataFile(iT, l2StressResult[0],  "l2StressErr");
  csvWriter.writeDataFile(iT, lInfStressResult[0], "lInfStressErr");

  std::vector<T> returnVec = {l2UResult[0], lInfUResult[0], l2StressResult[0], lInfStressResult[0]};
  clout << "N\t" << "L2 U Error\t" << "LInf U Error\t" << "L2 Stress Error\t" << "LInf Stress Error" << std::endl;
  clout << converter.getResolution() << "\t" << returnVec[0] << "\t" << returnVec[1] << "\t" << returnVec[2] << "\t\t" << returnVec[3] << std::endl;

  return returnVec;
}

void simulate(MyCase& myCase){
  using T = MyCase::value_t;
  auto& NClattice = myCase.getLattice(NavierCauchy{});
  const auto& converter = NClattice.getUnitConverter();
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  const T simTime = parameters.get<parameters::MAX_PHYS_T>();

  int maxIt = simTime / converter.getConversionFactorTime();

  OstreamManager clout(std::cout, "main");
  clout << "Awaiting " << maxIt << " Time steps. Starting simulation..." << std::endl;
  util::Timer<T> timer(converter.getLatticeTime(simTime),
                       geometry.getStatistics().getNvoxel());
  timer.start();

  int numDataPoints = 400.;
  std::vector<T> errors = {};
  std::ofstream fout;
  std::string dataFile = singleton::directories().getLogOutDir() + std::to_string(parameters.get<parameters::RESOLUTION>()) +  "_Err.dat";
  fout.open(dataFile.c_str(), std::ios::trunc);
  fout << "N;it;l2UErr;lInfUErr;l2StressErr;lInfStressErr" << std::endl;

  for (int iT = 0; iT <= maxIt; ++iT) {

    timer.update(iT);

    if (iT % (maxIt / numDataPoints) == 0) {
      errors = getResults(myCase, iT, maxIt);
      timer.printStep();
      fout << parameters.get<parameters::RESOLUTION>() << ";" << converter.getPhysTime(iT) << ";" << errors[0] << ";" << errors[1] << ";" <<errors[2] << ";" << errors[3] << std::endl;
    }
    NClattice.collideAndStream();
  }
  fout.close();
  timer.stop();
}

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  OstreamManager clout(std::cout, "main");
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION>(40);
    myCaseParameters.set<PHYS_CHAR_LENGTH>(1.0);
    myCaseParameters.set<DOMAIN_EXTENT>({1.5, 1.5});
    myCaseParameters.set<YOUNGS_MODULUS>(0.1);
    myCaseParameters.set<parameters::POISSON_RATIO>(0.7);
    myCaseParameters.set<PHYS_CHAR_DISPLACEMENT>(1.0);
    myCaseParameters.set<KAPPA>(1.0);
    myCaseParameters.set<MAX_PHYS_T>(60.);
  }
  myCaseParameters.fromCLI(argc, argv);

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  ellipse2D ellipseCase(myCaseParameters.get<parameters::PHYS_CHAR_LENGTH>());
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase, ellipseCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase, ellipseCase);

  /// === Step 7: Simulate ===
  simulate(myCase);
}
