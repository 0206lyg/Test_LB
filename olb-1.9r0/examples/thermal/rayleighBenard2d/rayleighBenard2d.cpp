/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2008 Orestis Malaspinas
 *                2025 Adrian Kummerlaender
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

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D2Q9<descriptors::FORCE>>,
  Temperature,  Lattice<double, descriptors::D2Q5<descriptors::VELOCITY>>
>;

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T physDeltaX = 0.1 / params.get<parameters::RESOLUTION>();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true,false});
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  //using T = MyCase::template value_t_of<NavierStokes>();
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  geometry.rename(0,2);
  geometry.rename(2,1,{0,1});

  // TODO Parameterization uncovers hidden duplications, deltaX should not be recomputed at use!
  const T physDeltaX = 0.1 / params.get<parameters::RESOLUTION>();

  std::vector<T> extend( 2, T(0) );
  extend[0] = params.get<parameters::DOMAIN_EXTENT>()[0];
  extend[1] = physDeltaX;
  std::vector<T> origin( 2, T(0) );
  IndicatorCuboid2D<T> bottom(extend, origin);

  origin[1] = params.get<parameters::DOMAIN_EXTENT>()[1]-physDeltaX;
  IndicatorCuboid2D<T> top(extend, origin);

  origin[0] = params.get<parameters::DOMAIN_EXTENT>()[0]/2.;
  origin[1] = physDeltaX;
  extend[0] = physDeltaX;
  extend[1] = physDeltaX;
  IndicatorCuboid2D<T> perturbation(extend, origin);

  /// Set material numbers for bottom, top and pertubation
  geometry.rename(2,2,1,bottom);
  geometry.rename(2,3,1,top);
  geometry.rename(1,4,perturbation);

  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase) {
  OstreamManager clout(std::cout,"prepareLattice");

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& NSlattice = myCase.getLattice(NavierStokes{});
  auto& ADlattice = myCase.getLattice(Temperature{});

  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using TDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;

  // TODO for now, to be combined with unit converter refactor
  const T physDeltaX = 0.1 / params.get<parameters::RESOLUTION>();
  const T Ra = params.get<parameters::RAYLEIGH>();
  const T Pr = params.get<parameters::PRANDTL>();
  const int N = params.get<parameters::RESOLUTION>();
  const T Tcold = params.get<parameters::T_COLD>();
  const T Thot = params.get<parameters::T_HOT>();

  NSlattice.setUnitConverter<ThermalUnitConverter<T,NSDESCRIPTOR,TDESCRIPTOR>>(
    (T) physDeltaX, // physDeltaX
    (T) 0.1 / (1e-5 / 0.1 * util::sqrt( Ra / Pr)) * 0.1 / N, // physDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX
    (T) 0.1,  // charPhysLength
    (T) 1e-5 / 0.1 * util::sqrt( Ra / Pr ), // charPhysVelocity
    (T) 1e-5,  // physViscosity
    (T) 1.0, // physDensity
    (T) 0.03, // physThermalConductivity
    (T) Pr * 0.03 / 1e-5 / 1.0,    // physSpecificHeatCapacity
    (T) Ra * 1e-5 * 1e-5 / Pr / 9.81 / (Thot - Tcold) / util::pow(0.1, 3), // physThermalExpansionCoefficient
    (T) Tcold, // charPhysLowTemperature
    (T) Thot // charPhysHighTemperature
  );

  const auto& converter = NSlattice.getUnitConverter();
  converter.print();

  ADlattice.setUnitConverter(converter);

  dynamics::set<ForcedBGKdynamics>(NSlattice, geometry.getMaterialIndicator({1,4}));
  boundary::set<boundary::BounceBack>(NSlattice, geometry.getMaterialIndicator({2,3}));

  dynamics::set<AdvectionDiffusionBGKdynamics>(ADlattice, geometry.getMaterialIndicator({1, 2, 3, 4}));

  boundary::set<boundary::AdvectionDiffusionDirichlet>(ADlattice, geometry, 2);
  boundary::set<boundary::AdvectionDiffusionDirichlet>(ADlattice, geometry, 3);

  NSlattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  ADlattice.setParameter<descriptors::OMEGA>(converter.getLatticeThermalRelaxationFrequency());

  clout << "Setup NSE-ADE coupling" << std::endl;

  const T boussinesqForcePrefactor = 9.81 / converter.getConversionFactorVelocity()
                                   * converter.getConversionFactorTime()
                                   * converter.getCharPhysTemperatureDifference()
                                   * converter.getPhysThermalExpansionCoefficient();

  auto& coupling = myCase.setCouplingOperator(
    "Boussinesq",
    NavierStokesAdvectionDiffusionCoupling{},
    names::NavierStokes{}, NSlattice,
    names::Temperature{},  ADlattice);
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::T0>(
    converter.getLatticeTemperature(Tcold));
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::FORCE_PREFACTOR>(
    boussinesqForcePrefactor * Vector<T,2>{0.0,1.0});

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& NSlattice = myCase.getLattice(NavierStokes{});
  auto& ADlattice = myCase.getLattice(Temperature{});

  const T Tcold = params.get<parameters::T_COLD>();
  const T Thot = params.get<parameters::T_HOT>();
  const T Tperturb = 1./5. * Tcold + 4./5. * Thot; // temperature of the perturbation

  /// for each material set Rho, U and the Equilibrium
  auto indicator = geometry.getMaterialIndicator({1,2,3,4});
  const auto& converter = NSlattice.getUnitConverter();

  /// define initial conditions
  AnalyticalConst2D<T,T> u0(0.0, 0.0);
  AnalyticalConst2D<T,T> T_cold(converter.getLatticeTemperature(Tcold));
  AnalyticalConst2D<T,T> T_hot(converter.getLatticeTemperature(Thot));
  AnalyticalConst2D<T,T> T_perturb(converter.getLatticeTemperature(Tperturb));

  momenta::setTemperature(ADlattice, geometry.getMaterialIndicator({1, 3}), Tcold);
  momenta::setTemperature(ADlattice, geometry.getMaterialIndicator({2}), Thot);
  momenta::setTemperature(ADlattice, geometry.getMaterialIndicator({4}), Tperturb);

  ADlattice.iniEquilibrium(geometry, 1, T_cold, u0);
  ADlattice.iniEquilibrium(geometry, 2, T_hot, u0);
  ADlattice.iniEquilibrium(geometry, 3, T_cold, u0);
  ADlattice.iniEquilibrium(geometry, 4, T_perturb, u0);

  NSlattice.initialize();
  ADlattice.initialize();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{ }

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  using T = MyCase::value_t;
  auto& NSlattice = myCase.getLattice(NavierStokes{});
  auto& ADlattice = myCase.getLattice(Temperature{});
  const auto& converter = NSlattice.getUnitConverter();

  SuperVTMwriter2D<T> vtkWriter("rayleighBenard2d");

  // TODO such times should also be parameters
  const int statIter = converter.getLatticeTime(10.0);
  const int saveIter = converter.getLatticeTime(10.0);

  if (iT == 0) {
    SuperLatticeCuboid2D cuboid(NSlattice);
    SuperLatticeRank2D rank(NSlattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();
  }

  if (iT % statIter == 0) {
    timer.update(iT);
    timer.printStep();
    NSlattice.getStatistics().print(iT, converter.getPhysTime(iT));
  }

  if (iT % saveIter == 0) {
    NSlattice.setProcessingContext(ProcessingContext::Evaluation);
    ADlattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticePhysVelocity2D velocity(NSlattice, converter);
    vtkWriter.addFunctor(velocity);
    SuperLatticePhysPressure2D pressure(NSlattice, converter);
    vtkWriter.addFunctor(pressure);

    using V = MyCase::value_t_of<NavierStokes>;
    using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
    using TDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;
    SuperLatticePhysTemperature2D temperature(
      ADlattice,
      static_cast<const ThermalUnitConverter<V,NSDESCRIPTOR,TDESCRIPTOR>&>(converter));
    vtkWriter.addFunctor(temperature);

    vtkWriter.write(iT);
  }
}

void simulate(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();

  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(
    parameters.get<parameters::MAX_PHYS_T>());

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    myCase.getLattice(NavierStokes{}).collideAndStream();
    myCase.getLattice(Temperature{}).collideAndStream();

    myCase.getOperator("Boussinesq").apply();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT);
  }

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<DOMAIN_EXTENT>({2, 1});
    myCaseParameters.set<RESOLUTION   >(10);
    myCaseParameters.set<RAYLEIGH     >(1e4);
    myCaseParameters.set<PRANDTL      >(0.71);
    myCaseParameters.set<MAX_PHYS_T   >(1000);
    myCaseParameters.set<T_HOT        >(274.15);
    myCaseParameters.set<T_COLD       >(273.15);
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
