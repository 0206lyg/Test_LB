/*  Lattice Boltzmann sample, written in C++, using the OpenLB library
 *
 *  Copyright (C) 2006-2025 Fedor Bukreev, Shota Ito,
 *  Jonas Latt, Mathias J. Krause, Vojtech Cvrcek,
 *  Peter Weisbrod, Adrian Kummerländer
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

 /* cylinder2d.h:
  * This example examines a steady flow past a cylinder placed in a channel.
  * The cylinder is offset somewhat from the center of the flow to make the
  * steady-state symmetrical flow unstable. At the inlet, a Poiseuille profile is
  * imposed on the velocity, whereas the outlet implements a Dirichlet pressure
  * condition set by p = 0.
  * Inspired by "Benchmark Computations of Laminar Flow Around
  * a Cylinder" by M.Schäfer and S.Turek.
  * An unsteady flow with Karman vortex street can be created by changing the
  * Reynolds number to Re=100.
  */

#include <olb.h>

using namespace olb;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<NavierStokes, Lattice<FLOATING_POINT_TYPE, descriptors::D2Q9<>>>;

namespace olb::parameters {

struct RADIUS_CYLINDER : public descriptors::FIELD_BASE<1> {};
struct CENTER_CYLINDER : public descriptors::FIELD_BASE<0, 1> {};

}

Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  // Data type for simulation
  using T = MyCase::value_t;

  // Get Domain parameters
  const T physDeltaX  = 2 * parameters.get<parameters::RADIUS_CYLINDER>() / parameters.get<parameters::RESOLUTION>();
  const T physLengthX = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physLengthY = parameters.get<parameters::DOMAIN_EXTENT>()[1] + physDeltaX;

  const Vector extent {physLengthX, physLengthY};
  const Vector origin {0, 0};

  // Create Simulation domain
  IndicatorCuboid2D<T> cuboid(extent, origin);

  // Create mesh from parameters
  Mesh<T, MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

// === Step 2: Prepare Geometry ===
void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  // Get data from case
  using T          = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();

  // Get cylinder data
  const T            radiusCylinder = parameters.get<parameters::RADIUS_CYLINDER>();
  const Vector<T, 2> center         = parameters.get<parameters::CENTER_CYLINDER>();

  // Get domain parameters
  const T lengthX = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T lengthY = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T L       = 0.1 / parameters.get<parameters::RESOLUTION>();

  // Create cylinder indicator
  IndicatorCircle2D<T> circle(center, radiusCylinder);

  geometry.rename(0, 2);
  geometry.rename(2, 1, {1, 1});

  Vector<T, 2> extend(lengthX, lengthY);
  Vector<T, 2> origin;
  //Set material number for inflow
  extend[0] = 2. * L;
  origin[0] = -L;
  IndicatorCuboid2D<T> inflow(extend, origin);
  geometry.rename(2, 3, 1, inflow);

  //Set material number for outflow
  origin[0] = lengthX - L;
  IndicatorCuboid2D<T> outflow(extend, origin);
  geometry.rename(2, 4, 1, outflow);

  // Set material number for cylinder
  geometry.rename(1, 5, circle);

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  geometry.checkForErrors();
  geometry.print();
}

/// @brief Set lattice dynamics
/// @param myCase The Case instance which keeps the simulation data
void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  // Get data from case
  using T          = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(NavierStokes {});

  // Get cylinder data
  const T            radiusCylinder = parameters.get<parameters::RADIUS_CYLINDER>();
  const Vector<T, 2> center         = parameters.get<parameters::CENTER_CYLINDER>();

  // Get simulation parameters
  const T L   = 0.1 / parameters.get<parameters::RESOLUTION>();
  const T CFL = parameters.get<parameters::CFL>();
  const T Re  = parameters.get<parameters::REYNOLDS>();

  // Set up the unit converter
  lattice.setUnitConverter(
      (T)L,                    // physDeltaX: spacing between two lattice cells in [m]
      (T)CFL * L / 0.2,        // physDeltaT: time step in [s]
      (T)2.0 * radiusCylinder, // charPhysLength: reference length of simulation geometry in [m]
      (T)0.2, // charPhysVelocity: highest expected velocity during simulation in [m/s]
      (T)0.2 * 2. * radiusCylinder / Re, // physViscosity: physical kinematic viscosity in [m^2/s]
      (T)1.0                             // physDensity: physical density in [kg/m^3]
  );

  // Print unit converter info
  lattice.getUnitConverter().print();

  // Material=1 -->bulk dynamics
  lattice.defineDynamics<BGKdynamics>(geometry, 1);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(lattice, geometry, 2);

  // Material=3 -->fixed velocity
  boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);

  // Material=4 -->fixed pressure
  boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 4);
  IndicatorCircle2D<T> circle(center, radiusCylinder);

  // Material=5 -->bouzidi
  setBouzidiBoundary(lattice, geometry, 5, circle);
}

/// Set initial condition for primal variables (velocity and density)
/// @param myCase The Case instance which keeps the simulation data
/// @note Be careful: initial values have to be set using lattice units
void setInitialValues(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");

  // Get data from case
  using T        = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& lattice  = myCase.getLattice(NavierStokes {});

  // Initial conditions
  AnalyticalConst2D<T, T> rhoF(1);
  AnalyticalConst2D<T, T> uF(0, 0);

  // Initialize all values of distribution functions to their local equilibrium
  lattice.defineRhoU(geometry, 1, rhoF, uF);
  lattice.iniEquilibrium(geometry, 1, rhoF, uF);
  lattice.setParameter<descriptors::OMEGA>(lattice.getUnitConverter().getLatticeRelaxationFrequency());

  // Initialize Lattice for simulation
  lattice.initialize();
  clout << "Prepare Lattice OK" << std::endl;
}

/// Update boundary values at times (and external fields, if they exist)
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
/// @note Be careful: boundary values have to be set using lattice units
void setTemporalValues(MyCase& myCase, std::size_t iT)
{
  using T          = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(NavierStokes {});
  auto& converter  = lattice.getUnitConverter();

  // Get global parameters
  const T maxPhysT = parameters.get<parameters::MAX_PHYS_T>();
  const T L        = 0.1 / parameters.get<parameters::RESOLUTION>();

  // Number of time steps for smooth start-up
  const std::size_t iTmaxStart = converter.getLatticeTime(maxPhysT * 0.4);
  const std::size_t iTupdate   = iTmaxStart / 1000;

  if (iT % iTupdate == 0 && iT <= iTmaxStart) {
    // Smooth start curve, polynomial
    PolynomialStartScale<T, T> StartScale(iTmaxStart, T(1));
    // Creates and sets the Poiseuille inflow profile using functors
    T iTvec[1] = {T(iT)};
    T frac[1]  = {};
    StartScale(frac, iTvec);
    T               maxVelocity   = converter.getCharLatticeVelocity() * 3. / 2. * frac[0];
    T               distance2Wall = L / 2.;
    Poiseuille2D<T> poiseuilleU(geometry, 3, maxVelocity, distance2Wall);
    lattice.defineU(geometry, 3, poiseuilleU);
    // Update velocity on GPU
    lattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }
}

/// Gets results and writes them to file
/// @param myCase The Case instance which keeps the simulation data
/// @param timer The timer
/// @param iT The time step
void getResults(MyCase& myCase, util::Timer<MyCase::value_t>& timer, std::size_t iT)
{
  OstreamManager clout(std::cout, "getResults");

  // Get data from case
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes {});
  auto& converter  = lattice.getUnitConverter();

  // Get cylinder data
  const T            radiusCylinder  = parameters.get<parameters::RADIUS_CYLINDER>();
  const Vector<T, 2> center          = parameters.get<parameters::CENTER_CYLINDER>();
  const T            centerCylinderX = center[0];
  const T            centerCylinderY = center[1];

  // Set writing intervals
  const std::size_t vtkIter  = converter.getLatticeTime(0.3);
  const std::size_t statIter = converter.getLatticeTime(0.8);

  // Create vtmWriter and add functors
  SuperVTMwriter2D<T>                       vtmWriter("cylinder2d");
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(lattice, converter);
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(lattice, converter);
  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);

  // Create pvd file on first step
  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  // Writes the vtk files
  if (iT % vtkIter == 0 && iT > 0) {
    // Send values from GPU to CPU for evaluation
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }

  // Writes the console log
  if (iT % statIter == 0) {
    // Send values from GPU to CPU for evaluation
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    // Timer console output
    timer.update(iT);
    timer.printStep();
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
    // Pressure drop
    AnalyticalFfromSuperF2D<T> intpolatePressure(pressure, true);
    T point1[2] = {};
    T point2[2] = {};
    point1[0]   = centerCylinderX - radiusCylinder;
    point1[1]   = centerCylinderY;
    point2[0]   = centerCylinderX + radiusCylinder;
    point2[1]   = centerCylinderY;
    T p1, p2;
    intpolatePressure(&p1, point1);
    intpolatePressure(&p2, point2);
    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;
    T pressureDrop = p1 - p2;
    clout << "; pressureDrop=" << pressureDrop << std::endl;
  }
}

/// Executes the simulation
/// @param myCase The Case instance that keeps the simulation data
void simulate(MyCase& myCase)
{
  // Get data from case
  using T          = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(NavierStokes {});
  auto& converter  = lattice.getUnitConverter();

  // Get max simulation time for timer
  const T maxPhysT = parameters.get<parameters::MAX_PHYS_T>();

  // Create timer
  std::size_t    iTmax = converter.getLatticeTime(maxPhysT);
  util::Timer<T> timer(iTmax, geometry.getStatistics().getNvoxel());

  timer.start();

  for (std::size_t iT = 0; iT < iTmax; ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setTemporalValues(myCase, iT);
    // === 6th Step: Collide and Stream Execution ===
    lattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults(myCase, timer, iT);
  }

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    using T = MyCase::value_t;
    myCaseParameters.set<RESOLUTION>(10);
    myCaseParameters.set<REYNOLDS>((T)20.);
    myCaseParameters.set<MAX_PHYS_T>((T)16.);
    myCaseParameters.set<RADIUS_CYLINDER>((T)0.05);
    myCaseParameters.set<CFL>((T)0.05);
    myCaseParameters.set<CENTER_CYLINDER>([&] {
      return Vector {(T)0.2, (T)0.2 + ((T)0.1 / myCaseParameters.get<RESOLUTION>()) / (T)2.};
    });
    myCaseParameters.set<DOMAIN_EXTENT>([&] {
      return Vector {(T)2.2, (T).41 + (T)0.1 / myCaseParameters.get<RESOLUTION>()};
    });
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
