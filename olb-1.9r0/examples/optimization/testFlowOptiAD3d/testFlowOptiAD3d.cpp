/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006-2025 Mathias J. Krause, Fabian Klemens,
 *  Julius Jessberger, Shota Ito
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

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include <olb.h>
#include "analyticalSolutionTestFlow3D.h"

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

using ADf = util::ADf<double, 3>;

/// @brief Step 1: Declare simulation structure.
/// State the type(s) of the used simulation cases and give a name to each case
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<descriptors::FORCE>>
>;
using MyADfCase = Case<
  NavierStokes, Lattice<ADf, descriptors::D3Q19<descriptors::FORCE>>
>;
using MyOptiCase = OptiCaseADf<
  Controlled, MyCase,
  Derivatives, MyADfCase,
  Reference, MyCase
>;

using OBJECTIVE = functors::L2DistanceF<functors::VelocityF>;
using CONTROLS = descriptors::FORCE;

/// @brief Step 3: Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
template<typename PARAMETERS>
auto createMesh(PARAMETERS& parameters) {
  using T = PARAMETERS::value_t;
  const int resolution = parameters.template get<parameters::RESOLUTION>();
  const T physLength = parameters.template get<parameters::PHYS_CHAR_LENGTH>();
  const Vector<T,3> origin{0, 0, 0};
  const Vector<int,3> extend{resolution+1, resolution+1, resolution+1};

  /// @li Create the mesh decomposed into `singleton::mpi().getSize()` cuboids
  const T physDeltaX = physLength / parameters.template get<parameters::RESOLUTION>();
  Mesh<T,3> mesh (origin, physDeltaX, extend, singleton::mpi().getSize());
  mesh.setOverlap(parameters.template get<parameters::OVERLAP>());
  return mesh;
}

/// @brief Step 5: Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers are used to assign physics to lattice nodes
template<typename CASE>
void prepareGeometry(CASE& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  /// @par Contents
  /// @li Store a reference for simplified access
  auto& geometry = myCase.getGeometry();

  /// @li Set material numbers
  geometry.rename(0, 2);
  geometry.rename(2, 1, {1,1,1});

  /// @li Print some information on geometry
  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

/// @brief Step 6: Set lattice dynamics and boundary conditions
/// @param myCase The Case instance which keeps the simulation data
template<typename CASE>
void prepareLattice(CASE& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T = CASE::value_t;
  using DESCRIPTOR = CASE::descriptor_t;
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& parameters = myCase.getParameters();

  /// @li Set up a unit converter with the characteristic physical units
  lattice.template setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    parameters.template get<parameters::RESOLUTION>(),               // resolution
    parameters.template get<parameters::LATTICE_RELAXATION_TIME>(),  // relaxation time
    parameters.template get<parameters::PHYS_CHAR_LENGTH>(),         // charPhysLength: reference length of simulation geometry in [m]
    parameters.template get<parameters::PHYS_CHAR_VELOCITY>(),       // charPhysVelocity: highest expected velocity during simulation in [m/s]
    parameters.template get<parameters::PHYS_CHAR_VISCOSITY>(),      // physViscosity: physical kinematic viscosity in [m^2/s]
    parameters.template get<parameters::PHYS_CHAR_DENSITY>()         // physDensity: physical density [kg/m^3]
  );
  auto& converter = lattice.getUnitConverter();
  converter.print();

  /// @li Material=1 --> bulk dynamics
  dynamics::set<ForcedBGKdynamics>(lattice, geometry.getMaterialIndicator({1}));
  /// @li Material=2,3 --> velocity boundary
  boundary::set<boundary::LocalVelocity>(lattice, geometry.getMaterialIndicator({2}));
  /// @li Set lattice relaxation frequency
  lattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
}

// Step 7: Set initial values for primal variables (e.g. velocity, density) and additional fields
/// @param myCase The Case instance which keeps the simulation data
/// @note Initial values have to be set using lattice units
template<typename CASE>
void setInitialValues(CASE& myCase)
{
  using T = CASE::value_t;
  using DESCRIPTOR = CASE::descriptor_t;
  /// @li Store references for simplified access
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();

  /// @li Set force field
  ForceTestFlow3D<T,T,DESCRIPTOR> forceF(converter);
  const T latticeScaling(converter.getConversionFactorMass() / converter.getConversionFactorForce());
  AnalyticalScaled3D<T,T> scaledForceF(forceF, latticeScaling);  // conversion to lattice units
  fields::set<descriptors::FORCE>(lattice, geometry.getMaterialIndicator({1}), scaledForceF);

  /// @li Initialize lattice
  lattice.initialize();
}

/// Step 8.1: Update boundary values at times (and additional fields, if needed)
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
/// In order to start the simulation slowly, we start with zero boundary velocity and gradually increase it
/// @note Boundary values have to be set using lattice units
template<typename CASE>
void setTemporalValues(CASE& myCase,
                       std::size_t iT)
{
  using T = CASE::value_t;
  using DESCRIPTOR = CASE::descriptor_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();
  const T physStartT = myCase.getParameters().template get<parameters::MAX_PHYS_T>() * (2.0 / 3.0);

  const std::size_t itStart = converter.getLatticeTime(physStartT);
  if (iT <= itStart) {
    /// @li Compute scaling factor
    PolynomialStartScale<T,T> startScaleF(itStart, T(1));
    const T iTvec[1] = {T(iT)};
    T frac[1] = {};
    startScaleF(frac, iTvec);

    /// @li Take analytical velocity solution, scale it to lattice units, set the boundary data
    VelocityTestFlow3D<T,T,DESCRIPTOR> velocityF(converter);
    AnalyticalScaled3D<T,T> uBoundaryStartF(velocityF, frac[0]);
    momenta::setVelocity(lattice, geometry.getMaterialIndicator({2}), uBoundaryStartF);

    /// @li Communicate the new boundary velocity to GPU (if needed)
    lattice.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }
}

/// Step 8.3: Compute simulation results at times
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
template<typename CASE>
void getResults(CASE& myCase,
                util::Timer<typename CASE::value_t>& timer,
                std::size_t iT)
{
  using T = CASE::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();
  const T physMaxT = myCase.getParameters().template get<parameters::MAX_PHYS_T>();

  /// @li Write vtk plots from time to time
  const std::size_t iTvtk = converter.getLatticeTime(physMaxT/5);
  const std::size_t iTlog = converter.getLatticeTime(physMaxT/5);

  SuperVTMwriter3D<T> vtmWriter("testFlow3d");
  SuperLatticePhysVelocity3D velocityF(lattice, converter);
  SuperLatticePhysPressure3D pressureF(lattice, converter);
  vtmWriter.addFunctor(velocityF);
  vtmWriter.addFunctor(pressureF);

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  // Writes the VTK files
  if (iT%iTvtk == 0 && iT > 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }

  /// @li Print some (numerical and computational) statistics
  if (iT%iTlog == 0) {
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.print(iT);
  }
}

/// Step 8: Execute simulation
/// @param myCase The Case instance which keeps the simulation data
/// Run time loop
template<typename CASE>
void simulate(CASE& myCase)
{
  using T = CASE::value_t;
  const T physMaxT = myCase.getParameters().template get<parameters::MAX_PHYS_T>();
  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(physMaxT);
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT)
  {
    /// @li Step 8.1: Update the Boundary Values and Fields at Times
    setTemporalValues(myCase, iT);

    /// @li Step 8.2: Collide and Stream Execution
    myCase.getLattice(NavierStokes{}).collideAndStream();

    /// @li Stripe off density offset due to Dirichlet boundary conditions
    myCase.getLattice(NavierStokes{}).stripeOffDensityOffset(
      myCase.getLattice(NavierStokes{}).getStatistics().getAverageRho() - T{1});

    /// @li Step 8.3: Computation and Output of the Results
    getResults(myCase, timer, iT);
  }

  /// @li Evaluate timer
  timer.stop();
  timer.printSummary();
}

void setInitialControl(MyOptiCase& optiCase) {
  // Intialize the three components of control
  optiCase.getController().set({0, 0, 0});
}

template<typename CASE>
void applyControl(const Controller<typename CASE::value_t>& controller, CASE& controlledCase) {
  using T = CASE::value_t;
  using DESCRIPTOR = CASE::descriptor_t;
  auto& lattice = controlledCase.getLattice(NavierStokes{});
  auto designDomain = controlledCase.getGeometry().getMaterialIndicator({1});
  auto& converter = lattice.getUnitConverter();

  /// Scale the force field component-wise by the control
  std::vector<T> control = controller.get();
  std::shared_ptr<AnalyticalF<3,T,T>> controlF
   = std::make_shared<AnalyticalConst3D<T,T>>(control);
  ForceTestFlow3D<T,T,DESCRIPTOR> forceF(converter);
  const T latticeScaling(converter.getConversionFactorMass() / converter.getConversionFactorForce());
  std::shared_ptr<AnalyticalF<3,T,T>> latticeForceF
   = std::make_shared<AnalyticalScaled3D<T,T>>(forceF, latticeScaling);
  latticeForceF = latticeForceF * controlF;
  fields::set<descriptors::FORCE>(lattice, designDomain, *latticeForceF);
  lattice.template setProcessingContext<Array<descriptors::FORCE>>(ProcessingContext::Simulation);
}

template<typename CASE>
CASE::value_t objectiveF(MyOptiCase& optiCase) {
  // decide whether we solve for value or derivatives
  using T = CASE::value_t;
  auto& controlledCase = optiCase.getCaseByType<T>();
  auto& controlledLattice = controlledCase.getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();
  auto objectiveDomain = controlledCase.getGeometry().getMaterialIndicator({1});

  auto& referenceLattice = optiCase.getCase(Reference{}).getLattice(NavierStokes{});
  auto& refConverter = referenceLattice.getUnitConverter();
  auto refObjectiveDomain = optiCase.getCase(Reference{}).getGeometry().getMaterialIndicator({1});

  // Evaluate functor for objective computation
  auto objectiveO = makeWriteFunctorO<OBJECTIVE,opti::J>(controlledLattice);
  objectiveO->restrictTo(objectiveDomain);

  // Write velocity field in the reference case
  writePhysFunctorTo<functors::VelocityF,OBJECTIVE::Reference>(referenceLattice,
                                                               refObjectiveDomain,
                                                               refConverter.getConversionFactorVelocity());

  // Get solution from the reference simulation for the inverse problem
  copyFields<OBJECTIVE::Reference,OBJECTIVE::Reference>(referenceLattice, controlledLattice);
  objectiveO->template setParameter<descriptors::CONVERSION>(converter.getConversionFactorVelocity());
  objectiveO->template setParameter<descriptors::NORMALIZE>(computeL2Norm<OBJECTIVE::Reference>(referenceLattice, refObjectiveDomain, refConverter.getPhysDeltaX()));
  objectiveO->apply();

  return integrateField<opti::J>(controlledLattice, objectiveDomain, converter.getPhysDeltaX())[0];
}

template<typename CASE>
CASE::value_t computeObjective(MyOptiCase& optiCase) {
  using T = CASE::value_t;
  // decide whether we solve for value or derivatives
  auto& controlledCase = optiCase.getCaseByType<T>();
  // Reset prior simulation lattice
  controlledCase.resetLattices();

  // Prepare new lattices
  prepareLattice(controlledCase);
  setInitialValues(controlledCase);

  // Set updated controls in the simulation
  applyControl(optiCase.getController<T>(), controlledCase);

  // Execute simulation
  simulate(controlledCase);

  // Compute Objective value from simulation results
  return objectiveF<CASE>(optiCase);
}

/// Steps 2-8: Set up and run a simulation
int main(int argc, char* argv[])
{
  /// @par Contents
  /// @li Step 2: Initialization
  initialize(&argc, &argv);
  MyCase::ParametersD myCaseParametersD;
  {
    using namespace parameters;
    myCaseParametersD.template set<PHYS_CHAR_VELOCITY     >(        1.0);
    myCaseParametersD.template set<PHYS_CHAR_LENGTH       >(        1.0);
    myCaseParametersD.template set<PHYS_CHAR_VISCOSITY    >(        0.1);
    myCaseParametersD.template set<PHYS_CHAR_DENSITY      >(        1.0);
    myCaseParametersD.template set<MAX_PHYS_T             >(        6.0);
    myCaseParametersD.template set<RESOLUTION             >(         11);
    myCaseParametersD.template set<LATTICE_RELAXATION_TIME>(        0.8);
  }
  MyADfCase::ParametersD myADfCaseParametersD;
  {
    using namespace parameters;
    myADfCaseParametersD.template set<PHYS_CHAR_VELOCITY     >( ADf{1.0});
    myADfCaseParametersD.template set<PHYS_CHAR_LENGTH       >( ADf{1.0});
    myADfCaseParametersD.template set<PHYS_CHAR_VISCOSITY    >( ADf{0.1});
    myADfCaseParametersD.template set<PHYS_CHAR_DENSITY      >( ADf{1.0});
    myADfCaseParametersD.template set<MAX_PHYS_T             >( ADf{6.0});
    myADfCaseParametersD.template set<RESOLUTION             >(       11);
    myADfCaseParametersD.template set<LATTICE_RELAXATION_TIME>(ADf{ 0.8});
  }

  /// @li Step 3: Create Mesh
  Mesh mesh = createMesh(myCaseParametersD);
  Mesh aDfMesh = createMesh(myADfCaseParametersD);

  // ==== BELOW HERE OPTIMIZATION SPECIFIC ====
  /// @li Step B: Create Cases and prepare
  // Prepare controlled case
  MyCase myCase(myCaseParametersD, mesh);
  prepareGeometry(myCase);

  // Create ADf-typed case for gradient computation
  MyADfCase aDfCase(myADfCaseParametersD, aDfMesh);
  prepareGeometry(aDfCase);

  // Compute solution for the objective functional
  MyCase referenceCase(myCaseParametersD, mesh);
  prepareGeometry(referenceCase);
  prepareLattice(referenceCase);
  setInitialValues(referenceCase);
  simulate(referenceCase);

  /// @li Step C: Create OptiCase and set cases
  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);
  optiCase.setCase<Derivatives>(aDfCase);
  optiCase.setCase<Reference>(referenceCase);

  /// @li Step D: Set initial control
  setInitialControl(optiCase);

  /// @li Step E: Define objective routine
  optiCase.setObjective(computeObjective<MyCase>, computeObjective<MyADfCase>);

  /// @li Step G: Create an Optimizer
  OptimizerLBFGS<MyCase::value_t,std::vector<MyCase::value_t>> optimizer(
    optiCase.getController().size(), 1.e-5, 10, 1., 20, "Wolfe", 20, 1.e-4);

  /// @li Step H: Optimize
  optimizer.optimize(optiCase);
}
