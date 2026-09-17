/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Shota Ito, Felix Schuhmann, Julius Jeßberger
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

/** \file
 * In this example, the pipebend optimization benchmark case is implemented,
 * cf. e.g. https://doi.org/10.1016/j.buildenv.2024.112508.
 * For given in- and outlets at the left and at the bottom side of a cavity,
 * the (curved) shape of a pipe has to be found, s.t. either the overall
 * pressure drop or the overall dissipation is minimized.
 *
 * We use adjoint optimization for this. Hence, we implement
 * - a solver class, which implements a (primal or dual) simulation
 * - parameter structures, which keep simulation- and otimization-related parameters
 * - two objective classes, which implement the objective computation (one for
 *   dissipation minimization, one for pressure drop minimization)
 *
 * Since, for low Reynolds numbers, leaving the complete design domain open leads to
 * a low dissipation, we penalize the volume of the pipe in order to prevent this.
 * So, we can establish a geometrically meaningful geometry.
 *
 * The values of simulation- and optimization-related parameters are set in the
 * corresponding file parameters.xml. Feel free to change the Reynolds number and
 * see how the optimal pipe shape changes. This may of course necessitate adjustments
 * of e.g. the resolution.
 */

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

using MyCase = Case<
  NavierStokes, Lattice<double,descriptors::PorousD2Q9Descriptor>
>;
using MyOptiCase = OptiCaseAdjoint<
  Controlled, MyCase,
  Adjoint, MyCase
>;
using DYNAMICS = PorousBGKdynamics<MyCase::value_t,MyCase::descriptor_t>;
using ControlledField = descriptors::POROSITY;
using ObjectiveF = functors::AddF<functors::AddF<functors::DissipationF,
                                                 functors::PorousDissipationF>,
                                  functors::TikhonovRegularizationF<ControlledField>>;

auto createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const T L = parameters.get<parameters::PHYS_CHAR_LENGTH>() / parameters.get<parameters::RESOLUTION>();
  const T lengthX = parameters.get<parameters::PHYS_CHAR_LENGTH>() + 2*L;
  const T lengthY = lengthX;
  const T inflowLength = 0.5*lengthX;
  const T outflowLength = inflowLength;
  const Vector<T,2> origin(-inflowLength, -outflowLength);
  const Vector<T,2> extend(lengthX + inflowLength, lengthY + outflowLength);
  IndicatorCuboid2D<T> domain(extend, origin);
  Mesh<T,MyCase::d> mesh{domain, L, singleton::mpi().getSize()};
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& superGeometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  const T L = parameters.get<parameters::PHYS_CHAR_LENGTH>() / parameters.get<parameters::RESOLUTION>();
  const T lengthX = parameters.get<parameters::PHYS_CHAR_LENGTH>() + 2*L;
  const T lengthY = lengthX;
  const T inflowY = 0.8*parameters.get<parameters::PHYS_CHAR_LENGTH>() + L;
  const T inflowRadius = 0.1*parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const T inflowLength = 0.5*lengthX;
  const T outflowX = inflowY;
  const T outflowRadius = inflowRadius;
  const T outflowLength = inflowLength;

  // Flow domain
  Vector<T,2> origin(0., 0.);
  Vector<T,2> extend(lengthX, lengthY);
  IndicatorCuboid2D<T> pipeDomain(extend, origin);

  // Inflow regions
  origin = Vector<T,2>(-inflowLength - 0.5*L, inflowY - inflowRadius);
  extend = Vector<T,2>(L, 2*inflowRadius);
  IndicatorCuboid2D<T> inflow(extend, origin);
  origin[0] = -0.5*L;
  IndicatorCuboid2D<T> objectiveDomainInflow(extend, origin);
  origin = Vector<T,2>(-inflowLength - 0.5*L, inflowY - inflowRadius - L);
  extend = Vector<T,2>(inflowLength + L, 2*inflowRadius + 2*L);
  IndicatorCuboid2D<T> extendedInflow(extend, origin);

  // Outflow regions
  origin = Vector<T,2>(outflowX - outflowRadius, -outflowLength - 0.5*L);
  extend = Vector<T,2>(2*outflowRadius, L);
  IndicatorCuboid2D<T> outflow(extend, origin);
  origin[1] = -0.5*L;
  IndicatorCuboid2D<T> objectiveDomainOutflow(extend, origin);
  origin = Vector<T,2>(outflowX - outflowRadius - L, -outflowLength - 0.5*L);
  extend = Vector<T,2>(2*outflowRadius + 2*L, outflowLength + L);
  IndicatorCuboid2D<T> extendedOutflow(extend, origin);

  // Material numbers
  // 1: fluid, 2: wall, 3: inflow, 4: outflow, 6: designDomain
  // 7,8: objectiveDomain (pressure)
  superGeometry.rename(0, 2, pipeDomain);
  superGeometry.rename(0, 2, extendedInflow);
  superGeometry.rename(0, 2, extendedOutflow);
  superGeometry.rename(2, 1, {1, 1});
  superGeometry.rename(2, 3, inflow);
  superGeometry.rename(2, 4, outflow);
  superGeometry.rename(1, 6, pipeDomain);
  superGeometry.rename(6, 7, objectiveDomainInflow);
  superGeometry.rename(6, 8, objectiveDomainOutflow);
}

void prepareLattice(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& superGeometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& parameters = myCase.getParameters();

  // define UnitConverter
  sLattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,MyCase::descriptor_t>>(
    (int) parameters.get<parameters::RESOLUTION>(),
    ( T ) parameters.get<parameters::LATTICE_RELAXATION_TIME>(),
    ( T ) parameters.get<parameters::PHYS_CHAR_LENGTH>(),
    ( T ) parameters.get<parameters::PHYS_CHAR_VELOCITY>(),
    ( T ) parameters.get<parameters::PHYS_CHAR_VISCOSITY>(),
    ( T ) parameters.get<parameters::PHYS_CHAR_DENSITY>()
  );
  auto& converter = sLattice.getUnitConverter();

  // define dynamics and bc
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,3,4,6,7,8});
  dynamics::set<DYNAMICS>(sLattice, bulkIndicator);
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);

  // set parameters
  sLattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& superGeometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes{});

  // Set porosity field
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,3,4,6,7,8});
  auto solidIndicator = superGeometry.getMaterialIndicator({0,2});
  AnalyticalConst2D<T,T> zero(0.);
  AnalyticalConst2D<T,T> one(1.);

  fields::set<descriptors::POROSITY>(sLattice, bulkIndicator, one);
  fields::set<descriptors::POROSITY>(sLattice, solidIndicator, zero);

  sLattice.initialize();
}

void setTemporalValues(MyCase& myCase, std::size_t iT) {
  using T = MyCase::value_t;
  auto& superGeometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& converter = sLattice.getUnitConverter();
  auto& parameters = myCase.getParameters();
  const T L = parameters.get<parameters::PHYS_CHAR_LENGTH>() / parameters.get<parameters::RESOLUTION>();
  const T physStartUpT = parameters.get<parameters::MAX_PHYS_T>() / 4;
  const T physBoundaryValueUpdateT = parameters.get<parameters::PHYS_BOUNDARY_VALUE_UPDATE_T>();

  std::size_t iTmaxStart = converter.getLatticeTime(physStartUpT);
  std::size_t iTupdate = converter.getLatticeTime(physBoundaryValueUpdateT);
  if ( iT%iTupdate==0 && iT<= iTmaxStart ) {
    PolynomialStartScale<T,T> StartScale( iTmaxStart, T( 1 ) );
    T iTvec[1] = {T( iT )}; T frac[1] = {};
    StartScale( frac,iTvec );
    T maxVelocity = converter.getCharPhysVelocity()*3./2.*frac[0];
    T distance2Wall = L/2.;
    Poiseuille2D<T> poiseuilleU(superGeometry, 3, maxVelocity, distance2Wall);
    momenta::setVelocity(sLattice, superGeometry.getMaterialIndicator({3}), poiseuilleU);
    sLattice.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

void getResults(MyCase& myCase, util::Timer<MyCase::value_t>& timer, std::size_t iT) {
  using T = MyCase::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& converter = sLattice.getUnitConverter();
  auto& parameters = myCase.getParameters();
  const T physMaxT = parameters.get<parameters::MAX_PHYS_T>();

  const std::size_t iTlog = converter.getLatticeTime(physMaxT/5.);

  if (iT%iTlog == 0) {
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.print(iT);
  }
}

void simulate(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& superGeometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& converter = sLattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  const std::size_t iTmax = converter.getLatticeTime(parameters.get<parameters::MAX_PHYS_T>());
  util::Timer<T> timer(iTmax, superGeometry.getStatistics().getNvoxel());
  timer.start();

  // main loop
  for (std::size_t iT=0; iT < iTmax; ++iT)
  {
    if (myCase.hasName(Controlled{})) {
      setTemporalValues(myCase, iT);
    }
    myCase.getLattice(NavierStokes{}).collideAndStream();
    getResults(myCase, timer, iT);
  }

  timer.stop();
  timer.printSummary();
  myCase.getLattice(NavierStokes{}).setProcessingContext(ProcessingContext::Evaluation);
}

void setInitialControl(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierStokes{});
  auto& geometry = controlledCase.getGeometry();
  auto& converter = lattice.getUnitConverter();
  auto& control = optiCase.getController();
  auto& parameters = controlledCase.getParameters();

  // Intialize controlled field in superLattice
  T porosity = projection::permeabilityToPorosity(parameters.get<parameters::INITIAL_CONTROL_SCALAR>(), converter);
  std::cout << "POROSITY: " << porosity << std::endl;
  AnalyticalConst2D<T,T> controls(porosity);
  fields::set<ControlledField>(lattice, geometry.getMaterialIndicator({6}), controls);
  control.template setProjection<projection::Sigmoid<T>>();
  control.template set<ControlledField>(geometry, 6, lattice);
}

void prepareAdjointLattice(MyOptiCase& optiCase) {
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& geometry = optiCase.getCase(Adjoint{}).getGeometry();

  adjointLattice.setUnitConverter(controlledLattice.getUnitConverter());

  // Define dual dynamics
  dynamics::set<DualPorousBGKDynamics>(adjointLattice, geometry.getMaterialIndicator({1,6,7,8}));
  boundary::set<boundary::BounceBack>(adjointLattice, geometry.getMaterialIndicator({2,3,4}));
  adjointLattice.template setParameter<descriptors::OMEGA>(adjointLattice.getUnitConverter().getLatticeRelaxationFrequency());
}

void setAdjointInitialValues(MyOptiCase& optiCase) {
  auto& adjointCase = optiCase.getCase(Adjoint{});
  auto& adjointLattice = adjointCase.getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& converter = adjointLattice.getUnitConverter();
  // Dissipation (6) or pressure drop (7,8)
  auto objectiveDomain = optiCase.getCase(Adjoint{}).getGeometry().getMaterialIndicator({6});
  //auto objectiveDomain = optiCase.getCase(Adjoint{}).getGeometry().getMaterialIndicator({7,8});

  setInitialValues(adjointCase);

  auto objectiveDerivativeO = makeWriteFunctorO<functors::DerivativeF<ObjectiveF, descriptors::POPULATION,DYNAMICS>,
                                                opti::DJDF>(controlledLattice);
  objectiveDerivativeO->restrictTo(objectiveDomain);
  objectiveDerivativeO->template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  objectiveDerivativeO->template setParameter<descriptors::DT>(converter.getPhysDeltaT());
  objectiveDerivativeO->template setParameter<descriptors::PHYS_VISCOSITY>(converter.getPhysViscosity());
  objectiveDerivativeO->template setParameter<descriptors::VISCOSITY>(converter.getLatticeViscosity());
  objectiveDerivativeO->template setParameter<descriptors::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
  objectiveDerivativeO->template setParameter<descriptors::DX>(converter.getPhysDeltaX());
  objectiveDerivativeO->apply();

  // Write primal population in dedicated fields
  writeFunctorTo<functors::PopulationF,opti::F>(controlledLattice,
                                                optiCase.getCase(Controlled{}).getGeometry().getMaterialIndicator({6,7,8}));

  // Provide fields required by the dual collision operator
  copyFields<ControlledField,ControlledField>(controlledLattice, adjointLattice);
  copyFields<opti::F,opti::F>(controlledLattice, adjointLattice);
  copyFields<opti::DJDF,opti::DJDF>(controlledLattice, adjointLattice);
}

void applyControl(MyOptiCase& optiCase) {
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierStokes{});

  optiCase.getController().template setUpdatedControlsOnField<ControlledField>(lattice);
  lattice.template setProcessingContext<Array<ControlledField>>(ProcessingContext::Simulation);
}

MyCase::value_t objectiveF(MyOptiCase& optiCase) {
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();
  // Dissipation (6) or pressure drop (7,8)
  auto objectiveDomain = optiCase.getCase(Controlled{}).getGeometry().getMaterialIndicator({6});
  //auto objectiveDomain = optiCase.getCase(Controlled{}).getGeometry().getMaterialIndicator({7,8});
  auto& parameters = optiCase.getCase(Controlled{}).getParameters();

  auto objectiveO = makeWriteFunctorO<ObjectiveF, opti::J>(controlledLattice);
  objectiveO->restrictTo(objectiveDomain);
  objectiveO->template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  objectiveO->template setParameter<descriptors::DT>(converter.getPhysDeltaT());
  objectiveO->template setParameter<descriptors::PHYS_VISCOSITY>(converter.getPhysViscosity());
  objectiveO->template setParameter<descriptors::DX>(converter.getPhysDeltaX());
  objectiveO->template setParameter<descriptors::VISCOSITY>(converter.getLatticeViscosity());
  objectiveO->template setParameter<descriptors::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
  objectiveO->template setParameter<opti::REG_ALPHA>(parameters.get<parameters::REGULARIZATION_FACTOR>());
  objectiveO->apply();

  auto volume_fraction = makeWriteFunctorO<functors::TikhonovRegularizationF<ControlledField>, opti::REGULARIZATION>(controlledLattice);
  volume_fraction->restrictTo(objectiveDomain);
  volume_fraction->template setParameter<opti::REG_ALPHA>(1. / 0.25);
  volume_fraction->apply();

  std::cout << "volume-fraction: " << integrateField<opti::REGULARIZATION>(controlledLattice, objectiveDomain, converter.getPhysDeltaX())[0] << std::endl;
  return integrateField<opti::J>(controlledLattice, objectiveDomain, converter.getPhysDeltaX())[0];
}

std::vector<MyCase::value_t> derivativeF(MyOptiCase& optiCase) {
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();
  auto objectiveDomain = optiCase.getCase(Controlled{}).getGeometry().getMaterialIndicator({6});
  auto& parameters = optiCase.getCase(Adjoint{}).getParameters();

  // Evaluate optimality condition
  auto optimalityO = makeWriteFunctorO<functors::OptimalityF<DYNAMICS,ControlledField>,
                                       opti::SENSITIVITY<ControlledField>>(adjointLattice);
  optiCase.getController().template setUpdatedProjectionDerivativesOnField<opti::DPROJECTIONDALPHA<ControlledField>>(adjointLattice);
  adjointLattice.template setProcessingContext<Array<opti::DPROJECTIONDALPHA<ControlledField>>>(ProcessingContext::Simulation);
  // Compute jacobian of collision operator regarding control variable
  auto dCDalphaO = makeWriteFunctorO<functors::DerivativeF<functors::CollisionF<DYNAMICS>,ControlledField,DYNAMICS>,
                                     opti::DCDALPHA<ControlledField>>(controlledLattice);
  dCDalphaO->template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  dCDalphaO->template setParameter<descriptors::DX>(1.0);
  dCDalphaO->apply();
  // Compute objective derivative regarding controlls
  auto objectiveDerivativeO = makeWriteFunctorO<functors::DerivativeF<ObjectiveF,ControlledField,DYNAMICS>,
                                                opti::DJDALPHA<ControlledField>>(controlledLattice);
  objectiveDerivativeO->restrictTo(objectiveDomain);
  objectiveDerivativeO->template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  objectiveDerivativeO->template setParameter<descriptors::PHYS_VISCOSITY>(converter.getPhysViscosity());
  objectiveDerivativeO->template setParameter<descriptors::VISCOSITY>(converter.getLatticeViscosity());
  objectiveDerivativeO->template setParameter<descriptors::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
  objectiveDerivativeO->template setParameter<descriptors::DX>(converter.getPhysDeltaX());
  objectiveDerivativeO->template setParameter<descriptors::DT>(converter.getPhysDeltaT());
  objectiveDerivativeO->template setParameter<opti::REG_ALPHA>(parameters.get<parameters::REGULARIZATION_FACTOR>());
  objectiveDerivativeO->apply();
  // Jacobian is computed on primal lattice as jacobian is evaluated for primal populations
  copyFields<opti::DCDALPHA<ControlledField>,opti::DCDALPHA<ControlledField>>(controlledLattice, adjointLattice);
  copyFields<opti::DJDALPHA<ControlledField>,opti::DJDALPHA<ControlledField>>(controlledLattice, adjointLattice);
  optimalityO->template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  optimalityO->apply();

  // Return serial vector containing total derivatives of objective regarding controls
  adjointLattice.setProcessingContext(ProcessingContext::Evaluation);
  return getSerializedFromField<opti::SENSITIVITY<ControlledField>>(adjointLattice, optiCase.getController().getDesignDomain<MyCase::descriptor_t>());
}

void getOptiResults(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();
  std::size_t iT = optiCase.getOptimizationStep();

  SuperVTMwriter2D<T> vtmWriter("pipeTopology2d");
  SuperLatticePhysVelocity2D velocityF(controlledLattice, converter);
  SuperLatticePhysPressure2D pressureF(controlledLattice, converter);
  SuperLatticeField2D<T,DESCRIPTOR,descriptors::POROSITY> porosityF(controlledLattice);
  SuperLatticeField2D<T,DESCRIPTOR,opti::J> jF(controlledLattice);
  SuperLatticeField2D<T,DESCRIPTOR,opti::DJDF> djdfF(adjointLattice);
  vtmWriter.addFunctor(velocityF);
  vtmWriter.addFunctor(pressureF);
  vtmWriter.addFunctor(porosityF);
  vtmWriter.addFunctor(jF);
  vtmWriter.addFunctor(djdfF);

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  // Writes the VTK files
  if (iT >= 0) {
    controlledLattice.setProcessingContext(ProcessingContext::Evaluation);
    adjointLattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }
}

MyCase::value_t computeObjective(MyOptiCase& optiCase) {
  auto& controlledCase = optiCase.getCase(Controlled{});
  // Reset prior simulation lattice
  controlledCase.resetLattices();

  // Prepare new lattices
  prepareLattice(controlledCase);
  setInitialValues(controlledCase);

  // Set updated controls in the simulation
  applyControl(optiCase);

  // Execute simulation
  simulate(controlledCase);

  // Compute Objective value from simulation results
  return objectiveF(optiCase);
}

std::vector<MyCase::value_t> computeDerivative(MyOptiCase& optiCase) {
  auto& adjointCase = optiCase.getCase(Adjoint{});
  // Reset prior simulation lattice
  adjointCase.resetLattices();

  // Prepare new lattices
  prepareAdjointLattice(optiCase);
  setAdjointInitialValues(optiCase);

  // Execute simulation
  simulate(adjointCase);

  // Optimization related IO
  getOptiResults(optiCase);

  // Compute Derivative value from simulation results
  return derivativeF(optiCase);
}

int main(int argc, char **argv) {
  initialize(&argc, &argv);

  MyCase::ParametersD myCaseParametersD;
  {
    using namespace parameters;
    myCaseParametersD.set<RESOLUTION                  >(    25);
    myCaseParametersD.set<LATTICE_RELAXATION_TIME     >(  0.51);
    myCaseParametersD.set<PHYS_CHAR_LENGTH            >(   0.5);
    myCaseParametersD.set<PHYS_CHAR_VELOCITY          >(   0.1);
    myCaseParametersD.set<PHYS_CHAR_VISCOSITY         >( 0.001);
    myCaseParametersD.set<PHYS_CHAR_DENSITY           >(   1.0);
    myCaseParametersD.set<MAX_PHYS_T                  >(   30.);
    myCaseParametersD.set<PHYS_BOUNDARY_VALUE_UPDATE_T>(   0.1);
    myCaseParametersD.set<INITIAL_CONTROL_SCALAR      >(1.5e-5);
    myCaseParametersD.set<REGULARIZATION_FACTOR       >( 0.003);
  }

  Mesh mesh = createMesh(myCaseParametersD);
  MyCase myCase(myCaseParametersD, mesh);
  prepareGeometry(myCase);
  prepareLattice(myCase);

  MyCase adjointCase(myCaseParametersD, mesh);
  prepareGeometry(adjointCase);

  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);
  optiCase.setCase<Adjoint>(adjointCase);

  setInitialControl(optiCase);
  optiCase.setObjective(computeObjective);
  optiCase.setDerivative(computeDerivative);

  OptimizerLBFGS<MyCase::value_t,std::vector<MyCase::value_t>> optimizer(
    optiCase.getController().size(), 1.e-10, 25
    , 1., 20, "Wolfe", 20, 1.e-4);
  optimizer.optimize(optiCase);
}
