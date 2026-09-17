/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2007, 2012 Jonas Latt, Mathias J. Krause
 *  Vojtech Cvrcek, Peter Weisbrod
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

/* poiseuille2d.cpp:
 * This example examines a 2D Poseuille flow
 * It illustrates the usage of different flow setups, boundary conditions
 * and computation of error norms.
 * The following simulation options can be combined freely:
 * - if the compiler flag ENABLE_MRT is set, mrt collision operators are used
 * - forced/ nonForced flow
 * - different boundary conditions
 * - simulation only or eoc convergence analysis
 */

// the main code of the simulation is in poiseuille2d.h as it is also used by the
// example ../../laminar/poiseuille2d

// set flag in order to use mrt collision operators instead of bgk
//#define ENABLE_MRT

#include "../../laminar/poiseuille2d/case.h"

/// Initialize gnuplot
static Gnuplot<T> gplot(
  "Velocity_and_StrainRate_eoc",
  false,
  "set terminal png size 720, 720 font 'Arial,10'",
  Gnuplot<T>::LOGLOG,
  Gnuplot<T>::LINREG);

namespace olb::parameters {

struct START_RESOLUTION : public descriptors::TYPED_FIELD_BASE<int,1> { };
struct NUM_SIMULATIONS : public descriptors::TYPED_FIELD_BASE<int,1> { };

}

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

    MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<START_RESOLUTION>(50);
    myCaseParameters.set<NUM_SIMULATIONS>(4);
    myCaseParameters.set<DOMAIN_EXTENT>({2., 1.});
    myCaseParameters.set<PHYS_CHAR_VELOCITY>(1.0);
    myCaseParameters.set<REYNOLDS>(10);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.8);
    myCaseParameters.set<MAX_PHYS_T>(30.);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1.0);
    myCaseParameters.set<TUNER_PARAM>(0.5);

    myCaseParameters.set<FLOW_TYPE>(FlowType::nonForced);
    myCaseParameters.set<BOUNDARY_TYPE>(BoundaryType::interpolated);

    myCaseParameters.set<CONVERGENCE_PRECISION>(1e-9);
    myCaseParameters.set<CONVERGED>(false);
    myCaseParameters.set<CONV_ITER>(0.25);
    myCaseParameters.set<VTK_ENABLED>(false);
    myCaseParameters.set<GNUPLOT_ENABLED>(false);
    myCaseParameters.set<COMPUTE_ERROR>(true);

    myCaseParameters.set<PHYS_VTK_ITER_T>([&] {
      return myCaseParameters.get<MAX_PHYS_T>()/20.;
    });

    myCaseParameters.set<PHYS_STAT_ITER_T>([&] {
      return myCaseParameters.get<MAX_PHYS_T>()/20.;
    });

  }
  myCaseParameters.fromCLI(argc, argv);

  switch(myCaseParameters.get<parameters::BOUNDARY_TYPE>()){
    case BoundaryType::freeSlip:
      throw std::invalid_argument(
        "eoc computation is currently not supported for free slip boundary conditions");
      break;
    case BoundaryType::partialSlip:
      throw std::invalid_argument(
        "eoc computation is currently not supported for partial slip boundary conditions");
  };

  // set the labels for the plot
  gplot.setLabel("Resolution test", "average Error");

  // loop over the different simulations
  for(int i = 0; i < myCaseParameters.get<parameters::NUM_SIMULATIONS>(); i++){

    myCaseParameters.set<parameters::RESOLUTION>(myCaseParameters.get<parameters::START_RESOLUTION>() + i*10);

    myCaseParameters.set<parameters::CONVERGED>(false);
    /// Run the simulations
    clout << "Starting next simulation with N = " << myCaseParameters.get<parameters::RESOLUTION>() << std::endl;

    /// === Step 3: Create Mesh ===
    Mesh mesh = createMesh(myCaseParameters);

    /// === Step 4: Create Case ===
    MyCase myCase(myCaseParameters, mesh);

    /// === Step 5: Prepare Geometry ===
    prepareGeometry(myCase);

    /// === Step 6: Prepare Lattice ===
    prepareLattice(myCase);

    /// === Step 7: Set Initial Conditions ===
    setInitialValues(myCase);

    /// === Step 8: Simulate ===
    simulate(myCase);

  }

  gplot.writePNG();

  return 0;
}