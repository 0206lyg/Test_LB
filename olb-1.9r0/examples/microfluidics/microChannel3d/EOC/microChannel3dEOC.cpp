/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Fedor Bukreev, Adrian Kummerländer,
 *  Mathias J. Krause
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
// This is an 2D example for a laminar flow at higher Knudsen numbers
// with Knudsen slip condition and thermal creep along the wall

#include "../microChannel3d.h"

using T = MyCase::value_t;

static Gnuplot<T> gplot(
  "Velocity_eoc",
  false,
  "set terminal png size 720, 720 font 'Arial,10'",
  Gnuplot<T>::LOGLOG,
  Gnuplot<T>::LINREG);

namespace olb::parameters {

struct START_RESOLUTION : public descriptors::TYPED_FIELD_BASE<int,1> { };
struct END_RESOLUTION : public descriptors::TYPED_FIELD_BASE<int,1> { };

}

int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<START_RESOLUTION>(21);
    myCaseParameters.set<END_RESOLUTION>(41);
    myCaseParameters.set<OVERLAP>(5);
    myCaseParameters.set<KNUDSEN>(0.45);
    myCaseParameters.set<B_COEFF>(-1);
    myCaseParameters.set<ACCOMODATION_COEFF>(1);
    myCaseParameters.set<RELAXATION_TIME>(1.1);
    myCaseParameters.set<MAX_PHYS_T>(2.e-6);
    myCaseParameters.set<START_TIME>([&] {
      return 0.05*myCaseParameters.get<parameters::MAX_PHYS_T>();
    });
    myCaseParameters.set<KINETIC_DIAMETER>(364e-12); //kinetic diamter of nitrogen [m]
    myCaseParameters.set<AVERAGE_PRESSURE>(1000);
    myCaseParameters.set<PRESSURE_DIFFERENCE>(0.5);
    myCaseParameters.set<MOLAR_MASS>(28.96e-3);
    myCaseParameters.set<R_SPECIFIC>(296.8);
    myCaseParameters.set<parameters::TEMPERATURE>(300.);
    myCaseParameters.set<parameters::TEMPERATURE_LEFT>(298.);
    myCaseParameters.set<parameters::TEMPERATURE_RIGHT>(302.);
    myCaseParameters.set<PHYS_CHAR_DENSITY>([&] {
      return util::idealGasDensity( myCaseParameters.get<MOLAR_MASS>(), myCaseParameters.get<AVERAGE_PRESSURE>(), myCaseParameters.get<parameters::TEMPERATURE>());
    });
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>([&] {
      return util::gasDynamicViscosity( myCaseParameters.get<MOLAR_MASS>(), myCaseParameters.get<AVERAGE_PRESSURE>(), myCaseParameters.get<parameters::TEMPERATURE>(), myCaseParameters.get<KINETIC_DIAMETER>())/myCaseParameters.get<PHYS_CHAR_DENSITY>();
    });
    myCaseParameters.set<MFP>([&] {
      return util::meanFreePath(myCaseParameters.get<MOLAR_MASS>(), myCaseParameters.get<AVERAGE_PRESSURE>(), myCaseParameters.get<parameters::TEMPERATURE>(), myCaseParameters.get<PHYS_CHAR_VISCOSITY>()*myCaseParameters.get<PHYS_CHAR_DENSITY>());
    });
    myCaseParameters.set<DOMAIN_EXTENT>([&] {
      return Vector{2.*myCaseParameters.get<MFP>()/myCaseParameters.get<KNUDSEN>(), myCaseParameters.get<MFP>()/myCaseParameters.get<KNUDSEN>(), myCaseParameters.get<MFP>()/myCaseParameters.get<KNUDSEN>()};
    });
    myCaseParameters.set<SLIPCOEFF>([&] {
      return (2. - myCaseParameters.get<ACCOMODATION_COEFF>())/myCaseParameters.get<ACCOMODATION_COEFF>()/(1. - myCaseParameters.get<B_COEFF>()*myCaseParameters.get<KNUDSEN>())*myCaseParameters.get<MFP>();
    });
    myCaseParameters.set<CREEPCOEFF>([&] {
      return 3./4.*myCaseParameters.get<PHYS_CHAR_VISCOSITY>()*myCaseParameters.get<PHYS_CHAR_DENSITY>()*myCaseParameters.get<R_SPECIFIC>()/myCaseParameters.get<AVERAGE_PRESSURE>();
    });
    myCaseParameters.set<PHYS_CHAR_VELOCITY>([&] {
      return myCaseParameters.get<DOMAIN_EXTENT>()[1]*myCaseParameters.get<DOMAIN_EXTENT>()[1]/2./myCaseParameters.get<PHYS_CHAR_VISCOSITY>()/myCaseParameters.get<PHYS_CHAR_DENSITY>()
      * (-myCaseParameters.get<PRESSURE_DIFFERENCE>())/myCaseParameters.get<DOMAIN_EXTENT>()[0]
      * (0.25 - 0.5 - myCaseParameters.get<SLIPCOEFF>()/myCaseParameters.get<DOMAIN_EXTENT>()[1])
      + myCaseParameters.get<CREEPCOEFF>()*(myCaseParameters.get<parameters::TEMPERATURE_RIGHT>()
      - myCaseParameters.get<parameters::TEMPERATURE_LEFT>())/myCaseParameters.get<DOMAIN_EXTENT>()[0];
    });
    myCaseParameters.set<PHYS_INTERVAL>(2.e-6/50.);
    myCaseParameters.set<RESIDUUM>(1e-6);
    myCaseParameters.set<HAS_CONVERGED>(false);
    myCaseParameters.set<ERROR_NORMS>({0,0,0,0,0,0});
  }
  myCaseParameters.fromCLI(argc, argv);

  gplot.setLabel("Resolution test", "average Error");

  for(int simuN =  myCaseParameters.get<parameters::START_RESOLUTION>();
      simuN <=  myCaseParameters.get<parameters::END_RESOLUTION>();
      simuN = simuN*2 - 1) {
    myCaseParameters.set<parameters::RESOLUTION>(simuN);
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

    auto errors = myCaseParameters.get<parameters::ERROR_NORMS>();
    gplot.setData (
      T(simuN),
      { errors[1], errors[3], errors[5]},
      { "velocity L1 rel Error","velocity L2 rel Error",
        "velocity Linf rel error"},
      "top right",
      { 'p','p','p' }
    );
  }

  gplot.writePNG();
}
