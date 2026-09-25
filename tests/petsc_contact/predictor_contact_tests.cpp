// Production failure regressions: geometry-safe Newton initialization and
// coupled normal/friction transitions, continued to the original LB endpoint.
#include "particleSubsteps.h"
#ifndef SLURRY_USE_PETSC
#error "This regression requires the real PETSc backend"
#endif
#include <filesystem>
#include <iomanip>
#include <iostream>

namespace g=slurry::gr_re2::graphite;
namespace d=g::particle_detail;

void require(bool value,const std::string& message){if(!value)throw std::runtime_error(message);}

void checkRejectedPredictor(const std::filesystem::path& file){
  auto x=d::readParticleReplay(file.string());const auto original=x.bodies;
  const std::size_t n=x.bodies.size();double L=0.;for(const auto& body:x.bodies)L=std::max(L,d::radius(body));
  std::vector<int> slots(x.contacts.size(),-1);d::Vector q(6*n);
  for(std::size_t i=0;i<n;++i)for(int k=0;k<3;++k){q[6*i+k]=x.bodies[i].velocity[k]*x.dt/L;q[6*i+k+3]=x.bodies[i].omega[k]*x.dt;}
  d::Residual residual{x.bodies,x.force,x.torque,x.settings,x.cache,x.contacts,slots,
                        x.dt,x.time,L,x.bodies.front().mass*L/(x.dt*x.dt),0};
  const auto predicted=q;std::string error;
  require(!d::contactFeasiblePredictor(residual,q,false,error),"Severe velocity predictor must be rejected");
  require(q==predicted&&residual.evaluations==0,"Geometry check must not mutate the guess or evaluate forces");
  d::Evaluation bad;
  require(residual(q,bad,error),"Old positive-gap guard admitted this predictor");
  require(bad.diagnostic.maxForce>1.,"Fixture must reproduce the >1 N repulsive trial force");
  for(double& value:q)value=0.;
  require(d::contactFeasiblePredictor(residual,q,true,error),"Saved state must provide a feasible restored guess");
  require(d::admissibleTrialGap(x.settings.rough.gap-2.*x.settings.contactGapTolerance,x.settings),
          "Jacobian probes retain the smooth trial extension outside contact feasibility");
  for(std::size_t i=0;i<n;++i)
    require(x.bodies[i].position==original[i].position&&x.bodies[i].rotation==original[i].rotation
            &&x.bodies[i].velocity==original[i].velocity,"Initial-guess restoration changed committed bodies");
  std::cout<<"PASS: rejects the singular predictor before force evaluation; immutable committed state\n";
}

void checkAdvancedShearPhase(){
  g::ParticleStepSettings settings;settings.shearRate=1.;settings.maxLineSearch=14;
  settings.rough.enabled=true;
  settings.rough.gap=2.e-9;settings.contactGapTolerance=1.e-12;
  const double radius=1.e-6,dt=1.e-4,L=radius;
  std::vector<g::Body> old(2);
  for(auto& body:old){body.axes={radius,radius,radius};body.mass=2.e-11;body.inertiaBody={8.e-24,8.e-24,8.e-24};}
  const double component=(2.*radius+settings.rough.gap)/std::sqrt(2.);
  old[0].position={4.e-6,.2e-6,5.e-6};
  old[1].position={4.e-6+component,settings.box[1]+.2e-6-component,5.e-6};
  const auto original=old;
  std::vector<g::Vec3> force(2),torque(2);std::vector<g::GapCache> cache(1);
  g::PersistentContactState contacts(1);std::vector<int> slots(1,-1);d::Vector q(12,0.);
  d::Residual residual{old,force,torque,settings,cache,contacts,slots,dt,0.,L,old[0].mass*L/(dt*dt),0};
  const auto oldGap=g::closestEllipsoidGap(old[0],d::closestImage(old[0],old[1],0.,settings));
  require(oldGap.gap>=settings.rough.gap-settings.contactGapTolerance,"LE fixture starts feasible");
  std::string error;
  require(!d::contactFeasiblePredictor(residual,q,false,error),"Zero displacement must see advanced image shear");
  require(d::contactFeasiblePredictor(residual,q,true,error),"Geometry restoration must handle advanced image shear");
  require(d::contactFeasiblePredictor(residual,q,false,error),"Restored geometry must independently pass");
  require(residual.evaluations==0,"Restoration must not evaluate the singular potential");
  for(int k=0;k<3;++k){
    require(std::abs(q[k]+q[6+k])<1.e-15,"Mass-weighted restoration preserves the center of mass");
    require(old[0].position[k]==original[0].position[k]&&old[1].position[k]==original[1].position[k],
            "LE restoration may not move committed positions");
  }
  std::cout<<"PASS: geometry-only restoration at advanced Lees--Edwards phase\n";
}

void continueFixture(const std::filesystem::path& file,bool expectFirstRestart=false){
  auto x=d::readParticleReplay(file.string());x.settings.solverBackend="petsc";
  d::ParticleDiagnosticScope diagnostics(x.settings);
  int newton=0,krylov=0;double maxForce=0.,maxTorque=0.,maxGap=0.;
  for(int k=x.substep;k<x.count;++k){
    std::vector<g::Body> next;g::ParticleStepDiagnostics info;std::string error;
    const int restartsBefore=diagnostics.normalGuessRestarts;
    const bool success=d::dispatchImplicitStep(x.bodies,x.force,x.torque,x.dt,x.outerTime+k*x.dt,x.settings,
                                               x.cache,x.contacts,next,info,error,x.outerTime,x.outerDt,x.count,k);
    require(success,file.filename().string()+" substep "+std::to_string(k)+": "+error);
    require(info.maxForceResidualRatio<=1.&&info.maxTorqueResidualRatio<=1.
            &&info.contactGapViolation<=x.settings.contactGapTolerance,"Original physical acceptance criteria failed");
    require(info.newtonIterations<=x.settings.maxNewtonIterations,"Original Newton budget exceeded");
    const int restarts=diagnostics.normalGuessRestarts-restartsBefore;
    require(restarts<=1,"A substep may restart its normal guess only once");
    if(expectFirstRestart&&k==x.substep)
      require(restarts==1,"Saved cycling case must exercise the normal-guess restart");
    x.bodies=std::move(next);newton+=info.newtonIterations;krylov+=info.krylovIterations;
    maxForce=std::max(maxForce,info.maxForceResidualRatio);maxTorque=std::max(maxTorque,info.maxTorqueResidualRatio);
    maxGap=std::max(maxGap,info.contactGapViolation);
  }
  std::cout<<std::setprecision(12)<<"PASS: "<<file.filename().string()<<" remaining="<<x.count-x.substep
           <<" Newton="<<newton<<" Krylov="<<krylov<<" force_ratio="<<maxForce
           <<" torque_ratio="<<maxTorque<<" gap_violation_m="<<maxGap
           <<" normal_guess_restarts="<<diagnostics.normalGuessRestarts<<std::endl;
}

int main(int argc,char** argv){
  if(PetscInitialize(&argc,&argv,nullptr,nullptr))return 1;
  int result=0;
  try{
    const auto directory=std::filesystem::path(__FILE__).parent_path()/"fixtures";
    checkRejectedPredictor(directory/"predictor_contact_108.dat");
    checkAdvancedShearPhase();
    continueFixture(directory/"predictor_contact_108.dat");
    continueFixture(directory/"adhesive_network_108.dat");
    continueFixture(directory/"normal_friction_corner_108.dat");
    continueFixture(directory/"normal_friction_corner_followup_108.dat");
    continueFixture(directory/"normal_release_metric_108.dat");
    continueFixture(directory/"normal_release_followup_108.dat");
    continueFixture(directory/"normal_mode_cycle_108.dat",true);
  }catch(const std::exception& error){std::cerr<<"FAIL: "<<error.what()<<'\n';result=1;}
  const auto finalize=PetscFinalize();return finalize?1:result;
}
