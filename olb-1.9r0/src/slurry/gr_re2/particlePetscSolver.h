#ifndef SLURRY_GR_RE2_PARTICLE_PETSC_SOLVER_H
#define SLURRY_GR_RE2_PARTICLE_PETSC_SOLVER_H
// Included inside graphite::particle_detail after the common particle kernels.
// PETSc owns globalization and Krylov iteration. The physical force laws and
// immutable beginning-of-substep contact history remain in Residual.

struct ParticleAttemptTrace {
  int iteration=0,kspIterations=0,activeContacts=0,slidingContacts=0,rollingContacts=0,activated=0,released=0,domainErrors=0,kspReason=0,candidateExpansion=0;
  int frictionBranchAttempts=0,frictionBranchCorrections=0;
  int contactStateUpdates=0,contactActivations=0,contactReleases=0;
  int nonlinearIteration=0,npcReason=0,totalKspIterations=0,residualEvaluations=0;
  double forceRatio=0.,torqueRatio=0.,gapViolation=0.,complementarityRatio=0.;
  double residualNorm=0.,fraction=1.,kspResidual=0.;
};
struct ParticleAttemptBundle {
  std::deque<ParticleAttemptTrace> trace;double outerTime=0.,outerDt=0.,subdt=0.;
  int count=0,substep=0,reason=0;std::string error;
};
struct ParticleDiagnosticScope;
inline ParticleDiagnosticScope*& pendingParticleDiagnostics() {
  static thread_local ParticleDiagnosticScope* scope=nullptr;return scope;
}
struct ParticleDiagnosticScope {
  const ParticleStepSettings& settings;ParticleDiagnosticScope* previous;
  std::deque<ParticleAttemptBundle> attempts;std::size_t rows=0;
  explicit ParticleDiagnosticScope(const ParticleStepSettings& s):settings(s),previous(pendingParticleDiagnostics()) {pendingParticleDiagnostics()=this;}
  ~ParticleDiagnosticScope(){pendingParticleDiagnostics()=previous;}
  void add(ParticleAttemptBundle attempt) {
    rows+=attempt.trace.size()+1;attempts.push_back(std::move(attempt));
    while(rows>2048 && attempts.size()>1){rows-=attempts.front().trace.size()+1;attempts.pop_front();}
  }
  void flushFailure();
};
inline std::string csvQuoted(const std::string& value) {
  std::string out="\"";for(char c:value){if(c=='\"')out+='\"';out+=c=='\n'?' ':c;}return out+'\"';
}
inline void writeAttemptTrace(const ParticleStepSettings& s,
    const std::deque<ParticleAttemptTrace>& trace,double outerTime,double outerDt,
    int count,int substep,double subdt,int reason,const std::string& error) {
  if(s.solverDiagnosticsPrefix.empty())return;
  if(auto* scope=pendingParticleDiagnostics()){scope->add({trace,outerTime,outerDt,subdt,count,substep,reason,error});return;}
  const std::string path=s.solverDiagnosticsPrefix+"_trace.csv";
  std::ifstream existing(path);const bool header=!existing.good()||existing.peek()==std::ifstream::traits_type::eof();
  std::ofstream file(path,std::ios::app);if(!file)return;
  if(header)file<<"outer_time_s,outer_dt_s,subdivision_count,substep_index,subdt_s,newton_iteration,force_ratio,torque_ratio,gap_violation_m,complementarity_ratio,residual_norm,step_fraction,ksp_iterations,ksp_residual_norm,snes_reason,status,message,active_contacts,sliding_contacts,rolling_contacts,activated_contacts,released_contacts,domain_errors,ksp_reason,candidate_expansion,friction_branch_attempts,friction_branch_corrections,contact_state_updates,contact_activations,contact_releases,nonlinear_iteration,npc_snes_reason,total_ksp_iterations,residual_evaluations,newton_step_fraction\n";
  file<<std::setprecision(17);
  for(const auto& r:trace)file<<outerTime<<','<<outerDt<<','<<count<<','<<substep<<','<<subdt<<','
    <<r.iteration<<','<<r.forceRatio<<','<<r.torqueRatio<<','<<r.gapViolation<<','<<r.complementarityRatio<<','
    <<r.residualNorm<<','<<r.fraction<<','<<r.kspIterations<<','<<r.kspResidual<<','<<reason<<",failed_iteration,"<<csvQuoted(error)<<','<<r.activeContacts<<','<<r.slidingContacts<<','<<r.rollingContacts
    <<','<<r.activated<<','<<r.released<<','<<r.domainErrors<<','<<r.kspReason<<','<<r.candidateExpansion
    <<','<<r.frictionBranchAttempts<<','<<r.frictionBranchCorrections
    <<','<<r.contactStateUpdates<<','<<r.contactActivations<<','<<r.contactReleases
    <<','<<r.nonlinearIteration<<','<<r.npcReason<<','<<r.totalKspIterations<<','<<r.residualEvaluations<<','<<r.fraction<<'\n';
  file<<outerTime<<','<<outerDt<<','<<count<<','<<substep<<','<<subdt
      <<",-1,0,0,0,0,0,0,0,0,"<<reason<<",failed_attempt,"<<csvQuoted(error)<<",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\n";
}
inline void ParticleDiagnosticScope::flushFailure() {
  auto* saved=pendingParticleDiagnostics();pendingParticleDiagnostics()=nullptr;
  for(const auto& a:attempts)writeAttemptTrace(settings,a.trace,a.outerTime,a.outerDt,a.count,a.substep,a.subdt,a.reason,a.error);
  pendingParticleDiagnostics()=saved;
}
inline void writeOuterOutcome(const ParticleStepSettings& s,double time,double dt,int count,bool success) {
  if(s.solverDiagnosticsPrefix.empty())return;
  std::ofstream f(s.solverDiagnosticsPrefix+"_outcomes.csv",std::ios::app);if(!f)return;
  if(f.tellp()==0)f<<"outer_time_s,outer_dt_s,subdivision_count,status\n";
  f<<std::setprecision(17)<<time<<','<<dt<<','<<count<<','<<(success?"recovered":"failed_outer_step")<<'\n';
}

// A versioned, dependency-free replay file contains the inputs of the final
// failed substep, including exact committed contact history and gap cache.
struct ParticleReplayInput {
  ParticleStepSettings settings;
  double dt=0.,time=0.,outerTime=0.,outerDt=0.;int count=1,substep=0;
  std::vector<Body> bodies;std::vector<Vec3> force,torque;std::vector<GapCache> cache;
  PersistentContactState contacts;
};
inline void writeParticleReplay(const ParticleReplayInput& x,const std::string& path) {
  std::ofstream f(path);if(!f)return;f<<std::setprecision(17)<<"GR_PARTICLE_REPLAY 1\n";
  const auto& s=x.settings;
  f<<x.dt<<' '<<x.time<<' '<<x.outerTime<<' '<<x.outerDt<<' '<<x.count<<' '<<x.substep<<'\n';
  f<<s.shearRate<<' '<<s.maxSubsteps<<' '<<s.maxNewtonIterations<<' '<<s.maxKrylovIterations<<' '<<s.maxLineSearch<<' '
    <<s.relativeTolerance<<' '<<s.forceAbsoluteTolerance<<' '<<s.torqueAbsoluteTolerance<<' '<<s.contactGapTolerance<<' '<<s.finiteDifferenceStep<<'\n';
  for(double a:s.box){f<<a<<' ';}f<<'\n';
  f<<s.pair.hamaker<<' '<<s.pair.sigma<<' '<<s.pair.switchGap<<' '<<s.pair.cutoffGap<<'\n';
  f<<s.nearField.viscosity<<' '<<s.nearField.matchingGap<<' '<<s.nearField.enabled<<' '<<s.nearField.tangential<<'\n';
  f<<s.rough.enabled<<' '<<s.rough.gap<<' '<<s.rough.friction<<' '<<s.rough.tangentialStiffness<<' '<<s.rough.rollingLength<<' '<<s.rough.rollingYieldAngle<<'\n';
  f<<x.bodies.size()<<'\n';
  for(std::size_t i=0;i<x.bodies.size();++i){const auto& b=x.bodies[i];
    for(const Vec3* a:{&b.position,&b.velocity,&b.omega,&b.axes})for(double v:*a)f<<v<<' ';
    for(double v:b.rotation){f<<v<<' ';}f<<b.mass<<' ';for(double v:b.inertiaBody){f<<v<<' ';}
    for(double v:x.force[i]){f<<v<<' ';}for(double v:x.torque[i]){f<<v<<' ';}f<<'\n';}
  f<<x.contacts.size()<<'\n';
  for(std::size_t i=0;i<x.contacts.size();++i){const auto& c=x.contacts[i];
    f<<c.active<<' '<<c.sliding<<' '<<c.rolling<<' ';
    for(const Vec3* a:{&c.normal,&c.elasticSlip,&c.elasticRoll})for(double v:*a)f<<v<<' ';
    for(double v:{c.rollingCap,c.rollingStiffness,c.normalLoad,c.elasticEnergy,c.plasticSlipWork,c.plasticRollWork,c.releasedEnergy,c.stepDuration})f<<v<<' ';
    for(const Vec3* a:{&c.tangentForce,&c.rollingTorque,&c.leverI,&c.leverJ})for(double v:*a)f<<v<<' ';
    const GapCache cache=i<x.cache.size()?x.cache[i]:GapCache{};
    f<<cache.valid<<' ';for(double v:cache.normal)f<<v<<' ';f<<'\n';}
}
inline ParticleReplayInput readParticleReplay(const std::string& path) {
  ParticleReplayInput x;std::ifstream f(path);std::string magic;int version=0;f>>magic>>version;
  if(magic!="GR_PARTICLE_REPLAY"||version!=1)throw std::runtime_error("Invalid particle replay header");
  auto& s=x.settings;f>>x.dt>>x.time>>x.outerTime>>x.outerDt>>x.count>>x.substep;
  f>>s.shearRate>>s.maxSubsteps>>s.maxNewtonIterations>>s.maxKrylovIterations>>s.maxLineSearch
    >>s.relativeTolerance>>s.forceAbsoluteTolerance>>s.torqueAbsoluteTolerance>>s.contactGapTolerance>>s.finiteDifferenceStep;
  for(double& a:s.box)f>>a;
  f>>s.pair.hamaker>>s.pair.sigma>>s.pair.switchGap>>s.pair.cutoffGap;
  f>>s.nearField.viscosity>>s.nearField.matchingGap>>s.nearField.enabled>>s.nearField.tangential;
  f>>s.rough.enabled>>s.rough.gap>>s.rough.friction>>s.rough.tangentialStiffness>>s.rough.rollingLength>>s.rough.rollingYieldAngle;
  std::size_t n=0;f>>n;if(n==0||n>100000)throw std::runtime_error("Unreasonable particle replay count");
  x.bodies.resize(n);x.force.resize(n);x.torque.resize(n);
  for(std::size_t i=0;i<n;++i){auto& b=x.bodies[i];
    for(Vec3* a:{&b.position,&b.velocity,&b.omega,&b.axes})for(double& v:*a)f>>v;
    for(double& v:b.rotation){f>>v;}f>>b.mass;for(double& v:b.inertiaBody){f>>v;}
    for(double& v:x.force[i]){f>>v;}for(double& v:x.torque[i]){f>>v;}}
  f>>n;if(n!=x.bodies.size()*(x.bodies.size()-1)/2)throw std::runtime_error("Invalid particle replay pair count");
  x.contacts.resize(n);x.cache.resize(n);
  for(std::size_t i=0;i<n;++i){auto& c=x.contacts[i];
    f>>c.active>>c.sliding>>c.rolling;
    for(Vec3* a:{&c.normal,&c.elasticSlip,&c.elasticRoll})for(double& v:*a)f>>v;
    f>>c.rollingCap>>c.rollingStiffness>>c.normalLoad>>c.elasticEnergy>>c.plasticSlipWork>>c.plasticRollWork>>c.releasedEnergy>>c.stepDuration;
    for(Vec3* a:{&c.tangentForce,&c.rollingTorque,&c.leverI,&c.leverJ})for(double& v:*a)f>>v;
    f>>x.cache[i].valid;for(double& v:x.cache[i].normal)f>>v;}
  if(!f)throw std::runtime_error("Truncated particle replay input");
  if(!(x.dt>0.)||!std::isfinite(x.dt)||!std::isfinite(x.time))throw std::runtime_error("Invalid particle replay time");
  for(std::size_t i=0;i<x.bodies.size();++i){const auto& b=x.bodies[i];
    if(!(b.mass>0.)||!std::isfinite(b.mass)||!finite(b.position)||!finite(b.velocity)||!finite(b.omega)||!finite(x.force[i])||!finite(x.torque[i]))
      throw std::runtime_error("Invalid particle replay body or load");
    for(double a:b.axes)if(!(a>0.)||!std::isfinite(a))throw std::runtime_error("Invalid particle replay semiaxis");
    for(double a:b.inertiaBody)if(!(a>0.)||!std::isfinite(a))throw std::runtime_error("Invalid particle replay inertia");}
  return x;
}

#ifdef SLURRY_USE_PETSC
static_assert(sizeof(PetscReal)==sizeof(double),"Graphite PETSc backend requires real double precision PETSc");
#ifdef PETSC_USE_COMPLEX
#error "Graphite particle solver requires a real-scalar PETSc build"
#endif
inline double fischerBurmeister(double a,double b) {
  const double r=std::hypot(a,b);
  if(a>=0.&&b>=0.)return r+a+b>0.?-2.*a*b/(r+a+b):0.;
  return r-a-b;
}
inline std::pair<double,double> fischerDerivative(double a,double b) {
  const double r=std::hypot(a,b);
  if(r>0.)return {a/r-1.,b/r-1.};
  const double corner=1./std::sqrt(2.)-1.;return {corner,corner};
}
struct NcpPreconditioner {
  BlockPreconditioner blocks;Vector bodyWeight;std::vector<double> gapFactor,loadFactor,schur;
  bool build(const Evaluation& e,const Residual& r,const Vector& weights) {
    if(!blocks.build(e,r))return false;
    bodyWeight=weights;
    gapFactor.clear();loadFactor.clear();schur.clear();
    for(const auto& c:blocks.contact) {
      const auto p=std::find_if(e.pairs.begin(),e.pairs.end(),[&](const PairLinearization& v){return v.slot==c.slot;});
      if(p==e.pairs.end())return false;
      const auto d=fischerDerivative((p->gap-r.settings.rough.gap)/r.settings.contactGapTolerance,
                                     e.normalLoads[p->index]/r.settings.forceAbsoluteTolerance);
      const double a=d.first*r.lengthScale/r.settings.contactGapTolerance;
      const double b=d.second*r.reactionScale/r.settings.forceAbsoluteTolerance;
      const double denominator=a*c.schur-b;
      if(!std::isfinite(denominator)||std::abs(denominator)<1.e-30)return false;
      gapFactor.push_back(a);loadFactor.push_back(b);schur.push_back(denominator);
    }return true;
  }
  Vector operator()(const Vector& r)const {
    Vector out(r.size(),0.),physical=r;
    for(std::size_t k=0;k<6*blocks.n;++k)physical[k]/=bodyWeight[k];
    for(std::size_t i=0;i<blocks.n;++i){Six x{};for(int k=0;k<6;++k)x[k]=physical[6*i+k];const auto y=multiply6(blocks.inverse[i],x);for(int k=0;k<6;++k)out[6*i+k]=y[k];}
    Vector corrected=physical;
    for(std::size_t j=0;j<blocks.contact.size();++j){const auto& c=blocks.contact[j];Six zi{},zj{};
      for(int k=0;k<6;++k){zi[k]=out[6*c.i+k];zj[k]=out[6*c.j+k];}
      const double lambda=(gapFactor[j]*(dot6(c.ci,zi)+dot6(c.cj,zj))-r[6*blocks.n+c.slot])/schur[j];
      out[6*blocks.n+c.slot]=lambda;
      for(int k=0;k<6;++k){corrected[6*c.i+k]-=c.bi[k]*lambda;corrected[6*c.j+k]-=c.bj[k]*lambda;}}
    for(std::size_t i=0;i<blocks.n;++i){Six x{};for(int k=0;k<6;++k)x[k]=corrected[6*i+k];const auto y=multiply6(blocks.inverse[i],x);for(int k=0;k<6;++k)out[6*i+k]=y[k];}
    return out;
  }
};
struct PetscParticleContext {
  Residual& residual;Vector weights,baseQ;Evaluation base;NcpPreconditioner pc;
  std::vector<unsigned char> engagement;std::deque<ParticleAttemptTrace> trace;
  std::string error;int totalKrylov=0,newtonAttempts=0,domainErrors=0,candidateExpansion=0,iterationBudget=0;
  int krylovOffset=0,evaluationOffset=0;bool missingCandidate=false,kspPending=false;
  SNES candidateSolver=nullptr;
  int frictionBranchAttempts=0,frictionBranchCorrections=0,frictionBranchKrylov=0;
  int frictionAttemptOffset=0,frictionCorrectionOffset=0;
  int iterationOffset=0,contactStateUpdates=0,contactActivations=0,contactReleases=0;
  std::vector<unsigned char> previousTraceContacts;
  std::vector<std::size_t> newCandidates;
  PetscParticleContext(Residual& r):residual(r),iterationBudget(r.settings.maxNewtonIterations){}
  void transform(const Vector& q,const Evaluation& e,Vector& out)const {
    out=e.residual;const std::size_t n=residual.old.size();
    for(std::size_t k=0;k<6*n;++k)out[k]*=weights[k];
    for(std::size_t p=0;p<residual.activeSlot.size();++p)if(residual.activeSlot[p]>=0)
      out[6*n+residual.activeSlot[p]]=fischerBurmeister((e.gaps[p]-residual.settings.rough.gap)/residual.settings.contactGapTolerance,
                  q[6*n+residual.activeSlot[p]]*residual.reactionScale/residual.settings.forceAbsoluteTolerance);
  }
  double complementarity(const Evaluation& e)const {
    double worst=0.;for(std::size_t p=0;p<residual.activeSlot.size();++p)if(residual.activeSlot[p]>=0)
      worst=std::max(worst,std::abs(fischerBurmeister((e.gaps[p]-residual.settings.rough.gap)/residual.settings.contactGapTolerance,
                                                   e.normalLoads[p]/residual.settings.forceAbsoluteTolerance)));
    return worst;
  }
  bool convergedInner(const Evaluation& e)const {
    // A retained contact is allowed to open in this fixed-regime subproblem.
    // Its finite rolling torque is removed only by the outer state update,
    // never from a residual evaluation midway through a line search.
    const auto& s=residual.settings;
    if(e.diagnostic.maxForceResidualRatio>1.||e.diagnostic.maxTorqueResidualRatio>1.
        ||e.diagnostic.contactGapViolation>s.contactGapTolerance||complementarity(e)>1.)return false;
    for(std::size_t p=0;p<residual.activeSlot.size();++p)if(residual.activeSlot[p]>=0) {
      const double load=e.normalLoads[p],gap=e.gaps[p]-s.rough.gap;
      if(load<0.||gap < -s.contactGapTolerance)return false;
      if(load>s.forceAbsoluteTolerance&&std::abs(gap)>s.contactGapTolerance)return false;
    }
    return true;
  }
  bool convergedPhysical(const Evaluation& e)const {
    if(!physicallyConverged(e,residual.settings)||complementarity(e)>1.)return false;
    const auto& s=residual.settings;
    for(std::size_t p=0;p<residual.activeSlot.size();++p)if(residual.activeSlot[p]>=0) {
      const double load=e.normalLoads[p],gap=e.gaps[p]-s.rough.gap;
      if(load<0. || gap < -s.contactGapTolerance)return false;
      if(load>s.forceAbsoluteTolerance && std::abs(gap)>s.contactGapTolerance)return false;
      if(gap>s.contactGapTolerance && e.contacts[p].active)return false;
    }return true;
  }
  bool evaluate(const Vector& q,Evaluation& e) {
    if(!residual(q,e,error)){++domainErrors;return false;}
    // Grow the candidate set only; complementarity handles release internally.
    // Expansion restarts PETSc with the same trial state and immutable history.
    if(residual.settings.rough.enabled)for(std::size_t p=0;p<residual.activeSlot.size();++p)
      if(residual.activeSlot[p]<0 && e.gaps[p]<residual.settings.rough.gap+residual.settings.contactGapTolerance) {
        missingCandidate=true;
        if(std::find(newCandidates.begin(),newCandidates.end(),p)==newCandidates.end())newCandidates.push_back(p);
      }
    return true;
  }
};
inline Vector petscRead(Vec v) {
  PetscInt n;VecGetLocalSize(v,&n);const PetscScalar* a=nullptr;VecGetArrayRead(v,&a);
  Vector x(a,a+n);VecRestoreArrayRead(v,&a);return x;
}
inline PetscErrorCode petscWrite(Vec v,const Vector& x) {
  PetscScalar* a=nullptr;PetscErrorCode ierr=VecGetArray(v,&a);if(ierr)return ierr;
  std::copy(x.begin(),x.end(),a);return VecRestoreArray(v,&a);
}
inline PetscErrorCode particlePetscFunction(SNES snes,Vec x,Vec f,void* pointer) {
  auto& c=*static_cast<PetscParticleContext*>(pointer);
  try {const Vector q=petscRead(x);Evaluation e;
    if(!c.evaluate(q,e)){SNESSetFunctionDomainError(snes);VecSet(f,0.);return 0;}
    Vector y;c.transform(q,e,y);return petscWrite(f,y);
  }catch(const std::exception& ex){c.error=ex.what();SNESSetFunctionDomainError(snes);VecSet(f,0.);return 0;}
}
inline bool particleNcpJacobian(PetscParticleContext& c,const Vector& direction,Vector& product) {
  auto& residual=c.residual;const auto& settings=residual.settings;const auto& q=c.baseQ;const auto& base=c.base;
  const auto& old=residual.old;const std::size_t n=old.size();const double dt=residual.dt,L=residual.lengthScale;
  double magnitude=0.;for(std::size_t k=0;k<6*n;++k)magnitude+=direction[k]*direction[k];magnitude=std::sqrt(magnitude);
  product.assign(q.size(),0.);
  if(magnitude>0.) {
    double bodyMagnitude=0.;for(std::size_t k=0;k<6*n;++k)bodyMagnitude+=q[k]*q[k];
    double epsilon=settings.finiteDifferenceStep*(1.+std::sqrt(bodyMagnitude))/magnitude;
    Vector trial=q;Evaluation displaced;bool valid=false;
    struct RestoreEngagement {
      Residual& residual;const std::vector<unsigned char>* previous;
      RestoreEngagement(Residual& value,const std::vector<unsigned char>* fixed):residual(value),previous(value.engagementOverride) {
        residual.engagementOverride=fixed;
      }
      ~RestoreEngagement(){residual.engagementOverride=previous;}
    } restoreEngagement(residual,&c.engagement);
    for(int attempt=0;attempt<8;++attempt) {
      for(std::size_t k=0;k<6*n;++k)trial[k]=q[k]+epsilon*direction[k];
      if(residual(trial,displaced,c.error)){valid=true;break;}epsilon*=-.5;
    }
    if(!valid)return false;
    for(std::size_t k=0;k<q.size();++k)product[k]=(displaced.residual[k]-base.residual[k])/epsilon;
    auto derivative=[](const Vec3& trial,const Vec3& change,double stiffness,double cap,double dk,double dc,bool yielded)->Vec3 {
      const double r=norm(trial);if(cap<=0.)return r>0.?scale(trial,-dc/r):Vec3{};
      // Use the branch selected by the actual return map.  Recomputing k*|s|
      // here can round differently from |k*s| exactly at a yield boundary.
      if(!yielded)return add(scale(change,-stiffness),scale(trial,-dk));
      const Vec3 v=scale(trial,1./r);return add(scale(sub(change,scale(v,dot(v,change))),-cap/r),scale(v,-dc));
    };
    for(const auto& p:base.pairs)if(p.slot>=0 && p.contact.active) {
      const auto found=std::find_if(displaced.pairs.begin(),displaced.pairs.end(),[&](const PairLinearization& a){return a.index==p.index;});
      if(found==displaced.pairs.end()){c.error="Candidate missing from geometric derivative";return false;}
      const auto& contact=p.contact;const auto& next=found->contact;
      const Vec3 ds=scale(sub(next.trialSlip,contact.trialSlip),1./epsilon),dr=scale(sub(next.trialRoll,contact.trialRoll),1./epsilon);
      const Vec3 dn=scale(sub(found->normal,p.normal),1./epsilon);
      const Vec3 dli=scale(sub(found->leverI,p.leverI),1./epsilon),dlj=scale(sub(found->leverJ,p.leverJ),1./epsilon);
      const auto& state=contact.candidateState;const auto& ns=next.candidateState;
      const Vec3 dft=derivative(contact.trialSlip,ds,settings.rough.tangentialStiffness,settings.rough.friction*std::max(0.,state.normalLoad),0.,0.,contact.sliding);
      const Vec3 dmr=derivative(contact.trialRoll,dr,state.rollingStiffness,state.rollingCap,
            (ns.rollingStiffness-state.rollingStiffness)/epsilon,(ns.rollingCap-state.rollingCap)/epsilon,contact.rolling);
      const Vec3 df=add(scale(dn,-state.normalLoad),dft);
      const Vec3 dti=add(add(cross(dli,contact.forceI),cross(p.leverI,df)),dmr);
      const Vec3 dtj=scale(add(add(cross(dlj,contact.forceI),cross(p.leverJ,df)),dmr),-1.);
      const Vec3 cf=sub(df,scale(sub(next.forceI,contact.forceI),1./epsilon));
      const Vec3 cti=sub(dti,scale(sub(next.torqueI,contact.torqueI),1./epsilon));
      const Vec3 ctj=sub(dtj,scale(sub(next.torqueJ,contact.torqueJ),1./epsilon));
      const double fi=old[p.i].mass*L/(dt*dt),fj=old[p.j].mass*L/(dt*dt);
      const double ti=*std::max_element(old[p.i].inertiaBody.begin(),old[p.i].inertiaBody.end())/(dt*dt);
      const double tj=*std::max_element(old[p.j].inertiaBody.begin(),old[p.j].inertiaBody.end())/(dt*dt);
      for(int k=0;k<3;++k){product[6*p.i+k]-=cf[k]/fi;product[6*p.j+k]+=cf[k]/fj;
        product[6*p.i+k+3]-=cti[k]/ti;product[6*p.j+k+3]-=ctj[k]/tj;}
    }
  }
  for(const auto& block:c.pc.blocks.contact) {
    const double d=direction[6*n+block.slot];
    for(int k=0;k<6;++k){product[6*block.i+k]+=block.bi[k]*d;product[6*block.j+k]+=block.bj[k]*d;}
  }
  for(std::size_t k=0;k<6*n;++k)product[k]*=c.weights[k];
  for(std::size_t j=0;j<c.pc.blocks.contact.size();++j){const auto& b=c.pc.blocks.contact[j];
    product[6*n+b.slot]=c.pc.gapFactor[j]*product[6*n+b.slot]+c.pc.loadFactor[j]*direction[6*n+b.slot];}
  return true;
}
inline PetscErrorCode particlePetscMatMult(Mat a,Vec x,Vec y) {
  void* pointer=nullptr;MatShellGetContext(a,&pointer);auto& c=*static_cast<PetscParticleContext*>(pointer);
  try {Vector product;if(!particleNcpJacobian(c,petscRead(x),product))return PETSC_ERR_USER;
    return petscWrite(y,product);
  }catch(const std::exception& ex){c.error=ex.what();return PETSC_ERR_USER;}
}
inline PetscErrorCode particlePetscPcApply(PC pc,Vec x,Vec y) {
  void* pointer=nullptr;PCShellGetContext(pc,&pointer);auto& c=*static_cast<PetscParticleContext*>(pointer);
  try{return petscWrite(y,c.pc(petscRead(x)));}catch(const std::exception& ex){c.error=ex.what();return PETSC_ERR_USER;}
}
inline PetscErrorCode particlePetscJacobian(SNES,Vec x,Mat,Mat,void* pointer) {
  auto& c=*static_cast<PetscParticleContext*>(pointer);
  try {c.baseQ=petscRead(x);if(!c.evaluate(c.baseQ,c.base))return PETSC_ERR_USER;
    c.engagement.assign(c.residual.activeSlot.size(),0);
    for(std::size_t p=0;p<c.engagement.size();++p)c.engagement[p]=c.base.contacts[p].active;
    if(!c.pc.build(c.base,c.residual,c.weights)){c.error="PETSc contact block preconditioner is singular";return PETSC_ERR_USER;}
    return 0;
  }catch(const std::exception& ex){c.error=ex.what();return PETSC_ERR_USER;}
}

// NGMRES accelerates the one-Newton candidate map using previously evaluated
// residuals. No extra yield-boundary predictor is used: every attempted Newton
// direction, including a failed KSP solve, consumes the same shared budget.
inline PetscErrorCode particlePetscKspPreSolve(KSP,Vec,Vec,void* pointer) {
  auto& c=*static_cast<PetscParticleContext*>(pointer);
  if(c.newtonAttempts>=c.iterationBudget){c.error="Particle Newton budget exhausted before KSP";return PETSC_ERR_NOT_CONVERGED;}
  ++c.newtonAttempts;c.kspPending=true;return 0;
}
inline PetscErrorCode particlePetscKspPostSolve(KSP ksp,Vec,Vec,void* pointer) {
  auto& c=*static_cast<PetscParticleContext*>(pointer);
  if(c.kspPending){PetscInt iterations=0;KSPGetIterationNumber(ksp,&iterations);c.totalKrylov+=iterations;c.kspPending=false;}
  return 0;
}
inline PetscErrorCode particlePetscCandidateConverged(SNES,PetscInt iteration,PetscReal,PetscReal,
                                                     PetscReal,SNESConvergedReason* reason,void* pointer) {
  auto& c=*static_cast<PetscParticleContext*>(pointer);
  // MAX_IT is the normal return of the one-step nonlinear preconditioner.
  // Only the outer callback and final physical checks may accept a substep.
  *reason=c.missingCandidate?SNES_DIVERGED_FUNCTION_DOMAIN:
          (iteration>=1?SNES_DIVERGED_MAX_IT:SNES_CONVERGED_ITERATING);
  return 0;
}
inline PetscErrorCode particlePetscNgmresUpdate(SNES snes,PetscInt) {
  void* pointer=nullptr;PetscErrorCode ierr=SNESGetApplicationContext(snes,&pointer);if(ierr)return ierr;
  auto& c=*static_cast<PetscParticleContext*>(pointer);SNES npc=nullptr;
  ierr=SNESGetNPC(snes,&npc);if(ierr)return ierr;
  const char* type=nullptr;ierr=SNESGetType(npc,&type);if(ierr)return ierr;
  PCSide side;ierr=SNESGetNPCSide(snes,&side);if(ierr)return ierr;
  PetscBool nested=PETSC_FALSE;ierr=SNESHasNPC(npc,&nested);if(ierr)return ierr;
  if(!type||std::string(type)!=SNESNEWTONLS||side!=PC_RIGHT||nested){
    c.error="Graphite NGMRES requires one right NEWTONLS nonlinear preconditioner without nesting";return PETSC_ERR_SUP;}
  // SNESSetUp applies NPC options again. This hook runs afterwards, before
  // EVERY candidate solve, so user max_it options cannot multiply the budget.
  ierr=SNESSetTolerances(npc,0.,0.,0.,1,100000);if(ierr)return ierr;
  ierr=SNESSetConvergenceTest(npc,particlePetscCandidateConverged,&c,nullptr);if(ierr)return ierr;
  ierr=SNESSetNormSchedule(npc,SNES_NORM_ALWAYS);if(ierr)return ierr;
  ierr=SNESSetFunctionType(npc,SNES_FUNCTION_UNPRECONDITIONED);if(ierr)return ierr;
  SNESLineSearch line;PetscInt lineMax;
  ierr=SNESGetLineSearch(npc,&line);if(ierr)return ierr;
  ierr=SNESLineSearchGetTolerances(line,nullptr,nullptr,nullptr,nullptr,nullptr,&lineMax);if(ierr)return ierr;
  ierr=SNESLineSearchSetTolerances(line,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,std::min<PetscInt>(lineMax,c.residual.settings.maxLineSearch));if(ierr)return ierr;
  KSP ksp;PetscReal rtol,atol,dtol;PetscInt maxit;
  ierr=SNESGetKSP(npc,&ksp);if(ierr)return ierr;
  ierr=KSPGetTolerances(ksp,&rtol,&atol,&dtol,&maxit);if(ierr)return ierr;
  ierr=KSPSetTolerances(ksp,rtol,atol,dtol,std::min<PetscInt>(maxit,c.residual.settings.maxKrylovIterations));if(ierr)return ierr;
  ierr=KSPSetPreSolve(ksp,particlePetscKspPreSolve,&c);if(ierr)return ierr;
  return KSPSetPostSolve(ksp,particlePetscKspPostSolve,&c);
}
inline PetscErrorCode particlePetscConverged(SNES snes,PetscInt iteration,PetscReal,PetscReal,
                                            PetscReal fnorm,SNESConvergedReason* reason,void* pointer) {
  auto& c=*static_cast<PetscParticleContext*>(pointer);*reason=SNES_CONVERGED_ITERATING;
  if(!std::isfinite(fnorm)){
#if PETSC_VERSION_GE(3,25,0)
    *reason=SNES_DIVERGED_FUNCTION_NANORINF;
#else
    *reason=SNES_DIVERGED_FNORM_NAN;
#endif
    return 0;}
  Vec solution;SNESGetSolution(snes,&solution);Vector q=petscRead(solution);Evaluation e;
  if(!c.evaluate(q,e)){*reason=SNES_DIVERGED_FUNCTION_DOMAIN;return 0;}
  if(c.missingCandidate){*reason=SNES_DIVERGED_FUNCTION_DOMAIN;return 0;}
  // At final acceptance only, remove roundoff-sized negative multipliers and
  // re-evaluate every physical residual. This cannot bypass a force criterion.
  bool changed=false;const std::size_t n=c.residual.old.size();
  for(std::size_t k=6*n;k<q.size();++k)if(q[k]<0. && -q[k]*c.residual.reactionScale<=c.residual.settings.forceAbsoluteTolerance){q[k]=0.;changed=true;}
  if(changed){Evaluation projected;if(c.evaluate(q,projected)&&c.convergedInner(projected)){
      petscWrite(solution,q);c.baseQ=q;c.base=std::move(projected);*reason=SNES_CONVERGED_FNORM_ABS;return 0;}}
  if(!changed && c.convergedInner(e)){c.baseQ=q;c.base=std::move(e);*reason=SNES_CONVERGED_FNORM_ABS;return 0;}
  if(c.newtonAttempts>=c.iterationBudget)*reason=SNES_DIVERGED_MAX_IT;
  return 0;
}
inline PetscErrorCode particlePetscMonitor(SNES snes,PetscInt iteration,PetscReal fnorm,void* pointer) {
  auto& c=*static_cast<PetscParticleContext*>(pointer);Vec solution;SNESGetSolution(snes,&solution);Evaluation e;
  if(!c.evaluate(petscRead(solution),e))return 0;
  ParticleAttemptTrace t;t.iteration=c.iterationOffset+c.newtonAttempts;
  t.nonlinearIteration=iteration;t.totalKspIterations=c.krylovOffset+c.totalKrylov;
  t.residualEvaluations=c.evaluationOffset+c.residual.evaluations;
  if(c.candidateSolver){SNESConvergedReason npcReason;SNESGetConvergedReason(c.candidateSolver,&npcReason);t.npcReason=static_cast<int>(npcReason);}
  t.forceRatio=e.diagnostic.maxForceResidualRatio;t.torqueRatio=e.diagnostic.maxTorqueResidualRatio;
  t.gapViolation=e.diagnostic.contactGapViolation;t.complementarityRatio=c.complementarity(e);Vector transformed;c.transform(petscRead(solution),e,transformed);t.residualNorm=length(transformed);
  t.activeContacts=e.diagnostic.contacts;t.slidingContacts=e.diagnostic.slidingContacts;t.rollingContacts=e.diagnostic.rollingContacts;t.domainErrors=c.domainErrors;t.candidateExpansion=c.candidateExpansion;
  if(c.previousTraceContacts.empty()){c.previousTraceContacts.resize(e.contacts.size());
    for(std::size_t p=0;p<e.contacts.size();++p)c.previousTraceContacts[p]=c.residual.startingContacts[p].active;}
  for(std::size_t p=0;p<e.contacts.size();++p){const bool on=e.contacts[p].active;
    t.activated+=on&&!c.previousTraceContacts[p];t.released+=!on&&c.previousTraceContacts[p];c.previousTraceContacts[p]=on;}
  SNES candidate=c.candidateSolver?c.candidateSolver:snes;
  // step_fraction is retained for old readers; both fraction fields describe
  // the Newton candidate search, not NGMRES's additive combination search.
  SNESLineSearch line;SNESGetLineSearch(candidate,&line);SNESLineSearchGetLambda(line,&t.fraction);
  KSP ksp;SNESGetKSP(candidate,&ksp);PetscInt count=0;KSPGetIterationNumber(ksp,&count);t.kspIterations=count;KSPGetResidualNorm(ksp,&t.kspResidual);
  KSPConvergedReason linearReason;KSPGetConvergedReason(ksp,&linearReason);t.kspReason=static_cast<int>(linearReason);
  t.frictionBranchAttempts=c.frictionAttemptOffset+c.frictionBranchAttempts;
  t.frictionBranchCorrections=c.frictionCorrectionOffset+c.frictionBranchCorrections;
  t.contactStateUpdates=c.contactStateUpdates;t.contactActivations=c.contactActivations;t.contactReleases=c.contactReleases;
  c.trace.push_back(t);if(c.trace.size()>256)c.trace.pop_front();return 0;
}
struct PetscParticleObjects {
  SNES snes=nullptr;Vec x=nullptr,f=nullptr;Mat jacobian=nullptr;
  ~PetscParticleObjects(){if(snes)SNESDestroy(&snes);if(jacobian)MatDestroy(&jacobian);if(f)VecDestroy(&f);if(x)VecDestroy(&x);}
};
inline bool implicitStepPetsc(const std::vector<Body>& old,const std::vector<Vec3>& force,
    const std::vector<Vec3>& torque,double dt,double time,const ParticleStepSettings& settings,
    std::vector<GapCache>& cache,PersistentContactState& contacts,std::vector<Body>& output,
    ParticleStepDiagnostics& diagnostic,std::string& error,double outerTime,double outerDt,int count,int substep) {
  error.clear();
  PetscBool initialized=PETSC_FALSE;PetscInitialized(&initialized);
  if(!initialized){error="PETSc particle solver requires PetscInitialize before advanceParticles";
    writeAttemptTrace(settings,{},outerTime,outerDt,count,substep,dt,0,error);return false;}
  if(settings.maxLineSearch<1){error="PETSc particle solver requires maxLineSearch >= 1";
    writeAttemptTrace(settings,{},outerTime,outerDt,count,substep,dt,0,error);return false;}
  const std::size_t n=old.size(),np=n*(n-1)/2;double L=0.;for(const auto& b:old)L=std::max(L,radius(b));
  const double reactionScale=old.front().mass*L/(dt*dt);
  std::vector<int> slots(np,-1);int activeCount=0;Vector q(6*n,0.);
  for(std::size_t i=0;i<n;++i)for(int k=0;k<3;++k){q[6*i+k]=old[i].velocity[k]*dt/L;q[6*i+k+3]=old[i].omega[k]*dt;}
  if(settings.rough.enabled){std::size_t p=0;
    for(std::size_t i=0;i<n;++i)for(std::size_t j=i+1;j<n;++j,++p){
      const Body image=closestImage(old[i],old[j],time,settings);
      const double range=std::max(settings.pair.cutoffGap,settings.nearField.matchingGap);
      if(norm(sub(image.position,old[i].position))>radius(old[i])+radius(image)+range&&!contacts[p].active)continue;
      const auto gap=closestEllipsoidGap(old[i],image,&cache[p]);
      if(gap.gap<settings.rough.gap-settings.contactGapTolerance){error="Initial accepted particle state violates rough contact gap";
        writeAttemptTrace(settings,{},outerTime,outerDt,count,substep,dt,0,error);return false;}
      const Vec3 vi=add(old[i].velocity,cross(old[i].omega,gap.leverI)),vj=add(image.velocity,cross(image.omega,gap.leverJ));
      const double predicted=gap.gap+dt*dot(sub(vj,vi),gap.normal);
      if(contacts[p].active || gap.gap<=std::max(settings.nearField.matchingGap,4.*settings.rough.gap)
            || predicted<=settings.rough.gap+settings.contactGapTolerance)slots[p]=activeCount++;
    }}
  q.resize(6*n+activeCount,0.);
  for(std::size_t p=0;p<np;++p)if(slots[p]>=0&&contacts[p].active)q[6*n+slots[p]]=std::max(0.,contacts[p].normalLoad)/reactionScale;
  std::vector<unsigned char> engagement(np,0);
  for(std::size_t p=0;p<np;++p)engagement[p]=slots[p]>=0&&contacts[p].active;
  auto previousTraceEngagement=engagement;
  std::deque<ParticleAttemptTrace> allTrace;int totalNewton=0,totalKrylov=0,totalEvaluations=0;int finalReason=0;
  int totalFrictionAttempts=0,totalFrictionCorrections=0;
  int stateUpdates=0,activations=0,releases=0,candidateExpansions=0;
  // Each inner solve owns one immutable engagement regime. Opening/closing
  // restarts the merit function only outside SNES, with the same beginning-of-
  // substep history. Neither a line search nor a derivative may release history.
  // Candidate growth and state restarts are bounded independently, while all
  // Newton work shares the original configured iteration budget.
  for(std::size_t outer=0;outer<=np+static_cast<std::size_t>(settings.maxNewtonIterations);++outer) {
    if(totalNewton>=settings.maxNewtonIterations) {
      std::ostringstream message;message<<"PETSc contact solve exhausted shared Newton budget across contact states: Newton="
        <<totalNewton<<", state updates="<<stateUpdates<<", activations="<<activations<<", releases="<<releases;
      error=message.str();finalReason=SNES_DIVERGED_MAX_IT;break;
    }
    Residual residual{old,force,torque,settings,cache,contacts,slots,dt,time,L,reactionScale,activeCount};residual.complementarity=true;
    residual.engagementOverride=&engagement;
    PetscParticleContext context(residual);context.candidateExpansion=candidateExpansions;
    context.iterationBudget=settings.maxNewtonIterations-totalNewton;context.iterationOffset=totalNewton;
    context.krylovOffset=totalKrylov;context.evaluationOffset=totalEvaluations;
    context.contactStateUpdates=stateUpdates;context.contactActivations=activations;context.contactReleases=releases;
    context.frictionAttemptOffset=totalFrictionAttempts;context.frictionCorrectionOffset=totalFrictionCorrections;
    context.previousTraceContacts=previousTraceEngagement;
    Evaluation initial;bool valid=false;
    for(int attempt=0;attempt<12;++attempt){if(context.evaluate(q,initial)){valid=true;break;}for(std::size_t k=0;k<6*n;++k)q[k]*=.5;}
    if(!valid){for(std::size_t k=0;k<6*n;++k)q[k]=0.;valid=context.evaluate(q,initial);}
    if(!valid){error=context.error;break;}
    if(context.missingCandidate){for(auto p:context.newCandidates)if(slots[p]<0)slots[p]=activeCount++;
      q.resize(6*n+activeCount,0.);++candidateExpansions;continue;}
    // Freeze physical tolerance row scales for this solve so the line search
    // compares one merit function. Final acceptance still uses current loads.
    context.weights.resize(6*n);
    for(std::size_t i=0;i<n;++i){
      const double fs=old[i].mass*L/(dt*dt),ts=*std::max_element(old[i].inertiaBody.begin(),old[i].inertiaBody.end())/(dt*dt);
      const double fw=fs/(settings.forceAbsoluteTolerance+settings.relativeTolerance*initial.forceReference[i]);
      const double tw=ts/(settings.torqueAbsoluteTolerance+settings.relativeTolerance*initial.torqueReference[i]);
      for(int k=0;k<3;++k){context.weights[6*i+k]=fw;context.weights[6*i+k+3]=tw;}}
    PetscParticleObjects objects;PetscErrorCode ierr=0;
    auto checked=[&](PetscErrorCode code){if(code){std::ostringstream msg;msg<<"PETSc setup failed with error "<<code;throw std::runtime_error(msg.str());}};
    checked(VecCreateSeq(PETSC_COMM_SELF,static_cast<PetscInt>(q.size()),&objects.x));
    if(ierr){error="PETSc failed to allocate particle vector";break;}
    checked(VecDuplicate(objects.x,&objects.f));checked(petscWrite(objects.x,q));
    checked(SNESCreate(PETSC_COMM_SELF,&objects.snes));checked(SNESSetOptionsPrefix(objects.snes,"gr_"));
    checked(SNESSetType(objects.snes,SNESNGMRES));
    checked(SNESSetNPCSide(objects.snes,PC_RIGHT));
    checked(SNESSetFunctionType(objects.snes,SNES_FUNCTION_UNPRECONDITIONED));
    checked(SNESNGMRESSetSelectType(objects.snes,SNES_NGMRES_SELECT_LINESEARCH));
    checked(SNESSetApplicationContext(objects.snes,&context));
    checked(SNESSetUpdate(objects.snes,particlePetscNgmresUpdate));
    SNES npc;checked(SNESGetNPC(objects.snes,&npc));context.candidateSolver=npc;
    checked(SNESSetType(npc,SNESNEWTONLS));
    checked(SNESSetTolerances(npc,0.,0.,0.,1,100000));
    checked(SNESSetConvergenceTest(npc,particlePetscCandidateConverged,&context,nullptr));
    checked(SNESSetFunction(objects.snes,objects.f,particlePetscFunction,&context));
    checked(MatCreateShell(PETSC_COMM_SELF,q.size(),q.size(),q.size(),q.size(),&context,&objects.jacobian));
    checked(MatShellSetOperation(objects.jacobian,MATOP_MULT,reinterpret_cast<void(*)(void)>(particlePetscMatMult)));
    checked(SNESSetJacobian(objects.snes,objects.jacobian,objects.jacobian,particlePetscJacobian,&context));
    checked(SNESSetTolerances(objects.snes,0.,0.,0.,context.iterationBudget,100000));
    checked(SNESSetConvergenceTest(objects.snes,particlePetscConverged,&context,nullptr));
    checked(SNESSetNormSchedule(objects.snes,SNES_NORM_ALWAYS));
    if(!settings.solverDiagnosticsPrefix.empty())checked(SNESMonitorSet(objects.snes,particlePetscMonitor,&context,nullptr));
    SNESLineSearch line;checked(SNESGetLineSearch(npc,&line));
    // At friction corners, monotone backtracking can shrink a locally correct
    // Newton direction to repeated microscopic steps. PETSc's secant search
    // samples the residual norm along that direction. Use its standard single
    // secant iteration, within the existing line-search work cap, and bound the
    // search to a full Newton step. Final physical acceptance is unchanged.
#if PETSC_VERSION_GE(3,24,0)
    checked(SNESLineSearchSetType(line,SNESLINESEARCHSECANT));
#else
    checked(SNESLineSearchSetType(line,SNESLINESEARCHL2));
#endif
    checked(SNESLineSearchSetTolerances(line,PETSC_DEFAULT,1.,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,std::min(1,settings.maxLineSearch)));
    KSP ksp;PC pc;checked(SNESGetKSP(npc,&ksp));checked(KSPSetType(ksp,KSPGMRES));
    checked(KSPGMRESSetRestart(ksp,settings.maxKrylovIterations));checked(KSPSetTolerances(ksp,.05,1.e-14,PETSC_DEFAULT,settings.maxKrylovIterations));
    checked(KSPSetPCSide(ksp,PC_RIGHT));checked(KSPSetNormType(ksp,KSP_NORM_UNPRECONDITIONED));
    checked(KSPGetPC(ksp,&pc));checked(PCSetType(pc,PCSHELL));checked(PCShellSetContext(pc,&context));checked(PCShellSetApply(pc,particlePetscPcApply));
    checked(SNESSetFromOptions(objects.snes));
    // Options may tune the candidate line search and linear solve, but may
    // not replace this architecture or approximate the physical residual.
    const char* type=nullptr;checked(SNESGetType(objects.snes,&type));
    if(!type||std::string(type)!=SNESNGMRES){error="Graphite contact backend requires -gr_snes_type ngmres";break;}
    checked(SNESGetType(npc,&type));
    PCSide npcSide;checked(SNESGetNPCSide(objects.snes,&npcSide));
    if(!type||std::string(type)!=SNESNEWTONLS||npcSide!=PC_RIGHT){error="Graphite NGMRES requires a right NEWTONLS nonlinear preconditioner";break;}
    PetscBool approximate=PETSC_FALSE,selected=PETSC_FALSE;char selector[32]={};
    checked(PetscOptionsGetBool(nullptr,"gr_","-snes_ngmres_approxfunc",&approximate,nullptr));
    checked(PetscOptionsGetString(nullptr,"gr_","-snes_ngmres_select_type",selector,sizeof(selector),&selected));
    if(approximate||(selected&&std::string(selector)!="linesearch")){
      error="Graphite NGMRES requires exact residuals and -gr_snes_ngmres_select_type linesearch";break;}
    checked(SNESSetFunctionType(objects.snes,SNES_FUNCTION_UNPRECONDITIONED));
    checked(SNESSetConvergenceTest(objects.snes,particlePetscConverged,&context,nullptr));
    checked(SNESSetNormSchedule(objects.snes,SNES_NORM_ALWAYS));
    checked(SNESSetTolerances(objects.snes,0.,0.,0.,context.iterationBudget,100000));
    checked(SNESSetUpdate(objects.snes,particlePetscNgmresUpdate));
    checked(SNESGetLineSearch(npc,&line));
    const char* lineType=nullptr;PetscInt lineMaxIterations=0;
    checked(SNESLineSearchGetType(line,&lineType));
    checked(SNESLineSearchGetTolerances(line,nullptr,nullptr,nullptr,nullptr,nullptr,&lineMaxIterations));
    checked(SNESLineSearchSetPreCheck(line,nullptr,nullptr));
    if(!ierr){PetscPushErrorHandler(PetscReturnErrorHandler,nullptr);ierr=SNESSolve(objects.snes,nullptr,objects.x);PetscPopErrorHandler();}
    SNESConvergedReason reason=SNES_CONVERGED_ITERATING;SNESGetConvergedReason(objects.snes,&reason);finalReason=static_cast<int>(reason);
    PetscInt iterations=0;SNESGetIterationNumber(objects.snes,&iterations);
    checked(SNESLineSearchGetType(line,&lineType));
    checked(SNESLineSearchGetTolerances(line,nullptr,nullptr,nullptr,nullptr,nullptr,&lineMaxIterations));
    // An error can return before KSP's post-solve hook; charge its work once.
    checked(particlePetscKspPostSolve(ksp,nullptr,nullptr,&context));
    if(!settings.solverDiagnosticsPrefix.empty()
        &&(ierr||reason<0||context.trace.empty()||context.trace.back().iteration!=context.iterationOffset+context.newtonAttempts))
      checked(particlePetscMonitor(objects.snes,iterations,0.,&context));
    totalNewton+=context.newtonAttempts;totalKrylov+=context.totalKrylov;totalEvaluations+=residual.evaluations;
    totalFrictionAttempts+=context.frictionBranchAttempts;totalFrictionCorrections+=context.frictionBranchCorrections;
    previousTraceEngagement=engagement;
    q=petscRead(objects.x);
    for(const auto& row:context.trace){allTrace.push_back(row);if(allTrace.size()>256)allTrace.pop_front();}
    if(context.missingCandidate){for(auto p:context.newCandidates)if(slots[p]<0)slots[p]=activeCount++;
      q.resize(6*n+activeCount,0.);++candidateExpansions;continue;}
    Evaluation final;
    if(!ierr && reason>0 && context.evaluate(q,final) && context.convergedInner(final)) {
      // An OFF candidate's tolerance-sized normal reaction is unresolved, not
      // evidence for recreating a friction/rolling contact. Project it before
      // deciding the regime, and verify momentum and complementarity again.
      // This projection never removes a resolved compressive reaction.
      bool projectedOpenReaction=false;
      for(std::size_t p=0;p<np;++p)if(slots[p]>=0&&!engagement[p]
          &&q[6*n+slots[p]]!=0.&&std::abs(final.normalLoads[p])<=settings.forceAbsoluteTolerance) {
        q[6*n+slots[p]]=0.;projectedOpenReaction=true;
      }
      if(projectedOpenReaction) {
        if(!context.evaluate(q,final)) {error=context.error;finalReason=SNES_DIVERGED_FUNCTION_DOMAIN;break;}
        if(!context.convergedInner(final))continue;
      }
      auto desired=engagement;int added=0,removed=0;
      for(std::size_t p=0;p<np;++p) {
        // Retain the CURRENT trial regime in the existing gap/force tolerance
        // band. Using committed history here resurrects an already released
        // spring at N=0 and can create an ON/OFF cycle with two valid endpoints.
        // Committed history remains immutable and is still the constitutive
        // input if a resolved compressive reaction requires reactivation.
        desired[p]=slots[p]>=0 && final.gaps[p]<=settings.rough.gap+settings.contactGapTolerance
                   && (engagement[p]||final.normalLoads[p]>settings.forceAbsoluteTolerance);
        added+=desired[p]&&!engagement[p];removed+=!desired[p]&&engagement[p];
      }
      if(added||removed) {
        ++stateUpdates;activations+=added;releases+=removed;
        if(stateUpdates>=settings.maxNewtonIterations) {
          std::ostringstream message;message<<"PETSc contact-state updates exhausted the configured iteration budget: updates="
            <<stateUpdates<<", activations="<<activations<<", releases="<<releases<<", Newton="<<totalNewton;
          error=message.str();finalReason=SNES_DIVERGED_MAX_IT;break;
        }
        // A repeated mask alone is not a repeated nonlinear state. Continue
        // under the shared Newton/state-update bounds rather than aborting the
        // first time a contact set is revisited with different positions/loads.
        engagement=std::move(desired);continue;
      }
      // An open candidate has exactly zero committed normal reaction. Project
      // only multipliers within the existing force tolerance, then re-evaluate
      // every physical residual; never erase a finite reaction after acceptance.
      bool projected=false,largeOpenReaction=false;
      for(std::size_t p=0;p<np;++p)if(slots[p]>=0&&!engagement[p]&&q[6*n+slots[p]]!=0.) {
        if(std::abs(final.normalLoads[p])>settings.forceAbsoluteTolerance){largeOpenReaction=true;break;}
        q[6*n+slots[p]]=0.;projected=true;
      }
      if(!largeOpenReaction&&(!projected||context.evaluate(q,final))&&context.convergedPhysical(final)) {
        output=std::move(final.bodies);cache=std::move(final.cache);contacts=std::move(final.contacts);diagnostic=final.diagnostic;
        diagnostic.newtonIterations=totalNewton;diagnostic.krylovIterations=totalKrylov;diagnostic.residualEvaluations=totalEvaluations;
        diagnostic.frictionBranchAttempts=totalFrictionAttempts;diagnostic.frictionBranchCorrections=totalFrictionCorrections;
        diagnostic.contactStateUpdates=stateUpdates;diagnostic.contactActivations=activations;diagnostic.contactReleases=releases;
        for(auto& b:output){wrap(b,time+dt,settings);}
        return true;
      }
      if(projected&&!largeOpenReaction)continue;
    }
    KSPConvergedReason lastLinearReason=KSP_CONVERGED_ITERATING;KSPGetConvergedReason(ksp,&lastLinearReason);
    SNESConvergedReason npcReason=SNES_CONVERGED_ITERATING;SNESGetConvergedReason(npc,&npcReason);
    std::ostringstream message;message<<"PETSc SNES failed: solver=ngmres, candidate=newtonls, selector=linesearch, reason="<<finalReason<<" ("<<SNESConvergedReasons[reason]<<")"
      <<", NPC="<<static_cast<int>(npcReason)<<" ("<<SNESConvergedReasons[npcReason]<<")"
      <<", KSP="<<static_cast<int>(lastLinearReason)<<" ("<<KSPConvergedReasons[lastLinearReason]<<")"
      <<", line search="<<(lineType?lineType:"unknown")<<", line search max_it="<<lineMaxIterations
      <<", ierr="<<ierr<<", Newton="<<totalNewton
      <<", Krylov="<<totalKrylov<<", residual evaluations="<<totalEvaluations
      <<", friction branch predictors=disabled"
      <<", contact state updates="<<stateUpdates<<", activations="<<activations<<", releases="<<releases;
    if(context.evaluate(q,final))message<<", force residual ratio="<<final.diagnostic.maxForceResidualRatio
      <<", torque residual ratio="<<final.diagnostic.maxTorqueResidualRatio<<", gap violation="<<final.diagnostic.contactGapViolation
      <<" m, complementarity ratio="<<context.complementarity(final);
    if(!context.error.empty())message<<", last_trial_error="<<context.error;
    error=message.str();break;
  }
  if(finalReason>=0) {
    finalReason=SNES_DIVERGED_MAX_IT;
    if(error.empty())error="PETSc contact solve exhausted bounded contact-state restarts without physical convergence";
  }
  writeAttemptTrace(settings,allTrace,outerTime,outerDt,count,substep,dt,finalReason,error);
  return false;
}
#endif // SLURRY_USE_PETSC

inline bool dispatchImplicitStep(const std::vector<Body>& old,const std::vector<Vec3>& force,
    const std::vector<Vec3>& torque,double dt,double time,const ParticleStepSettings& settings,
    std::vector<GapCache>& cache,PersistentContactState& contacts,std::vector<Body>& output,
    ParticleStepDiagnostics& diagnostic,std::string& error,double outerTime,double outerDt,int count,int substep) {
  bool okay=false;
  const bool saveReplay=count>settings.maxSubsteps/2 && !settings.solverDiagnosticsPrefix.empty();
  const auto startingCache=saveReplay?cache:std::vector<GapCache>{};
  const auto startingContacts=saveReplay?contacts:PersistentContactState{};
  try {
  if(settings.solverBackend=="legacy")okay=implicitStep(old,force,torque,dt,time,settings,cache,contacts,output,diagnostic,error);
  else if(settings.solverBackend=="petsc") {
#ifdef SLURRY_USE_PETSC
    okay=implicitStepPetsc(old,force,torque,dt,time,settings,cache,contacts,output,diagnostic,error,outerTime,outerDt,count,substep);
#else
    throw std::runtime_error("particle_solver=petsc requires a build with SLURRY_USE_PETSC; no legacy fallback is permitted");
#endif
  }else throw std::invalid_argument("Unknown particle_solver: "+settings.solverBackend);
  }catch(const std::exception& ex){error=ex.what();
    writeAttemptTrace(settings,{},outerTime,outerDt,count,substep,dt,0,error);okay=false;}
  if(!okay && saveReplay) {
    writeParticleReplay({settings,dt,time,outerTime,outerDt,count,substep,old,force,torque,startingCache,startingContacts},settings.solverDiagnosticsPrefix+"_failure.dat");
  }
  return okay;
}
#endif
