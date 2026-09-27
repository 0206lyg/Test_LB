#ifndef SLURRY_GR_RE2_PARTICLE_PAIR_JACOBIAN_H
#define SLURRY_GR_RE2_PARTICLE_PAIR_JACOBIAN_H

// Included inside graphite::particle_detail, after Residual and Block6.
// Assemble local interaction derivatives once per Newton iterate.  Each local
// residual is the same Residual implementation used by the global solve, with
// two particles and the corresponding immutable contact history.  Krylov and
// trust-region products then require no new force or closest-gap evaluations.
struct PairAssembledJacobian {
  struct PairBlock {
    std::array<std::size_t,13> indices{};
    std::array<double,169> values{};
    int size=12;
  };
  std::size_t dimension=0;
  std::vector<Block6> inertiaBlocks,bodyBlocks;
  std::vector<PairBlock> pairs;
  int pairEvaluations=0;

  static Six inertiaResidual(const Body& old,const Six& q,double dt,double lengthScale) {
    Body next=old;
    const Vec3 angle{q[3],q[4],q[5]};
    next.rotation=rotateLaboratory(old.rotation,angle);
    next.omega=scale(angle,1./dt);
    const Vec3 change=scale(sub(worldMomentum(next),worldMomentum(old)),1./dt);
    const double torqueScale=*std::max_element(old.inertiaBody.begin(),old.inertiaBody.end())/(dt*dt);
    Six result{};
    for(int k=0;k<3;++k) {
      result[k]=q[k]-old.velocity[k]*dt/lengthScale;
      result[k+3]=change[k]/torqueScale;
    }
    return result;
  }

  // Keep the selected slip/rolling return-map derivative at a yield corner,
  // as in the existing semismooth matrix-free Newton product.
  static void correctContactColumn(const Residual& r,const Evaluation& base,
                                   const Evaluation& moved,double epsilon,Vector& column) {
    if(base.pairs.empty()||!base.pairs.front().contact.active)return;
    const auto& p=base.pairs.front();
    const auto& next=moved.pairs.front();
    const auto& c=p.contact;const auto& nc=next.contact;
    auto tangent=[](const Vec3& trial,const Vec3& change,double stiffness,double cap,
                    double dk,double dc,bool yielded)->Vec3 {
      const double magnitude=norm(trial);
      if(cap<=0.)return magnitude>0.?scale(trial,-dc/magnitude):Vec3{};
      if(!yielded)return add(scale(change,-stiffness),scale(trial,-dk));
      const Vec3 unit=scale(trial,1./magnitude);
      return add(scale(sub(change,scale(unit,dot(unit,change))),-cap/magnitude),scale(unit,-dc));
    };
    const Vec3 ds=scale(sub(nc.trialSlip,c.trialSlip),1./epsilon);
    const Vec3 dr=scale(sub(nc.trialRoll,c.trialRoll),1./epsilon);
    const Vec3 dn=scale(sub(next.normal,p.normal),1./epsilon);
    const Vec3 dli=scale(sub(next.leverI,p.leverI),1./epsilon);
    const Vec3 dlj=scale(sub(next.leverJ,p.leverJ),1./epsilon);
    const auto& state=c.candidateState;const auto& ns=nc.candidateState;
    const Vec3 dft=tangent(c.trialSlip,ds,r.settings.rough.tangentialStiffness,
        r.settings.rough.friction*std::max(0.,state.normalLoad),0.,0.,c.sliding);
    const Vec3 dmr=tangent(c.trialRoll,dr,state.rollingStiffness,state.rollingCap,
        (ns.rollingStiffness-state.rollingStiffness)/epsilon,
        (ns.rollingCap-state.rollingCap)/epsilon,c.rolling);
    const Vec3 df=add(scale(dn,-state.normalLoad),dft);
    const Vec3 dti=add(add(cross(dli,c.forceI),cross(p.leverI,df)),dmr);
    const Vec3 dtj=scale(add(add(cross(dlj,c.forceI),cross(p.leverJ,df)),dmr),-1.);
    const Vec3 cf=sub(df,scale(sub(nc.forceI,c.forceI),1./epsilon));
    const Vec3 cti=sub(dti,scale(sub(nc.torqueI,c.torqueI),1./epsilon));
    const Vec3 ctj=sub(dtj,scale(sub(nc.torqueJ,c.torqueJ),1./epsilon));
    const double fi=r.old[0].mass*r.lengthScale/(r.dt*r.dt);
    const double fj=r.old[1].mass*r.lengthScale/(r.dt*r.dt);
    const double ti=*std::max_element(r.old[0].inertiaBody.begin(),r.old[0].inertiaBody.end())/(r.dt*r.dt);
    const double tj=*std::max_element(r.old[1].inertiaBody.begin(),r.old[1].inertiaBody.end())/(r.dt*r.dt);
    for(int k=0;k<3;++k) {
      column[k]-=cf[k]/fi;column[6+k]+=cf[k]/fj;
      column[3+k]-=cti[k]/ti;column[9+k]-=ctj[k]/tj;
    }
  }

  bool build(const Residual& r,const Vector& q,const Evaluation& global,
             const Vector& weights,std::string& error) {
    dimension=q.size();pairEvaluations=0;pairs.clear();
    const std::size_t n=r.old.size();
    inertiaBlocks.assign(n,Block6{});bodyBlocks.assign(n,Block6{});
    for(std::size_t i=0;i<n;++i) {
      Six localQ{};for(int k=0;k<6;++k)localQ[k]=q[6*i+k];
      const Six initial=inertiaResidual(r.old[i],localQ,r.dt,r.lengthScale);
      for(int col=0;col<6;++col) {
        if(col<3) {inertiaBlocks[i][6*col+col]=weights[6*i+col];continue;}
        const double epsilon=r.settings.finiteDifferenceStep*(1.+std::abs(localQ[col]));
        Six changed=localQ;changed[col]+=epsilon;
        const Six displaced=inertiaResidual(r.old[i],changed,r.dt,r.lengthScale);
        for(int row=3;row<6;++row)
          inertiaBlocks[i][6*row+col]=(displaced[row]-initial[row])/epsilon*weights[6*i+row];
      }
    }
    bodyBlocks=inertiaBlocks;
    const double range=std::max(r.settings.pair.cutoffGap,std::max(r.settings.nearField.matchingGap,
        r.settings.rough.enabled?r.settings.rough.gap:0.));
    struct PairTask {std::size_t i,j,index;};
    std::vector<PairTask> tasks;
    std::size_t index=0;
    for(std::size_t i=0;i<n;++i)for(std::size_t j=i+1;j<n;++j,++index) {
      if(r.activeSlot[index]>=0 || global.gaps[index]<range)tasks.push_back({i,j,index});
    }
    pairs.resize(tasks.size());
    std::vector<int> evaluations(tasks.size(),0);
    std::vector<unsigned char> success(tasks.size(),0);
    std::vector<std::string> errors(tasks.size());
    // Every worker owns its local Residual, cache, contact history and output
    // block.  No MPI/PETSc calls occur here.  Reduce blocks in the original
    // lexicographic pair order below, making thread counts bitwise equivalent.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(tasks.size()>=128)
#endif
    for(std::ptrdiff_t task=0;task<static_cast<std::ptrdiff_t>(tasks.size());++task) {
      auto assemblePair=[&]()->bool {
      const auto i=tasks[task].i,j=tasks[task].j,index=tasks[task].index;
      const int slot=r.activeSlot[index];
      int& pairEvaluations=evaluations[task];
      std::string& error=errors[task];
      const bool candidate=slot>=0;
      const std::vector<Body> old{r.old[i],r.old[j]};
      const std::vector<Vec3> zero(2);
      const std::vector<GapCache> cache{r.startingCache[index]};
      const PersistentContactState contacts{r.startingContacts[index]};
      const std::vector<int> slots{candidate?0:-1};
      const std::vector<unsigned char> engaged{static_cast<unsigned char>(global.contacts[index].active)};
      Residual local{old,zero,zero,r.settings,cache,contacts,slots,r.dt,r.time,r.lengthScale,r.reactionScale,candidate?1:0};
      local.complementarity=r.complementarity;local.engagementOverride=&engaged;
      Vector localQ(candidate?13:12,0.);
      for(int k=0;k<6;++k){localQ[k]=q[6*i+k];localQ[6+k]=q[6*j+k];}
      if(candidate)localQ[12]=q[6*n+slot];
      Evaluation base;
      if(!local(localQ,base,error))return false;
      ++pairEvaluations;
      PairBlock block;block.size=candidate?13:12;
      for(int k=0;k<6;++k){block.indices[k]=6*i+k;block.indices[6+k]=6*j+k;}
      if(candidate)block.indices[12]=6*n+slot;
      double gapFactor=0.,loadFactor=0.;
      if(candidate) {
        const double a=(base.gaps[0]-r.settings.rough.gap)/r.settings.contactGapTolerance;
        const double b=localQ[12]*r.lengthScale/r.settings.contactGapTolerance;
        if(a<=b)gapFactor=r.lengthScale/r.settings.contactGapTolerance;
        else loadFactor=r.lengthScale/r.settings.contactGapTolerance;
      }
      for(int col=0;col<12;++col) {
        // Pair forces depend on relative translation/velocity. Their second
        // translational block is the negative of the first, including LE images.
        if(col>=6&&col<9) {
          for(int row=0;row<block.size;++row)block.values[13*row+col]=-block.values[13*row+col-6];
          continue;
        }
        double epsilon=r.settings.finiteDifferenceStep*(1.+std::abs(localQ[col]));
        Vector trial=localQ;Evaluation moved;bool valid=false;
        const double adhesionEdge=r.settings.pair.roughnessGap+r.settings.pair.adhesionRange;
        for(int attempt=0;attempt<12;++attempt) {
          trial[col]=localQ[col]+epsilon;
          ++pairEvaluations;
          if(local(trial,moved,error)) {
            const bool changesAdhesionBranch=r.settings.pair.surfaceAdhesion
                && ((base.gaps[0]<adhesionEdge)!=(moved.gaps[0]<adhesionEdge));
            if(!changesAdhesionBranch){valid=true;break;}
          }
          epsilon*=attempt%2==0?-1.:-.25;
        }
        if(!valid){error="Cannot differentiate particle pair within its constitutive branch: "+error;return false;}
        Vector column(block.size,0.);
        for(int row=0;row<block.size;++row)column[row]=(moved.residual[row]-base.residual[row])/epsilon;
        if(candidate)correctContactColumn(local,base,moved,epsilon,column);
        // Remove the two-particle residual's inertia; global inertia is inserted
        // once above, independently of each particle's number of neighbours.
        const int side=col/6;
        Six beforeQ{},afterQ{};
        for(int k=0;k<6;++k){beforeQ[k]=localQ[6*side+k];afterQ[k]=trial[6*side+k];}
        const Six before=inertiaResidual(old[side],beforeQ,r.dt,r.lengthScale);
        const Six after=inertiaResidual(old[side],afterQ,r.dt,r.lengthScale);
        for(int k=0;k<6;++k)column[6*side+k]-=(after[k]-before[k])/epsilon;
        for(int row=0;row<12;++row)block.values[13*row+col]=column[row]*weights[block.indices[row]];
        if(candidate)block.values[13*12+col]=gapFactor*column[12];
      }
      if(candidate) {
        const auto& p=base.pairs.front();
        const double fi=old[0].mass*r.lengthScale/(r.dt*r.dt),fj=old[1].mass*r.lengthScale/(r.dt*r.dt);
        const double ti=*std::max_element(old[0].inertiaBody.begin(),old[0].inertiaBody.end())/(r.dt*r.dt);
        const double tj=*std::max_element(old[1].inertiaBody.begin(),old[1].inertiaBody.end())/(r.dt*r.dt);
        for(int k=0;k<3;++k) {
          block.values[13*k+12]=-p.forcePerNormalLoad[k]*r.reactionScale/fi*weights[6*i+k];
          block.values[13*(6+k)+12]=p.forcePerNormalLoad[k]*r.reactionScale/fj*weights[6*j+k];
          block.values[13*(3+k)+12]=-p.torqueIPerNormalLoad[k]*r.reactionScale/ti*weights[6*i+3+k];
          block.values[13*(9+k)+12]=-p.torqueJPerNormalLoad[k]*r.reactionScale/tj*weights[6*j+3+k];
        }
        block.values[13*12+12]=loadFactor;
      }
      pairs[task]=std::move(block);
      return true;
      };
      // Exceptions must not cross an OpenMP region; report the first failed
      // pair in deterministic order after all workers finish.
      try {success[task]=assemblePair();}
      catch(const std::exception& ex) {errors[task]=ex.what();}
      catch(...) {errors[task]="Unknown exception during pair Jacobian assembly";}
    }
    for(std::size_t task=0;task<tasks.size();++task)pairEvaluations+=evaluations[task];
    for(std::size_t task=0;task<tasks.size();++task) {
      if(!success[task]){error=errors[task];return false;}
      const auto i=tasks[task].i,j=tasks[task].j;
      const auto& block=pairs[task];
      for(int row=0;row<6;++row)for(int col=0;col<6;++col) {
        bodyBlocks[i][6*row+col]+=block.values[13*row+col];
        bodyBlocks[j][6*row+col]+=block.values[13*(row+6)+col+6];
      }
    }
    for(const auto& block:inertiaBlocks)for(double value:block)
      if(!std::isfinite(value)){error="Non-finite inertial Jacobian";return false;}
    for(const auto& block:pairs)for(double value:block.values)
      if(!std::isfinite(value)){error="Non-finite pair Jacobian";return false;}
    return true;
  }

  // Replace only the block approximation, retaining the exact operator and
  // the NCP derivative selected by NcpPreconditioner::build.  Unlike the old
  // resistance-only blocks these also include conservative adhesive stiffness
  // and the complete rotating-inertia derivative.  Failure leaves the existing
  // preconditioner unchanged, including every cached Schur complement.
  template<class Preconditioner>
  bool updatePreconditioner(Preconditioner& pc)const {
    const std::size_t n=bodyBlocks.size();
    if(pc.blocks.n!=n||pc.bodyWeight.size()<6*n
        ||pc.gapFactor.size()!=pc.blocks.contact.size()
        ||pc.loadFactor.size()!=pc.blocks.contact.size())return false;
    std::vector<Block6> inverse(n);
    for(std::size_t i=0;i<n;++i) {
      Block6 physical=bodyBlocks[i];
      for(int row=0;row<6;++row) {
        const double weight=pc.bodyWeight[6*i+row];
        if(!(weight>0.)||!std::isfinite(weight))return false;
        for(int col=0;col<6;++col)physical[6*row+col]/=weight;
      }
      if(!inverse6(physical,inverse[i]))return false;
      for(double value:inverse[i])if(!std::isfinite(value))return false;
      // Guard a numerically singular inverse as well as a literal zero pivot.
      // This is a preconditioner fallback, not a change to the Newton Jacobian.
      double identityError=0.;
      for(int row=0;row<6;++row)for(int col=0;col<6;++col) {
        double product=0.;for(int k=0;k<6;++k)product+=physical[6*row+k]*inverse[i][6*k+col];
        identityError=std::max(identityError,std::abs(product-(row==col?1.:0.)));
      }
      if(!std::isfinite(identityError)||identityError>1.e-5)return false;
    }
    std::vector<double> contactSchur(pc.blocks.contact.size()),ncpSchur(pc.blocks.contact.size());
    for(std::size_t p=0;p<pc.blocks.contact.size();++p) {
      const auto& c=pc.blocks.contact[p];
      if(c.i>=n||c.j>=n)return false;
      contactSchur[p]=dot6(c.ci,multiply6(inverse[c.i],c.bi))
                    +dot6(c.cj,multiply6(inverse[c.j],c.bj));
      ncpSchur[p]=pc.gapFactor[p]*contactSchur[p]-pc.loadFactor[p];
      if(!std::isfinite(contactSchur[p])||!std::isfinite(ncpSchur[p])
          ||std::abs(ncpSchur[p])<1.e-30)return false;
    }
    pc.blocks.inverse=std::move(inverse);
    for(std::size_t p=0;p<pc.blocks.contact.size();++p)pc.blocks.contact[p].schur=contactSchur[p];
    pc.schur=std::move(ncpSchur);
    return true;
  }

  void multiply(const Vector& x,Vector& y)const {
    y.assign(dimension,0.);
    for(std::size_t i=0;i<inertiaBlocks.size();++i)for(int row=0;row<6;++row)
      for(int col=0;col<6;++col)y[6*i+row]+=inertiaBlocks[i][6*row+col]*x[6*i+col];
    for(const auto& block:pairs)for(int row=0;row<block.size;++row)
      for(int col=0;col<block.size;++col)y[block.indices[row]]+=block.values[13*row+col]*x[block.indices[col]];
  }
  void multiplyTranspose(const Vector& x,Vector& y)const {
    y.assign(dimension,0.);
    for(std::size_t i=0;i<inertiaBlocks.size();++i)for(int row=0;row<6;++row)
      for(int col=0;col<6;++col)y[6*i+col]+=inertiaBlocks[i][6*row+col]*x[6*i+row];
    for(const auto& block:pairs)for(int row=0;row<block.size;++row)
      for(int col=0;col<block.size;++col)y[block.indices[col]]+=block.values[13*row+col]*x[block.indices[row]];
  }
};

#endif
