#!/usr/bin/env python3
"""Prepare graphite-water RE-squared plus rough-contact shear; Python standard library only.

Example: python3 run_graphite.py --output runs/quick100 --ranks 16
Use --dry-run to inspect the Mach-based time mapping without running OpenLB.
"""
import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda: source.read(1048576), b''):
            value.update(block)
    return value.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def positive(value, name, zero=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or (value < 0 if zero else value <= 0):
        raise ValueError(name+' must be finite and '+('nonnegative' if zero else 'positive'))
    return value


def integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(name+' must be an integer >= '+str(minimum))
    return value


def resolve(config, shear_rate=None, max_steps=0, target_mach=None, time_step=None, end_strain=None):
    cfg = copy.deepcopy(config)
    if cfg.get('schema_version') != 2:
        raise ValueError('Only schema_version=2 is supported by this RE-squared runner')
    cfg['flow'].setdefault('end_strain', 10.0)
    defaults = {'particle_tolerance':1e-4,
                'particle_force_absolute_tolerance_N':1e-15,
                'particle_torque_absolute_tolerance_N_m':1.65e-21,
                'contact_gap_tolerance_m':1e-12}
    for key,value in defaults.items():
        cfg['numerics'].setdefault(key,value)
    cfg['numerics'].setdefault('particle_solver', 'petsc')
    cfg['numerics'].setdefault('particle_max_krylov_iterations', 120)
    cfg['numerics'].setdefault('solver_diagnostics', True)
    contact = cfg.setdefault('rough_contact',{})
    for key,value in {'checkpoint_every_steps':0,'checkpoint_every_seconds':21000.0,
                      'checkpoint_keep':2}.items():
        cfg['output'].setdefault(key,value)
    for key,value in {'enabled':True,'roughness_gap_m':2e-9,'sliding_friction':.5,
                      'tangential_stiffness_N_m':9.0,'rolling_length_m':100e-9,
                      'rolling_yield_angle_rad':.01}.items():
        contact.setdefault(key,value)
    # Missing local_gap_fraction means the original RE2 interaction, including
    # when running older JSON cases with the new executable.
    for key,value in {'local_gap_m':.3e-9,'local_gap_fraction':0.0,
                      'local_switch_excess_gap_m':2e-9,
                      'local_cutoff_excess_gap_m':10e-9}.items():
        cfg['interaction'].setdefault(key,value)
    if shear_rate is not None:
        cfg['flow']['shear_rate_s_inv'] = shear_rate
    if end_strain is not None:
        cfg['flow']['end_strain'] = end_strain
    if target_mach is not None:
        cfg['numerics']['target_mach'] = target_mach
    if time_step is not None:
        cfg['numerics']['time_step_s'] = time_step
    f, g, p, flow, n, interaction, out = (cfg[k] for k in ('fluid','geometry','particles','flow','numerics','interaction','output'))
    for group, names in ((f, ('density_kg_m3','dynamic_viscosity_Pa_s','temperature_K')),
                         (g, ('box_length_m','dx_m')),
                         (p, ('diameter_m','thickness_m','density_kg_m3')),
                         (flow, ('shear_rate_s_inv','end_strain')),
                         (n, ('nu_lattice','target_mach','epsilon_cells')),
                         (interaction, ('sigma_lj_m','switch_gap_m','cutoff_gap_m'))):
        for name in names:
            positive(group[name], name)
    positive(interaction['hamaker_J'], 'hamaker_J', zero=True)
    positive(interaction['local_gap_m'], 'local_gap_m')
    positive(interaction['local_gap_fraction'], 'local_gap_fraction', zero=True)
    if interaction['local_gap_fraction'] > 1:
        raise ValueError('local_gap_fraction must be at most one')
    positive(interaction['local_switch_excess_gap_m'], 'local_switch_excess_gap_m', zero=True)
    positive(interaction['local_cutoff_excess_gap_m'], 'local_cutoff_excess_gap_m')
    if not interaction['local_switch_excess_gap_m'] < interaction['local_cutoff_excess_gap_m']:
        raise ValueError('Require local_switch_excess_gap_m < local_cutoff_excess_gap_m')
    positive(n['time_step_s'], 'time_step_s', zero=True)
    positive(p['minimum_gap_m'], 'minimum_gap_m', zero=True)
    integer(p['count'], 'particles.count', 1)
    integer(p['seed'], 'particles.seed')
    integer(n['particle_max_substeps'], 'particle_max_substeps', 1)
    integer(n['particle_max_iterations'], 'particle_max_iterations', 1)
    integer(n['particle_max_krylov_iterations'], 'particle_max_krylov_iterations', 1)
    if n['particle_solver'] not in ('petsc', 'legacy'):
        raise ValueError('particle_solver must be petsc or legacy')
    if not isinstance(n['solver_diagnostics'], bool):
        raise ValueError('solver_diagnostics must be boolean')
    positive(n['particle_tolerance'], 'particle_tolerance')
    if n['particle_tolerance'] >= 1:
        raise ValueError('particle_tolerance must be below one')
    for name in ('particle_force_absolute_tolerance_N','particle_torque_absolute_tolerance_N_m','contact_gap_tolerance_m'):
        positive(n[name],name)
    positive(n['lubrication_cutoff_cells'], 'lubrication_cutoff_cells', zero=True)
    if not isinstance(contact['enabled'],bool):
        raise ValueError('rough_contact.enabled must be a JSON boolean')
    for name in ('roughness_gap_m','tangential_stiffness_N_m','rolling_yield_angle_rad'):
        positive(contact[name],'rough_contact.'+name)
    for name in ('sliding_friction','rolling_length_m'):
        positive(contact[name],'rough_contact.'+name,zero=True)
    if not n['contact_gap_tolerance_m'] < contact['roughness_gap_m'] < interaction['cutoff_gap_m']:
        raise ValueError('Require contact_gap_tolerance_m < roughness_gap_m < cutoff_gap_m')
    if contact['enabled'] and p['minimum_gap_m'] < contact['roughness_gap_m']:
        raise ValueError('particles.minimum_gap_m must be at least rough_contact.roughness_gap_m')
    if interaction['local_gap_fraction'] > 0:
        if not contact['enabled']:
            raise ValueError('Local-gap adhesion requires rough_contact.enabled=true')
        if interaction['local_gap_m'] > contact['roughness_gap_m']:
            raise ValueError('Local-gap adhesion requires local_gap_m <= roughness_gap_m')
        if contact['roughness_gap_m'] + interaction['local_cutoff_excess_gap_m'] > interaction['switch_gap_m']:
            raise ValueError('Local-gap adhesion requires roughness_gap_m + local_cutoff_excess_gap_m <= switch_gap_m')
    integer(max_steps, 'max_steps')
    integer(out['sample_every_steps'], 'sample_every_steps', 1)
    integer(out['vtk_every_steps'], 'vtk_every_steps')
    integer(out['checkpoint_every_steps'], 'checkpoint_every_steps')
    positive(out['checkpoint_every_seconds'], 'checkpoint_every_seconds', zero=True)
    integer(out['checkpoint_keep'], 'checkpoint_keep', 1)
    length, dx, diameter, thickness = g['box_length_m'],g['dx_m'],p['diameter_m'],p['thickness_m']
    rate, eta, nu_lb = flow['shear_rate_s_inv'],f['dynamic_viscosity_Pa_s'],n['nu_lattice']
    cells = round(length/dx)
    if cells < 8 or not math.isclose(length/dx, cells, rel_tol=0, abs_tol=1e-7):
        raise ValueError('box_length_m / dx_m must be an integer >= 8')
    if not 0 < thickness <= diameter < length:
        raise ValueError('Require 0 < thickness <= diameter < box length')
    if not interaction['sigma_lj_m'] < interaction['switch_gap_m'] < interaction['cutoff_gap_m']:
        raise ValueError('Require sigma_lj_m < switch_gap_m < cutoff_gap_m')
    if 2*(diameter+interaction['cutoff_gap_m']) >= length:
        raise ValueError('Pair cutoff requires 2*(diameter + cutoff_gap) < box length')
    volume = math.pi*diameter**2*thickness/6
    phi = p['count']*volume/length**3
    if phi >= 1:
        raise ValueError('Particle volume fraction must be below one')
    # Centered affine flow ux = rate * (y-Ly/2), with LB sound speed 1/sqrt(3).
    affine_speed = 0.5*rate*length
    dt = n['time_step_s'] or n['target_mach']*dx/(math.sqrt(3)*affine_speed)
    rho_num = eta*dt/(nu_lb*dx**2)
    density_scale = rho_num/f['density_kg_m3']
    re_phys = f['density_kg_m3']*rate*diameter**2/eta
    re_num = rho_num*rate*diameter**2/eta
    steps = math.ceil(flow['end_strain']/(rate*dt))
    stop = min(steps,max_steps) if max_steps else steps
    metadata = {
        'interaction':copy.deepcopy(interaction), 'rough_contact':copy.deepcopy(contact),
        'shear_rate_s_inv':rate, 'dt_s':dt,
        'time_mapping':'manual_time_step' if n['time_step_s'] else 'centered_affine_target_mach',
        'target_mach':n['target_mach'], 'affine_mach':math.sqrt(3)*affine_speed*dt/dx,
        'particle_solver':n['particle_solver'],
        'particle_tolerance':n['particle_tolerance'],
        'particle_force_absolute_tolerance_N':n['particle_force_absolute_tolerance_N'],
        'particle_torque_absolute_tolerance_N_m':n['particle_torque_absolute_tolerance_N_m'],
        'contact_gap_tolerance_m':n['contact_gap_tolerance_m'],
        'particle_max_substeps':n['particle_max_substeps'],
        'particle_max_iterations':n['particle_max_iterations'],
        'particle_max_krylov_iterations':n['particle_max_krylov_iterations'],
        'solver_diagnostics':n['solver_diagnostics'],
        'shear_increment_per_step':rate*dt,
        'nu_lattice':nu_lb, 'tau_lattice':0.5+3*nu_lb,
        'physical_fluid_re_D':re_phys, 'numerical_fluid_re_D':re_num,
        'physical_particle_St_D':re_phys*p['density_kg_m3']/f['density_kg_m3'],
        'numerical_particle_St_D':re_num*p['density_kg_m3']/f['density_kg_m3'],
        'fluid_inertial_density_kg_m3':rho_num,
        'particle_inertial_density_kg_m3':p['density_kg_m3']*density_scale,
        'inertial_density_scale':density_scale,
        'particle_volume_m3':volume, 'actual_volume_fraction':phi,
        'actual_mass_fraction':phi*p['density_kg_m3']/(phi*p['density_kg_m3']+(1-phi)*f['density_kg_m3']),
        'grid_cells_per_direction':cells, 'nominal_bulk_cells':cells**3,
        'd3q19_single_population_bytes':cells**3*19*8,
        'thickness_cells':thickness/dx, 'diameter_cells':diameter/dx,
        'steps_to_requested_end_strain':steps, 'this_run_max_step':stop,
        'this_run_target_strain':stop*rate*dt,
        'notes':[
            'Mach target refers to the centered imposed affine profile; measured maximum Mach includes particle-induced flow.',
            'Dynamic viscosity is preserved; fluid and particle inertial densities are scaled together.',
            'At the default target Mach, numerical inertia is deliberately increased for a quick magnitude check; no 10-percent error bound is assumed.',
            'Re_D uses diameter and St_D is Re_D times the particle/fluid density ratio; these are diagnostics, not run gates.',
            'The LB time step stays fixed throughout a run. End strain specifies its length, not rheological convergence.'
        ]
    }
    return cfg,metadata


def solver_values(cfg, output, particles, max_steps):
    f,g,p,flow,n,interaction,out = (cfg[k] for k in ('fluid','geometry','particles','flow','numerics','interaction','output'))
    contact = cfg['rough_contact']
    return {
        'shear_rate':flow['shear_rate_s_inv'],
        'box_x':g['box_length_m'], 'box_y':g['box_length_m'], 'box_z':g['box_length_m'], 'dx':g['dx_m'],
        'diameter':p['diameter_m'], 'thickness':p['thickness_m'],
        'rho_particle':p['density_kg_m3'], 'rho_fluid':f['density_kg_m3'],
        'dynamic_viscosity':f['dynamic_viscosity_Pa_s'],
        'nu_lattice':n['nu_lattice'], 'target_mach':n['target_mach'], 'time_step_s':n['time_step_s'],
        'epsilon_cells':n['epsilon_cells'],
        'particle_max_substeps':n['particle_max_substeps'],
        'particle_max_iterations':n['particle_max_iterations'],
        'particle_max_krylov_iterations':n['particle_max_krylov_iterations'],
        'particle_solver':n['particle_solver'],
        'solver_diagnostics':int(n['solver_diagnostics']),
        'particle_tolerance':n['particle_tolerance'],
        'particle_force_absolute_tolerance':n['particle_force_absolute_tolerance_N'],
        'particle_torque_absolute_tolerance':n['particle_torque_absolute_tolerance_N_m'],
        'contact_gap_tolerance':n['contact_gap_tolerance_m'],
        'lubrication_cutoff_cells':n['lubrication_cutoff_cells'],
        'rough_contact_enabled':int(contact['enabled']),
        'roughness_gap':contact['roughness_gap_m'],
        'sliding_friction':contact['sliding_friction'],
        'tangential_stiffness':contact['tangential_stiffness_N_m'],
        'rolling_length':contact['rolling_length_m'],
        'rolling_yield_angle':contact['rolling_yield_angle_rad'],
        'hamaker':interaction['hamaker_J'], 'sigma_lj':interaction['sigma_lj_m'],
        'switch_gap':interaction['switch_gap_m'], 'cutoff_gap':interaction['cutoff_gap_m'],
        'local_gap':interaction['local_gap_m'], 'local_gap_fraction':interaction['local_gap_fraction'],
        'local_switch_excess_gap':interaction['local_switch_excess_gap_m'],
        'local_cutoff_excess_gap':interaction['local_cutoff_excess_gap_m'],
        'end_strain':flow['end_strain'], 'max_steps':max_steps,
        'sample_every':out['sample_every_steps'], 'vtk_every':out['vtk_every_steps'],
        'checkpoint_every':out['checkpoint_every_steps'], 'checkpoint_seconds':out['checkpoint_every_seconds'],
        'checkpoint_keep':out['checkpoint_keep'],
        'output_dir':str(output), 'particles_csv':str(particles)
    }


def require_local_adhesion_build(cfg, build_info):
    if cfg['interaction']['local_gap_fraction'] > 0 and build_info.get('local_gap_adhesion') is not True:
        raise ValueError('This configuration requires local-gap adhesion. Rebuild with build_slurry_cpu.sbatch before running.')


def read_checkpoint(directory):
    directory=Path(directory).expanduser().resolve()
    data=json.loads((directory/'checkpoint.json').read_text(encoding='utf-8'))
    if data.get('engine')!='pure_gr' or data.get('format_version')!=1 or data.get('complete') is not True:
        raise ValueError('Not a completed pure_gr checkpoint: '+str(directory))
    integer(data['step'],'checkpoint step')
    integer(data['ranks'],'checkpoint ranks',1)
    positive(data['dt_s'],'checkpoint dt_s')
    positive(data['shear_rate_s_inv'],'checkpoint shear rate')
    expected={'state.bin','initial_particles.csv'}|{'lattice_rank_%d.bin'%rank for rank in range(data['ranks'])}
    if set(data['files'])!=expected:
        raise ValueError('Checkpoint file inventory is incomplete: '+str(directory))
    for name,size in data['files'].items():
        integer(size,'checkpoint file size',1)
        path=directory/name
        if not path.is_file() or path.stat().st_size!=size:
            raise ValueError('Missing or truncated checkpoint file: '+str(path))
    return directory,data


def latest_checkpoint(run):
    run=Path(run).expanduser().resolve()
    if (run/'checkpoint.json').is_file():
        return read_checkpoint(run)
    candidates=[]
    for directory in (run/'checkpoints').glob('checkpoint_*'):
        if directory.is_dir() and directory.name[11:].isdigit() and (directory/'checkpoint.json').is_file():
            candidates.append((int(directory.name[11:]),directory))
    if not candidates:
        raise ValueError('No completed pure_gr checkpoint in '+str(run)+
                         '. Old history.csv/particle failure files do not contain the fluid state.')
    return read_checkpoint(max(candidates)[1])


def restart_targets(folder):
    folder=Path(folder).expanduser().resolve()
    if (folder/'checkpoint.json').is_file() or (folder/'checkpoints').is_dir():
        return [latest_checkpoint(folder)]
    scope=folder/'pure_gr' if (folder/'pure_gr').is_dir() else folder
    runs=sorted({path.parent for path in scope.rglob('checkpoints') if path.is_dir()})
    if not runs:
        raise ValueError('No pure_gr checkpoint under '+str(folder)+
                         '. Checkpoints are available only for runs made with the restart update.')
    return [latest_checkpoint(run) for run in runs]


def validate_restart(cfg,meta,checkpoint,ranks,max_steps=0,allow_complete=False):
    values=solver_values(cfg,Path('run'),Path('particles.csv'),max_steps)
    values.update(dt_s=meta['dt_s'],particle_count=cfg['particles']['count'],ranks=ranks)
    changed=[]
    for key,saved in checkpoint['immutable_config'].items():
        actual=values.get(key)
        if actual is None or not math.isclose(float(actual),float(saved),rel_tol=1e-13,abs_tol=0.0):
            changed.append(key)
    if changed:
        raise ValueError('Restart physical settings/MPI layout differ: '+', '.join(changed)+
                         '. Keep the saved physical settings; end_strain, solver limits/tolerances and output intervals may change.')
    endpoint=math.ceil(cfg['flow']['end_strain']/(cfg['flow']['shear_rate_s_inv']*meta['dt_s']))
    if max_steps:endpoint=min(endpoint,max_steps)
    if endpoint<=checkpoint['step'] and not allow_complete:
        raise ValueError('Restart endpoint must exceed saved step %d; increase flow.end_strain or --max-steps.'%checkpoint['step'])
    return endpoint>checkpoint['step']


def copy_history_prefix(checkpoint,data,output):
    # Keep the source run intact. Copy exactly the flushed bytes represented by
    # the checkpoint, excluding any later samples from an interrupted attempt.
    source=checkpoint.parent.parent if checkpoint.parent.name=='checkpoints' else None
    if source is None:return False
    sizes=data.get('output_bytes',{})
    if not all((source/name).is_file() for name in ('history.csv','particles.csv')):
        return False
    for name in ('history.csv','particles.csv'):
        size=sizes.get(name,0)
        integer(size,'checkpoint output size',1)
        if (source/name).stat().st_size<size:
            raise ValueError('Source output is shorter than its checkpoint: '+str(source/name))
    for name in ('history.csv','particles.csv'):
        left=sizes[name]
        with (source/name).open('rb') as src,(output/name).open('wb') as dst:
            while left:
                block=src.read(min(left,1048576))
                if not block:raise ValueError('Source output was truncated while copying '+name)
                dst.write(block);left-=len(block)
    return True


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=Path(__file__).with_name('config.json'))
    parser.add_argument('--shear-rate',type=float)
    parser.add_argument('--target-mach',type=float)
    parser.add_argument('--time-step',type=float,help='Explicit fixed LB step in seconds; 0 uses target Mach')
    parser.add_argument('--end-strain',type=float,help='Accumulated strain to run; default config value is 10')
    parser.add_argument('--max-steps',type=int,default=0,help='Explicit step limit; 0 uses end strain')
    parser.add_argument('--dry-run',action='store_true')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--executable',type=Path,default=Path(__file__).with_name('build')/'current'/'graphiteCouette3d')
    parser.add_argument('--generator',type=Path,default=Path(__file__).with_name('generate_particles.py'))
    parser.add_argument('--ranks',type=int,default=1)
    parser.add_argument('--restart',type=Path,help='Completed pure_gr checkpoint or previous case directory')
    args=parser.parse_args()
    try:
        restart,checkpoint=latest_checkpoint(args.restart) if args.restart else (None,None)
        rate=args.shear_rate
        if checkpoint:
            if rate is not None and not math.isclose(rate,checkpoint['shear_rate_s_inv'],rel_tol=1e-13):
                raise ValueError('Restart retains the checkpoint shear rate')
            rate=checkpoint['shear_rate_s_inv']
        cfg,meta=resolve(json.loads(args.config.read_text(encoding='utf-8')),rate,args.max_steps,
                         args.target_mach,args.time_step,args.end_strain)
        integer(args.ranks,'ranks',1)
        if checkpoint:validate_restart(cfg,meta,checkpoint,args.ranks,args.max_steps)
        if args.dry_run:
            print(json.dumps({'config':cfg,'derived':meta},indent=2,allow_nan=False))
            return 0
        if args.output is None:
            raise ValueError('--output is required')
        output=args.output.expanduser().resolve()
        executable=args.executable.expanduser().resolve()
        generator=args.generator.expanduser().resolve()
        if not executable.is_file() or not os.access(executable,os.X_OK):
            raise ValueError('Executable not found or not executable: '+str(executable))
        # Python 3.6 compatibility: capture_output and text require Python 3.7.
        build_query=subprocess.run([str(executable),'--build-info'],
                                   stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                                   universal_newlines=True)
        try:
            build_info=json.loads(build_query.stdout)
        except (ValueError,TypeError):
            raise ValueError('Executable does not provide valid --build-info; rebuild this RE-squared application')
        if build_query.returncode or not isinstance(build_info,dict) or not isinstance(build_info.get('mpi_enabled'),bool):
            raise ValueError('Executable did not report its MPI build mode; rebuild this application')
        if build_info.get('rough_contact') is not True:
            raise ValueError('Executable predates rough contact and the corrected particle solver; run build_slurry_cpu.sbatch again')
        require_local_adhesion_build(cfg, build_info)
        if build_info.get('pure_gr_checkpoint_version')!=1:
            raise ValueError('This driver requires pure_gr checkpoint support. Rebuild with build_slurry_cpu.sbatch.')
        if cfg['numerics']['particle_solver'] == 'petsc' and build_info.get('particle_solver') != 'petsc':
            raise ValueError('This configuration requires the PETSc contact solver. Rebuild with build_slurry_cpu.sbatch before running.')
        if args.ranks>1 and not build_info['mpi_enabled']:
            raise ValueError('--ranks > 1 requires an MPI-enabled executable; use build_slurry_cpu.sbatch')
        if not restart and not generator.is_file():
            raise ValueError('Particle generator not found: '+str(generator))
        if any((output/name).exists() for name in ('manifest.json','history.csv','resolved_run.cfg','initial_particles.csv')):
            raise ValueError('Output already contains a run; specify a new --output directory')
        mpirun=shutil.which('mpirun') if args.ranks>1 else None
        if args.ranks>1 and not mpirun:
            raise ValueError('mpirun not found')
    except (ValueError,OSError,KeyError,TypeError) as error:
        parser.error(str(error))
    output.mkdir(parents=True,exist_ok=True)
    effective=output/'effective_config.json'
    write_json(effective,cfg)
    particles=output/'initial_particles.csv'
    if restart:
        shutil.copy2(str(restart/'initial_particles.csv'),str(particles))
        copy_history_prefix(restart,checkpoint,output)
    else:
        generation=subprocess.run([sys.executable,str(generator),'--config',str(effective),'--output',str(particles)])
        if generation.returncode:return generation.returncode
    config_path=output/'resolved_run.cfg'
    values=solver_values(cfg,output,particles,args.max_steps)
    values['restart_dir']=str(restart) if restart else ''
    with config_path.open('w',encoding='utf-8') as stream:
        for name,value in values.items():
            if '\n' in str(value) or '\r' in str(value):
                raise ValueError('Newline in solver parameter '+name)
            stream.write(str(name)+'='+str(value)+'\n')
    argv=([mpirun,'-np',str(args.ranks)] if mpirun else [])+[str(executable),'--config',str(config_path)]
    write_json(output/'manifest.json',{
        'created_utc':datetime.now(timezone.utc).isoformat(), 'slurm_job_id':os.environ.get('SLURM_JOB_ID'),
        'config':cfg,'derived':meta,'ranks':args.ranks,'argv':argv,'executable_build':build_info,
        'restart_directory':str(restart) if restart else None,
        'restart_step':checkpoint['step'] if checkpoint else None,
        'petsc_options':os.environ.get('PETSC_OPTIONS',''),
        'sha256':{'executable':digest(executable),'driver':digest(__file__),
                  'generator':digest(generator) if generator.is_file() else None,'particles':digest(particles)}
    })
    print(json.dumps(meta,indent=2),flush=True)
    print('Run directory: '+str(output),flush=True)
    start=time.monotonic()
    requested=[]
    def request_stop(number,frame):
        requested.append(number)
        (output/'STOP_REQUEST').touch()
        print('Stop requested; saving a collective pure_gr checkpoint at the next completed LB step.',flush=True)
    for number in (signal.SIGUSR1,signal.SIGTERM,signal.SIGINT):signal.signal(number,request_stop)
    with (output/'solver.log').open('w',encoding='utf-8') as logfile:
        process=subprocess.Popen(argv,cwd=str(output),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                                 universal_newlines=True,bufsize=1,start_new_session=True)
        for line in process.stdout:
            print(line,end='',flush=True)
            logfile.write(line)
            logfile.flush()
        rc=process.wait()
    status_path=output/'status.json'
    solver_status=json.loads(status_path.read_text(encoding='utf-8')) if rc==0 and status_path.is_file() else {}
    if rc==0 and solver_status.get('status') not in ('COMPLETED','MAX_STEPS','CHECKPOINTED'):
        print('Solver did not write a valid completion/checkpoint status.',file=sys.stderr);rc=1
    status={'exit_code':rc,'wall_seconds':time.monotonic()-start,
            'status':solver_status.get('status','FAILED') if rc==0 else 'FAILED',
            'stop_requested':bool(requested),'solver_status':solver_status}
    write_json(output/'driver_status.json',status)
    print(json.dumps(status,indent=2),flush=True)
    summarizer=output/'summarize_particle_solver.py'
    if rc!=0 and summarizer.is_file():
        report=subprocess.run([sys.executable,str(summarizer)],cwd=str(output),
                              stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                              universal_newlines=True)
        (output/'particle_solver_summary.txt').write_text(report.stdout,encoding='utf-8')
        print(report.stdout,end='',flush=True)
    return rc if rc>=0 else 128-rc


if __name__=='__main__':
    raise SystemExit(main())
