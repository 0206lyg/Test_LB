# PETSc 입자 solver 업데이트

기준 소스: `0206lyg/Test_LB`, commit
`ad0befc0f75b40e16bbd81fc3f29c53e76a8b551`.
이 패키지는 해당 저장소에 적용하는 소스 업데이트입니다. OpenLB 전체나 실행 파일은 포함하지 않습니다.

## 설치와 재빌드

다운로드한 ZIP을 MIT의 **`/home/lyjania/OpenLB`**에 놓고 아래 명령을 실행합니다.
이 폴더가 `build_slurry_cpu.sbatch`, `run_slurry_cpu.sbatch`, `olb-1.9r0`, `slurry`가 있는 작업 루트입니다.

```bash
cd /home/lyjania/OpenLB
unzip -o petsc_solver_update.zip
python3 petsc_solver_update/install.py
sbatch build_slurry_cpu.sbatch
```

**빌드는 기존 `build_slurry_cpu.sbatch` 하나입니다.** 이 작업 안에서 PETSc 확인/설치와
전체 공통 실행 파일 빌드를 순서대로 처리합니다. CMC와 graphite 모두 같은
`build/slurry/current/slurry`를 사용합니다. 추가 빌드 작업을 제출할 필요가 없습니다.

앞서 잘못 추가한 별도 빌드/설치 shell 스크립트가 이미 있으면 백업 후 제거합니다.
설치 스크립트는 바뀌는 기존 파일을 `solver_update_backups/<시각>/`에 백업합니다.
`slurry/cases/pure_gr.json`은 현재 파일을 읽어 `particle_solver=petsc`,
`solver_diagnostics=true`를 설정하고, 없으면 `particle_max_krylov_iterations=120`을 추가합니다.
기존 물성·입자 배치 설정·마찰·굴림 저항·Mach 수·수렴 허용오차·substep/Newton 상한은 보존합니다.
직접 수정했던 C++/실행 스크립트는 백업본과 비교해 필요한 변경을 병합하십시오.

`build_slurry_<jobid>.out`에서 `BUILD COMPLETE`를 확인한 뒤 실행합니다.

```bash
sbatch run_slurry_cpu.sbatch --cases pure_gr --shear-rates 100
```

첫 빌드는 호환되는 PETSc가 없으면 PETSc 3.25.5와 BLAS/LAPACK을 사용자 경로에 설치합니다.
이때 다운로드가 가능해야 합니다. 이후에는 설치한 의존성을 재사용합니다.
기존 GCC/OpenMPI 모듈과 Slurm 설정을 유지합니다. 다른 클러스터나 다운로드가 막힌 노드에서의
설치는 `README_PETSC_BUILD.md`를 참고하십시오. 이전 실행 파일은 재사용할 수 없습니다.

기존 `build_slurry_22944246.out`에서 PETSc 설치가 완료된 경우, 같은 MPI 모듈로
다시 빌드하면 설치본을 재사용합니다. MPI 사전검사는 `mpirun -np 1`로 실행하도록 수정했습니다.

## 변경 내용

- 입자 비선형 풀이를 PETSc `SNESNEWTONLS`와 backtracking line search, `KSPGMRES`로 변경했습니다.
  법선 접촉은 Fischer–Burmeister 상보성 잔차로 풀며, 행렬 없는 Jacobian 연산과 접촉 블록
  preconditioner를 사용합니다. PETSc SNES VI를 직접 호출하는 구현은 아닙니다.
- `roughContact.h`의 접촉·Coulomb 미끄러짐·접착 기반 굴림 저항, RE2와 lubrication 식은 유지했습니다.
  접촉 이력은 substep 시작 상태에서 평가하고, 수렴한 해만 확정합니다.
- PETSc의 종료 판정 뒤에도 힘·토크·접촉 조건을 검사합니다. 기존 허용오차를 완화하거나
  실패 시 자동으로 기존 solver로 돌아가지 않습니다.
- 설정과 실행 파일의 solver 종류를 확인하고 PETSc/MPI 의존성 및 빌드 정보를 기록합니다.
  실패한 실행에서는 요약 스크립트를 자동 호출합니다.

## 실패 기록과 재현

입자 단계가 모든 subdivision 시도 뒤 최종 실패하면 결과 폴더에 rank별 파일을 저장합니다.

| 파일 | 내용 |
|---|---|
| `particle_solver_rankN_trace.csv` | 반복별 힘·토크 오차 비율, gap/상보성 위반, line-search 계수, SNES/KSP 종료 이유, 접촉·미끄러짐·굴림 개수와 전환 |
| `particle_solver_rankN_outcomes.csv` | 실패한 외부 시간 단계와 subdivision 수 |
| `particle_solver_rankN_failure.dat` | 마지막 실패 substep의 입자 상태·외력·접촉 이력·설정 |
| `particle_solver_summary.txt` | 실행 드라이버가 생성한 진단 요약 |

정상 실행 및 subdivision으로 회복한 단계는 진단 파일을 쓰지 않습니다.
메모리에 보관하는 반복 기록은 제한되어 있어 최종 실패에 가까운 시도부터 남습니다.
기록은 입자 solver의 수치 실패용입니다. Slurm 시간 제한, 프로세스 강제 종료, 파일시스템 장애에는
실패 파일이 생성되지 않을 수 있습니다. 전체 유체장 checkpoint/restart 기능은 아닙니다.

결과 폴더에서 요약을 다시 볼 수 있습니다.

```bash
python3 summarize_particle_solver.py
```

저장소 루트에서 재현 도구를 한 번 빌드합니다. 빌드 작업과 같은 GCC/OpenMPI 모듈을 사용하십시오.

```bash
source build/petsc/env.sh
python3 slurry/tools/compile_petsc_test.py tests/petsc_contact/replay_particle_step.cpp
build/petsc-tests/replay_particle_step /결과폴더/particle_solver_rank0_failure.dat
```

이 도구는 저장된 마지막 입자 substep만 풉니다. 같은 실패를 재현한 다음 solver 변경 효과를
비교할 수 있습니다. `--max-newton 80` 등은 비교 실험용이며 기본 설정을 바꾸지는 않습니다.

## 검증 범위

실제 PETSc 3.25.5를 빌드해 solver 및 재현 도구를 검증했습니다.
OpenLB 전체 실행 파일의 직렬 빌드와 축소된 입자/격자의 짧은 실행으로 연결을 확인했습니다.
물리 파라미터를 사용한 oblate 입자의 접착·미끄러짐·굴림 저항 테스트도 포함합니다.
단위 및 통합 검증의 구체적인 결과는 `VALIDATION.md`에 기록합니다.

이 환경에서는 사용자 클러스터의 MPI 다중 rank 실행과 108개 입자의 전체 strain 계산을
완료하지 않았습니다. 기존 실패 시점까지 안정적으로 진행하는지는 재빌드 후 실제 조건에서
확인해야 합니다. 새 solver의 채택만으로 전체 실행의 수렴을 보장하지 않습니다.

추가 solver 검사는 저장소 루트에서 실행할 수 있습니다.

```bash
source build/petsc/env.sh
python3 slurry/tools/compile_petsc_test.py tests/petsc_contact/contact_solver_tests.cpp --run
python3 -m unittest discover -s slurry/tests -p 'test_particle_diagnostics.py'
```
