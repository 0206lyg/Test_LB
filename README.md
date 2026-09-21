# OpenLB 배터리 슬러리 — 통합 사용 설명서

OpenLB 1.9r0과 CMC/graphite 슬러리 확장 모듈의 통합 저장소입니다.
**빌드는 `build_slurry_cpu.sbatch` 하나, 실행은 `run_slurry_cpu.sbatch` 하나**를 사용합니다.
프로젝트의 설치·설정·재시작·모델·검증 설명은 이 README에서 관리합니다.
OpenLB 원본과 외부 라이브러리의 문서·LICENSE는 각 소스 트리에 보존합니다.

| 모델 | 내용 | 설정 파일 |
| --- | --- | --- |
| `pure_cmc` | CMC Cross 유체 | `slurry/cases/run.json`의 `cmc`, `cmc250k_cross_parameters.csv` |
| `pure_gr` | RE² graphite, 국소 부착, lubrication, rough contact, Lees–Edwards, checkpoint/restart | `slurry/cases/pure_gr.json` |
| `gr_baseline` | 기존 graphite Couette 및 기존 checkpoint/restart | `slurry/cases/gr_baseline.json` |

## 1. 설치, 빌드, 실행

MIT 작업 루트는 `/home/lyjania/OpenLB`입니다. 이 안에 `olb-1.9r0/`, `slurry/`,
두 sbatch 파일이 있어야 합니다. 소스 업데이트를 적용한 뒤 다음 명령 하나로 빌드합니다.

```bash
cd /home/lyjania/OpenLB
sbatch build_slurry_cpu.sbatch
```

기존의 분산된 프로젝트 설명서는 이 빌드에서 자동으로 정리합니다. 지정된 이전 문서만
`build/documentation_backup/obsolete_documents_<시각>_<hash>.zip` 하나로 백업하고, 원본 내용과
백업을 대조한 뒤 개별 파일을 제거합니다. 사용자 JSON과 OpenLB·외부 라이브러리 문서는
보존하며, 별도 설치·정리 명령은 필요 없습니다.

`build_slurry_<jobid>.out`에서 `BUILD COMPLETE` 또는 `BUILD REUSED`를 확인한 뒤 실행합니다.

```bash
sbatch run_slurry_cpu.sbatch --cases pure_gr --shear-rates 100
```

여러 모델은 쉼표로 지정합니다. 공통 실행 파일은 `build/slurry/current/slurry`이고,
결과는 `runs/` 아래 실행별 폴더에 저장됩니다.

```bash
sbatch run_slurry_cpu.sbatch --cases pure_cmc,pure_gr --shear-rates 100
```

소스 변경에는 재빌드가 필요합니다. JSON의 물성·출력 설정만 바꿀 때는 새 실행부터 적용되며,
재빌드는 필요하지 않습니다. 기존 실행과 생성된 결과에는 소급 적용되지 않습니다.

### PETSc와 클러스터 환경

기본 모듈은 `gcc/12.2.0`, `openmpi/4.1.4`, Slurm partition은 `mit_normal`입니다.
같은 컴파일러/MPI 환경으로 의존성 설치, 빌드, 실행을 진행합니다.
다른 클러스터에서는 실제 모듈명을 `SLURRY_GCC_MODULE`, `SLURRY_MPI_MODULE`로 지정하고,
partition은 `sbatch --partition=...`로 지정합니다.

PETSc 준비는 **동일한 빌드 작업 내부**에서 수행합니다. 호환 설치본이 없으면 공식 PETSc 3.25.5와
BLAS/LAPACK을 내려받아 `$HOME/.local/slurry-petsc/3.25.5/<toolchain-id>/install`에 설치하고
이후 재사용합니다. 첫 설치에는 다운로드가 가능해야 하며 sudo나 별도 MPI 설치는 필요 없습니다.
다운로드가 막힌 환경에서는 호환 PETSc를 준비하고 `PETSC_DIR`을 지정한 뒤 같은 빌드 명령을 씁니다.
in-place PETSc 빌드에는 `PETSC_ARCH`도 지정합니다. 저장된 환경은 `build/petsc/env.sh`입니다.

실수형 배정밀도 PETSc 3.18 이상이 필요합니다. MPI 실행 파일에는 MPIUNI PETSc를 사용할 수 없습니다.
빌드 전 컴파일·링크·PETSc 초기화·SNES 생성과 MPI 실행을 검사하며, MPI 검사는 선택된 MPI의
`mpirun -np 1` 또는 `mpiexec`로 수행합니다. 앞선 빌드가 이 사전검사에서만 실패한 경우에도
완료된 PETSc 설치는 재사용합니다. 설치·컴파일러·라이브러리 정보는 빌드 fingerprint와
`build_manifest.json`에 기록됩니다. 실행 전에도 의존성과 지원 기능을 검사합니다.

개발용 비교가 필요한 경우에만 같은 빌드 진입점에 `--particle-solver legacy`를 전달할 수 있습니다.
현재 기본은 PETSc입니다. 공식 설치 정보:
[PETSc installation](https://petsc.org/release/install/install/),
[PETSC_DIR/PETSC_ARCH](https://petsc.org/release/install/multibuild/).

## 2. pure_gr 설정과 부착력 조절

설정 단위는 SI이며 입력은 `slurry/cases/pure_gr.json`입니다. 아래는 저장소의 기준값입니다.
사용자가 바꾼 설정은 실제 실행에 사용한 JSON과 결과 폴더의 설정 기록을 기준으로 확인합니다.

| JSON 키 | 기준값 | 의미 |
| --- | ---: | --- |
| `interaction.hamaker_J` | `9.9e-20` J | Hamaker 상수 |
| `interaction.sigma_lj_m` | `4.197e-10` m | LJ 반발 길이, 0.4197 nm |
| `interaction.local_gap_m` | `3e-10` m | 기하학적 접촉에서의 국소 간격 D₀ |
| `interaction.local_gap_fraction` | `1.0` | 국소 보정의 유효 기여율 α |
| `interaction.local_switch_excess_gap_m` | `2e-9` m | h−h₀에 대한 감쇠 시작 |
| `interaction.local_cutoff_excess_gap_m` | `1e-8` m | h−h₀에 대한 보정 종료 |
| `rough_contact.roughness_gap_m` | `2e-9` m | 기하학적 접촉 간격 h₀ |
| `rough_contact.sliding_friction` | `1.0` | 마찰계수 μ |
| `rough_contact.tangential_stiffness_N_m` | `80.0` N/m | 접선 탄성 강성 kₜ |
| `rough_contact.rolling_length_m` | `1e-7` m | rolling 길이 |
| `rough_contact.rolling_yield_angle_rad` | `0.01` rad | rolling 항복각 |

**새 부착 보정을 낮추려면 `interaction.local_gap_fraction`을 낮춥니다.**
예를 들어 `1.0 → 0.8`은 동일한 입자 간격과 배향에서 추가 보정 항을 20% 줄입니다.
D₀, σ, h₀, μ를 함께 바꿀 필요는 없습니다.

| α | 동일한 FF 접촉의 총 순인력 |
| ---: | ---: |
| 1.0 | 936.215 nN |
| 0.8 | 754.586 nN |
| 0.7 | 663.772 nN |
| 0.5 | 482.144 nN |

이 접촉에서는 총 순인력이 `28.073 + α × (936.215 − 28.073)` nN입니다.
벌크 응력은 접촉망과 배향 변화에도 의존하므로 같은 비율로 줄어든다고 가정하지 않습니다.
국소 간격 항이 없는 이전 설정은 α=0으로 읽어 기존 RE²를 사용합니다.

## 3. pure_gr 재시작

이 기능으로 생성한 **완료 checkpoint가 있는 실행부터** 재시작할 수 있습니다.
예전 `history.csv`, `particles.csv`, 입자 실패 파일만으로는 전체 유체 상태를 복원할 수 없습니다.

```bash
sbatch run_slurry_cpu.sbatch --restart 폴더명
```

`runs/폴더명`에서 최신 완료 checkpoint를 선택합니다. `runs/`를 포함한 상대 경로,
절대 경로, 개별 전단율 폴더, 특정 `checkpoints/checkpoint_...`도 지정할 수 있습니다.
`--restart`만 쓰면 모델은 자동으로 `pure_gr`이고, 각 checkpoint의 기존 전단율을 유지합니다.
여러 전단율이 있으면 현재 종료 지점에 도달하지 않은 각 케이스를 이어갑니다.
아직 시작하지 않아 checkpoint가 없는 전단율은 별도 신규 실행 대상입니다.

현재 `/home/lyjania/OpenLB/slurry/cases/pure_gr.json`을 읽습니다.
`--config` 또는 공통 `run.json`에 다른 파일을 지정했다면 그 파일을 사용합니다.

| 설정 | 재시작 시 처리 |
| --- | --- |
| `flow.end_strain` | 현재 값. 처음부터의 **총 누적 strain** 목표 |
| `--max-steps` | 처음부터의 **총 LB step 수** 상한 |
| 출력·checkpoint 간격, 보관 개수 | 현재 값 |
| solver 반복 상한·수렴 허용오차·진단 옵션 | 현재 값 |
| 물리값, 입자 수·크기, 격자, 시간 간격, MPI rank 수 | 저장 당시와 같아야 함 |
| 초기 배치 seed·minimum gap | 배치를 다시 만들지 않고 저장 상태 사용 |

strain 4.2에서 중단하고 `end_strain=10`이면 4.2부터 10까지 진행합니다.
완료된 계산을 더 연장하려면 JSON의 종료 strain을 높이거나 다음처럼 지정합니다.

```bash
sbatch run_slurry_cpu.sbatch --restart 폴더명 --end-strain 20
```

α, D₀, σ, 마찰, 접촉 강성 등을 바꾸면 호환성 검사에서 중단합니다.
부착력 변경 비교는 새 계산으로 시작합니다. 이전 실행에 `--ranks`를 지정했다면
재시작에서도 동일한 값을 사용합니다.

결과는 새 `runs/slurry_...` 폴더에 저장합니다. 원래 결과를 보존하고 `history.csv`와
`particles.csv`는 저장 지점까지 복사한 뒤 중복 행 없이 이어 씁니다.
저장 지점 이후 원래 실행이 기록한 행은 복사하지 않습니다. 이전 VTK는 원래 폴더에 남습니다.
CSV 없이 checkpoint만 옮긴 경우 새 CSV는 저장 시각부터 시작합니다.

### 저장 주기와 파일

기본 주기는 실제 실행 시간(wall time) **350분**이며 step 기준 주기 저장은 꺼져 있습니다.
현재 `pure_gr.json`에 아래 키가 없으면 자동으로 이 기본값을 사용합니다.
다른 주기를 원하거나 기존에 명시한 주기가 있다면 `output` 객체에서 수정합니다.

```json
"checkpoint_every_steps": 0,
"checkpoint_every_seconds": 21000.0,
"checkpoint_keep": 2
```

첫 주기는 신규/재시작 실행의 계산 루프 시작부터, 이후 주기는 마지막 저장 완료부터
21,000초(350분)를 셉니다. 시뮬레이션상의 물리 시간이 아니라 실제 경과 시간이며,
Slurm 대기 시간과 재시작 전 중단 시간은 포함하지 않습니다.
저장은 항상 **완료된 LB step 경계**에서 수행하므로 한 step이 오래 걸리면 350분을 넘을 수 있습니다.
**0 step과 재시작 직후에는 저장하지 않습니다.** 정상 종료, `--max-steps` 도달,
종료 요청 시 저장은 유지하고, 최근 완료본 두 개를 보관합니다.
각 주기를 0으로 설정하면 해당 주기적 저장만 끕니다. step 기준을 다시 활성화하려는 경우에만
`checkpoint_every_steps`를 양수로 지정합니다. 0 step의 CSV 기록은 계속 남습니다.

sbatch의 `#SBATCH --signal=B:USR1@180` 예고를 받으면 현재 step을 마친 뒤 저장하고
`CHECKPOINTED`로 종료합니다. 강제 종료가 먼저 발생하거나 입자 풀이가 실패하면 직전 완료
checkpoint를 사용합니다. 오래 걸리는 step과 파일 쓰기가 180초 안에 항상 끝나는 것은 아닙니다.

| 파일 | 보존 내용 |
| --- | --- |
| `lattice_rank_N.bin` | rank별 전체 유체 분포함수와 직렬화 대상 필드 |
| `state.bin` | 입자 위치·배향·속도·각속도, gap cache, 접촉·미끄럼·rolling 이력, 각가속도, 진단, step·누적 시간 |
| `checkpoint.json` | 호환 설정, 파일 길이, CSV 저장 경계 |
| `initial_particles.csv` | 원래 초기 입자 목록 |

LB step으로 `checkpoints/checkpoint_00000000000000000200`처럼 이름을 붙입니다.
모든 rank가 파일을 닫은 뒤 `.partial`을 완료 폴더로 바꾸고 `latest_checkpoint.txt`를 갱신합니다.
복원 시 바이너리 checksum을 검사합니다. 형식 1은 같은 OpenLB 자료형 배치와 MPI rank 수를
전제로 합니다. 유체 전체를 저장하므로 용량은 입자 수보다 격자 크기에 좌우됩니다.

## 4. 입자 solver, 실패 기록, 재현

PETSc 비선형 풀이, Newton line search, GMRES 선형 풀이와 접촉 블록 preconditioner를 사용합니다.
법선 접촉은 Fischer–Burmeister 상보성 잔차로 표현합니다. PETSc 종료 판정 후에도
물리적인 힘·토크·접촉 조건을 검사하고, 수렴한 substep만 접촉 이력에 반영합니다.
실패했다고 물리 허용오차를 자동 완화하거나 legacy solver로 돌아가지 않습니다.

모든 subdivision 시도 후 입자 단계가 실패하면 결과 폴더에 다음 진단을 남깁니다.

| 파일 | 내용 |
| --- | --- |
| `particle_solver_rankN_trace.csv` | 반복별 힘·토크·gap·상보성 잔차, line search, SNES/KSP 종료 이유, 접촉·미끄럼·rolling 전환 |
| `particle_solver_rankN_outcomes.csv` | 실패한 외부 시간 단계와 subdivision 수 |
| `particle_solver_rankN_failure.dat` | 마지막 실패 substep의 입자 상태·외력·접촉 이력·설정 |
| `particle_solver_summary.txt` | 실행 드라이버의 진단 요약 |

정상 실행과 subdivision으로 회복한 단계는 이 실패 파일을 쓰지 않습니다.
반복 trace는 제한된 메모리에 최근 기록부터 보관합니다. 강제 종료나 파일시스템 장애에는
진단이 남지 않을 수 있습니다. 이 파일은 **입자 substep 재현용**이며 유체 checkpoint와 다릅니다.
결과 폴더에서 `python3 summarize_particle_solver.py`로 요약을 다시 생성할 수 있습니다.

저장소 루트에서, 빌드 때와 같은 모듈로 실패 상태를 재현합니다.

```bash
source build/petsc/env.sh
python3 slurry/tools/compile_petsc_test.py tests/petsc_contact/replay_particle_step.cpp
build/petsc-tests/replay_particle_step /결과폴더/particle_solver_rank0_failure.dat
```

`--max-newton 80` 등의 replay 옵션은 원인 분리용 비교 설정이며 생산 JSON을 바꾸지 않습니다.

### MIT job 23255422: min-map의 힘·거리 스케일 수정

이전 min-map 수정만으로 실제 실행의 수렴 문제가 해결된 것은 아니었습니다.
이번 실행은 checkpoint 33129에서 재시작하여 33316까지 완료한 뒤, step 33317의
subdivision 256 / index 31에서 실패했습니다. 이전 실행보다 239 step 이릅니다.
현재 min-map 코드로 같은 실패를 재현했습니다: Newton 60, Krylov 842,
force/torque ratio 49.8724/37.6167. Gap 위반은 6.71e-16 m,
기존 FB complementarity ratio는 0.00127로, 법선 조건은 이미 충분히 작았습니다.

문제는 법선식의 형태뿐 아니라 **힘과 거리를 비교하는 수치 스케일**에도 있었습니다.
기존 solver 식 `min((h-h0)/eps_g, N/eps_F)`는 `eps_g/eps_F=10 m/N`으로
두 항을 환산합니다. 해당 substep의 코드상 관성 기준은
`C_ref = L/reactionScale = subdt^2/m_ref = 4.5071e-6 m/N`이므로 약 222만 배
차이가 났습니다. 여기서 `m_ref`는 기존 `reactionScale`에 쓰이는 첫 입자의 관성 질량입니다.

현재 solver 식과 Jacobian/preconditioner는 같은 관성 기준을 사용합니다.

```text
g = h - h0
C_ref = L / reactionScale = subdt^2 / m_ref
solver normal residual = min(g, C_ref*N) / contact_gap_tolerance
```

`C_ref>0`이므로 영점은 여전히 `g>=0, N>=0, g*N=0`입니다.
이 환산은 Newton 반복에만 사용하며 접촉 spring이나 새로운 힘을 추가하지 않습니다.
**최종 판정은 이전의 FB 함수와 force/torque/gap 허용오차, 음의 N 금지 검사 그대로**입니다.
물리 파라미터, dt, 최대 subdivision 256, Newton 60, Krylov 120도 유지합니다.
`residual_norm`의 수치 의미는 바뀌지만 `complementarity_ratio`는 기존 FB 기준입니다.

원인 분리 시험에서 NGMRES는 잔차 74.3을 56.0으로 낮춘 Newton 후보를 다시 361.5로
악화시키기도 했습니다. 그러나 plain Newton/backtracking도 마찰 경계에서 정체했고,
trust-region 비교안은 index 31을 통과한 뒤 index 32에서 실패했습니다.
따라서 solver 이름이나 line-search 옵션만 바꾸는 안은 포함하지 않았습니다.
현재 수정은 기존 NGMRES/SECANT 구조를 유지하며 위 정규화와 일관된 미분만 바꿉니다.

PETSc 3.25.5 serial 재현 결과(원래 물리 합격 기준, 각 LB endpoint까지 연속 계산):

| 입력 | 남은 substep | 최대 force ratio | 최대 torque ratio |
| --- | ---: | ---: | ---: |
| 이번 실패 index 31 | 225 | 0.9044 | 0.9562 |
| trust-region 비교안의 후속 실패 index 32 | 224 | 0.9276 | 0.8998 |
| 이전 네 실패 fixture | 517 | 0.9701 | 0.9649 |

총 966개 substep을 통과했습니다. 새 실패 단독은 Newton 12 / Krylov 190,
후속 실패 단독은 Newton 4 / Krylov 58입니다. 새 입력은
`normal_release_metric_108.dat`, `normal_release_followup_108.dat`이며
기존 `predictor_contact_tests.cpp`에 포함했습니다. 기존 접촉 회귀 25개도 통과했습니다.

이 검증은 **입자 실패 상태를 원래 LB step 끝까지 재현한 범위**입니다.
첨부 dump에는 재시작에 필요한 전체 LB 분포장이 없으므로, 유체를 갱신하는 장시간
완전 결합 계산의 안정성을 확인한 결과로 해석하면 안 됩니다.
이번 history의 마지막 fluid Ma는 0.0592, density drift는 11.71%였으며,
수렴 수정이 유체 정확도까지 보장하지는 않습니다.
설정 파일 변경 없이 `sbatch build_slurry_cpu.sbatch`로 빌드합니다.

### 이전 수정 기록 — MIT job 23181002: 법선 상보조건의 Newton 선형화와 마찰 전환

100/s, target Ma=0.02에서 step 33556을 계산하던 중 subdivision 256의 index 104에서
정체했습니다. 힘/토크 잔차 비율은 3.07228/1.01268이었고, Newton 60회 동안 수렴하지
않았습니다. SECANT를 full Newton으로 바꾼 비교 실험은 그 지점을 통과했지만 index 115에서
다시 실패했습니다. 따라서 line-search 종류 변경은 수정에 포함하지 않았습니다.

원인은 약한 법선 반력과 작은 gap 위반을 함께 푸는 Fischer–Burmeister(FB) 식의
유한 Newton 보정에 있었습니다. index 115의 입자 17–64에서 gap 위반은 1.0579 pm,
N은 0.01070 nN이었습니다. FB 선형화는 gap을 0.0196 pm만 복구하면서 N을 2.1338 nN
늘리면 상보성 잔차가 거의 사라진다고 예측했습니다. 실제 FB 잔차는 1.0631에서 1.0383으로
거의 줄지 않았고, 커진 Coulomb 상한이 sliding→sticking 전환을 일으켜 힘 잔차도 악화했습니다.
현재 분기 내부의 Jacobian이 맞아도 유한 보정의 예측은 이처럼 크게 틀릴 수 있습니다.

당시 수정은 `particlePetscSolver.h`의 solver용 법선 잔차와 그 미분이었습니다.
아래 당시 식의 힘 스케일은 위 job 23255422 수정으로 대체되었습니다.

```text
a = (h - h0) / contact_gap_tolerance
b = N / force_absolute_tolerance
solver normal residual = min(a, b)
derivative = (1, 0) if a <= b, otherwise (0, 1)
```

`min(a,b)=0`과 기존 FB=0은 모두 `h>=h0, N>=0, (h-h0)*N=0`을 뜻합니다.
음의 gap·양의 반력 분기에서는 잔차가 gap 자체여서, 반력만 늘려 gap 오차가 해결될 것처럼
예측하는 경로가 없어집니다. Jacobian과 block preconditioner에도 같은 미분을 사용합니다.
**최종 합격 판정과 trace의 complementarity ratio는 기존 FB 검사 그대로**입니다.
힘·토크·gap·음의 반력 검사, 물리 파라미터, 시간 간격, subdivision 수와 반복 상한도 같습니다.
`residual_norm`은 새 solver 잔차의 norm이므로 이전 FB solver norm과 직접 비교하지 않습니다.

독립 재현 결과(PETSc 3.25.5, serial, 기존 합격 기준):

| 입력 | Newton 횟수 | force ratio | torque ratio | gap 위반 |
| --- | ---: | ---: | ---: | ---: |
| 원본 실패 index 104 | 3 | 0.01399 | 0.009885 | 1.81e-17 m |
| full Newton 비교 실험의 후속 실패 index 115 | 11 | 0.001131 | 0.000963 | 4.24e-18 m |
| index 104–255 전체 152개 substep | 합계 231 | 최대 0.9253 | 최대 0.8877 | 최대 1.95e-15 m |
| 후속 실패 index 115–255 전체 141개 substep | 합계 232 | 최대 0.9619 | 최대 0.9964 | 최대 1.77e-15 m |

회귀 입력은 `normal_friction_corner_108.dat`와
`normal_friction_corner_followup_108.dat`이며, 기존 `predictor_contact_tests.cpp`가
이 두 상태와 이전 두 실패 상태를 각각 원래 LB endpoint까지 연속 계산합니다.
네 fixture의 총 517개 substep과 기존 접촉 회귀 25개가 통과했습니다.
통합 OpenLB 실행 파일의 serial/PETSc 빌드도 통과했습니다.

```bash
python3 slurry/tools/compile_petsc_test.py tests/petsc_contact/predictor_contact_tests.cpp --run
```

이 검사는 해당 LB step 동안 고정된 유체 하중을 사용합니다. 새 유체장을 계산하는 장시간
coupled restart의 검증은 별도이며, 이전 저전단 유체 Mach 폭주나 다른 domain 오류까지
해결했다는 뜻은 아닙니다. 이번 로그에는 step 33129 checkpoint 저장이 기록되어 있습니다.
기존 JSON 그대로 재빌드 후 해당 run의 checkpoint에서 이어 실행할 수 있습니다.

### 확인된 과거 오류: MIT job 23075858

10/s, LB step 3264, subdivision 256의 index 112에서 실패한 108입자 상태를 그대로 재현했습니다.
rank 0과 rank 10의 업로드 파일은 byte 단위로 같았습니다. GMRES 상대 기준 0.05가
첫 Newton 방향을 Krylov 2회 만에 받아들이면서 부정확한 방향을 사용했고, 비선형 잔차가
force/torque 비율 3.1189/4.9400에서 정체했습니다.

선형 상대 기준을 0.001로 강화한 뒤 같은 상태는 Newton 1회/Krylov 6회에 통과했습니다.
α=1, 힘 법칙, 물리 수렴 기준, substep 길이, 반복 상한은 그대로였습니다.
force/torque 비율은 0.2413/0.2025, gap 위반은 8.40e-13 m였습니다.
같은 frozen 유체 하중으로 남은 144개 입자 substep도 통과했습니다.
**이 결과는 해당 실패 상태의 해결이며, 이후 다른 접촉망 상태에서의 수렴까지 입증하지 않습니다.**

원본 fixture는 `tests/petsc_contact/fixtures/adhesive_network_108.dat`입니다.

```bash
python3 slurry/tools/compile_petsc_test.py tests/petsc_contact/replay_particle_step.cpp --run -- tests/petsc_contact/fixtures/adhesive_network_108.dat
```

성공 기준은 `REPLAY ACCEPTED`와 물리 잔차 통과입니다. 반복 횟수는 컴파일러/PETSc에 따라
달라질 수 있습니다. 이 명령은 저장된 입자 substep 하나를 풀며 앞선 전체 LB 경로를 재실행하지 않습니다.

### MIT job 23100417: 비침투 조건을 검사하지 않은 Newton 초기 예측값

이번 1/s, α=0.8 실행은 네 번째 LB step에서 실패했습니다. 최종 subdivision은 256,
실패 index는 176, subdt는 `4.5105489780439514e-7` s입니다. 첫 Newton 갱신 전에
GMRES가 120회 상한을 소진했고, 실패 시 sliding/rolling 접촉은 모두 0이었습니다.

직접 재현에서 입자 48–98은 저장 상태에서 아직 접촉 전이었습니다. 초기 속도 외삽이
이 쌍을 약 315 pm 가까이 옮겨, 기하학적 gap을 2.033 nm에서 1.718 nm로 줄였습니다.
`d = D0 + h - h0`는 0.018185 nm가 되어 순수 초기 추정의 pair 힘이 약 64.06 nN에서
123.879 N으로 증가했습니다. 이는 물리적으로 승인한 상태의 힘이 아니라
**비선형 풀이의 출발점이 국소 반발 특이점 가까이 들어간 결과**입니다.

기존 코드는 초기 예측값에서 d>0만 만족하면 이 상태를 허용했습니다. 그 위치에서
방향별 유한차분 Jacobian의 선형성이 나빠졌습니다. 원래 설정 그대로 재현한
GMRES의 120회 종료 시 내부 추정 잔차는 54.38이었으나, 실제 matrix-vector 연산으로
다시 계산한 잔차는 약 2.769e6이었습니다. 반복 상한을 늘리는 것으로 해결할 문제는 아닙니다.
이전 0.05→0.001 변경은 특정 선형 방향의 정확도 조정이었고 이 초기화 결함을 고치지 못했습니다.

수정은 `particlePetscSolver.h`의 초기값 생성 절차입니다.

1. 힘을 계산하기 전에 예측 위치·배향의 실제 gap을 검사합니다.
2. 기존 접촉 허용오차를 넘는 초기 변위는 축소하여 비침투 조건을 만족시키고,
   통과한 초기값에서 힘·토크 residual과 row scaling을 만듭니다.
3. Lees–Edwards 위상 진행 때문에 변위 0도 불가능한 경우, 초기 추정에만
   질량 가중 병진 복원을 수행하고 모든 pair의 gap을 다시 검사합니다.
   제한된 복원이 실패하면 해당 시도를 거부하고 기존 subdivision 절차로 돌아갑니다.

승인된 입자 상태와 접촉 이력을 투영하지 않습니다. A_H, σ, D0, α, μ, k_t, rolling 법칙,
시간 간격, 물리 허용오차, Newton/Krylov/substep 상한은 바꾸지 않았습니다.
미분용 작은 perturbation에는 기존의 매끄러운 trial-domain 확장을 유지합니다.
모든 Jacobian probe에 엄격한 비침투를 강제하면 다중 접촉 모서리에서 차분 자체가
불가능해질 수 있기 때문입니다.

새 fixture는 `tests/petsc_contact/fixtures/predictor_contact_108.dat`입니다.
동일 설정·동일 frozen 유체 하중으로 새 실패 이후 80개 substep과 이전 실패 이후
144개 substep을 모두 이어서 검사합니다.

```bash
python3 slurry/tools/compile_petsc_test.py tests/petsc_contact/predictor_contact_tests.cpp --run
```

| 입력 | 남은 substep | 최대 force ratio | 최대 torque ratio | 최대 gap 위반 |
| --- | ---: | ---: | ---: | ---: |
| job 23100417 | 80 | 0.6212 | 0.9271 | 0.6851 pm |
| job 23075858 | 144 | 0.6272 | 0.7577 | 0.8404 pm |

합격 기준은 각각 ratio≤1, gap 위반≤1 pm이며 그대로입니다. PETSc 접촉 회귀 25개,
국소 gap/replay 검사, 새 초기값·Lees–Edwards 복원 검사, Python 검사 42개도 통과했습니다.
의도적인 실패를 검사하는 회귀 fixture 한 개는 새 초기값에서 Newton 두 번에 수렴하므로
시험용 예산만 한 번으로 낮춰 실제 예산 소진과 rollback을 계속 검사합니다.
생산 실행의 60회 상한은 바뀌지 않았습니다.

전체 유체–입자 결합도 원래 108입자·200³ 격자, 1/s, α=0.8, μ=1, k_t=80,
같은 시간 간격과 초기 생성 seed로 8 LB step까지 실행했습니다. 초기 pair 수·gap·힘·응력은
업로드 실행과 약 1e-11 상대 오차 이내에서 일치했습니다. 기존 실패 지점인 네 번째 step을
지나 여덟 번째 step과 checkpoint 저장까지 정상 종료했습니다. 최대 force/torque ratio는
0.990935/0.952935, 최대 gap 위반은 5.408e-14 m로 기존 기준을 만족했습니다.
이 검사는 PETSc 3.25.5의 serial 실행이며 MIT의 32-rank 장시간 실행 검사는 아닙니다.

같은 짧은 실행에서도 `fluid_density_drift`는 8번째 step에서 0.2083,
`max_fluid_mach`는 실행 중 최대 0.2416이었습니다. 따라서 이번 통과는 보고된
입자 solver 실패의 회복을 확인한 결과입니다. 유체의 시간·속도 해상도와 장시간 물리 정확도가
충분하다는 판정으로 확대하지 않습니다.

## 5. 국소 부착 모델의 정의와 범위

기하학적 접촉 간격 h₀는 유지하면서 국소 표면이 D₀까지 접근하는 효과를 유효 포텐셜로
표현합니다. 조도 형상을 직접 해상하거나 해당 graphite 분말의 측정 부착력에 맞춘 모델은 아닙니다.
α=1은 다음 혼합식에서 국소 상호작용을 전부 사용하는 기준 조건이며, 실제 접촉 면적률 ψ가 아닙니다.

q를 배향과 나머지 기하 정보라 하면:

```text
s = h − h₀
d = D₀ + s
U(h,q) = U_RE²(h,q) + α S(s) [U_RE²(d,q) − U_RE²(h,q)]
```

두 RE² 항은 같은 Hamaker 상수와 σ, 배향 인자를 쓰며 인력과 반발을 모두 포함합니다.
α=1, S=1에서는 기존 근접 기여를 국소 기여로 교체합니다.
기하학적 구속조건은 계속 h≥h₀이고 입자 모양이나 모멘트 팔은 줄이지 않습니다.

S=1은 s≤2 nm, S=0은 s≥10 nm에 적용하며 중간은
`t=(s−2 nm)/(8 nm)`, `S=1−10t³+15t⁴−6t⁵`입니다.
즉 실제 h=4 nm부터 감쇠해 h=12 nm에서 보정이 끝납니다.
기존 장거리 switch/cutoff 400/500 nm는 유지합니다.
힘과 양쪽 입자의 토크는 switch와 배향을 포함한 전체 에너지의 미분으로 구합니다.
출력의 인력·반발 성분, 에너지와 virial도 같은 보정 평가를 사용합니다.

적분된 LJ 12-6 평판 에너지의 평형 간격은 `D₀=(2/15)^(1/6) σ`입니다.
D₀=0.30 nm에 대응하는 σ=0.419725 nm를 입력에서는 0.4197 nm로 반올림했습니다.
반발을 포함한 기준 부착일은 `w=A_H/(16πD₀²)≈21.9 mJ/m²`입니다.
이는 가정한 분자 간격에서 유도한 기준이며 독립적인 분말 측정값은 아닙니다.

반축 (1.65, 1.65, 0.20) µm의 정렬된 접촉에서 계산한 순인력은 다음과 같습니다.

| 배향 | 이전 σ=3 nm, α=0 | 현재 σ, α=0 | 현재 σ, α=1 |
| --- | ---: | ---: | ---: |
| FF | 17.469 nN | 28.073 nN | 936.215 nN |
| EF | 0.48275 nN | 0.77837 nN | 25.951 nN |
| EE | 0.25504 nN | 0.41247 nN | 13.748 nN |

FF 값은 유효 곡률 반경 약 6.806 µm와 `2πR_eff w`에서 얻는 크기와 일치합니다.
이를 벌크 항복응력의 같은 배율 증가나 500 Pa 재현으로 해석하지 않습니다.

Sliding 상한은 μN이며 새 pair 힘이 법선 평형의 반력 N에 들어갑니다.
Rolling 상한은 `rolling_length × adhesiveBirthForce`이며, 생성 시 보정된 총 순인력을 사용하고
생성 후 cap과 강성을 고정하는 기존 방식은 유지합니다.
kₜ=80 N/m, N=936 nN, μ=1인 FF 접촉의 항복 접선변위는 약 11.7 nm입니다.
kₜ는 항복 전 탄성 변위를 정하며 **μN 자체를 높이지 않습니다.**
Newton 중간 시도의 d≤0은 거부하고 수락 상태에는 h≥h₀를 요구합니다.

물리식의 출처는 [Everaers–Ejtehadi RE²](https://doi.org/10.1103/PhysRevE.67.041710),
[Li et al.의 graphite–water Hamaker 상수](https://doi.org/10.1103/PhysRevB.71.235412)입니다.
국소 간격 혼합식과 α=1의 선택은 이 코드에서 채택한 유효 접촉 가정입니다.

## 6. 검증 명령과 확인 범위

일반 실행에는 아래 개발용 검사를 따로 빌드할 필요가 없습니다.

```bash
python3 -m unittest discover -s slurry/tests
g++ -std=c++17 -O2 -I olb-1.9r0/src/slurry/gr_re2 tests/graphite_adhesion/adhesion_tests.cpp -o /tmp/graphite_adhesion_tests
/tmp/graphite_adhesion_tests
source build/petsc/env.sh
python3 slurry/tools/compile_petsc_test.py tests/petsc_contact/contact_solver_tests.cpp --run
python3 slurry/tools/compile_petsc_test.py tests/petsc_contact/local_gap_solver_tests.cpp --run
```

의도적으로 직렬 MPIUNI PETSc를 쓰는 개발 환경에서는 C++ 검사 명령에 `--compiler g++ --serial`을
추가합니다. MPI 생산 빌드에서 이 옵션을 사용할 필요는 없습니다.

재시작 검사는 접촉하는 두 입자의 연속 16 step과 `8 step → 저장 → 16 step`을 비교합니다.

```bash
python3 tests/checkpoint/check_restart.py --work /tmp/gr_restart_check
python3 tests/checkpoint/check_stop.py --source /tmp/gr_restart_check/part --config /tmp/gr_restart_check/config.json
```

MPI 바이너리에서는 첫 명령에 `--ranks 2`를 추가합니다. 두 번째 명령은 source checkpoint의 rank 수를 읽습니다. 첫 검사는 활성 접촉과 접촉 탄성 에너지를 확인하고,
시간 측정 열을 제외한 물리 CSV와 최종 rank별 lattice 파일을 비교합니다.
기본 350분 주기에서는 짧은 실행이 종료될 때만 저장되는지, 0 step과 재시작 직후 저장이
없는지도 확인합니다. 첫 명령에 `--checkpoint-seconds 0.001`과 별도 `--work` 경로를 주면
검사용 시간 간격만 줄여 시간 기준 저장 경로를 실행할 수 있습니다.
두 번째 검사는 주기 저장이 없는 상태에서 실제 종료 신호 후 저장·재시작 및 손상 파일 거부를 확인합니다.

### 기록된 검증 결과

아래는 각 변경 당시 실제 수행한 결과입니다. 현재 모든 MPI 생산 조건을 검증했다는 뜻은 아닙니다.
직렬 환경은 Linux/G++ 13.3, real double PETSc 3.25.5 MPIUNI이며,
재시작 MPI 검사는 별도 OpenMPI 4.1.4와 legacy 입자 backend를 사용했습니다.

| 변경 | 수행한 검증 |
| --- | --- |
| PETSc 도입 | 전체 OpenLB 직렬 빌드·링크, 4입자/80³/7 step 통합 실행, 물리 fixture 14개, 진단 10개, 설치 보존 10개 통과; 저장 실패 재현 |
| 빌드 진입점 통합 | 의존성 준비·환경 보존·캐시 재사용·명령 전달 등 mock 기반 검사 10개 통과 |
| MIT MPI 사전검사 수정 | launcher/오류 구분/설치 재사용 등 mock 검사 8개, 실제 직렬 PETSc 초기화 통과 |
| 국소 부착 | 에너지 미분·각운동량·switch·한계 검증 6그룹, PETSc 접촉 25개, Python 26개, 전체 직렬 빌드 통과 |
| 국소 부착 통합 | 두 FF 입자/80³/100/s/7 step 통과; 최종 pair 힘 936.2185 nN, gap≈2 nm, 최대 gap 위반 9.94e-14 m |
| job 23075858 수정 | 업로드 실패 substep 및 남은 144개 substep 통과; PETSc 25개/Python 26개와 전체 직렬 재빌드 통과 |
| checkpoint/restart | Python 36개(재시작 10개 포함), 전체 직렬 PETSc 및 2-rank OpenMPI legacy 빌드 통과 |
| 재시작 동등성 | 직렬 PETSc 및 2-rank legacy 각각 물리 CSV 전체 일치, 최종 rank별 lattice byte 일치 |
| 종료·손상 처리 | 위 두 환경에서 SIGUSR1 전달·CHECKPOINTED 종료·복구, 동일 길이 바이너리 손상 거부 통과 |
| 350분 저장 주기 | 직렬 PETSc 재빌드·Python 42개 통과. 기본 주기에서 시작 저장 없음·종료 저장·정확한 재시작 확인. 검사용 0.001초 주기로 시간 기준 저장 확인. SIGUSR1 저장·재시작·손상 거부 통과 |

재시작 단위검사는 미완료 `.partial`, 누락 rank 파일, checkpoint 없는 이전 실행,
물리값/dt/rank 호환성, 누적 종료점, 저장 CSV 경계, 혼합 전단율 batch의 완료 케이스 건너뛰기도
확인합니다. 단위·소규모 검사는 전체 108입자 MPI 경로의 장시간 안정성이나 실험 유변학 검증과는
구분합니다. MIT 32-rank 생산 실행과 실제 Slurm timeout 전달은 위 검증에 포함되지 않습니다.

job 23075858의 제공된 history는 10/s 계산으로 마지막 표본 응력이 243.99 Pa였고,
별도로 언급된 100/s 결과는 포함하지 않았습니다. 이 history에는 최대 fluid Mach 0.622와
밀도 편차 39.0%도 기록되어 있습니다. 입자 substep 재현 통과만으로 이러한 유체 시간/속도
척도의 문제까지 해결됐다고 판단하지 않습니다.

## 7. CMC 계수의 출처와 선택적 재피팅

`slurry/fitting/cmc250k_digitized_black_lines.csv`는 논문 SI CMC 250k 패널의 검은 피팅 선에서
추출한 좌표입니다. LBM 결과가 아닙니다. 같은 폴더의 `cmc250k_digitization_metadata.json`에
원본 SI SHA-256, 축 보정과 농도 배정이 있고, `cmc250k_cross_parameters.json`에 당시 피팅
계수·설정·Python 환경이 기록되어 있습니다.

LBM이 읽는 물성 파일은 **`slurry/cases/cmc250k_cross_parameters.csv`**입니다.
피팅 폴더의 JSON을 바꿔도 실행 물성은 바뀌지 않으며 보통 실행에는 재피팅이 필요 없습니다.
계수를 다시 추정할 때만 Python 3.12의 별도 환경에서 실행합니다.
MIT에서는 환경 준비와 피팅을 할당된 계산 노드에서 진행합니다.

```bash
cd slurry/fitting
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements-fit.txt
python3 fit_cmc250k.py
```

결과는 `reproduced_cmc250k/`에 저장됩니다. seed는 41이고 각 농도에서 여러 초기값을 사용한
log-viscosity least squares입니다. 원본 환경은 Python 3.12.14, NumPy 2.3.5,
SciPy 1.17.0, Matplotlib 3.10.8입니다. LBM 실행 driver는 Python 3.6 이상 표준 라이브러리만 씁니다.
원본 SI에서 추출부터 반복하려면 해당 PDF와 PyMuPDF가 추가로 필요합니다.

```bash
python3 fit_cmc250k.py --pdf /path/to/nn6c10201_si_001.pdf --outdir reproduced_from_pdf
python3 check_eta_inf_sensitivity.py
```

두 번째 명령은 `sensitivity_reproduced/`에 무한전단 점도 민감도 결과를 저장합니다.
계수는 인쇄된 피팅 선을 재현하며 저자의 원 계수나 실험 불확실도를 복원한 값은 아닙니다.
일부 `eta_inf=0`은 최적화 하한에 도달한 결과이고 농도 사이의 연속 보간 함수는 정의하지 않았습니다.

## 8. 소스 위치

| 경로 | 내용 |
| --- | --- |
| `olb-1.9r0/src/slurry/` | 물리 모듈 |
| `olb-1.9r0/examples/slurry/` | 통합 C++ 실행 진입점 |
| `slurry/cases/` | 실행 설정 |
| `slurry/drivers/`, `slurry/tools/` | 공통 빌드·실행 및 모델 driver |
| `slurry/fitting/` | CMC 피팅 입력·출처·재현 도구 |
| `tests/`, `slurry/tests/` | 물리·solver·재시작 및 Python 검증 |
| `runs/` | 실행 결과, Git 제외 |

과거 파일 이관 경로와 원본 hash는 `slurry/source_manifest.json`에 기록되어 있습니다.
그 기록에 있는 이전 설명서 경로는 역사적 출처이며 현재 설명서는 이 README로 통합했습니다.
