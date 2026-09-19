# pure_gr checkpoint / restart

## 설치와 실행

이 업데이트 ZIP을 `/home/lyjania/OpenLB` 안에서 압축 해제한 뒤,
평소처럼 한 번 다시 빌드합니다. ZIP은 기존 `slurry/cases/pure_gr.json`을
덮어쓰지 않습니다. 아래 저장 옵션은 JSON에 없어도 기본값으로 적용됩니다.

```bash
cd /home/lyjania/OpenLB
unzip -o OpenLB_pure_gr_restart.zip
sbatch build_slurry_cpu.sbatch
```

빌드 로그의 `BUILD COMPLETE` 또는 `BUILD REUSED`를 확인한 뒤 처음 실행:

```bash
sbatch run_slurry_cpu.sbatch --cases pure_gr --shear-rates 100
```

예를 들어 이 실행 결과가 `runs/slurry_123456_20260919T120000Z`라면:

```bash
sbatch run_slurry_cpu.sbatch --restart slurry_123456_20260919T120000Z
```

`runs/`가 포함된 상대 경로, 절대 경로, 개별 전단율 폴더,
특정 `checkpoints/checkpoint_...` 폴더도 지정할 수 있습니다.
`--restart`만 사용하면 대상 모델은 자동으로 `pure_gr`입니다.
CMC 계산은 시작하지 않습니다.

재시작 결과는 새 `runs/slurry_...` 폴더에 저장됩니다. 원래 결과는 보존하고,
`history.csv`와 `particles.csv`는 선택된 저장 지점까지 복사한 뒤 이어 씁니다.
저장 이후 실패한 실행에서 남긴 행은 복사하지 않으며, 저장 지점의 행을
중복 출력하지 않습니다. 이전 VTK 파일은 원래 폴더에 남습니다.
checkpoint만 별도로 옮겨 CSV가 없는 경우에는 저장 시각부터 새 CSV를 만듭니다.

여러 전단율 폴더가 있으면 각 폴더의 최신 완료 checkpoint를 선택합니다.
현재 설정의 종료 지점에 이미 도달한 케이스는 건너뜁니다. 아직 한 번도
실행하지 않아 checkpoint가 없는 전단율은 별도의 신규 실행 대상입니다.

## 현재 pure_gr.json이 적용되는 범위

기본 입력은 `/home/lyjania/OpenLB/slurry/cases/pure_gr.json`입니다.
`--config` 또는 공통 `run.json`에서 다른 파일을 지정하면 그 파일을 사용합니다.
전단율은 각 checkpoint에 저장된 값을 유지합니다.

| 항목 | 재시작 시 처리 |
| --- | --- |
| `flow.end_strain` | 현재 JSON 사용. 시작 이후 **총 누적 strain**의 목표 |
| `--max-steps` | 시작 이후 **총 LB step 수**의 상한 |
| 출력·checkpoint 간격과 보관 개수 | 현재 JSON 사용 |
| 입자 solver 반복 상한·수렴 허용오차·진단 옵션 | 현재 JSON 사용 |
| 물리값, 입자 수·크기, 격자, 시간 간격, MPI rank 수 | 저장 당시와 같아야 함 |
| 입자 배치 seed·초기 minimum gap | 배치를 새로 생성하지 않고 저장 상태 사용 |

예를 들어 strain 4.2에서 중단했고 `end_strain=10`이면 4.2부터 10까지
진행합니다. strain 10에서 완료한 계산을 20까지 연장하려면 현재 JSON의
`flow.end_strain`을 20으로 바꾸거나 다음 명령을 사용합니다.

```bash
sbatch run_slurry_cpu.sbatch --restart 폴더명 --end-strain 20
```

`local_gap_fraction`(α), `sigma_lj_m`, `local_gap_m`, 마찰·접촉 강성 등을
바꾸면 같은 물리계의 연속 계산이 아니므로 재시작 전에 차이를 보고하고
중단합니다. 변경된 인력으로 비교하려면 해당 JSON으로 신규 실행합니다.
이 업데이트는 기존 인력 크기나 PETSc 수렴 판정을 바꾸지 않습니다.

동일한 rank 수를 유지해야 rank별 격자 배치가 일치합니다. 보통 기존과 같은
`run_slurry_cpu.sbatch`를 쓰면 됩니다. 이전 실행에 `--ranks`를 사용했다면
재시작에도 같은 값을 전달합니다.

## 자동 저장

아래 키를 기존 JSON의 `output` 객체에 추가하면 저장 주기를 바꿀 수 있습니다.

```json
"checkpoint_every_steps": 200,
"checkpoint_every_seconds": 900.0,
"checkpoint_keep": 2
```

- 200 LB step 간격 또는 직전 저장 후 wall time 900초 중 먼저 해당하는 시점에 저장합니다.
- 모든 저장은 **완료된 LB step의 경계**에서 실행합니다. 한 step이 오래 걸리면 900초를 넘을 수 있습니다.
- 시작 시점, 정상 종료, `--max-steps` 도달 때도 저장합니다.
- 기존 sbatch의 `#SBATCH --signal=B:USR1@180` 종료 예고를 받으면 현재 step을 끝내고 저장한 뒤 `CHECKPOINTED`로 종료합니다.
- 최근 완료 checkpoint 2개를 보관합니다. 각 주기를 0으로 설정하면 해당 주기적 저장만 끕니다. 시작·종료·신호 저장은 유지됩니다.

강제 종료가 저장 완료보다 빠르면 직전 완료 checkpoint로 돌아갑니다.
입자 solver 실패 시에도 진행 중인 불완전한 LB 상태를 저장하지 않고,
이전 완료 checkpoint를 사용합니다. 오래 걸리는 한 step이나 파일 기록을
Slurm의 180초 예고 시간 안에 항상 끝낼 수 있다는 보장은 없습니다.

checkpoint는 입자 수보다 **전체 유체 격자 크기**가 용량을 좌우합니다.
rank별 격자 전체를 저장하므로 CSV보다 큽니다. 저장 간격과 보관 수는
실행 폴더의 실제 `checkpoints/` 용량을 보고 조절할 수 있습니다.

## 보존하는 상태와 파일

- 유체 lattice의 분포함수와 직렬화 대상 필드: `lattice_rank_N.bin`
- 입자 위치·배향·속도·각속도, pair gap cache, 활성 접촉·미끄럼·rolling 이력,
  각가속도와 마지막 입자 solver 진단: `state.bin`
- LB step, 누적 시간·strain, Lees–Edwards 위상, 누적 실행 시간
- 초기 입자 목록과 호환성·파일 길이·CSV 경계를 기록한 `checkpoint.json`

`checkpoints/checkpoint_00000000000000000200`처럼 LB step으로 이름을 붙입니다.
모든 rank의 파일 쓰기가 끝난 뒤에만 `.partial` 폴더를 완료 폴더로 바꾸고
`latest_checkpoint.txt`를 갱신합니다. 바이너리 파일은 checksum을 검사한 뒤
복원하므로 같은 길이의 손상도 검출합니다. 현재 형식은 같은 OpenLB 버전,
자료형 배치와 MPI rank 수를 전제로 합니다.

## 이전 실행의 한계

**이 업데이트로 생성한 checkpoint가 있는 실행부터 재시작할 수 있습니다.**
기존 `history.csv`, `particles.csv`, 입자 실패 replay 파일에는 전체 유체 상태가
없으므로 이전 실행을 정확한 연속 상태로 변환할 수 없습니다.

## 검증 재현

빌드 후 작은 접촉 케이스로 연속 16 step과 `8 step → 저장 → 16 step`을 비교:

```bash
python3 tests/checkpoint/check_restart.py --work /tmp/gr_restart_check
```

MPI 빌드라면 `--ranks 2`를 추가합니다. 이 fixture는 해당 바이너리에 맞는
입자 solver를 사용합니다. 물리 CSV 값과 최종 rank별 lattice를 정확히 비교하고,
활성 접촉과 저장된 접촉 탄성 에너지가 실제로 존재하는지도 검사합니다.

실제 종료 신호와 이후 복구, 의도적으로 손상시킨 상태의 거부를 검사:

```bash
python3 tests/checkpoint/check_stop.py \
  --source /tmp/gr_restart_check/part \
  --config /tmp/gr_restart_check/config.json
python3 -m unittest discover -s slurry/tests -p 'test_graphite_restart.py'
```

MIT Slurm의 32-rank 생산 실행은 이 소규모 검증에 포함되지 않습니다.
