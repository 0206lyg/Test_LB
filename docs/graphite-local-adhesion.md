# Graphite 국소 접근거리와 부착 보정

`pure_gr`는 기하학적 접촉 간격을 유지하면서 국소 graphite 표면이 더 가까이
접근하는 경우를 유효 포텐셜로 표현한다. 이 변경은 조도 형상을 직접 해상하는
모델이나 측정한 graphite 입자의 부착력 보정 결과가 아니다. `alpha=1`은 아래
혼합 모델에서 국소 상호작용을 전부 사용하는 강한 기준 조건이다.

## 실행 설정

`slurry/cases/pure_gr.json`의 단위는 SI이다.

| JSON 설정 | 값 | 의미 |
| --- | ---: | --- |
| `interaction.hamaker_J` | `9.9e-20` J | 기존 Hamaker 상수 유지 |
| `interaction.sigma_lj_m` | `4.197e-10` m | LJ 반발 길이, 약 0.4197 nm |
| `interaction.local_gap_m` | `3e-10` m | 기하학적 접촉에서의 국소 간격 D0 |
| `interaction.local_gap_fraction` | `1.0` | 국소 상호작용의 유효 기여율 alpha |
| `interaction.local_switch_excess_gap_m` | `2e-9` m | h-h0에 대한 보정 감쇠 시작 |
| `interaction.local_cutoff_excess_gap_m` | `1e-8` m | h-h0에 대한 보정 종료 |
| `rough_contact.roughness_gap_m` | `2e-9` m | 기존 기하학적 접촉 간격 h0 유지 |
| `rough_contact.sliding_friction` | `1.0` | 기존 마찰계수 유지 |
| `rough_contact.tangential_stiffness_N_m` | `80.0` N/m | 접선 탄성 강성 |
| `rough_contact.rolling_length_m` | `1e-7` m | 기존 rolling 길이 유지 |
| `rough_contact.rolling_yield_angle_rad` | `0.01` rad | 기존 rolling 항복각 유지 |

기존 400/500 nm 장거리 switch/cutoff, 물의 점도, 입자 크기, 입자 수,
초기 seed, 시간 간격 및 수렴 허용오차는 변경하지 않는다. 새 국소 간격 항이
없는 이전 설정은 `local_gap_fraction=0`으로 처리되어 기존 RE²를 사용한다.

## 포텐셜

`q`를 두 입자의 배향 및 나머지 기하 정보라 하면

```text
s = h - h0
d = D0 + s
U(h,q) = U_RE2(h,q) + alpha S(s) [U_RE2(d,q) - U_RE2(h,q)]
```

두 RE² 항은 동일한 Hamaker 상수와 sigma, 배향 인자를 사용하며 인력과 반발을
모두 포함한다. `alpha`는 실제 접촉 면적률 psi가 아니다. `alpha=1, S=1`에서는
기존 기여를 더하는 대신 국소 기여로 교체하므로 같은 상호작용을 두 번 세지
않는다. 기하학적 구속조건은 계속 `h >= h0`이고 입자 모양이나 모멘트 팔을
줄이지 않는다.

`S=1`은 `s <= 2 nm`, `S=0`은 `s >= 10 nm`에서 적용한다. 중간 구간은
`t=(s-2 nm)/(8 nm)`에 대해 `S=1-10t^3+15t^4-6t^5`이다. 따라서 보정은
실제 간격 h가 4 nm일 때 감쇠를 시작하고 12 nm에서 끝난다. 기존 장거리
switch는 그대로 사용한다.

힘과 두 입자의 토크는 switch 및 모든 배향 인자를 포함한 전체 에너지의
미분으로 구한다. 접촉 여부만으로 상수 힘을 켜거나 힘을 임의로 제한하지
않는다. 출력의 attractive/repulsive force, energy 및 stress는 각각 보정된
인력/반발 성분이며 합이 전체 pair 상호작용이다.

## 길이와 힘의 기준

적분된 LJ 12-6 평판 에너지에서 국소 평형 간격은
`D0=(2/15)^(1/6) sigma`이다. D0=0.30 nm는 sigma=0.419725 nm에 대응하며,
입력은 합의한 반올림 값 0.4197 nm를 사용한다. 반발을 포함한 부착일 기준은
`w=AH/(16 pi D0^2)`, 약 21.9 mJ/m²이다. 이는 가정한 분자 간격에서 유도한
기준값이며, 해당 분말의 독립적인 측정값은 아니다.

현재 축 길이 (1.65, 1.65, 0.20) µm의 정렬된 접촉에서 코드 계산값은 다음과 같다.

| 배향 | 기존 sigma=3 nm, alpha=0 | 새 sigma, alpha=0 | 새 sigma, alpha=1 |
| --- | ---: | ---: | ---: |
| FF | 17.469 nN | 28.073 nN | 936.215 nN |
| EF | 0.48275 nN | 0.77837 nN | 25.951 nN |
| EE | 0.25504 nN | 0.41247 nN | 13.748 nN |

FF 최대 인력은 곡률 반경 약 6.806 µm에서 `2 pi R_eff w`와 일치하는 크기다.
이 강도 증가가 suspension 응력의 같은 배율 증가나 500 Pa 재현을 보장하지는 않는다.

## 접촉 법칙 및 solver 연결

- Sliding 한계는 기존 `mu N`이다. 새 pair 힘이 법선 평형의 반력 N에 들어간다.
- Rolling 한계는 기존 `rolling_length * adhesiveBirthForce`이고, 접촉 생성 시
  보정된 총 pair 인력을 사용한다. 생성 후 cap과 강성의 동결 방식은 유지한다.
- 접선 강성은 80 N/m이다. N=936 nN, mu=1인 FF 기준으로 항복 접선변위는
  약 11.7 nm이며, 강성 자체가 `mu N`의 상한을 바꾸지는 않는다.
- 에너지와 virial은 같은 보정된 pair 평가 결과에서 가져온다.
- Newton/PETSc의 중간 시도에서도 `d > 0`을 지킨다. 허용되지 않는 간격의
  시도는 solver에서 거부하며, 실제 수락된 접촉 조건 `h >= h0`는 유지한다.
- 실패 재현 파일은 새 매개변수를 보존한다. 이전 replay 파일은 국소 보정이
  꺼진 설정으로 읽는다.

## 빌드 및 실행

저장소 루트에서 기존 명령을 사용한다. 새 JSON을 이전 바이너리에 적용하는
실수는 build-info의 국소 부착 지원 여부를 검사하여 차단한다.

```bash
sbatch build_slurry_cpu.sbatch
```

빌드 로그에서 `BUILD COMPLETE` 또는 `BUILD REUSED`를 확인한 뒤:

```bash
sbatch run_slurry_cpu.sbatch --cases pure_gr --shear-rates 100
```

## 구현 확인

PETSc 없이 포텐셜과 접촉 법칙을 확인할 수 있다.

```bash
g++ -std=c++17 -O2 -I olb-1.9r0/src/slurry/gr_re2 tests/graphite_adhesion/adhesion_tests.cpp -o /tmp/graphite_adhesion_tests
/tmp/graphite_adhesion_tests
python3 -m unittest discover -s slurry/tests
```

PETSc 설치 환경에서 기존 solver 회귀 검증을 실행한다.

```bash
source build/petsc/env.sh
python3 slurry/tools/compile_petsc_test.py tests/petsc_contact/contact_solver_tests.cpp --run
python3 slurry/tools/compile_petsc_test.py tests/petsc_contact/local_gap_solver_tests.cpp --run
```

물리 근거: [Everaers and Ejtehadi, RE²](https://doi.org/10.1103/PhysRevE.67.041710),
[Li et al., graphite-water Hamaker constant](https://doi.org/10.1103/PhysRevB.71.235412).
위 국소 간격 혼합식과 alpha=1의 선택은 본 모델의 명시적인 유효 접촉 가정이다.
