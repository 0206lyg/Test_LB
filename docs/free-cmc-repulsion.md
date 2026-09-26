# Free CMC에 따른 유효 반발 포텐셜

90k CMC–graphite 모델에 free CMC 농도에 비례하는 보존 반발항을 추가한다.
흡착 CMC에 따른 접착력 감소는 기존 식을 유지하며, free CMC 농도는 별도로 입력한다.
이 문서는 구현한 최소 모델과 단위 환산을 설명한다. 여기 제시한 시험 매개변수로
slurry의 yielding–Newtonian 전이를 재현했다는 뜻은 아니다.

## 입력과 실행

`slurry/cases/gr_CMC.json`의 `cmc` 설정은 다음과 같다. 이 객체 외의 graphite,
물성, 입자, 접촉 및 수치 설정은 해당 case 파일에 계속 필요하다.

```json
"cmc": {
  "adsorbed_g_L": 3.8,
  "free_g_L": 0.0,
  "adsorbed_saturation_g_L": 2.8984214285714285,
  "q_sat": 0.33,
  "free_repulsion": {
    "enabled": true,
    "strength": 1.0,
    "decay_length_m": 5e-9,
    "degree_of_substitution": 0.7,
    "repeat_unit_molar_mass_kg_mol": 0.218,
    "osmotic_coefficient": 0.5
  }
}
```

`adsorbed_g_L`과 `free_g_L`의 분모는 모두 수상 1 L이다. 흡착량에만
`adsorbed_saturation_g_L` 상한을 적용하며 초과분을 free 농도로 옮기지 않는다.
기본 입력 3.8 g/L의 실제 적용 흡착량은 2.8984214 g/L이고 free 농도는 0이다.
Free 농도에는 별도의 포화 상한을 두지 않는다.

| 비교 기준 | 총농도 4 g/L에 대응하는 free 입력 | 총농도 16 g/L에 대응하는 free 입력 |
| --- | ---: | ---: |
| 흡착량을 2.8 g/L로 반올림한 논의 | 1.2 | 13.2 |
| 현재 case의 포화량 2.8984214 g/L 사용 | 1.1015786 | 13.1015786 |

이 표는 사용자가 총농도를 독립 입력으로 나눌 때 참고하는 계산이다. 코드가
`free_g_L`에서 흡착량을 다시 빼지는 않는다. 용액 질량분율, 전체 slurry 기준 농도
또는 총 free CMC 질량을 사용하는 실험 데이터라면 먼저 수상 기준 g/L로 바꿔 입력한다.

이번 소스 변경 후 한 번 재빌드한다. 이후 JSON 농도·반발 매개변수만 바꾸면
새 실행부터 적용되므로 다시 컴파일할 필요가 없다.

```bash
sbatch build_slurry_cpu.sbatch
# BUILD COMPLETE / BUILD REUSED를 확인한 뒤 실행
sbatch run_slurry_cpu.sbatch --cases gr_cmc --shear-rates 100
```

Solver 실행 전에 변환된 설정을 확인하려면 다음 명령을 사용한다.

```bash
python3 slurry/drivers/gr_re2/run_graphite.py --config slurry/cases/gr_CMC.json --dry-run
```

## 농도와 압력의 단위 환산

`free_g_L`은 질량이 아니라 이미 질량농도이다.

\[
1\ {\rm g/L}=\frac{10^{-3}\ {\rm kg}}{10^{-3}\ {\rm m^3}}
=1\ {\rm kg/m^3}.
\]

따라서 코드에서 SI 질량농도 \(c_f\)의 수치는 `free_g_L`과 같다. 입자
질량분율, 수상 부피분율, 밀도 등을 추가로 곱하지 않는다. 반복단위와 명목 전하의
몰농도는 각각 다음과 같다.

\[
c_{\rm repeat}=\frac{c_f}{M_0}\quad[\mathrm{mol/m^3}],\qquad
c_{\rm charge}^{\rm nominal}=DS\frac{c_f}{M_0}\quad[\mathrm{mol/m^3}].
\]

명목 전하 농도는 모든 전하가 독립적인 자유 이온처럼 작용한다는 주장이 아니다.
유효 osmotic coefficient를 곱하여 기준 압력을 정의한다.

\[
\Pi_b(c_f)=\phi_{\rm osm}\frac{DS}{M_0}RTc_f,\qquad
P_R(c_f)=\alpha_R\Pi_b(c_f).
\]

\(R=8.31446261815324\ \mathrm{J/(mol\,K)}\)를 사용한다. 기본값
\(DS=0.7\), \(M_0=0.218\ \mathrm{kg/mol}\), \(T=298.15\ \mathrm K\),
\(\phi_{\rm osm}=0.5\)에서 \(\Pi_b/c_f\)는 약
**3979.98 Pa/(g/L)**이다. `strength`는 무차원 \(\alpha_R\)이다.
온도는 기존 `fluid.temperature_K`를 사용하므로 별도의 반발 온도 입력을 두지 않는다.

| 입력/변환값 | 단위 |
| --- | --- |
| `free_g_L` | g/L aqueous phase |
| \(c_f\) | kg/m³ aqueous phase |
| \(M_0\) | kg/mol repeat unit |
| \(c_{\rm repeat}\), \(c_{\rm charge}^{\rm nominal}\) | mol/m³ |
| \(\Pi_b\), \(P_R\) | Pa |
| \(\ell\), \(\lambda_{ij}^{D}\) | m |
| \(W_f\) | J/m² |
| \(\Phi_f\) | J |

## 보존 포텐셜과 힘

탄소 표면 간격을 \(H\), 기존 rough contact 기준 간격을 \(h_0\)라 두고
물리적인 비중첩 구간에서 \(s=H-h_0\ge0\)를 정의한다. 유효 초과압력,
평판 상호작용 에너지 및 입자쌍 에너지는 다음과 같다.

\[
\Pi_f(s,c_f)=P_R(c_f)e^{-s/\ell},
\]
\[
W_f(s,c_f)=\int_s^\infty\Pi_f(u,c_f)\,du
=P_R(c_f)\ell e^{-s/\ell},
\]
\[
\Phi_f(s,c_f)=\pi\lambda_{ij}^{D}\int_s^\infty W_f(u,c_f)\,du
=\pi\lambda_{ij}^{D}P_R(c_f)\ell^2e^{-s/\ell}.
\]

기존 원거리 switch를 \(\psi_{\rm out}\), switch 전의 기존 배경·흡착 접착
포텐셜을 \(\Phi_{\rm old}^{\rm inner}\)라 하면

\[
\Phi_{ij}=\psi_{\rm out}(H)
\left[\Phi_{\rm old}^{\rm inner}(H,x_a)+\Phi_f(H,c_f)\right].
\]

고정된 곡률, switch가 1인 구간의 정상 반발력은
\(F_f=\pi\lambda_{ij}^{D}P_R\ell e^{-s/\ell}\)이다. 실제 구현은 이 스칼라
힘만 별도로 더하지 않고 전체 포텐셜을 자동 미분한다. 간격, 국소 접촉 곡률 및
원거리 switch의 미분을 포함하므로 힘·토크가 같은 에너지에서 나온다.

기존 흡착 법칙은 그대로이다.

\[
\theta=\min(x_a/x_{a,\rm sat},1),\qquad
q=\left[1-(1-\sqrt{q_{\rm sat}})\theta\right]^2,
\]
\[
W_{\rm eff}=W_{\rm bg}+q(W_{\rm bare}-W_{\rm bg}).
\]

`interaction.adhesion_work_J_m2`에는 계속 bare graphite의 기준 부착일을 입력한다.
Free 반발을 보상하기 위해 이 부착일을 자동으로 다시 조정하지 않는다. 기존 접촉
간격, 수직 접촉 법칙, sliding friction coefficient, 접선 강성 및 rolling 법칙은
유지한다. 새 항이 입자 배치와 접촉 하중을 바꾸면 결과적인 마찰 응력은 달라질 수 있다.
연속상 및 lubrication 점도는 `free_g_L`에 따라 자동 변경하지 않는다.

## 시험 매개변수와 물리적 해석

`decay_length_m=5e-9`와 `strength=1`은 합의한 시험 출발값이다. 90k CMC–graphite의
직접 표면력 측정으로 얻은 값이 아니다. \(\phi_{\rm osm}=0.5\)도 유효 근사이며
측정한 반발압과 벌크 삼투압을 동일시하지 않는다. \(\alpha_R\)는 벌크 압력 척도가
실제 접근 반발에 연결되는 정도를 대표한다.

흡착 포화 이후의 안정화와 charged polymer의 표면력 연구가 이 유효 반발항의
물리적 동기이다. 단일 지수형은 그 효과를 최소한의 매개변수로 표현하는 선택이다.
사슬 overlap 농도에서 반발을 켜는 스위치는 없다. 농도에 따라 반발은 매끄럽게
변하지만 기존 인력과 경쟁하여 접근 장벽이 생기면 접촉망은 급격하게 달라질 수 있다.

이 항의 주된 역할은 접촉 생성·재형성을 억제하는 것이다. 이미 형성된 강한 접착
접촉망이 농도 증가만으로 즉시 분리되지는 않으며 초기 구조와 preshear 이력의
영향을 평가해야 한다. 4 g/L yielding과 16 g/L Newtonian-like 거동은 실제 다입자
계산에서 검증할 대상이다. \(G'\), \(G''\)를 이 항의 보정 대상으로 사용하지 않는다.

## 기록과 재시작

`effective_config.json`에 독립 입력과 반발 매개변수를 보존한다. `manifest.json`의
`derived.cmc`에는 적용 흡착량과 잔존 접착력을 기록한다. 그 아래 `free_repulsion`
객체에는 다음 환산 결과가 포함된다.

| Manifest 키 | 의미 |
| --- | --- |
| `temperature_K` | 기존 fluid 설정에서 읽은 계산 온도 |
| `free_concentration_kg_m3` | SI free 질량농도 |
| `repeat_unit_concentration_mol_m3` | 반복단위 몰농도 |
| `nominal_charge_concentration_mol_m3` | DS를 곱한 명목 전하 몰농도 |
| `bulk_osmotic_pressure_Pa` | \(\Pi_b\) |
| `effective_repulsion_pressure_Pa` | 활성 설정을 적용한 \(P_R\) |
| `surface_energy_at_contact_J_m2` | \(P_R\ell\) |
| `pair_energy_per_derjaguin_length_J_m` | \(\pi P_R\ell^2\) |
| `active` | 실제 반발항의 활성 여부 |
| `model_version` | free 반발 모델 버전 |

활성 반발의 압력과 길이는 `resolved_run.cfg`의 `free_cmc_repulsion_pressure`와
`free_cmc_repulsion_length`로 전달하므로 C++ 계산에서 농도 환산을 다시 하지 않는다.

`free_g_L=0`, `strength=0` 또는 `enabled=false`이면 반발항이 0이며 기존 pair
법칙을 유지한다. CMC 객체가 없는 `pure_gr`에도 이 반발항이 적용되지 않는다.
기존 free 반발 0인 checkpoint는 새 항도 0인 설정으로 재시작할 수 있다.
활성 반발압이나 감쇠 길이 또는 흡착 접착력이 달라지는 설정으로 기존 checkpoint를
이어 돌리는 것은 물리값 호환성 검사에서 거부한다. 이러한 농도 비교는 새 계산으로
수행한다. 재시작은 기존과 같이 `--cases gr_cmc --restart 폴더명`으로 지정한다.

## 관련 근거

- Gwag et al., ACS Nano, [DOI: 10.1021/acsnano.6c10201](https://doi.org/10.1021/acsnano.6c10201):
  흡착·free polymer의 역할과 90k slurry의 유변학적 비교. 현재 흡착 포화량 환산의 원자료.
- Moazzami-Gudarzi et al., Physical Review Letters 117, 088001 (2016),
  [DOI: 10.1103/PhysRevLett.117.088001](https://doi.org/10.1103/PhysRevLett.117.088001):
  비흡착성 charged polymer 용액의 표면력. 해당 계의 함수·매개변수를 CMC로 직접 이전하지 않는다.
- Schön and von Klitzing, Beilstein Journal of Nanotechnology 9, 1095–1107 (2018),
  [DOI: 10.3762/bjnano.9.101](https://doi.org/10.3762/bjnano.9.101):
  다른 분산계에서 추가 지수형 유효 반발을 사용하는 실험적 선례.
