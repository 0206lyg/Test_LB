# 흡착 CMC와 접촉을 유지하는 free-CMC 반발

90k CMC–graphite 모델은 흡착 CMC에 따른 추가 접착력 감소, 기존 외측 지수형 반발,
접촉 근처에만 작용하는 내측 반발을 사용한다. 기계적 접촉은 탄소 표면 간격
\(h_0=2\ \mathrm{nm}\)에서 유지한다. 피복으로 접촉면을 4 nm까지 옮기거나 free CMC가
전체 인력을 곱셈 계수로 차폐하는 구성식은 사용하지 않는다.

## 입력과 농도

`slurry/cases/gr_CMC.json`은 `pure_gr.json`을 자동 상속하지 않는 독립 설정이다.
기본 CMC 설정은 다음과 같다.

```json
"cmc": {
  "adsorbed_g_L": 2.8984214285714285,
  "free_g_L": 13.10157857142857,
  "adsorbed_saturation_g_L": 2.8984214285714285,
  "q_sat": 0.33,
  "free_repulsion": {
    "enabled": true,
    "strength": 1.0,
    "decay_length_m": 5e-9,
    "degree_of_substitution": 0.7,
    "repeat_unit_molar_mass_kg_mol": 0.218,
    "osmotic_coefficient": 0.5,
    "inner_work_per_g_L_J_m2": 0.000543218,
    "inner_range_m": 6.7e-10,
    "inner_exponent": 2.1
  }
}
```

흡착량과 free 농도는 **수상 1 L 기준의 독립 입력**이다. 흡착량은 포화값까지만
적용하며, 초과 입력을 free CMC로 옮기지 않는다. `free_g_L`에는 상한을 두지 않는다.
`1 g/L = 1 kg/m³`이므로 흡착량을 다시 빼거나 수상 부피분율을 곱하지 않는다.

| 총농도를 나누는 예 | 흡착량 (g/L) | free 농도 (g/L) |
| --- | ---: | ---: |
| 4 g/L | 2.8984214285714285 | 1.1015785714285715 |
| 16 g/L, 기본 설정 | 2.8984214285714285 | 13.10157857142857 |

이는 입력 예이며 드라이버가 질량 수지를 자동 재분배하지는 않는다. 포화값은
[Gwag et al., DOI 10.1021/acsnano.6c10201](https://doi.org/10.1021/acsnano.6c10201)의
Gr+CB 대비 겉보기 흡착량 약 0.37 ± 0.09 wt%를 모델 기준으로 환산한 값이다.

\[
x_{a,\mathrm{sat}}=0.0037\times997\times\frac{0.44}{0.56}
=2.8984214285714285\ \mathrm{g/L}.
\]

작은 CMC 질량을 총 고형분 환산에서 생략했으며, 논문이 직접 보고한 수상 농도는 아니다.
`gr_CMC.json`(16 g/L)과 `gr_CMC_4g_L.json`(4 g/L)의 초기 배치 최소 간격은 모두 10 nm이다.
이 초기 배치 조건은 기계적 접촉 간격 2 nm와 별개이다.

## 흡착 CMC: 같은 포화 잔류 접착력

\[
\theta_a=\min(x_a/x_{a,\mathrm{sat}},1),\qquad
q_a=\left[1-(1-\sqrt{q_{\mathrm{sat}}})\theta_a\right]^2,
\]
\[
W_{\mathrm{bg}}(h_0)=\frac{A_H}{12\pi h_0^2}
\left[1-\frac{(\sigma/h_0)^6}{30}\right],\qquad
\Delta W_{\mathrm{ads}}=q_a[W_{\mathrm{bare}}-W_{\mathrm{bg}}(h_0)].
\]

`interaction.adhesion_work_J_m2`에는 bare 기준값을 입력한다. 드라이버가 전달하는
\(W_{\mathrm{bg}}+\Delta W_{\mathrm{ads}}\)에서 C++이 배경 일을 한 번 빼서 추가
접착력을 계산한다. `q_sat=0.33`은 포화 시 **추가 접착력**의 잔존율이며,
RE² 배경 인력 전체의 잔존율이나 실측 벌크 응력비가 아니다. 포화 흡착량을 쓰는
4·16 g/L 조건의 이 항은 동일하다.

탄소 표면 간격을 \(H\), 접촉으로부터의 개방 거리를 \(s=H-h_0\)라 두면 기존
추가 접착 포텐셜은 다음과 같다. \([z]_+=\max(z,0)\)이다.

\[
U_{\mathrm{adh}}=-\frac{\pi\lambda_{ij}^{D}\Delta W_{\mathrm{ads}}\delta_a}{2}
\left[1-\frac{s}{\delta_a}\right]_+^2,
\qquad \delta_a=0.67\ \mathrm{nm}.
\]

## Free CMC: 외측 반발과 내측 반발

외측 반발은 `7004d4762b`의 압력 환산과 지수형 포텐셜을 유지한다.

\[
P_f=\texttt{strength}\;\chi RT\frac{\mathrm{DS}\,c_f}{M_0},\qquad
U_{\mathrm{out}}=\pi\lambda_{ij}^{D}P_f\ell^2\exp(-s/\ell).
\]

\(c_f\)는 kg/m³, \(M_0\)는 kg/mol이며, 기본값은
`strength=1`, \(\chi=0.5\), \(\mathrm{DS}=0.7\), \(M_0=0.218\ \mathrm{kg/mol}\),
\(\ell=5\ \mathrm{nm}\)이다. \(P_f\)의 단위는 Pa이다.
이는 명목 전하 농도를 이용한 유효 반발 압력 척도이며, 측정된 이온 농도나
미시적인 electrostatic force law를 그대로 재현한 식은 아니다.

내측 반발은 같은 기계적 접촉면을 기준으로 유한 범위에서만 작용한다.

\[
\boxed{U_{\mathrm{in}}=
\frac{\pi\lambda_{ij}^{D}W_s\delta_i}{p}
\left[1-\frac{s}{\delta_i}\right]_+^{p}},\qquad
W_s=\kappa c_f,
\]
\[
\delta_i=0.67\ \mathrm{nm},\qquad p=2.1,\qquad
\kappa=0.000543218\ \frac{\mathrm{J/m^2}}{\mathrm{g/L}}.
\]

여기서 \(c_f\)는 g/L 단위 입력이다. \(W_s\)는 J/m²이며 \(U_{\mathrm{in}}\)은 J이다.
\(p>2\)이므로 내측 cutoff에서 에너지와 1·2차 미분이 연속적으로 0이 된다.
\(\kappa\), 범위, 지수는 접촉 인력을 남기면서 근접 응집 에너지를 줄이기 위한
모델 선택값이며 측정된 CMC 물성값은 아니다.

전체 포텐셜은

\[
U_{ij}=U_{\mathrm{RE^2}}+U_{\mathrm{adh}}+U_{\mathrm{out}}+U_{\mathrm{in}}.
\]

기존 RE² 곡률 보정과 원거리 switch를 유지한다. 각 반발항도 스칼라 포텐셜에
더하므로 간격뿐 아니라 이동하는 접촉점의 곡률 길이 \(\lambda_{ij}^{D}\)까지
자동 미분하여 힘과 토크를 얻는다. 계산 후 인력을 잘라내거나 접촉 법선력만
상쇄하지 않는다. 전단률에 따라 계수나 흡착량을 강제로 바꾸는 항은 없다.

`strength`는 외측 반발만 조절한다. `enabled=false`는 두 반발항을 모두 끈다.
기존 JSON에서 `inner_work_per_g_L_J_m2`를 생략하면 기본값은 **0**이므로 기존 외측
반발만 작동한다. 위 기본 CMC 파일은 내측 계수를 명시하여 새 항을 활성화한다.

## 입자쌍 기준 확인값

다음은 포화 흡착, 정면 face-to-face 배향, \(T=298\ \mathrm K\), 위 계수의 입자쌍
기준값이다. 양의 힘은 분리 방향 반발, 음의 힘은 인력이다. 기본 실행의 free 농도
13.10157857 g/L 및 온도 298.15 K와 구별한다.

| 항목 | \(c_f=13.2\ \mathrm{g/L}\) |
| --- | ---: |
| 접촉점 순힘 | −10.000 nN |
| 분리 중 최대 인력의 크기 | 17.121 nN |
| 접촉 상태에서 장벽까지의 탈출 에너지 | \(12.224\times10^{-18}\) J |
| 외측 인력 우물에서 장벽까지의 진입 에너지 | \(13.067\times10^{-18}\) J |
| 접근 경로의 최대 반발력 | 1.921 nN |

같은 조건에서 \(c_f=1.1\ \mathrm{g/L}\)의 접촉점 순힘은 −301.38 nN이다.
고농도에서도 접촉점 인력과 유한한 탈출 장벽을 남기며, 내측 범위
\(0<s<\delta_i\)에 별도의 비접촉 안정 간격을 만들지 않는 선택이다.
이 입자쌍 확인은 다입자 구조, 4 g/L의 yield stress 또는 16 g/L의
약 1 Pa·s Newtonian plateau를 검증한 결과가 아니다.

## 접촉, 호환성, 실행

기계적 접촉은 기존 \(H\ge h_0\), \(N\ge0\), \(N(H-h_0)=0\)을 유지한다.
Sliding·rolling 법칙, 마찰계수, 접선 강성, 연속상 점도, lubrication 및 입자 solver는
변경하지 않는다. 기존 adhesive rolling 기준력은 기존 방식대로 생성 시의
순접착력에서 정해지며 새 반발항의 영향은 그 힘에 반영된다.

옛 `cmc.free_cohesion`과 `cmc.contact_offset_at_saturation_m` 입력은 드라이버가
거부한다. 해당 키를 제거하고 위 `free_repulsion` 객체를 사용한다.
CMC 객체가 없거나 흡착·free 농도가 모두 0이면 기존 pure-Gr 힘 법칙을 회복한다.
흡착량을 유지하고 free 농도만 0으로 만들면 반발항만 사라진다.

새 내측 C++ 입력은 `free_cmc_inner_repulsion_work`, `free_cmc_inner_repulsion_range`,
`free_cmc_inner_repulsion_power`이며 빌드 정보의
`free_cmc_inner_repulsion_version`은 1이어야 한다. 외측 입력은 기존
`free_cmc_repulsion_pressure`, `free_cmc_repulsion_length`를 유지한다.
`manifest.json`의 `derived.cmc.free_repulsion`에는 환산 압력·내측 표면 일 척도와
활성 여부를 기록한다. 활성 물리 계수는 checkpoint signature에도 반영한다.
**동일 물리 설정의 checkpoint만 재시작**하며 농도나 포텐셜을 바꾼 비교는 새 계산으로
수행한다. 비활성 CMC 항은 기존 pure-Gr signature를 바꾸지 않는다.

소스 변경 후 기존 명령으로 한 번 재빌드한 다음 실행한다.

```bash
sbatch build_slurry_cpu.sbatch
sbatch run_slurry_cpu.sbatch --cases gr_cmc --shear-rates 100
```

두 번째 명령은 빌드 완료 후 제출한다. 이후 JSON 값만 바꿀 때는 재빌드가 필요 없다.
