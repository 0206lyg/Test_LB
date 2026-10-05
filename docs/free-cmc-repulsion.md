# CMC 90k: 기존 반발 모델과 고농도 유효 순포텐셜

## 기본 설정과 모델의 의미

`slurry/cases/gr_CMC.json`은 16 g/L에서 전체 유효 순포텐셜을 교체한다.
`gr_CMC_4g_L.json`은 이전 설정을 유지한다. 기계적 접촉은 두 경우 모두 탄소 간격
h0=2 nm이다. 흡착량을 전단에 따라 바꾸거나 별도의 CMC 점성저항을 추가하지 않는다.

새 고농도 모델의 9 nm cutoff와 힘의 크기는 접촉 활성화·해제를 위한 모델 선택값이다.
**9 nm에서 실제 van der Waals 인력이 없어짐을 의미하지 않으며, 측정된 CMC 차폐식이 아니다.**
완전 적용 시 RE² 꼬리도 함께 제거한다. 기존 인력을 밑에 남기면 아래 표와 다른 모델이다.

## 농도와 보간

흡착량과 free 농도는 수상 1 L 기준의 독립 입력이다. `1 g/L = 1 kg/m³`이다.
흡착량이 포화값을 넘으면 흡착 계산에만 포화값을 쓰며 초과량을 free로 옮기지 않는다.

| 총농도 기준 | adsorbed CMC (g/L) | free CMC (g/L) |
| --- | ---: | ---: |
| 4 g/L | 2.8984214285714285 | 1.1015785714285715 |
| 16 g/L | 2.8984214285714285 | 13.101578571428572 |

포화값은 Gwag et al., DOI [10.1021/acsnano.6c10201](https://doi.org/10.1021/acsnano.6c10201)의
Gr+CB 대비 겉보기 흡착량 약 0.37 ± 0.09 wt%를 기준으로,
`0.0037 × 997 × 0.44 / 0.56`으로 옮긴 모델 기준값이다.
논문이 직접 보고한 수상 2.8984 g/L 측정값은 아니다.

`cmc.net_potential`을 생략하거나 `enabled=false`로 두면 새 에너지는 사용하지 않는다.
활성화한 경우 다음 농도 보간을 사용한다.

\[
w=\operatorname{clip}\left(\frac{c_f-c_{\rm start}}
 {c_{\rm full}-c_{\rm start}},0,1\right),\qquad
U=(1-w)U_{\rm previous}(H,c_f)+wU_{\rm net}(H).
\]

`U_previous`는 같은 입력 농도의 기존 RE²+흡착 접착+free 반발 모델이다.
4 g/L의 포텐셜을 모든 중간 농도에 고정해서 쓰는 것이 아니다.
`w=0`은 기존 계산 경로를 그대로 실행하며, `w=1`은 이전 에너지를 계산·상쇄하지 않고
새 에너지를 직접 계산한다. 따라서 완전 적용 시 반올림으로 남는 RE² 꼬리가 없다.
두 농도 사이의 선형 보간은 구성식 선택이다. 중간 농도에는 기존 꼬리와 안정점이
남을 수 있으며 아래의 단순한 힘 부호 구간은 `w=1`의 기준 입자쌍에 해당한다.

```json
"net_potential": {
  "enabled": true,
  "start_free_g_L": 1.1015785714285715,
  "full_free_g_L": 13.101578571428572,
  "contact_force_N": 1.5e-10,
  "barrier_force_N": 2e-11,
  "attraction_range_m": 1e-9,
  "repulsion_range_m": 6e-9
}
```

## 새 순포텐셜의 정의

La=1 nm, Lr=6 nm, Hs=h0+La=3 nm, Hr=Hs+Lr=9 nm이다.
FA=150 pN, FB=20 pN이며 기준 곡률 길이는 입력 입자 반경·반두께에서
lambda_ref=a²/c=13.6125 µm로 계산한다.

\[
U_{\rm net}=\frac{\lambda^D_{ij}}{\lambda_{\rm ref}}(U_A^0+U_R^0),\qquad
U_A^0=-\frac{F_A L_a}{3}\left[1-\frac{H-h_0}{L_a}\right]_+^3.
\]

반발 에너지는 다음과 같다. z=(Hr-H)/Lr이다.

\[
U_R^0=\begin{cases}
\frac{8}{15}F_B L_r,&H\le H_s,\\
16F_B L_r\left(\frac{z^3}{3}-\frac{z^4}{2}+\frac{z^5}{5}\right),&H_s<H<H_r,\\
0,&H\ge H_r.
\end{cases}
\]

정면 face-to-face에서는 lambda_D=lambda_ref이고, 분리 방향을 양으로 잡은 힘은
2–3 nm에서 `-FA(1-x)²` (x=(H-h0)/La), 3–9 nm에서
`16 FB y²(1-y)²` (y=(H-Hs)/Lr), 9 nm 이상에서 0이다.

3 nm와 9 nm에서 에너지 및 1·2차 미분이 연속이다. 3 nm는 장벽의 최대점이고,
그 안쪽에는 비접촉 안정점이 없다. 9 nm 밖은 상호작용이 없는 평탄한 영역이다.
양의 간격 `H<h0` solver trial에는 안쪽 다항식을 연장하되, 허용되는 접촉 상태는
기존 `H>=h0` 제약을 따른다.

임의 배향에서는 곡률 길이의 위치·회전 미분까지 포함한다. 힘만 사후에 잘라내거나
토크를 별도 법칙으로 주지 않는다. attractive/repulsive 출력은 각각 위 두 에너지의
미분이며, 임의 배향에서는 곡률 변화에 따른 성분도 포함한다.

| 정렬 입자쌍 | 최대 접근 반발 (pN) | 최대 분리 인력 (pN) | 접근 일 (J) | 접촉 탈출 일 (J) |
| --- | ---: | ---: | ---: | ---: |
| Face–face | 20.0 | 150.0 | 6.40e-20 | 5.00e-20 |
| Edge–face | 0.5545 | 4.1589 | 1.7745e-21 | 1.3863e-21 |
| Edge–edge | 0.2938 | 2.2039 | 9.4031e-22 | 7.3462e-22 |

이 표는 입자쌍의 보존력 검증값이다. 16 g/L의 1 Pa·s를 계산한 벌크 유변학 결과가 아니다.
접촉 에너지는 무한 분리 상태보다 FF 기준 1.4e-20 J 높다. 접촉은 장벽 안의 국소 상태이며
분산 상태보다 낮은 에너지의 영구 응집 상태로 설계하지 않았다.
기존 모델처럼 Brownian 힘은 추가하지 않았다. EF/EE 에너지 장벽은 298.15 K의 kBT보다
작으므로 이 표를 열적 결합 수명의 예측으로 해석하지 않는다.

## 보존하는 기존 모델 (w=0)

흡착 CMC는 추가 접착력을 줄인다.

\[
\theta_a=\min(c_a/c_{a,\rm sat},1),\quad
q_a=[1-(1-\sqrt{q_{\rm sat}})\theta_a]^2,\quad q_{\rm sat}=0.33.
\]

배경 표면 일은 `Wbg=AH[1-(sigma/h0)^6/30]/(12 pi h0²)`이고,
추가 접착력은 `DeltaW=qa(Wbare-Wbg)`를 쓴다.
`interaction.adhesion_work_J_m2`는 항상 bare 입력으로 유지하여 재입력 때
흡착 감소를 두 번 적용하지 않는다.

\[
U_{\rm adh}=-\frac{\pi\lambda^D_{ij}\Delta W\delta_a}{2}
 [1-(H-h_0)/\delta_a]_+^2,\qquad\delta_a=0.67\ {\rm nm}.
\]

이전 free CMC 구성식은 외측 반발과 내측 반발이다.

\[
P_f=\texttt{strength}\,\chi RT\frac{\mathrm{DS}\,c_f}{M_0},\qquad
U_{\rm out}=\pi\lambda^D_{ij}P_f\ell^2e^{-(H-h_0)/\ell},
\]
\[
W_s=\kappa c_f,\qquad
U_{\rm in}=\frac{\pi\lambda^D_{ij}W_s\delta_i}{p}
 [1-(H-h_0)/\delta_i]_+^p.
\]

기존 설정은 strength=1, ell=5 nm, chi=0.5, DS=0.7,
M0=0.218 kg/mol, kappa=0.000543218 (J/m²)/(g/L), delta_i=0.67 nm, p=2.1이다.
명목 전하 농도를 이용한 유효 압력 환산이며 측정한 미시적 표면력 그 자체가 아니다.
`free_repulsion.enabled`는 기존 두 반발항만 제어한다. 새 모델을 끄려면
**`net_potential.enabled=false`**를 사용한다.
기존 inner 계수를 생략하면 0이고, 새 net 객체를 생략하면 이전 모델 그대로다.

## 설정 전달, 수치 허용오차, 재시작

새 C++ 입력은 `cmc_net_blend`, `cmc_net_contact_force`, `cmc_net_barrier_force`,
`cmc_net_attraction_range`, `cmc_net_repulsion_range`, `cmc_net_reference_length`이다.
빌드 기능 키 `cmc_net_potential_version`은 1이다. rough_contact와 surface_adhesion을
사용하고 접촉면 이동/기존 cohesion-retention과 혼용하지 않는다.
두 새 범위의 합과 h0는 기존 far switch보다 작거나 같아야 한다.

활성 순포텐셜의 여섯 계수와 버전은 checkpoint signature에 포함된다.
순포텐셜 비활성 시 해당 키를 signature에 넣지 않아 이전 계산과의 호환성을 보존한다.
서로 다른 물리 설정의 checkpoint는 섞지 않는다. 이전 16 g/L 결과는 새 계산으로 시작한다.
Particle replay는 새 버전 7을 쓰며 기존 버전 1–6을 읽을 수 있다.

16 g/L에서 절대 힘 허용오차는 1e-15 N, 토크는 1.65e-21 N·m이다.
기존 1e-13 N은 약 0.3 pN인 새 EE 장벽의 약 1/3이므로 낮췄다.
상대 허용오차, solver 알고리즘, 초기 입자 배치, friction/rolling, lubrication,
연속상 점도는 유지한다. 접촉 수명 등의 새 시계열 로그는 추가하지 않는다.
4 g/L JSON은 이번 순포텐셜 변경에서 그대로 보존한다.

```bash
sbatch build_slurry_cpu.sbatch
# 빌드 완료 후
sbatch run_slurry_cpu.sbatch --cases gr_cmc --shear-rates 10,100
```

입력/manifest의 계수와 빌드 버전으로 실행 모델을 구별한다.
