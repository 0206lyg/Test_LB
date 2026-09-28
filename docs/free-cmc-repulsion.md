# CMC 피복 접촉과 free-CMC 유효 응집 인력 차폐

이 문서는 90k CMC–graphite 모델의 현재 구성식을 설명한다. 파일 이름은 기존 링크를
유지하기 위한 것이며, 현재 모델은 이전의 지수형 반발 포텐셜을 대체한다.
흡착 CMC는 접촉면과 추가 접착 에너지를 바꾸고, free CMC는 전체 응집 인력을
차폐한다. 고농도에서 피복층끼리 압축 접촉하여 기존 마찰을 전달하는 것이 시험할
기계적 그림이다. 목표 점도 또는 yielding–Newtonian 전이를 이미 재현했다는 뜻은 아니다.

## 입력

`slurry/cases/gr_CMC.json`은 pure-Gr 설정을 복사한 독립 파일이다. `cmc` 객체 외의
입자·유체·상호작용·수치 설정도 필요하며, `pure_gr.json`을 자동 상속하지 않는다.

```json
"cmc": {
  "adsorbed_g_L": 3.8,
  "free_g_L": 0.0,
  "adsorbed_saturation_g_L": 2.8984214285714285,
  "q_sat": 0.33,
  "contact_offset_at_saturation_m": 2e-9,
  "free_cohesion": {
    "enabled": true,
    "nonadhesive_concentration_g_L": 13.2
  }
}
```

두 농도는 수상 1 L 기준으로 독립 입력한다. 흡착량만 포화값에서 제한하며 초과량을
free CMC로 옮기지 않는다. 기본 흡착량 입력 3.8 g/L의 적용값은 2.8984214 g/L이고
free 입력은 0으로 남는다. `free_g_L` 자체에는 상한이 없다.

\[
1\ \mathrm{g/L}=1\ \mathrm{kg/m^3}.
\]

따라서 SI 농도는 수치가 같고, 흡착량을 다시 빼거나 수상 부피분율을 곱하지 않는다.
다른 분모로 측정한 농도라면 사용자가 먼저 수상 기준 g/L로 환산한다.

| 총농도를 나누는 기준 | 4 g/L의 free 입력 | 16 g/L의 free 입력 |
| --- | ---: | ---: |
| 흡착량을 2.8 g/L로 설정 | 1.2 | 13.2 |
| 현재 포화값 2.8984214 g/L 사용 | 1.1015786 | 13.1015786 |

이 표는 질량 수지의 예이며, 드라이버가 두 독립 입력을 자동 재분배하지 않는다.
특히 두 번째 행의 13.1015786 g/L는 기본 차폐 끝점 13.2 g/L보다 작으므로 엄밀히
같은 `q_f=0` 끝점은 아니다. 끝점 검증에서는 의도한 free 입력과
`nonadhesive_concentration_g_L`을 일치시킨다.

## 흡착량, 피복 접촉면, 접착 에너지

탄소 표면 간격을 \(H\), bare 접촉 기준을 \(h_0\), 피복의 기계적 접촉 기준을
\(H_c\)라 둔다.

\[
\theta_a=\min(x_a/x_{a,\mathrm{sat}},1),\qquad
H_c=h_0+\Delta H_{\mathrm{sat}}\theta_a.
\]

`contact_offset_at_saturation_m`은 두 피복층이 합쳐서 더하는 유효 접촉 거리
\(\Delta H_{\mathrm{sat}}\)이다. 단일 층 두께가 아니다. 기본값 2 nm는 시험 가정이며,
포화 CMC 층 두께를 직접 측정한 값이 아니다. 흡착 포화 후에는 \(H_c\)가 고정되고
free CMC는 이 거리 기준을 바꾸지 않는다.

현재 포화 흡착량은 Gwag et al.의 Gr+CB 대비 겉보기 흡착량 약 0.37 ± 0.09 wt%를
모델의 고형분 질량분율 0.44와 물 밀도 997 g/L로 환산한 것이다.

\[
x_{a,\mathrm{sat}}=0.0037\times997\times\frac{0.44}{0.56}
=2.8984214\ \mathrm{g/L}.
\]

논문이 직접 보고한 수상 농도가 아니며, 작은 CMC 질량을 총 고형분 환산에서 생략했다.
분말·고형분 기준을 바꾸면 이 입력도 재검토한다.

흡착에 따른 추가 접착력 잔존율은 기존 식을 유지한다.

\[
q_a=\left[1-(1-\sqrt{q_{\mathrm{sat}}})\theta_a\right]^2,
\]
\[
W_{\mathrm{bg}}(h_0)=\frac{A_H}{12\pi h_0^2}
\left[1-\frac{(\sigma/h_0)^6}{30}\right],\qquad
\Delta W_{\mathrm{ads}}=q_a[W_{\mathrm{bare}}-W_{\mathrm{bg}}(h_0)].
\]

`interaction.adhesion_work_J_m2`에는 bare 기준값을 입력한다. 드라이버는
\(W_{\mathrm{bg}}(h_0)+\Delta W_{\mathrm{ads}}\)를 C++에 전달하고, C++은 같은 bare
배경 일을 빼서 \(\Delta W_{\mathrm{ads}}\)를 얻는다. 접촉 위치가 변해도 이 에너지
기준은 그대로이다. 입력을 저장하고 다시 읽어도 흡착 감소를 두 번 적용하지 않는다.

기존 유한 범위 접착 포텐셜의 거리 기준만 \(s_c=H-H_c\)로 옮긴다. 곡률 근접 구간에서
기존 Derjaguin 길이 \(\lambda_{ij}^{D}\)를 쓰면 다음과 같다.

\[
\Phi_{\mathrm{adh}}(s_c;x_a)=
-\frac{\pi}{2}\lambda_{ij}^{D}\Delta W_{\mathrm{ads}}\delta_a
\left(1-\frac{s_c}{\delta_a}\right)^2,
\qquad 0\le s_c<\delta_a,
\]
\[
\Phi_{\mathrm{adh}}=0,\qquad s_c\ge\delta_a.
\]

RE² 배경은 실제 탄소 간격 \(H\)를 계속 사용한다. 접착 구간이 곡률 근접 구간을
벗어나지 않도록 \(H_c+\delta_a\le\texttt{curvature_switch_gap_m}\)를 검사한다.
초기 배치 최소 간격이 \(H_c\)보다 작으면 그 값까지만 높이고 유효 설정에 남긴다.
입자 축 길이·질량·관성 또는 유체 격자 형상을 피복 두께만큼 팽창시키지는 않는다.

`q_sat=0.33`은 추가 접착력의 포화 잔존율에 대한 초기 모델값이다. 직접 측정한 CMC
상수나 벌크 응력비를 뜻하지 않는다.

## Free CMC의 최소 차폐 구성식

\[
u=\min(x_f/x_f^\star,1),\qquad
q_f=1-10u^3+15u^4-6u^5=(1-u)^3(1+3u+6u^2).
\]

\(x_f^\star\)는 `nonadhesive_concentration_g_L`이며 기본값은 13.2 g/L이다.
활성 상태에서 \(x_f=0\)이면 \(q_f=1\), \(x_f\ge x_f^\star\)이면 \(q_f=0\)이다.
`free_cohesion.enabled=false`이면 농도와 무관하게 \(q_f=1\)이다. 양 끝에서 1·2차
미분이 0인 함수로 보간하며, 코드에서는 음수가 되는 반올림 오차를 피하기 쉬운
곱 형태를 사용한다. 농도 입력에 상한을 두는 것이 아니라 차폐 효과가 포화된다.

최종 입자쌍 포텐셜은

\[
\boxed{
\Phi_{ij}=q_f\left[\Phi_{\mathrm{RE^2,att}}(H)
+\Phi_{\mathrm{adh}}(H-H_c;x_a)\right]
+\Phi_{\mathrm{RE^2,rep}}(H).
}
\]

각 항은 기존 곡률 보정과 원거리 switch를 포함한다. 추가 접착력뿐 아니라 RE²의
인력도 같은 잔존율로 차폐하므로, 완전 차폐 끝점에 RE²의 잔여 접착을 남기지 않는다.
Bare Hamaker 입력은 보존하지만 CMC 내에서 작용하는 유효 분산 인력은 실제로 변한다.
\(q_f\)는 입자 좌표·배향과 무관한 상수이므로 기존 AD 구조에서 스칼라 포텐셜을
미분해 힘·토크를 일관되게 얻는다. 계산된 힘을 사후에 자르거나 법선력만 상쇄하지 않는다.

이 식은 **유효 응집 인력 차폐라는 구성적 가정**이다. CMC–graphite 표면력의 직접
측정값이나 electrostatic 이론으로 유도한 함수가 아니다. 기존 exponential을 크게
증폭하여 접근 장벽을 만드는 대신, 비접착 압축 접촉을 허용하는 최소 모델을 시험한다.
13.2 g/L는 현재 고농도 끝점 가정이며 보편 물성값이 아니다. 독립된 반발 길이·압력·
이온 농도 환산 매개변수는 이 구성식에 사용하지 않는다.

## 기계적 접촉과 마찰

\[
g=H-H_c\ge0,\qquad N\ge0,\qquad Ng=0,
\qquad |\mathbf F_t|\le\mu N.
\]

법선 반력 \(N\)은 기존 접촉 제약 solver가 계산한다. \(\mu=0.1\), 접선 강성과 기존
sliding·rolling 법칙을 유지하며 별도의 법선 스프링, rolling, twisting 모드를 추가하지
않는다. 새 접촉의 adhesive rolling 기준력은 기존과 같이 접촉 생성 시 순접착력에서
계산되므로 완전 차폐로 그 힘이 0이면 함께 0이 된다.

| 상태 | 접촉 하중의 기원 | 시험하려는 거동 |
| --- | --- | --- |
| 낮은 free 농도 | 잔존 접착 예압과 유동 압축 | 접착성 마찰을 통한 yielding |
| 완전 차폐 끝점 | 유동이 만든 압축 | 비접착성 마찰을 포함하는 점성 응력 |

연속상 점도와 lubrication 점도·기하학적 간격은 이번 변경에서 유지한다. Lubrication의
탄소 간격을 \(H-H_c\)로 바꾸어 인위적 발산을 만들지 않는다. 고농도에서 뉴턴 점성이
가능한 구조라고 해서 목표 1–2 Pa·s가 자동으로 따라오지는 않는다. 고농도 비접착 끝점의
접촉 빈도·접촉 응력·총 점도를 먼저 계산하고 저농도 끝점과 비교한다. 중간 농도는 그 뒤
검토하며 \(G'\), \(G''\)는 calibration에 쓰지 않는다.

## 기존 exponential 입력에서 전환

다음 옛 객체는 삭제한다. `enabled=false`라도 구 모델을 조용히 새 모델로 바꾸지 않도록
Python 입력 검증에서 거부한다.

```json
"free_repulsion": {
  "enabled": true,
  "strength": 1.0,
  "decay_length_m": 5e-9,
  "degree_of_substitution": 0.7,
  "repeat_unit_molar_mass_kg_mol": 0.218,
  "osmotic_coefficient": 0.5
}
```

위 객체 대신 `cmc`에 다음 설정을 둔다.

```json
"contact_offset_at_saturation_m": 2e-9,
"free_cohesion": {
  "enabled": true,
  "nonadhesive_concentration_g_L": 13.2
}
```

옛 압력·감쇠 길이에서 새 차폐 끝점 농도로의 자동 환산은 없다. 두 모델을 동시에
작동시키지 않는다. C++의 예전 exponential 계산은 과거 실패 dump/replay 호환에만
남겨 두며, 새 생산 JSON에서는 선택하지 않는다.

## 실행, 기록, 재시작

소스 변경 후 한 번 재빌드한다. 이후 JSON의 농도·모델값 변경은 새 실행에 적용한다.

```bash
sbatch build_slurry_cpu.sbatch
# BUILD COMPLETE / BUILD REUSED 확인 후:
sbatch run_slurry_cpu.sbatch --cases gr_cmc --shear-rates 100
```

배치나 solver 없이 환산값을 확인하려면 다음을 사용한다.

```bash
python3 slurry/drivers/gr_re2/run_graphite.py --config slurry/cases/gr_CMC.json --dry-run
```

`effective_config.json`에는 독립 입력과 bare 부착일을 남긴다.
`manifest.json`의 `derived.cmc`에는 다음 계산 결과를 기록한다.

| 키 | 의미 |
| --- | --- |
| `effective_adsorbed_g_L`, `adsorbed_clamped` | 적용 흡착량과 포화 제한 여부 |
| `theta`, `q` | 흡착률과 추가 접착력 잔존율 |
| `bare_adhesion_work_J_m2`, `background_work_J_m2`, `effective_adhesion_work_J_m2` | bare 기준·배경 일·흡착 감소를 적용한 일 |
| `bare_contact_gap_m`, `contact_gap_m` | \(h_0\), \(H_c\) |
| `contact_model_version` | 접촉·차폐 모델 버전 |
| `free_cohesion.free_concentration_kg_m3` | SI free 농도 |
| `free_cohesion.concentration_fraction` | 제한된 농도비 \(u\) |
| `free_cohesion.cohesion_retention` | \(q_f\) |
| `free_cohesion.active`, `free_cohesion.model_version` | 차폐 활성 여부와 모델 버전 |

활성 CMC 접촉/차폐의 C++ 입력은 `cmc_contact_gap`, `cmc_cohesion_retention`,
`cmc_contact_version`이다. CMC 객체가 없거나 \(x_a=x_f=0\)이면 기존 pure-Gr 힘
법칙을 정확히 회복한다. 흡착만 있는 경우 \(q_f=1\)이지만 접촉면은 이동한다.
`free_cohesion.enabled=false`는 흡착층 접촉 이동까지 끄는 설정이 아니다.

새 피복 접촉·차폐와 물리적으로 다른 옛 CMC checkpoint는 거부한다. 농도·접촉 간격
등 물리를 바꾼 비교는 새 계산으로 수행한다. 동일 물리 설정에서 재시작은 다음과 같다.

```bash
sbatch run_slurry_cpu.sbatch --cases gr_cmc --restart 폴더명
```

별도의 수치 정책 `numerics.pass_max=1`(기본)은 최대 반복에 도달한 유한·평가 가능한
마지막 상태를 허용하고 다음 substep으로 진행한다. 0이면 기존 실패 재시도·종료 정책을
사용한다. 이 값은 횟수 한도가 아닌 0/1 스위치이다. 허용된 횟수는 history의
`max_iteration_passes`(해당 LB step)와 `max_iteration_passes_total`(누적)에 기록하며,
누적값은 checkpoint에 저장된다. NaN·정의되지 않는 상태는 이 정책으로 통과시키지 않는다.
수렴 기준을 충족한 뒤 강제로 최소 비선형 반복 수를 채우는 기능은 현 코드에 없다.

## 실험과 모델 가정의 구분

- [Gwag et al., ACS Nano, DOI 10.1021/acsnano.6c10201](https://doi.org/10.1021/acsnano.6c10201):
  90k CMC의 농도별 유변학과 흡착량이 모델을 시험할 실험 기준이다. 유변 곡선만으로
  이 문서의 특정 차폐 함수나 접촉 기구가 직접 측정된 것은 아니다.
- [Wu et al., Soft Matter 11, 587 (2015), DOI 10.1039/C4SM02380C](https://doi.org/10.1039/C4SM02380C):
  CMC의 graphite 흡착 구조에 관한 근거이다. 희석 조건의 개별 사슬 AFM 높이를
  포화 층의 균일 두께로 간주하지 않으며, 현재의 두 층 합계 2 nm를 그 측정값이라고 쓰지 않는다.

이번 구현의 식, 접촉 간격 증가량, 차폐 끝점 농도는 위 실험 관찰을 설명하기 위한
최소 구성 가정이다. 장시간 다입자 계산이나 oscillatory calibration을 문서 작성 단계에서
수행한 것으로 보고하지 않는다.
