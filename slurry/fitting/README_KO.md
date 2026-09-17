# Figure S2 계수의 출처와 재피팅

`cmc250k_digitized_black_lines.csv`는 논문 SI의 CMC 250k 패널에서 추출한
검은 피팅 선의 좌표입니다. LBM 계산 결과가 아니라, 물성 계수를 얻기 위한 입력 데이터입니다.
`cmc250k_digitization_metadata.json`에는 원본 SI의 SHA-256, 축 보정, 농도 배정이 있습니다.
`cmc250k_cross_parameters.json`은 당시 피팅의 계수·설정·Python 환경을 보존한 기록입니다.

LBM 실행이 실제로 읽는 유일한 물성 파일은
`../cmc_cross_openlb/cmc250k_cross_parameters.csv`입니다. 이 폴더의 JSON을
수정해도 실행 입력은 바뀌지 않습니다. 정상적인 복원·LBM 실행에는 재피팅이 필요 없습니다.

계수를 다시 추정하려는 경우에만, Python 3.12 환경에서 아래를 실행합니다.
MIT에서는 Python 환경 준비와 피팅을 할당된 계산 노드에서 진행합니다.

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements-fit.txt
python3 fit_cmc250k.py
```

동봉 CSV를 읽어 `reproduced_cmc250k/`에 새 계수와 그래프를 만듭니다.
고정 난수 seed는 41이며, 각 농도에서 여러 초기값을 사용한 log-viscosity least squares입니다.
원본 피팅 환경은 Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8입니다.
실행용 LBM driver는 별도로 Python 3.6 이상 표준 라이브러리만 사용합니다.

원본 PDF에서 추출 단계부터 반복하려면 본인이 보관한 SI와 PyMuPDF가 추가로 필요합니다.

```bash
python3 fit_cmc250k.py --pdf /path/to/nn6c10201_si_001.pdf --outdir reproduced_from_pdf
```

무한전단 점도 민감도 계산은 `python3 check_eta_inf_sensitivity.py`로 재현하며,
`sensitivity_reproduced/`에 출력됩니다. 생성 결과 폴더는 Git에서 제외됩니다.

계수는 인쇄된 피팅 선을 재현한 값이며 저자의 원 계수나 실험 불확실도를 복원한 값은
아닙니다. 특히 일부 `eta_inf=0`은 최적화 하한에 도달한 결과입니다. 농도 사이의
연속 보간 함수는 아직 정의하지 않았습니다.
