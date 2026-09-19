# OpenLB 1.9 배터리 슬러리

OpenLB 1.9r0과 슬러리 확장 모듈을 포함한 통합 소스 저장소입니다.
하나의 실행 파일과 공통 build/run 배치를 사용합니다.

## 지원 모델
- pure_cmc: CMC Cross 모델
- pure_gr: RE² graphite, lubrication, rough contact, Lees–Edwards, checkpoint/restart
- gr_baseline: 기존 graphite Couette 및 checkpoint/restart

## MIT OnDemand 실행
저장소 최상위 폴더에서 실행합니다.
환경: gcc/12.2.0, openmpi/4.1.4, mit_normal.

1. `sbatch build_slurry_cpu.sbatch`
2. 빌드 로그에서 BUILD COMPLETE 또는 BUILD REUSED 확인
3. `sbatch run_slurry_cpu.sbatch --cases pure_gr --shear-rates 100`

CMC와 Gr을 순차 실행하려면:
`sbatch run_slurry_cpu.sbatch --cases pure_cmc,pure_gr --shear-rates 100`

## 개발 위치
- `olb-1.9r0/src/slurry/`: 물리 모듈
- `olb-1.9r0/examples/slurry/`: 통합 실행 진입점
- `slurry/`: 설정, 드라이버, 공통 도구, 피팅 자료

설정은 `slurry/cases/`에서 관리합니다.
결과는 `runs/`에 저장되며 Git 관리 대상에서 제외합니다.
OpenLB 원본의 LICENSE와 저작권 표시는 해당 소스 트리에 보존합니다.

`pure_gr`의 국소 부착 포텐셜과 새 매개변수는
[Graphite 국소 부착 설명](docs/graphite-local-adhesion.md)에 정리되어 있습니다.

`pure_gr` 재시작: `sbatch run_slurry_cpu.sbatch --restart 폴더명`
(`runs/` 아래 폴더). 저장·복구 범위와 설정은
[Graphite 재시작 설명](docs/graphite-restart.md)을 참고하세요.
