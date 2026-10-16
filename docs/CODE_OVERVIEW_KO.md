# 공개 코드 한눈에 보기 (code_release)

작성 2026-10-15. 세부 내용(스크립트별 입력·출력, 외부화 규칙, 테스트)은 [CODE_STRUCTURE_KO.md](CODE_STRUCTURE_KO.md)에 있다.

## 1. 세 층 구조

| 층 | 위치 | 하는 일 |
|---|---|---|
| 공통 라이브러리 | `src/fragpolar/` | config 로더, CpG call 정제와 marker 충분통계(`extract.py`), 한 번 읽기 엔진(`engine.py`), 12지표 참조 구현(`feature_definitions.py`), BAM 추출기 사본(`vendor/`), 테스트용 순수 수식(`metrics.py`) |
| 분석 파이프라인 | `pipelines/01–09` | 원고 결과를 만드는 9개 분석 단계 |
| 그림 | `figures/main`, `figures/supp` | 본문 Fig. 1–6, 보충 S1–S18(공통 스타일 `style_v8.py`) |

경로, 입력 파일, 열 이름, 표본 목록은 모두 `config/config.example.yaml`(자리표시자; paths/inputs/columns/sample_lists/params)에서 읽는다. 간염은 config 목록의 표본만 쓴다. 표본 ID, 기관명, 서버 경로, 사용자 ID, 원자료와 결과 파일은 포함하지 않았다.

## 2. 파이프라인과 원고의 대응

| 폴더 | 내용 | 원고 |
|---|---|---|
| 01_group_effects | marker 행렬로 59개 층의 HCC−건강인 효과, BFV 분해, K·길이 보정, 방향 margin (Test) | 결과 2, 결과 3(margin), Fig. 2, S3, S5 |
| 02_train_rebuild | Train 777개를 raw에서 같은 규칙으로 다시 집계하고 01을 Train에 적용 | 결과 2·3의 Train 수치 |
| 03_fpi_k | K 보정 FPI와 근사 오차 | 결과 2, S4 |
| 04_generators | 합성 생성기(G1–G3) 반응, 관측 변화와의 cosine, 엄격 양극화 판정 | 결과 1·3, Fig. 1d, 3f, S6 |
| 05_direct_fragments | 매칭 쌍의 fragment 재집계, sOUF 제거와 같은 수 대조, read 분할 | 결과 3, Fig. 3a,c–e |
| 06_tumour_hepatitis_age | 종양분율 용량–반응, 간염 세 군 비교, 연령 분석 | 결과 4–6, Fig. 4c, 5, 6, S10, S12 |
| 07_pbmc | PBMC BAM 추출과 혈장 대비 | 결과 4, S9 |
| 08_component_origin | 정상 세포 atlas(PAT)와 조직 혼합 기대값, 효과 벡터 cosine | 결과 4·5, Fig. 4d, S10 |
| 09_mixture_and_length | 혼합 기대 BFV, fragment 길이 지표와 보정 | 결과 2·3, Fig. 3b, S7, S8 |

실행 순서: 01 → 02 → 03, 04, 05, 06, 07, 08 → 09 → figures (자세한 순서와 명령은 README).

## 3. 검증

- 배포 전 문법 점검(2026-10-16): Python 파일 116개, Python 3.8.5에서 실패 0.
- 단위 테스트 43개 통과: 분산 분해(WFV + BFV = μ(1−μ)), FPI·FPI_K, sOUF, 제거 대조, read 분할 hash, 혼합 기대 BFV, 추출–엔진 일치, config 매핑.
- 실제 자료 대조: Fig. 1b, 4b–d, 5a, 5b, 6a–c, S10a–c 원표와 매칭 쌍 파일이 원본과 같다. Fig. 5c와 Supp. Fig. 10은 공개 코드로 다시 생성했다(CODE_STRUCTURE_KO.md 9절).
- 식별자 점검: 서버 경로, 사용자 ID, 표본 ID, 기관명, 내부 폴더·파일·열 이름 grep 결과 배포 대상에서 발견되지 않음. 내부용 경로·이름 대응표는 배포본에서 제거했다.

## 4. 배포 상태와 후속 확인

1. 내부용 경로·이름 대응표는 배포본에서 제거했다. 데이터와 실제 설정 파일은 Git에서 제외한다.
2. 원고 결과 5의 간염 μ CI 하한(1.1 → 1.0)과 U 기준 층 sOUF 증가 CI 층 수(7 → 5), Fig. 6d 간염 증가 층 수 변화 반영 여부(CODE_STRUCTURE_KO.md 9절).
3. 코드에는 MIT 라이선스를 적용했다. 저장소 URL은 README에 기록하며 논문 DOI는 출판 후 보완한다.
4. Methods의 Python 버전 표기.
