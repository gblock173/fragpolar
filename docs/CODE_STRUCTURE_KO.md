# 코드 구조와 실행 안내 (code_release)

- 작성 2026-10-15, 개정 2026-10-15(사용자 결정 반영: 간염 제외 분기 제거, 내부 이름 일반화, Fig. 5c·Supp. Fig. 10 공개 코드 재생성).
- 원본 코드는 이동·수정하지 않았다. 모든 파일은 복사본이며, 수정은 경로·식별자·이름의 외부화와 간염 표본 선택(목록만 사용)에 한정했다. 재구성 스크립트 2개는 7절.
- 배포본에서는 내부 경로·이름 대응표와 내부 실행 기록을 제외했다(2026-10-16).

## 1. 폴더 구조

```
code_release/
├── README.md, environment.yml, environment_py38.yml, requirements.txt, pyproject.toml, .gitignore
├── config/config.example.yaml    paths / inputs / columns / sample_lists / params (자리표시자만)
├── src/fragpolar/                config, extract, engine, feature_definitions, metrics, vendor/
├── pipelines/
│   ├── 01_group_effects/         ge1_aggregate, ge2_analysis, ge3_strict_subset_check, ge4_draft_figures
│   ├── 02_train_rebuild/         tr0_make_ids, tr1_extract, tr2_aggregate, tr3_analysis, tr4_strict_subset_check,
│   │                             tr5_compare_legacy_raw, tr6_parity, tr7_compare_strict_recompute
│   ├── 03_fpi_k/                 fk1_aggregate, fk2_raw_invK, fk3_khist, fk4_analysis, fk5_jensen_raw, fk6_draft_figures
│   ├── 04_generators/            onepass_engine(진입점), gen_lib, gen_panel_job, gen_observed, gen_observed_rebuilt,
│   │                             gen_analysis(주), gen_analysis_initial, gen_seed_consistency, gen_draft_figures
│   ├── 05_direct_fragments/      df0_initial_pairs, df0b_pairs_from_generators(재구성), df1_extract, df2_sample, df3_stratum,
│   │                             df0c_groups_full (--compare full), df4_stats, df5_draft_figures, df6_raster
│   ├── 06_tumour_hepatitis_age/  tha_common, tha0_strict_extract, tha1_extract, tha2_tumour_fraction, tha3_hepatitis,
│   │                             tha3b_hepatitis_contrasts(재구성), tha4_ageing, tha5_draft_figures
│   ├── 07_pbmc/                  pb0_manifest, pb1_extract, pb2_plasma_strict, pb3_qc_strict, pb4_summarize, pb5_analyze,
│   │                             pb6_draft_figure, pb7_parity_check
│   ├── 08_component_origin/      co0_pat_tasks, co1_marker_cpg_ranges, co2_pat_extract, co3_tissue_suff, co4_marker_obs,
│   │                             co5_analysis, co6_draft_figures
│   └── 09_mixture_and_length/    ml1_mixture_model, ml2_mixture_draft_figures, ml3_length_manifest, ml4_length,
│                                 ml5_length_analysis, ml6_length_draft_figures, ml7_length_tumour_fraction
├── figures/main/                 style_v8, common, fig456_common, prep_fig1, prep_fig1b_cards, prep_fig2, prep_fig3,
│                                 prep_fig3_raster, fig456_prep, prep_fig4a_pairs, prep_fig4c_stats, prep_fig6d, fig5a_auc_make, fig1 … fig6
├── figures/supp/                 sc(공통), s1 … s18
├── scheduler/job_template.sh
├── tests/                        conftest, test_metrics, test_extract_engine, test_config
└── docs/CODE_STRUCTURE_KO.md
```

출력: `<work_dir>/<패키지>/{results,data,work,figures,logs}`, 패키지 = group_effects, train_rebuild, fpi_k, generators,
direct_fragments, tumour_hepatitis_age, pbmc, component_origin, mixture_and_length. 그림은 `<figure_dir>`.

## 2. 설정과 이름 규칙

- `paths`: 입력 루트와 출력 루트. `inputs`: 모든 입력 파일·폴더를 일반 이름으로 지정(data_dir 기준 상대경로 또는 `{reference_dir}`, `{extraction_dir}` 확장). `columns`: 일반 열 이름 → 사용자 파일의 실제 열 이름(`healthy_direction`, `hof_feature`, `legacy_source_label`). 원래 자료의 파일·열 이름을 바꾸지 않고 config만으로 연결한다.
- `src/fragpolar/config.py`: `CFG.inp(name)`, `CFG.col(name)`, `CFG.featname(f)`(HOF 특징의 파일 stem/라벨), `CFG.unmap_features(df)`(이전 표의 HOF 라벨을 `hof`로), `CFG.path`, `CFG.ids`, `CFG.param`, `CFG.work`, `CFG.cpg_calls`, `CFG.pipeline_dir`.
- 출력 파일 접두사: 군 비교 `ge_`, FPI_K `fpik_`, 생성기 `gen_`, 직접 대조 `df_`, 종양분율 `tf_`, 간염 `hep_`, 연령 `age_`, 혼합 기대 `mix_`, 길이 `len_`. 학습군 재구성 반영 변형은 `_rebuilt`, 생성기 쌍 기준 결과 폴더는 `results_genpairs`.
- HOF 특징 이름은 코드에서 `hof`. 이전 분석 표와의 parity 비교는 `prior` 라벨(`priorset`, `priorown`, `*_vs_prior`).
- 계산식, seed, 반복 수는 바꾸지 않았다. read-split 해시의 salt·hash key 문자열은 계산 상수이므로 그대로 두었다.

## 3. 간염 표본

간염은 `sample_lists.hepatitis_included`의 표본만 쓴다. 이 목록 밖의 간염 표본을 계산·출력하는 분기는 공개본에서 모두 지웠다:
`tha1_extract`(표본×층 값 생성 단계에서 목록으로 제한), `tha3_hepatitis`(세 군 대비), `tha3b_hepatitis_contrasts`, `co4_marker_obs`, `co5_analysis`,
`ml3_length_manifest`, `ml5_length_analysis`(참조 간염군 비교 제거), `fig456_prep`, `fig5a_auc_make`, `s8`. config의 참조 간염군 키는 삭제했다.

## 4. 공통 라이브러리

| 모듈 | 역할 |
|---|---|
| config.py | 위 2절 |
| extract.py | 업스트림 추출 모듈의 aggregate/raw_extract/check_metrics 그대로, synthetic_checks |
| engine.py | one-pass 엔진(marker_table, clean_sort, aggregate, gold_summary, gold, trim_keep, main) |
| feature_definitions.py | 12지표 참조 구현, 예시 집합, 값 표기 |
| metrics.py | 파이프라인 인라인 수식의 테스트용 재진술(파이프라인은 import하지 않음; 수식 문자열 동기화 테스트) |
| vendor/ | extract_cpg_calls.py, methyl_features.py (BAM → CpG call) |

## 5. 단계별 역할과 원고 대응

| 단계 | 핵심 스크립트 | 주요 출력 | 원고 |
|---|---|---|---|
| 01 group_effects | ge1_aggregate [batch] → ge2_analysis | ge_strata_effects_all_metrics, ge_deltaB_decomposition, ge_regression_*, ge_binstandardized_*, ge_margin_* | Fig. 2b–d(Test), Supp. Fig. 5 |
| 02 train_rebuild | tr1_extract [batch], tr2_aggregate [batch], tr3_analysis; parity tr5–tr7 | 같은 ge_* (Train), legacy_vs_raw_*, parity_* | Fig. 2b–d(Train), Suppl. Methods 재집계 |
| 03 fpi_k | fk1–fk5 | fpik_* | Fig. 2d, Supp. Fig. 4 |
| 04 generators | onepass_engine, gen_panel_job [batch], gen_observed(_rebuilt), gen_analysis | gen_response_TrainSeed20_rebuilt, gen_cosine_by_stratum_rebuilt, gen_polarization_*_rebuilt, gen_matched_pairs_rebuilt | Fig. 1d, 2a, 3f, Supp. Fig. 3, 6 |
| 05 direct_fragments | df0b → df1, df2 [batch] → df3 → df4 | results_genpairs/df_native_reduction_by_stratum, df_matched_dose_by_stratum, df_readsplit_* | Fig. 3c–e |
| 06 tumour_hepatitis_age | tha0 [batch], tha1 → tha2 → tha3 → tha3b → tha4 | sample_strata_long, tf_dose_response, tf_f_samples, hep_contrasts, hep_contrasts_included, age_* | Fig. 2e, 4c, 5b–c, 6, Supp. Fig. 10, 12 |
| 07 pbmc | pb1 [batch, py3.8], pb2 [batch], pb3 → pb5 | group_means, effects, effect_correlation, paired_diffs | Supp. Fig. 9 |
| 08 component_origin | co2, co3 [batch], co4 → co5 | similarity_*, observed_vs_observed_* | Fig. 4d, Supp. Fig. 10 |
| 09 mixture_and_length | ml1 [batch], ml4 [batch], ml5 | mix_strata_point(_Test), mix_summary_boot, len_* | Fig. 3b, Supp. Fig. 7, 8 |
| figures | prep_* / fig456_prep / fig5a_auc_make → fig1–6, s1–s18 | source_data, Figure_n, Supplementary_Figure_Sn | 전 그림 |

사용자가 준비하는 표본 목록 파일(공개하지 않음): sample_lists 4종, `pbmc/data/plasma_missing_strict.tsv`, `component_origin/data/tissue_tasks.tsv`, tha0 대상 목록, fk2/fk5의 표본 ID 인자. 실행 순서는 README의 "Run order".

## 6. 포함/제외

- 포함: 9개 분석 패키지와 그림 코드 전부(아래 제외 제외), 상위 코드에서 필요한 함수(engine, extract, feature_definitions, vendor).
- 제외: 원고 조립 스크립트, 모든 .pbs(→ scheduler/job_template.sh), 특정 표본 ID에 고정된 엔진·패널 parity 스크립트 2개, 미사용 parity 참조 import, SUMMARY·결과·데이터.

## 7. 재구성 스크립트

- `df0b_pairs_from_generators.py`: 생성기 쌍 → 직접 대조 입력. 내부 원 결과 3개 파일과 완전 동일.
- `tha3b_hepatitis_contrasts.py`: 간염 세 군 대비(그림 5c, 6d, Supp. Fig. 10 입력). 점추정은 원 표와 동일, CI 끝값은 원 생성 코드가 없어 난수 순서가 달라 최대 0.0025 차이. 사용자 결정에 따라 이 공개 코드 결과로 Fig. 5·Supp. Fig. 10을 다시 생성했다(9절).

## 8. 검증

- 배포 전 확인(2026-10-16): Python 3.8.5 AST 구문 검사 116개 파일 통과, 합성 입력 pytest 43개 통과. 아래는 코드 정리 과정의 검증 기록이다.

- py_compile: 108개 .py, Python 3.11.14와 3.8.5 모두 실패 0. pyflakes 미정의 이름 0.
- pytest: 43 passed (합성 입력; 분산 분해, FPI/FPI_K, sOUF, 제거 대조, read-split hash, 혼합 기대 BFV, 수식 동기화 11건, extract/engine 일치, gold_summary vs 참조 구현, config의 inputs·columns·feature 매핑).
- 실데이터 스모크(이름 변경 후, 내부 config로 원래 파일·열 이름에 연결): Fig1b, Fig5a, Fig4b, Fig4c, Fig4d, Fig5b, Fig6a–c 원본과 동일; 생성기 쌍 직접 대조 입력 3개 동일; 간염 대비표 점추정 동일; S10a–c 원본과 완전 동일.

## 9. Fig. 5c·Supp. Fig. 10 재생성과 원고 수치

공개 코드(tha3b_hepatitis_contrasts → fig456_prep → fig5, s10)를 원 데이터로 실행해 Figure_5, Supplementary_Figure_S10, Fig5c.tsv, Fig5b_contrasts_AllMarkers.tsv, S10a–c.tsv를 다시 만들었다(이전 파일은 figures/_pre_regen/에 백업). 원고 결과 5 대조:

| 항목 | 원고 | 재생성 | 판정 |
|---|---|---|---|
| sOUF HCC−H / 간염−H / HCC−간염 (전체 marker) | +0.75 / −0.36 / +1.11%p, 모두 CI가 0 제외 | +0.75 (0.59–0.91) / −0.36 (−0.49 – −0.22) / +1.11 (0.91–1.32) | 일치 |
| strict-mixed | −3.0 / −2.7 / −0.3 (CI −0.8–0.2) | −3.04 / −2.73 / −0.31 (−0.83–0.18) | 일치 |
| 간염 μ | +1.4 (1.1–1.8) | +1.40 (1.01–1.76) → 1.0–1.8 | **하한 1.1 → 1.0** |
| 간염 sOUF·BFV 감소 층(CI) | 31·34 | 31·34 | 일치 |
| M 기준 층 sOUF 감소(CI) | 13/15 | 13/15 | 일치 |
| U 기준 층 sOUF 증가 / CI | 12/14, CI 7 | 12/14, CI **5** | **7 → 5** |
| strict-mixed 감소 층 간염·HCC | 45·47 | 45·47 | 일치 |

그림 5의 표시 값(CI 제외 층 수 31·34 등)은 변하지 않았다. Fig. 6d(세 원천 방향 요약)도 같은 표를 쓰며 간염 증가 층 수(CI)가 sOUF 10→8, BFV 6→5, strict-mixed 7→6으로 바뀐다. 사용자 지시 범위(Fig. 5, S10) 밖이라 Figure 6과 Fig6d.tsv는 다시 만들지 않았다(_pre_regen에 백업만). 원고는 고치지 않았다.

## 10. 남은 문제와 결정 사항

1. 원고 수정 필요 여부: 간염 μ CI 하한(1.1→1.0), U 기준 sOUF 증가 CI 층 수(7→5), Fig. 6d 간염 증가 층 수(위 표).
2. Methods의 Python 버전 표기(3.8.5)와 실제 주 분석 환경(3.11.14) 불일치.
3. 코드의 라이선스: MIT (`LICENSE` 참조). 외부 소프트웨어와 데이터에는 해당 제공자의 조건이 적용된다.
4. source 열로 묶는 점검 출력(ge3/tr4_strict_subset_check, ml5의 출처별 중앙값)은 공개하지 않는다.
5. 라이선스는 MIT. 저장소 URL은 README에 기록하며, 논문 DOI와 인용 문구는 출판 후 보완한다.

## 11. 식별자 점검

서버 경로, 사용자 ID, 표본 ID 형식, 기관·출처 이름, conda 경로, 내부 저장소·패키지·폴더·파일·열 이름 패턴으로 배포 대상 파일을 점검했다. 코드·config·README·tests·배포 문서에서 해당 식별자는 발견되지 않았다. 일반 PBS 템플릿의 문법 단어만 남는다.

## 12. 보충 그림 S14–S18, Fig. 6d 재설계, CpG island 표기 (2026-10-15)

- 새 스크립트: `figures/supp/s14.py`–`s18.py`, `figures/main/prep_fig6d.py`, `pipelines/09_mixture_and_length/ml7_length_tumour_fraction.py`(종양분율 46명에서 길이 지표·sOUF의 Spearman ρ와 HCC 표본 bootstrap CI 5,000회, seed 20261015; 출력 `len_tumour_fraction_spearman.tsv`, `len_tumour_fraction_samples.tsv`).
- 수정: `style_v8.py`(`ann_label`, `stratum_label`: 'Island' → 'CpG island', 'nonIsland' → 'non-CpG island' 표시; `strata_x` 층 배치), `sc.py`·`fig2.py`(heatmap 열 라벨), `s12.py`(d 라벨, M/mid 비교의 Train 행 누락 수정), `fig456_prep.py`(Fig. 6b/c 라벨, 옛 Fig. 6d 계산 삭제), `fig6.py`(d 재설계).
- `s14.py`는 생성기 표의 HOF 라벨을 `CFG.unmap_features`로 일반화한다.
- 검증: 원 결과를 새 이름으로 연결한 작업 폴더에서 공개본을 실행해 S12abc, S12d, S14–S18, Fig6d source data와 len_tumour_fraction_* 가 내부본과 바이트 단위로 같았다. pytest 43 passed.

## 13. 05_direct_fragments: 비교 방식 옵션 (matched | full)
- `df4_stats.py --compare matched`(기본): 매칭 쌍 짝 차이, 쌍 재추출 bootstrap.
- `df4_stats.py --compare full`: HCC 군 평균 − 건강인 군 평균(Train 104 대 673, Test 73 대 336), 군별 재추출 bootstrap. 입력은 `df0c_groups_full.py`가 만드는 `results/df_groups.tsv`(환경변수 `GROUPS`로 변경 가능).
- full 순서: `df0c_groups_full` → (추출 안 된 표본만) `df1_extract`, `df2_sample` (`SAMPLES=results/df_extra_samples.tsv`, 배치 array) → `df3_stratum` (`SAMPLES=results/df_union_samples.tsv`, `PKL=work/stratum_values_full.pkl`) → `df4_stats --compare full` (`PKL=...`, `OUTD=<results_full>`).
- full 결과표의 `n_pairs` 열은 HCC 쪽 유효 표본 수이다. 설정(seed, 반복 수, 대조 정의)은 matched와 같다.
