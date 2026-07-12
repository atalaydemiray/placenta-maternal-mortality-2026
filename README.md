# Long-Term Maternal Mortality After Pathologic Placental Separation

Analysis code and results for a study of long-term maternal mortality following placental abruption and placental retention, published in *PLOS ONE* (2026).

## Citation

> Jasani S, Demiray A, Stevenson J, Krawiec C. A comparison of long-term maternal mortality associated with pathologic placental separation: highlighting possible trends and mechanisms. *PLoS ONE*. 2026;21(4):e0338586. https://doi.org/10.1371/journal.pone.0338586

A machine-readable version of this citation is in [`CITATION.cff`](CITATION.cff).

## Authors

- **Sona Jasani** — Department of Obstetrics, Gynecology and Reproductive Sciences, Division of Obstetric Specialties and Midwifery, Yale School of Medicine, New Haven, CT, USA
- **Atalay Demiray** — Department of Health Policy and Management, Yale School of Public Health, New Haven, CT, USA
- **Julia Stevenson** — Yale School of Medicine, New Haven, CT, USA
- **Conrad Krawiec** — Department of Pediatrics, Division of Critical Care Medicine, Penn State Health, Hershey, PA, USA

## Summary

This study compares long-term all-cause maternal mortality among patients with two types of abnormal placental separation, abruption and retention, against patients with normal placental separation. The cohort included 638,911 vaginal deliveries drawn from a federated research network (625,890 normal, 5,435 abruption, 7,586 retention). After adjusting for age, race, and social determinants of health, placental retention carried a 95% higher long-term mortality risk (HR 1.95, 95% CI 1.59 to 2.40) and placental abruption a 59% higher risk (HR 1.59, 95% CI 1.21 to 2.08).

## Repository contents

```
.
├── README.md                               This file
├── CITATION.cff                            Machine-readable citation
├── S7_Analysis_Code.R                      Complete analysis script (S7 Appendix)
├── Figure1_Flow_Diagram.tiff               Cohort flow diagram
│
├── Table1_Baseline_Characteristics.csv     Baseline characteristics by group
├── Table1.docx
├── Table2_Hazard_Ratios.csv                Cox regression hazard ratios
├── Table2.docx
├── Table3_Health_Outcomes.csv              Health outcome associations
├── Table3.docx
│
├── S1_ICD_codes.docx                       ICD-10 codes used to define the cohort
├── S2_EPV_Analysis.csv                     Events-per-variable check
├── S3_PH_Assumption_Tests.csv              Proportional hazards tests
├── S4_Time_Specific_HR.csv                 Piecewise (time-specific) hazard ratios
├── S5_Multiple_Comparison_Corrections.csv  Bonferroni and FDR corrections
└── S6_Sensitivity_Analyses.csv             Sensitivity analyses
```

Every table and figure in this repository matches the published article and its supplementary materials.

## Data availability

The data underlying this study were accessed through a global federated research network of de-identified electronic health record data from participating healthcare organizations. The authors accessed these data under an institutional license governed by a data use agreement. Because the data are proprietary and subject to contractual licensing restrictions, the authors cannot publicly share, redistribute, or deposit them.

**No patient-level data are included in this repository.** Researchers who wish to reproduce the analysis can request equivalent access to the same network through their own institution.

## Reproducing the analysis

### Requirements

- R version 4.5.2 or later
- R packages: `dplyr`, `tidyr`, `survival`, `officer`, `flextable`

### Input data

The script reads a single analysis-ready cohort file from a local `data/` folder that is not distributed with this repository. Once you have obtained equivalent data access, prepare the cohort and save it as:

```
data/analysis_data.RData
```

The file should contain one data frame named `filtered_data`, with one row per patient and the delivery-group, follow-up-time, death-indicator, and covariate columns used by the script. The expected column names are documented inline in `S7_Analysis_Code.R`.

### Run

```r
source("S7_Analysis_Code.R")
```

The script regenerates every table and supplementary result reported in the paper.

## Key findings

1. Both placental abruption (HR 1.59, 95% CI 1.21 to 2.08) and placental retention (HR 1.95, 95% CI 1.59 to 2.40) were associated with significantly higher long-term maternal mortality than normal delivery.
2. Retention showed persistently elevated risk across every time horizon examined (HR 1.75 to 2.02). Abruption showed its highest risk in the first year (HR 1.88), attenuated around the 42-day mark (HR 1.31), and remained significant over the long term (HR 1.49).
3. The proportional hazards assumption held for the primary exposure (p = 0.62).

## How to cite

If you use this code, please cite the article above. GitHub reads `CITATION.cff` to generate a formatted citation from the "Cite this repository" link.

## License

This analysis code is released under the MIT License. See [LICENSE](LICENSE) for details. The license covers the code only. The underlying patient data remain restricted under the data use agreement described above. If you build on this work, please cite the published article.

## Contact

Atalay Demiray, atalay.demiray@yale.edu
