# Long-term Maternal Mortality in Pathologic Placental Separation

## Study Information

**Title:** A comparison of long-term maternal mortality associated with pathologic placental separation: highlighting possible trends and mechanisms

**Authors:** Jasani S, Demiray A, Stevenson J, Krawiec C

**Corresponding Author:** Atalay Demiray (atalay.demiray@yale.edu)

## Repository Structure

```
placenta/
├── README.md                           # This file
├── S7_Analysis_Code.R                  # Complete reproducible analysis code
├── Figure1_Flow_Diagram.tiff           # CONSORT flow diagram (PLOS format)
│
├── Table1_Baseline_Characteristics.csv # Baseline demographics
├── Table1.docx                         # Formatted Table 1
├── Table2_Hazard_Ratios.csv            # Cox regression results
├── Table2.docx                         # Formatted Table 2
├── Table3_Health_Outcomes.csv          # Health outcome associations
├── Table3.docx                         # Formatted Table 3
│
├── S1_ICD_codes.docx                   # ICD-10 codes for cohort definitions
├── S2_EPV_Analysis.csv                 # Events per variable analysis
├── S3_PH_Assumption_Tests.csv          # Proportional hazards assumption tests
├── S4_Time_Specific_HR.csv             # Piecewise Cox regression results
├── S5_Multiple_Comparison_Corrections.csv  # Bonferroni/FDR corrections
├── S6_Sensitivity_Analyses.csv         # Sensitivity analysis results
│
└── data (DO NOT SHARE)/                # Data files (not for redistribution)
    ├── filtered_patient_data.RData     # Analysis-ready cohort
    └── vaginal_patient_data.csv        # Raw data from TriNetX
```

## Data Availability

The data used in this study were obtained from TriNetX, a federated health research network. Per the TriNetX data use agreement, raw data cannot be redistributed. Researchers interested in replicating this analysis may request access to TriNetX through their institution.

## Reproducing the Analysis

### Requirements
- R version 4.5.2 or later
- Required R packages: dplyr, tidyr, survival, officer, flextable

### Instructions
1. Obtain data access through TriNetX
2. Place `filtered_patient_data.RData` in the working directory or `data/` folder
3. Run `S7_Analysis_Code.R`

```r
source("S7_Analysis_Code.R")
```

The script will generate all tables and supplementary materials.

## Key Findings

1. **Primary Result:** Both placental abruption (HR 1.59, 95% CI 1.21-2.08) and placental retention (HR 1.95, 95% CI 1.59-2.40) are associated with significantly elevated long-term maternal mortality compared to normal delivery.

2. **Temporal Patterns:**
   - Retention: Persistently elevated risk across all time horizons (HR 1.75-2.02)
   - Abruption: Highest risk short-term (HR 1.88), attenuated at 42 days (HR 1.31), significant long-term (HR 1.49)

3. **Methodological Note:** The proportional hazards assumption is satisfied for the primary exposure variable (p=0.62).

## Contact

For questions about this analysis, please contact:
- Atalay Demiray: atalay.demiray@yale.edu

## License

This code is provided for academic and research purposes.
