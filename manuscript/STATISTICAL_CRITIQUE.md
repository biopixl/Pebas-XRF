# Critical Review: Statistical and Methodological Concerns

## Reviewer Perspective on Hypothesis Testing Analysis

**Date:** 2026-03-05
**Purpose:** Self-critical evaluation before manuscript submission

---

## Major Concerns

### 1. Pseudo-replication and Non-Independence

**Problem:** The analysis treats each XRF measurement (n=1,938) as an independent observation, but measurements at 3mm spacing along continuous cores are **highly autocorrelated**.

| Issue | Impact |
|-------|--------|
| Effective sample size | Dramatically smaller than reported n |
| p-values | Severely inflated (too small) |
| Confidence intervals | Too narrow |
| Effect sizes | May be biased |

**Example:** Reporting "p < 0.001" with n=876 vs n=1,062 is misleading when the true effective sample size might be n~20-50 (number of independent stratigraphic intervals).

**Reviewer would ask:** "How did you account for spatial autocorrelation in your statistical tests?"

**Solution needed:**
1. Calculate autocorrelation length for each proxy
2. Use block bootstrap or mixed-effects models with section as random effect
3. Report effective degrees of freedom
4. Alternatively, aggregate to section-level means before comparison

---

### 2. Inappropriate Unit of Replication

**Problem:** The comparison is framed as "TAM vs SC" but the actual replication units are:
- TAM: 11 sections
- SC: 16 sections

**Correct approach:** Section-level statistics should be the unit of analysis

```
# What we did (wrong):
wilcox.test(all_TAM_measurements, all_SC_measurements)

# What we should do:
wilcox.test(TAM_section_means, SC_section_means)  # n=11 vs n=16
```

**Consequence:** A reviewer would likely reject the current statistical framework outright.

---

### 3. Multiple Comparisons Problem

**Problem:** We tested 10 proxies without correction for multiple comparisons.

| Proxies tested | Family-wise error rate (uncorrected) |
|----------------|-------------------------------------|
| 10 | 1 - (1-0.05)^10 = 40% |

**Solution needed:**
- Apply Bonferroni correction: α = 0.05/10 = 0.005
- Or use Benjamini-Hochberg FDR correction
- Or frame as exploratory rather than confirmatory

---

### 4. Circular Reasoning in Facies Classification

**Problem 1:** The Fe/Mn > 50 threshold for "reducing conditions" is applied, then we report "85% of TAM is reducing" as a finding. This is circular - we defined the threshold, then counted what exceeded it.

**Problem 2:** Thresholds (Ca/Ti > 2, Fe/Mn > 50, Zr/Rb > 1) are arbitrary and not validated against independent facies observations.

**Reviewer would ask:**
- "How were these thresholds determined?"
- "Were they calibrated against thin-section petrography, CT density, or other independent data?"
- "What is the sensitivity of your conclusions to threshold choice?"

**Solution needed:**
1. Justify thresholds from literature or calibrate against independent observations
2. Perform sensitivity analysis with different threshold values
3. Use cluster analysis or mixture models instead of arbitrary cutoffs

---

### 5. Spectral Analysis Deficiencies

**Problems with current approach:**

| Issue | Description |
|-------|-------------|
| No significance testing | Peaks not tested against red noise null |
| No prewhitening | AR(1) trend not removed |
| Irregular sampling | Data gaps not properly handled |
| No confidence intervals | Peak wavelengths reported without uncertainty |
| Depth-to-time conversion | Cannot interpret in orbital terms without age model |

**Reviewer would ask:**
- "How do you distinguish true periodic signals from red noise?"
- "What is the 95% confidence level for spectral peaks?"
- "How do you convert depth-domain cycles to time-domain without radiometric ages?"

**Solution needed:**
1. Use REDFIT or MTM spectral analysis with red noise confidence levels
2. Apply AR(1) prewhitening
3. Report peaks only if they exceed 95% CI against red noise
4. Acknowledge that orbital interpretation is speculative without age control

---

### 6. Sr Anomaly Detection is Tautological

**Problem:** Defining anomalies as ">95th percentile" guarantees that ~5% of data will be flagged as anomalous. This is not detection - it's definition.

**Result reported:** "5.0% of TAM and 5.1% of SC are anomalies" - this is expected by construction.

**Reviewer would ask:** "How is this different from simply selecting the top 5% of values?"

**Solution needed:**
1. Use independent criteria (e.g., Sr/Ca > 0.5 based on marine shell values)
2. Compare to literature values for marine vs freshwater carbonates
3. Require multiple proxy co-occurrence (elevated Sr + Ba + S)
4. Validate against ichnological or faunal evidence

---

### 7. Confounding Variables Not Addressed

**Problem:** TAM and SC differ in multiple ways:

| Factor | TAM | SC | Confounded? |
|--------|-----|-----|-------------|
| Location | Different | Different | Yes |
| Stratigraphic interval | May differ | May differ | Unknown |
| Number of sections | 11 | 16 | Yes |
| Section lengths | Variable | Variable | Yes |
| Preservation | Unknown | Unknown | Possibly |

**Cannot conclude:** "TAM represents a deeper sub-environment" without ruling out:
- Stratigraphic age differences
- Diagenetic overprinting
- Sampling bias
- Preservation differences

**Solution needed:**
1. Constrain stratigraphic ages of both localities
2. Compare time-equivalent intervals only
3. Assess diagenetic alteration indicators
4. Discuss alternative explanations

---

### 8. Correlation ≠ Mechanism

**Problem:** We interpret correlations mechanistically without validation:

| Claim | Evidence | Gap |
|-------|----------|-----|
| "High Fe/Mn = reducing" | Literature | Not calibrated to Pebas specifically |
| "High Zr/Rb = coarser" | Literature | Not validated against grain size measurements |
| "High Ca = carbonate" | Assumed | Could be detrital carbonate |

**Reviewer would ask:** "Did you validate your proxy interpretations with independent measurements?"

**Solution needed:**
1. Compare Fe/Mn to sediment color or organic matter content
2. Compare Zr/Rb to laser grain size or CT density
3. Compare Ca to carbonate content (LOI or XRD)

---

## Minor Concerns

### 9. Data Quality Assumptions

- MSE threshold of 3.0 is arbitrary
- CPS threshold of 20,000 is arbitrary
- Foam exclusion zones are manually defined (subjective)
- No inter-laboratory calibration or standards reported

### 10. Reproducibility

- Raw spectral data and intermediate calculations not archived
- Specific software versions not documented
- Random seeds not set for any stochastic procedures

### 11. Presentation Issues

- Log-scale plots may obscure important patterns
- Color choices may not be accessible to colorblind readers
- Effect sizes should be reported alongside p-values (partially done)

---

## Recommended Revisions

### High Priority (Must Address)

1. **Reframe statistics at section level:**
   - Calculate section-level summary statistics
   - Use section means as unit of analysis (n=11 vs n=16)
   - Apply mixed-effects models if using individual measurements

2. **Address autocorrelation:**
   - Calculate semivariograms or autocorrelation functions
   - Report effective sample sizes
   - Use appropriate spatial statistics

3. **Fix spectral analysis:**
   - Implement REDFIT or MTM with red noise confidence
   - Report only significant peaks
   - Remove orbital interpretation without age control

4. **Validate facies thresholds:**
   - Calibrate against independent observations
   - Perform sensitivity analysis
   - Consider data-driven classification (clustering)

### Medium Priority (Should Address)

5. Apply multiple comparison corrections
6. Validate proxy interpretations with independent data
7. Discuss confounding variables and alternative explanations
8. Revise Sr anomaly detection criteria

### Lower Priority (Could Address)

9. Improve color accessibility
10. Archive intermediate data products
11. Document software versions

---

## Revised Conclusions (More Defensible)

### Original Claim:
> "TAM is significantly more reducing than SC (85.3% vs 49.0%), indicating a deeper, more stratified lacustrine sub-environment."

### Revised Claim:
> "Section-level median Fe/Mn ratios are higher at TAM (median = X, IQR = Y, n = 11 sections) compared to SC (median = X, IQR = Y, n = 16 sections; Wilcoxon W = X, p = Y). This pattern is consistent with, but does not prove, more persistent water column stratification at TAM. Alternative explanations including stratigraphic age differences and diagenetic effects cannot be excluded without additional constraints."

### Original Claim:
> "Spectral analysis reveals dominant periodicities at 42.4 cm and 95.3 cm, potentially corresponding to obliquity and eccentricity cycles."

### Revised Claim:
> "Spectral analysis of the TAM Fe/Mn profile reveals variance concentration at wavelengths of approximately 40-45 cm and 90-100 cm. Without independent age control, these depth-domain cycles cannot be converted to time-domain periodicities. The interpretation as orbital cycles remains speculative pending radiometric or biostratigraphic age constraints."

---

## Summary

The current analysis provides a useful exploratory framework but requires substantial revision before the conclusions can be considered statistically defensible. The core issues are:

1. **Pseudo-replication** - Individual measurements treated as independent
2. **Circular definitions** - Thresholds applied then counted
3. **Unvalidated proxies** - Interpretations not calibrated
4. **Missing significance tests** - Spectral peaks not tested against null

A reviewer familiar with spatial statistics or sediment geochemistry would likely request major revisions addressing these concerns.

---

*Self-critique completed: 2026-03-05*
