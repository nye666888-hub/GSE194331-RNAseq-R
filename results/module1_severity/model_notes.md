# Severity model notes
- AP-only samples: MAP 57, MSAP 20, SAP 10. Input values are log2(existing QN expression + 1), explicitly exploratory due unknown source units.
- Primary trend design: ~ severity_score with MAP=1, MSAP=2, SAP=3; limma robust empirical Bayes moderated t on per-step slope.
- This tests association with an ordinal score under a constant per-step log-expression slope; it does not establish a linear biological dose response.
- Sensitivity model: ~ severity categorical with a global 2-df moderated F test. Comparisons B, C, D provide pairwise AP category contrasts.
- All 87 AP samples retained. SAP has only 10 samples; effect sizes/95% intervals accompany the trend and pairwise tables; nonsignificance is not absence.
- No clinical covariates or batches were available.
