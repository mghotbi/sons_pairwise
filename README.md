# sons_pairwise

Reproducible R implementation of the SONS (Shared OTUs and Similarity) metrics from Schloss & Handelsman (2006).

This repository implements:

- Classic incidence-based **Jaccard** and **Sørenson** (Eqs. 3–4)
- Abundance-based shared-OTU overlap estimators **Û** and **V̂** (Eqs. 5–6)
- Abundance-based **Jaccard** and **Sørenson** derived from Û and V̂ (Eqs. 7–8)
- Yue & Clayton **θYC** community structure similarity (Eq. 9)

## References

- Schloss PD & Handelsman J. 2006. *Introducing SONS: Shared OTUs and Similarity.* Appl Environ Microbiol.
- (Paper PDF link): https://www.schlosslab.org/assets/pdf/2006_schloss_a.pdf
- θYC (mothur wiki): https://mothur.org/wiki/thetayc/

---

## Quick start

Clone or download this repository, then source the R script:

```r
source("https://raw.githubusercontent.com/mghotbi/sons_pairwise/Rhizosphere-nitrogen-fate/R/sons.R")

```
### minimal test

```r

otu_table <- matrix(
  c(10, 0,
     5, 3,
     1, 1,
     0, 7,
     2, 1),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(
    paste0("OTU", 1:5),
    c("Community_A", "Community_B")
  )
)

res_none <- sons_from_otu_table(
  otu_table,
  "Community_A",
  "Community_B",
  f2_correction = "none"
)

res_add1 <- sons_from_otu_table(
  otu_table,
  "Community_A",
  "Community_B",
  f2_correction = "add1"
)

res_none
res_add1

```

Bootstrap example

```r

res_bs <- sons_from_otu_table(
  otu_table,
  "Community_A",
  "Community_B",
  f2_correction = "add1",
  bootstrap = TRUE,
  B = 499,
  conf = 0.95,
  seed = 1
)

res_bs$bootstrap$intervals
```
 
 
## Citation

If you use this code in academic work, please cite:

- Schloss PD & Handelsman J (2006) for the SONS method
- Cummins et al. (2025, bioRxiv; DOI: 10.1101/2025.11.24.690264) for an applied use case
