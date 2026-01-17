# sons_pairwise
Reproducible R implementation of the SONS (Shared OTUs and Similarity) metrics
from Schloss & Handelsman (2006).

The code above reproduces the SONS overlap estimators (Û, V̂, abundance-based Jaccard/Sørenson) exactly from the equations in Schloss & Handelsman (2006).https://www.schlosslab.org/assets/pdf/2006_schloss_a.pdf 
For “community structure” comparisons commonly discussed alongside SONS, thetaYC is widely used and is implemented in mothur (formula cited) and reproduced above.
https://mothur.org/wiki/thetayc/?utm_source=chatgpt.com

**NOTE**
The only “judgement call” above is what to do when f_.2 = 0 (no shared doubletons in B) or f_2. = 0 (no shared doubletons in A), because Eq. (5)–(6) then have a division by zero. The paper does not provide an explicit alternative in that case; the code above uses the conservative convention setting the correction term to 0 rather than returning Inf/NA.


## use

Just source the file:

```r
source("R/sons.R")
```
### minimal test

```r
otu_table <- matrix(
  c(10,0,
    5,3,
    1,1,
    0,7,
    2,1),
  nrow = 5, byrow = TRUE,
  dimnames = list(paste0("OTU", 1:5), c("Community_A","Community_B"))
)

sons_from_otu_table(otu_table, "Community_A", "Community_B", f2_correction = "none")
sons_from_otu_table(otu_table, "Community_A", "Community_B", f2_correction = "add1")
```
