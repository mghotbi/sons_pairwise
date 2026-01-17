# sons_pairwise
SONS overlap estimators

The code above reproduces the SONS overlap estimators (Û, V̂, abundance-based Jaccard/Sørenson) exactly from the equations in Schloss & Handelsman (2006).https://www.schlosslab.org/assets/pdf/2006_schloss_a.pdf 
For “community structure” comparisons commonly discussed alongside SONS, thetaYC is widely used and is implemented in mothur (formula cited) and reproduced above.
https://mothur.org/wiki/thetayc/?utm_source=chatgpt.com

**NOTE**
The only “judgement call” above is what to do when f_.2 = 0 (no shared doubletons in B) or f_2. = 0 (no shared doubletons in A), because Eq. (5)–(6) then have a division by zero. The paper does not provide an explicit alternative in that case; the code above uses the conservative convention setting the correction term to 0 rather than returning Inf/NA.
