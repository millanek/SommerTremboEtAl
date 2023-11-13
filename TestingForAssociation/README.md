#  Testing for association between exploratory behavior and genetic variants

Both the GLM and pGLS approaches are included here, together with a small test dataset with allele frequencies for the first ~100,000 SNPs of chromosome 1.

To run the code on the example data:
* GLM: `Rscript GLMrun.R NC_031969_AlleleFrequencies_fromProbabilities_n100000.txt.gz`
* pGLS: `Rscript PGLSrun.R NC_031969_AlleleFrequencies_fromProbabilities_n100000.txt.gz`<br>
the pGLS analysis requires the R [caper](https://cran.r-project.org/web/packages/caper/index.html) package: `install.packages("caper")`
 
 

