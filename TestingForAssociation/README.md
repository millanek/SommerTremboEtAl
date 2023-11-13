#  Testing for association between exploratory behavior and genetic variants

Both the GLM and pGLS approaches are included here, together with a small test dataset with allele frequencies for the first ~100,000 SNPs of chromosome 1.

## Running the code on example data

* GLM: `Rscript GLMrun.R NC_031969_AlleleFrequencies_fromProbabilities_n100000.txt.gz`
* pGLS: `Rscript PGLSrun.R NC_031969_AlleleFrequencies_fromProbabilities_n100000.txt.gz`<br>
the pGLS analysis requires the R [caper](https://cran.r-project.org/web/packages/caper/index.html) package: `install.packages("caper")`
 
 ## Preparing the full dataset

The full dataset cannot be stored on GitHub because of space limitations but VCF files can be obtained from the associated DataDryad [repository](https://datadryad.org/stash/share/kDGaxgviYannD9DNVDUQP2pOYwU-I1yGBj_1CWzlI90). To get the allele frequency estimates from these VCF files, you should use the [evo](https://github.com/millanek/evo/tree/e8ed8b1cd19aed36be568fbee05a2692ebcfa6e4) software as follows:

`evo alleleFreq -n TanganyikaGWAS_wOreTan_NC_031969_fromProbabilities --use-genotype-probabilities TanganyikaGWAS_wOreTan_NC_031969_SelectedVariants_variableSnpsBiallelicNoSingletons.vcf.gz medians_for_GWAS_withoutfilter_speciesCodesAndGenomeIDs_wOretan.txt`

