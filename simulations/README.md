#  Using simulations to control for incomplete lineage sorting and gene flow

## Requirements:

* [msprime v1.0.2](https://pypi.org/project/msprime/1.0.2/)<br> The simulation script is written for this version. It may or may not work with newer versions. In newer versions there are certainly easier ways to perform some of the tasks (e.g., use msprime to directly write out a VCF file).

* A number of other python modules listed in the `simulate.py` script.    

## Goal

* Simulate genetic data from the Lake Tanganyika radiation, taking into account incomplete lineage sorting and gene flow. 

## Details

The Tanganyika phylogeny and gene flow estimates come from: <br> 
Ronco, F. et al. (2021) Drivers and dynamics of a massive adaptive radiation in cichlid fishes. Nature 589, 76–81. doi: [https://doi.org/10.1038/s41586-020-2930-4](https://doi.org/10.1038/s41586-020-2930-4)

As in Ronco et al. (2021), we use 20 randomly subsampled gene-flow matrices (see that paper's supplement for details). 

Within each locus, recombination occurs at the rate of 2.2e-8 per bp per generation. The mutation rate is 3.5e-9 per bp per generation. Loci are unlinked.

In the example below, as in the paper, we set Ne to 80,000 throughout and we allow gene flow between any two species to persist for up to 1 million generations after split, but not longer. 

## Usage

Executing `./runTwoSampleWithRecomb.sh 0 80000 1000000` will run 200 simulations, each generating a 10kb segment of the genome.     

To get the entire set of 20,000 unlinked loci as for the paper, run `for i in {0..99}; do ./runTwoSampleWithRecomb.sh ${i} 80000 1000000; done`. You probably want to parallelise this on a compute cluster.  

This simulates the entire Lake Tanganyika radiation. To subset the data to the species of interest for this paper, you can use for example [bcftools](https://samtools.github.io/bcftools/bcftools.html) like this: <br>
`bcftools view -c 1:minor -O z -o simulate.0_1_twoSampleWithRecomb_Ne_80000_p_1000000_SpeciesSubset.vcf.gz --samples-file SimulationSubsetSpecies.txt simulate.0_1_twoSampleWithRecomb_Ne_80000_p_1000000.vcf`


