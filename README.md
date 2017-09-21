How to Select K in fastPHASE
================

Description
-----------

The repo contains the following files:

-   `make_test_files.R`: R script to generate fastPHASE input files with artificially missing genotypes

-   `helpers.R`: functions for `make_test.files.R`

Workflow to select the best K?
------------------------------

1.  Prepare a file in the format of fastPHASE input file (\*.inp).

2.  Generate test files with `make_test_files.R`

3.  Impute the missing genotypes in each test files with fastPHASE. Apply different K for each set of test files.

4.  Estimate the imputation quality with `EstimateErrors()` R function and chose K that minimizes the imputation error.

The usage of make\_test\_files.R?
---------------------------------

A. Download `make_test_files.R` and `halpers.R` files.

B. Let `make_test_files.R` to be executable. Run in command line:

    chmod +x make_test_files.R

C. Check that you have the following R packages installed:

    library("optparse")
    library("plyr")

D. Run `make_test_files.R` to generate `n` test files of fastPHASE format (\*.inp), each having `p` proportion of genotypes randomly masked.

Let `chr1.inp` to contain the genotypes of a population from the chromosome one.

The following command will result in 5 test files with 10% percent of genotypes missed.

    ./make_test_files.R chr1.inp 

The default parameters (-n 5 -p 0.1 -o test/test) are applied here. Five test files are named `test.m{1:5}.inp` and saved in `/test` subdirectory that script created.

The following command will produce 3 test files with 5% of genotypes masked. The test files `chr1.m{1:3}.inp` are saved in `/masked` subdirectory.

    ./make_test_files.R -n 3 -p 0.05 -o masked/chr1 chr1.inp 

If there are missing genotypes before running the script, the actual proportion of missing genotypes will be higher, since mask adds missing genotypes to those that exist in original data set.

The usage of fastPHASE tool
---------------------------

According to fastPHASE manual, we can adjust the following option arguments:

-   -T, the number of seeds to launch EM cycles
-   -C, the number of EM cycles
-   -K, the number of haplotype clusters

There are some flags I advice to apply:

-   -H-1, to turn off the phasing, since we need only imputation
-   -n, to tell that we have simplified input, since `make_test_files.R` produces that sort of files
-   -Z, to simplify the format of output files

I recommend to save the results of imputation in a such way that each K has its own folder (`/k10`, `/k15`, etc. )

An example of command for imputation:

    for i in $(seq 1 5); do fastPHASE -T10 -C25 -K10 -H-1 -n -Z -ok10/chr1.m$i masked/chr1.m$i.inp; done

Here we assume that 5 test files (`chr1.m{1:5}.inp`) are in folder `/masked`.

The code will produce 5 `chr1.m{1:5}_genotypes.out` files, where `chr1` is your identifier, `m{1:5}` is added by the above command, and `_genotypes.out` is given by fastPHASE tool. The imputed data sets are saved in `/k10` subdirectory.

The above code helps us to agree input/output from different stages.

Be aware that imputation is time consuming stage. Upon accomplishing, we can estimate the quality of imputation.

The usage of EstimateErrors()
-----------------------------

There is a lot of metrics to estimate the imputation quality of genotypes (Chan et al, 2016). So far I applied the proportion of correctly imputed genotypes. To compute it, use `EstimateErrors()` function. It returns a data frame with three columns: "alleles", "genotypes", and "K". The first two contains the errors, the third one - K (number of clusters).

The error is counted as 1 - accuracy, where accuracy is a proportion of correctly imputed genotypes (alleles). The function returns values for one set of test files.

    # Source functions
    source("helpers.R")

    # Count errors for one set of test files
    errors <- EstimateErrors(origin = "chr1.inp",
                             mask = "masks.RDS", 
                             imputed = imputed,
                             K = 15)

In this code we assume that `chr1.inp` has original genotypes, `masks.RDS` contains the list of all masks (matrices) generated, and `imputed` is a vector with the full paths to `*_genotypes.out` produced by fastPHASE. The order of files in `imputed` is the same as the order of masks applied upstream. All test files were treated with fastPHASE applying K = 15.

References
----------

1.  Scheet P, Stephens M. A Fast and Flexible Statistical Model for Large-Scale Population Genotype Data: Applications to Inferring Missing Genotypes and Haplotypic Phase. American Journal of Human Genetics. 2006;78(4):629-644. [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/16532393)

2.  fastPHASE software. [link](http://scheet.org/software.html)

3.  fastPHASE 1.4 manual. [download](http://scheet.org/code/fastphase_doc_1.4.pdf)

4.  Chan AW, Hamblin MT, Jannink J-L. Evaluating Imputation Algorithms for Low-Depth Genotyping-By-Sequencing (GBS) Data. Feltus FA, ed. PLoS ONE. 2016;11(8):e0160733. [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/27537694)

Author
------

Gennady Khvorykh, a bioinformatician, [inZilico.com](http://inzilico.com)

Interested in contributing to the project? Suggestions, questions, and comments are open! Feel free [to drop me the message](http://www.inzilico.com/contacts/).
