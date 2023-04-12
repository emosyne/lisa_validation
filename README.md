# Genetic effects of tissue-specific enhancers in schizophrenia and hypertrophic cardiomyopathy
## Doctoral thesis by Emanuele Felice Osimo for Imperial College London
The link to the thesis will be included once publicly deposited.
Contact me at eosimo at ic.ac.uk

### Chapter Chapter 4 - EP-WAS external validation on a PGC sample.



This repository contains the code for work contained in Chapter 4. This is a Nextflow pipeline, in which I validated the EP-WAS, utilising the same methodology as in the previous section -- this time on an external PGC cohort. 
Here I calculated partitioned PRSs for each partition -- the enhancer-based partition based on UKBB dominant and recessive EP-WASes -- and the residual and original GWAS partitions, both based on the xs234 leave-one-out (LOO) PGC original additive GWAS. Then, using the xs234 PGC cohort as the target population, I calculated the coefficients of determination (CoDs) for each partition, and for the multivariate models. 

Run this pipeline with:

```
nextflow run main.nf -profile servername
```


