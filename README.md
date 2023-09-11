A suite for automatically detecting the best MSA subsampling parameters for the purpose of predicting changes in the conformational landscapes of proteins in response to mutations with subsampled AlphaFold2

# Preface:

This is a work in progress. Below is an example of one possible workflow using this method in which no prior information about a protein system is necessary but its sequence and the sequences of its variants of interest.

https://github.com/GMdSilva/af2_subsampling_opt/blob/main/notebook_demo.ipynb

Steps in bold have been implemented and are part of this codebase, other steps are under development

**1**. Build Deep MSA with jackhmmer for target protein (implemented at https://github.com/GMdSilva/rel_state_pop_af2_raw_data/blob/main/rel_state_populations_af2_msa_generation.ipynb, to be added to main workflow)
   (Step 1 can be skipped if the user prefers to submit their own MSAs)
   
2. Determine the initial range of max_seq:extra_seq MSA subsampling parameter values based on MSA depth 

**3**. Systematically run AlphaFold2 through the colabfold_batch wrapper with chosen range of max_seq:extra_seq parameter values to make predictions for the target protein, 160 seeds per prediction

**4.** Analyze output ensembles using MDanalysis and identify patterns correlated to the accurate prediction of relative state populations

   **i**. Calculate the Root Mean Square Fluctuation (RMSF) of C-alpha atomic positions when aligned to the #1 ranked prediction (by average pLDDT) or to user-selected reference **for each prediction ensemble**
   
   **ii.** Using the find_peaks method from scikit-learn, identify contiguous residue ranges of significant variation across each prediction ensemble
   
   **iii.** For the N selected regions, calculate the Root Mean Square Deviations (RMSD) of C-alpha atomic positions within the selected range for each structure within a prediction ensemble with regards to the #1 ranked prediction (by average pLDDT) or to user-selected reference
   
   **iv.** For each RMSD distribution, apply a kernel density estimation function and use find_peaks on the processed data to get the modes (if any) of the distribution. Identified modes (peaks) might correlate to discrete conformational states
   
   **v.** For RMSD distributions with at least two modes, calculate peak height (correlates to state population), peak width (correlates to state flexibility), distance between peaks (correlates to magnitude of conformational change), among others
   
   **vi.** Score the statistics calculated above using a function whose output grows based on the proximity of the statistics to those observed for the benchmarking systems tested in https://www.biorxiv.org/content/10.1101/2022.10.17.512570v1.full
   
   **vii.** Select subsampling parameters that led to the largest score (i.e., parameters whose AF2 prediction results had observables of similar distributions to those observed in the benchmarking tests)
   
5. Re-run AlphaFold2 with the selected parameter set, now with 160*N seeds for increased statistical power

6. Build deep MSAs with jackhmmer for variants of target protein
   
**7**. For each variant, run AlphaFold2 with the previously selected parameter set using the corresponding MSA built from 6
   
**8**. For each prediction ensemble, repeat the steps in 4 pertaining to peak detection, and contrast peak statistics with wild-type results
    
9. Group results based on their effects (if any) to the relative populations of predicted conformations

# Installation Guide:

(work in progress)

We strongly recommend the use of our Google Collab notebook for generating MSAs. Installation for its usage can be found on the notebook itself.

# Demo:

(work in progress)

# Instructions for Reproduction:

(work in progress) 
