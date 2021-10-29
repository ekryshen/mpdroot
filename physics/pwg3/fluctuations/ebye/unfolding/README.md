## Intro

This code was developed to make use of the Unfolding technique [1] that de-convolute distributions "spoiled" by the detector by means of the detector response matrix proper inversion. The elements of this matrix are the probabilities to measure, for example, N_rec charged particles in the event, while in reality in this event N_sim was produced. Also this method allows one to account for lost (by the event selection procedure, for instance) and fake (by trigger malfunctioning, for example) events or tracks. There are a lot of other methods to correct the data, for example, bin-by-bin weighting of distributions, corrections at the level of moments, etc. But Unfolding is a powerful tool since it accounts for event migrations (and therefore can be applied when the response matrix is non-diagonal, which is impossible for bin-by-bin weighting method). Also, Unfolding procedure is easy to use since it has already been implemented in Root (see RooUnfold, [1]).

[1] https://gitlab.cern.ch/RooUnfold/RooUnfold

## Table of content:

1. create_tree
2. run_N_unf
3. run_NFNB_unf
4. run_PTN_unf
5. run_pt_unf
6. SMASH_txt

## The procedure

The analysis consists of two steps. 

At **first** step, one has to create the specific tree applying the scripts in "create tree" folder and using the dataset that contains simulated and reconstructed data. For example, in the folder "SMASH_txt" there are several file lists (as an example, 3 out of 200 are committed to git) that in sum are the full statistics of the prepared dataset (Request 6: PWG1 - SMASH, BiBi @ 9.46 GeV, min. bias, GEANT3, 4 млн событий,[https://mpdforum.jinr.ru/t/request-6-pwg1-smash-bibi-9-46-gev-min-bias-geant3/240]). 

Output at this step: an event tree in root file.

Before the second step, one has to install RooUnfold package from https://gitlab.cern.ch/RooUnfold/RooUnfold.

At the **second** step, one can perform unfolding of different quantities. To simplify the procedure and decrease the run-time the codes are presented separately:

`run_N_unf` - performs the one-dimensional unfolding of multiplicity distributions in different pseudorapidity regions

`run_NFNB_unf` - performs the two-dimensional unfolding of NF-NB distributions (the correlation between multiplicities in forward and backward pseudorapidity regions at different distances)

`run_PTN_unf` - performs the two-dimensional unfolding of PT-N distributions (the correlation between scalar sum of the transverse momentum of particles in the event and event multiplicity) in different pseudorapidity regions

`run_pt_unf` - performs the one-dimensional unfolding of particle pT spectrum in different pseudorapidity regions

Output at this step: root file with simulated (pure data), reconstructed (biased data) and unfolded (corrected data) histograms for different charges of particles in several pseudorapidity intervals.

Please note that at each step one needs to properly set paths (in particular, for mpdroot and RooUnfold and others), filenames, and other parameters inside the scripts.

## 1. Create_tree

The file `mpd_create_rm.cpp` creates the root file with event tree that keeps info about sim and rec data simultaneously. This is needed to build the response matrix to perform RooUnfold further. There are several parameter to set up. All sim and rec event and track cuts are presented in this script. If needed one has to modify them here.

The script can be compiled with the attached Makefile. To use `MpdMiniDst` format that is used in this analysis one has to compile in mpdroot as well the `MpdMiniDst` library to have the `libMpdMiniDst.so` file that should be copied to the working directory. After all, one can submit the executable file to cluster by `./launch.sh` -- there will be 200 jobs on files from `SMASH_txt` folder. When it is done (can be checked by "qstat" command) one has to merge 200 output root files in just one tree. This can be done by:

```
root -l
.L mergeScript.cc
mergeScript("path_to_200_output_files/*/*.root", "where_and_with_which_name_to_save_merged_file.root")
```


## 2. Do unfolding (there are 4 different separate possibilities). 

There are parameters to set up in each .cc file. Some parameters to establish before running on cluster are in `my_start.sh` file. The RooUnfold libraries are uploaded up in `runme.C` file. Number of unfolding iterations (since in this example we use Bayesian iterative method) and number of bootstrap samples should be set up in file `launch.sh`.

To submit the job on cluster (the new job for each combination of n_unf_iterations and n_subsamples) one has to type `./launch.sh` in the unfolding folder of interest.

In the output there will be a root file with plenty of histograms (sim, rec, unfolded; all for all charged, positively and negatively charged + other combinations; all at different distances between forward and backward pseudorapidity intervals or at different pseudorapidity distances) that can be used for any analysis.

