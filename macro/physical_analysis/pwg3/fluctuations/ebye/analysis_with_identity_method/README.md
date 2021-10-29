## Intro
Precise determination of the moments of multiplicity distributions of identified particles could be challenging due to the misidentification in detectors. The so-called Identity Method [1-3] allows one to solve this problem. Using the code placed in this folder, performance of the Identity Method was tested on the A–A events simulated in the conditions of the MPD experiment at NICA (using TPC dE/dx vs mometum histograms). Namely, moments within a single kinematic window as well as coefficients of forward-backward pseudorapidity correlations are extracted. 

[1] M. Gazdzicki et al., Phys.Rev. C 83, 054907 (2011)

[2] M. Gorenstein, Phys.Rev. C 84, 024902 (2011)

[3] A. Rustamov, M. Gorenstein, Phys.Rev. C 86, 044906 (2012)

## Table of content:
1. [Analysis of MPD simulated data](#sec1)
2. [Extraction of the special tree with tracks for identity method](#sec2)
3. [Calculation of "true" moments, extraction of dE/dx distributions](#sec3)
4. [Analysis with Identity Method](#sec4)
5. [Plotting final results](#sec5)


The analysis consists of several steps. Note that at each step one needs to properly set paths, filenames, and other parameters inside the scripts.

<a name="sec1"></a>
## 1. Analysis of MPD simulated data
As the first step, trees with event- and track-level information are extracted from the DST root-files (that contain MC and reconstructed data)
using the task `get_tree_MPD_task.cc` from the folder `analysis_task`.
For that, DST files from a chosen directory (e.g. on NICA cluster) are merged into a chain (file names are taken from `inputFilesByHand.txt`).
Execution of the analysis can be started using the `bashRunJobs.sh` script.

<a name="sec2"></a>
## 2. Extraction of the special tree with tracks for identity method
Obtained trees can be then processed using a script `get_info_and_identity_tree.cxx` in the `analysis_of_output_tree` folder.
The script does several things:
1) it extracts a special tree with tracks for Identity Method (called "treeIM") - to be used below,
2) it also performs calculations for Forward-Backward correlations analysis at SIM, REC_PRIMARY and REC_ALL levels,
3) a lot of supplementary information is extracted as well.

<a name="sec3"></a>
## 3. Calculation of "true" moments, extraction of dE/dx distributions
After that, the script `get_dEdx_vs_mom_and_true_pi_K_p_moments_from_tree.C` runs over obtained `an_results.root` file in order to:
1) extract the dEdx vs pTPC distributions and put them in `file_dEdx_in_bins.root` (**IMPORTANT**: this is done with utilizing the PID info from simulations! With real data, this is not possible to do, and sophisticated fits of dE/dx distributions in narrow momentum slices will be required.)
2) calculate the true (i.e. primary-particle level) moments (first, second, cross-moments) and put them in `file_truth_moments_IM.root`.

<a name="sec4"></a>
## 4. Analysis with Identity Method
In the identity_analysis folder, the script `fit_dEdx_pt_slices.C` takes the `file_dEdx_in_bins.root` file and fit narrow momentum slices for each particle species with a function (i.e. Gaus). 
The results for all slices are then stored into `fitted_dEdx_in_pt_bins.root`.

Now everything is ready to run the analysis of the moments with the Identity Method itself. The main script is `analyse_identity_dEdx.C`, which requires the `fitted_dEdx_in_pt_bins.root` and `file_dEdx_in_bins.root` files,
as well as the root-file `an_results.root` with the identity tree "treeIM" (which has a structure shown in the cartoon `structure_tree_for_Identity_Method_analysis_cartoon.png`).
This code is based on the Identity module created by A. Rustamov and M. Arslandok, NIM A 946 (2019) 162622, https://arxiv.org/abs/1807.06370, and also a more advanced script provided by A.Rustamov.
In order to load necessary libraries from the identity_analysis/src directory, the analysis should be started via the `run.C` script.
The main output is a `file_extracted_moments_IM.root` file.

<a name="sec5"></a>
## 5. Plotting final results
Finally, truth moments from `file_truth_moments_IM.root` can be compared with those reconstructed with the Identity Method (`file_extracted_moments_IM.root`) by running a script `draw_comparison_truth_IM.C` in the `draw_output` folder.

A set of results obtained with this code can be found in Altsybeev, Petrov, J. Phys.: Conf. Ser. 1690 012124 (2020), https://arxiv.org/abs/2010.08742.






