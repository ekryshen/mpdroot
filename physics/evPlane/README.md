# <b>evPlane wagon for event plane measurements </b>

This wagon performs measurements of the event plane using FHCal (first order event plane based on v1 signal) and TPC (second order event plane based on v2 signal).
Requires `evCentrality` wagon.

## Basic structure
The main classes are in the evPlane directory:

* `MpdEventPlaneAllParams`: stores input parameters from the config txt file, tools to parse input config txt files
* `MpdEventPlaneAll`: main class that performs event plane measurements

## Input file parameters
Input txt file (see `evPlane/macros/pEP.txt` for example) can be configured as follows:
```
#-------Parameters used for analysis------
# Event selection: 
mZvtxCut 130 #  cut on vertex z coordinate

# Track selection: 
mNofHitsCut   16  # minimal number of hits to accept track
mEtaCut      1.5  # maximal pseudorapidity accepted
mEtaGapCut   0.1  # minimal pseudorapidity accepted: abs(eta)>0.05 for mEtaGap=0.1
mPtminCut    0.1  # minimal pt used in analysis
mPtmaxCut    2.0  # maximal pt used in analysis
mDcaCut      2.0  # maximal DCA accepted
# Event plane corrections:
mInFileEpCorr   pEpQa.root  # input file with QA histograms and EP corrections profiles
``` 

One can comment lines or parts of the line by using `#`.

## Usage
The wagon performs event plane corrections aimed to reduce any bias coming from the non-uniform (azimuthal) detector acceptance.
In case of the ideal azimuthal coverage of TPC and FHCal (when all sections/modules of the detector are working properly) those corrections are optional.
But if there are some issues with azimuthal acceptance (inefficiencies, etc.), those corrections should be done.
The EP corrections are the following:
* Recentering;
* Shift (or flattening).
One can found more information in [this article](https://arxiv.org/pdf/nucl-ex/9805001.pdf).

In order to do this corrections, main analysis macro (`evPlane/macros/RunAnalyses.C`) should run 3 times since the corrections are applied iteratively: first run collects information for the recentering, second run applies recentering and collects information for the shift correction, the third run applies both recentering and shift corrections. For the first run one should put `ANY` for the `mInFileEpCorr` parameter. On the second and the third run `mInFileEpCorr` should take output (`pEP.root` by default) from the previous run. Event plane resolution as a function of centrality can be calculated then using `evPlane/macros/getResolution.C` which has 2 arguments: input file (`pEP.root` by default) and output file (`pEP_resolutions.root` by default).

If one wants to use event plane angle in their analysis wagon, it is stored in the `MpdAnalysisEvent *event`:
```
event.fMpdEP.GetPhiEP_FHCal_F_all() // event plane angle from FHCal N+S
event.fMpdEP.GetPhiEP_FHCal_N_all() // event plane angle from FHCal N (eta<0)
event.fMpdEP.GetPhiEP_FHCal_S_all() // event plane angle from FHCal S (eta>0)
event.fMpdEP.GetPhiEP_TPC_N_all()   // event plane angle from TPC N (eta<0)
event.fMpdEP.GetPhiEP_TPC_S_all()   // event plane angle from TPC S (eta>0)
```

similarly to `event.getCentrTPC()` to get event centrality.