# <b>Classes for physics analyses </b>
Directories contain classes for physics analyses and basic analysis framework


Framework **MpdAnalysisManager** allows simultaneous analysis of same dataset with several analyses tasks. This reduces time to access data or MC and useful for e.g. estimates of systematic uncertainties when different tasks use different sets of cuts etc.

Class **MpdAnalysisManager** requires list of input files, list of branches to be used for analysis and list of tasks to process these files. Finally MpdAnalysisManager takes care of writing output objects for each task with special list. See file photons/macros/ConversionPi0.C as an example for usecase. 

Taks which will be called by MpdAnalysisManager should be derived from **MpdAnalysisTask** and have several methods implemented:
  - void UserInit(); // Users should prepare objects to fill (histograms, trees etc)
  - void ProcessEvent(MpdAnalysisEvent &event) ; // method is called for each event and event data are provided by container MpdAnalysisEvent
  - void Finish() ; //method is called when scan in finished but class data are not written yet.
Example of implementation of such class one can find in photons/MpdConvPi0.h class  

Class **MpdAnalysisEvent** contains references to all branched containing data for this event.  MpdAnalysisManager can be configured to read only few branches reaslly necessary for analysis.

# <b>ebye </b>
# <b>femto </b>
# <b>photons </b>
Directory with classes for photon, pi0 and other neutral meson analyses