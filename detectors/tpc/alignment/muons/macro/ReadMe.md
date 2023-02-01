 **AFile.C** ---  macro generate the file  with a TPC alignment
~~~ 
void AFile(Double_t rcm,Double_t adg,Int_t Print, const char* FileName="adir/afile.root")
~~~           
 
[-rcm,+rcm] - the interval of the random position of the center of local sector coordinate system, cm  
[-adg,+adg] - the interval of the random Euler's angles of the of local sector coordinate system, angular degrees  
Print == 0  - no any data printed in stdout

*Examples* 
~~~
mkdir -p ~/alignment/data ~/alignment/log
root -b -l -q 'AFile.C(0.5,0.5,1,"~/alignment/data/0r5_0g5.root")' > ~/alignment/log/AFile.txt
~~~
zero alignment
~~~
root -b -l -q 'AFile.C(0.,0.,1,"~/alignment/data/0r0_0g0.root")' > ~/alignment/log/AFileZero.txt
~~~

-----------------------------------------------------------------------------  

-----------------------------------------------------------------------------  
  
**missMC.C** --- The macro generate a sample of miniMC events from laser rays. 
~~~
void missMC(const char *aReal = "0r5_0g5.root", const char *aUsed = "0r0_0g0.root", const char *outDir = "miniDST",
            const char *outFile = "lrc0_0r5_0g5_cl0_dt08evt10000")
~~~   
- aReal: file name with simulated alignment
- aUsed: file name with zero alignment
- outDir: directory where output root file with track hits is stored
- outFile: name of the output file without the .root suffix   

Detailed description of other parameters is in file TpcMissAlignment.h
 
 *Examples*
 ~~~
 mkdir -p ~/alignment/miniDST
 root -b -l -q 'missMC.C("~/alignment/data/0r5_0g5.root", "~/alignment/data/0r0_0g0.root", "~/alignment/miniDST", "lrc0_0r5_0g5_cl0_dt08evt10000")' > ~/alignment/log/lrc0_0r5_0g5_cl0_dt08evt10000.txt
 ~~~
 - generate file: ~/alignment/miniDST/lrc0_0r5_0g5_cl0_dt08evt10000.root with track hits from laser rays
  
-----------------------------------------------------------------------------   

-----------------------------------------------------------------------------  

 **missMF.C** --- The macro generate a sample of miniMC events in a magnetic field.
 ~~~
 void missMF(const char *aReal = "0r5_0g5.root", const char *aUsed = "0r0_0g0.root", const char *outDir = "miniDST",
             const char *outFile = "0r5_0g5_80R130_dt08_evt10000")
 ~~~
- aReal: file name with simulated alignment
- aUsed: file name with zero alignment
- outDir: directory where output root file with track hits is stored
- outFile: name of the output file without the .root suffix     

Detailed description of other parameters is in file TpcMissAlignment.h
 
 *Examples* 
 ~~~
 root -b -l -q 'missMF.C("~/alignment/data/0r5_0g5.root", "~/alignment/data/0r0_0g0.root", "~/alignment/miniDST", "0r5_0g5_80R130_dt08_evt10000")' > ~/alignment/log/0r5_0g5_80R130_dt08_evt10000.txt
 ~~~
 - generate file: ~/alignment/miniDST/0r5_0g5_80R130_dt08_evt10000.root with track hits from muons. Magnetic field is on. Muon momentums
 are in the interval when radii of the helix trajectory is  in the interval [80,130]cm 
     
-----------------------------------------------------------------------------   

-----------------------------------------------------------------------------   
  
**tpcAlignment.C** 
~~~
tpcAlignment.C(const char* InDataFile,const char* InAlignmentFile=" ")
~~~ 
- the macro uses track hits, found by a given alignment, to find the alignment which produces a minimum xi2 for the sum of all tracks.  
- InDataFile - the file with track hits and the information about the detector magnetic field.  
- InAlignmentFile==" " - InDataFile has the miniMC DST stucture, which includes real and used for reconstruction alignments and track hits.  

*Important*:  
- to be able to use tpcAlignment.C macro you must recompile tpcAlignment library with ALIGNMENT_XI2 option
~~~
cd ~/mpdroot/build && rm -rf lib/libtpcAlignment.*
cmake -DALIGNMENT_XI2=ON .. && make -j8 install
~~~  
- if you want to use again scripts missMC.C and missMF.C you have to recompile without ALIGNMENT_XI2 option
~~~
cd ~/mpdroot/build && rm -rf lib/libtpcAlignment.*
cmake .. && make -j8 install
~~~  
                             
 *Examples*   
 ~~~
 root -b -l -q 'tpcAlignment.C("~/alignment/miniDST/lrc0_0r5_0g5_cl0_dt08evt10000.root"," ")' > ~/alignment/log/tpcAlignment.txt
 ~~~
 The input file has track hits reconstructed using zero-alignment when the real one is from the file ~/alignment/data/0r5_0g5.root.
  
-----------------------------------------------------------------------------   

-----------------------------------------------------------------------------   
  