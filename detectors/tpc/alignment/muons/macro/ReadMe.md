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

