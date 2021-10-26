# <b>The MpdRoot Framework </b>
The main framework for the MPD experiment at NICA.  
http://mpdroot.jinr.ru/  

<img src="eventdisplay/evepic.png" width="600">

Based on:  
[FairSoft](https://github.com/FairRootGroup/FairSoft) 
and 
[FairRoot](https://github.com/FairRootGroup/FairRoot)

[Getting started](http://mpdroot.jinr.ru/mpdroot-start-guide/)

### Install Prerequisites 
Look for your OS in the link below to install required dependencies  
[Dependencies link](https://github.com/FairRootGroup/FairSoft/blob/master/legacy/dependencies.md)

Additionally development package for fftw3 library is required to install MpdRoot
(fftw3-devel for RedHat based OS, fftw3-dev for Debian based OS)

NOTE: If you are new user we strongly suggest you install FairSoft and FairRoot by copying the scripts/install_fairApr21 to your $HOMEDIR
and run:  
```
chmod +x install_fairApr21
./install_fairApr21 ~/fairApr21 8
```
- parameter ~/fairApr21 is the directory where your local installation of FairSoft/FairRoot will be located  
- parameter 8 is the number of cpu_threads used in compilation  
(you can of course, replace these by your own values)  
- cloned FairSoft/FairRoot repositories and their builds will be located in fair_build directory  

Proceed with installation of MpdRoot by copying scripts/install_mpdroot to your $HOMEDIR and run:  
```
chmod +x install_mpdroot
./install_mpdroot --developer 8
```
if you are a developer with ssh access or:  
```
chmod +x install_mpdroot
./install_mpdroot --regular 8
```
if you are a regular user.  
Again 8 is the number of cpu_threads used in compilation.  

Once the script succeeds, your installation will be located in mpdroot directory.  
If the script fails, please email the output to hnatics@jinr.ru   

#### NOTE: When using MPDRoot SetEnv.sh and config.sh must be invoked in each new terminal by
```
source ~/mpdroot/SetEnv.sh && source ~/mpdroot/build/config.sh
```
alternatively, you can  add this line to your ~/.bashrc file, if you don't want to type this each time you open a new terminal


## Manual Install
  (Advanced users)  
  
The FairSoft, FairRoot and MpdRoot installation steps below are explanation of the installation procedure from the above scripts  

### Install FairSoft 

Set installation directory  
```
export fair_dir=~/fairApr21  
```
Set number of cpu threads used in compilation  
```
export cpu_threads=8
```
Set working directory for builds  
```
export work_dir=~/fair_build
```
Create and enter working directory  
```
mkdir $work_dir && cd $work_dir
```
Clone FairSoft apr21p1 release into working dir, create & enter build directory, run cmake configuration  
```
git clone -b apr21p1 https://github.com/FairRootGroup/FairSoft.git FairSoft  
cd FairSoft && mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$fair_dir/FairSoft -DGEANT4MT=OFF ..

```
Build Fairsoft  
(note: build is couple times restarted as some files downloaded during the build are broken and must be fixed in between.
All questions about this "installation feature" should be directed to the FairSoft developers.)  
```
make -j $cpu_threads
sed -i '/#include "UserDefaults.h"/a #include <thread>' Source/dds/dds-info/src/main.cpp
make -j $cpu_threads
sed -i '/#include <vector>/a #include <thread>' Source/fairmq/fairmq/sdk/DDSSession.cxx
sed -i 's/  int dummy;/  int dummy = 0;/g' Source/fairmq/extern/googletest/googletest/src/gtest-death-test.cc
sed -i 's/  bool result;/  bool result = false;/g' Source/fairmq/extern/googletest/googletest/src/gtest-death-test.cc
make -j $cpu_threads
cd $fair_dir/FairSoft/lib
ln -s libpythia6.so libPythia6.so
export SIMPATH=$fair_dir/FairSoft
```

### Install FairRoot

Enter working directory  

```
cd $work_dir
```
   
Clone FairRoot v18.6.4 release into working dir, create & enter build directory, run cmake configuration 
```
git clone -b v18.6.4 https://github.com/FairRootGroup/FairRoot.git FairRoot
cd FairRoot && mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$fair_dir/FairRoot ..

```
Build and install FairRoot
```
make -j $cpu_threads
make install
```

### Install MpdRoot
dev - developer branch - currently the only supported branch   
  
HTTPS (read-only access, e.g. for regular users)
```        
git clone -b dev --recursive https://git.jinr.ru/nica/mpdroot.git  

```

SSH (for developers)  

``` 
git clone -b dev --recursive git@git.jinr.ru:nica/mpdroot.git 
```

Edit mpdroot/SetEnv.sh file and replace the lines
```
export SIMPATH=/opt/fairsoft/install
export FAIRROOTPATH=/opt/fairroot/install
```

with the ones having correct FairSoft/FairRoot environment variables

```  
export SIMPATH=~/fairApr21/FairSoft
export FAIRROOTPATH=~/fairApr21/FairRoot
```

Prevent the local  SetEnv.sh file from  being overwritten by the future git synchronizations
```
git rm --cached SetEnv.sh
```

Source the environment file
```
cd mpdroot
. SetEnv.sh
```

Create & enter local build dir, run cmake config
```
mkdir build && cd build && cmake .. 
```

Build 
```
make -j $cpu_threads
```

Source new configuration
```
. config.sh  
```


