# <b>The MpdRoot Framework </b>
The main framework for the MPD experiment at NICA.  
http://mpd.jinr.ru/  

<img src="eventdisplay/evepic.png" width="600">

Based on:  
[FairSoft](https://github.com/FairRootGroup/FairSoft) 
and 
[FairRoot](https://github.com/FairRootGroup/FairRoot)

[Getting started](http://mpd.jinr.ru/mpdroot-start-guide/)

### Install Prerequisites
for RedHat-based OS (eg, CentOS, Scientific Linux):
```   
sudo su 
yum install subversion git make cmake gcc-gfortran gcc-c++ \
binutils file patch redhat-lsb-core libX11-devel libXmu-devel \
libXpm-devel libXft-devel libXext-devel mesa-libGLU-devel \
libxml2-devel expat-devel zlib-devel postgresql-devel \
mysql-devel openssl-devel curl-devel automake libtool fftw3-devel \
```   
for Debian-based OS (eg, Ubuntu):
```    
sudo su 
apt install subversion git make cmake g++ gcc gfortran binutils \
patch lsb-release libx11-dev libxmu-dev libxpm-dev libxft-dev \
libxext-dev dpkg-dev xlibmesa-glu-dev libglew-dev libxml2-dev \
libexpat1-dev zlib1g-dev libpqxx-dev libmysqlclient-dev libssl-dev \
libcurl4-openssl-dev automake libtool fftw3-dev 
```   
### Install FairSoft 
After cloning FairSoft make sure to check fairsoft/DEPENDENCIES  
```
export INSTALLATION_PATH=/opt  
cd $INSTALLATION_PATH  
git clone https://github.com/FairRootGroup/FairSoft.git fairsoft  
cd fairsoft  
git checkout jun19_patches
./configure.sh

 1) GCC (on Linux)  
 1) No Debug Info / or Optimize and Debug if you prefer  
 2) No (do not install FairMQ Only)  
 1) Yes (install Simulation engines and event generators)  
 2) Internet (install G4 files from internet)  
 2) No (do not compile Geant4 in multihreaded mode)  
 2) No (do not install the python bindings)   
 path: $INSTALLATION_PATH/fairsoft/install  
```     
 #### Install FairRoot
```   
export INSTALLATION_PATH=/opt  
cd $INSTALLATION_PATH  
export SIMPATH=$INSTALLATION_PATH/fairsoft/install
export PATH=$SIMPATH/bin:$PATH
git clone https://github.com/FairRootGroup/FairRoot.git fairroot
cd fairroot 
git checkout v18.2_patches
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX="$INSTALLATION_PATH/fairroot/install" ..
make
make install
```       
#### Install MpdRoot
dev - developer branch - contains latest changes  
pro - master stable branch    
  
HTTPS (read-only access, e.g. for regular users)
```        
git clone -b dev --recursive https://git.jinr.ru/nica/mpdroot.git  
```  
SSH (for developers)  
``` 
git clone -b dev --recursive git@git.jinr.ru:nica/mpdroot.git 
```  

    cd mpdroot 
    mkdir build
    cp config/SetEnv.sh.in SetEnv.sh
By default, in the SetEnv.sh file SIMPATH points to /opt/fairsoft/mpd/new, and   
FAIRROOTPATH – /opt/fairroot/mpd/new directories.  If you installed FairSoft or FairRoot to another directory,  
please, change SIMPATH and FAIRROOTPATH variables in the file to correct install paths.  
e.g.  
```   
export SIMPATH=/opt/fairsoft/install  
export FAIRROOTPATH=/opt/fairroot/install  
```   
```   
    . SetEnv.sh  
    cd build  
    cmake ..  
    make  
    . config.sh  
```   
Make sure you use the config file in each new terminal.

If you are using Visual Studio Code as your IDE, you can use
provided configuration files. In the main directory run the
command
```
     cp -r config/vscode.in/ .vscode
```
and edit paths inside files .vscode/launch.json and
.vscode/settings.json depending on your configuration.
