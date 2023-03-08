**MPD TPC detector API** (Design by Contract)

- API contains abstract module interfaces, abstract primitives, base class invariants for TPC detector encapsulated in library libtpc.so
- all MPD TPC modules must implement this API. Implementations of specific ModuleName are encapsulated in library libtpcModuleName.so.
- module performance is subject to testing by Acceptance TDD paradigm. Tests access only API entities (they do not access implementation details) and are by definition module requirements translated into computer language.

---

*STATUS*  
*Abstract module interfaces*  
AbstractTpcClusterHitFinder

*Abstract primitives*  
AbstractTpcDigit  
AbstractTpc2dCluster  
AbstractTpcHit  

*Base class invariants*  
BaseTpcSectorGeo

*IMPLEMENTATIONS*  
alignment - alignment of misaligned data module  
clusterHitFinder - cluster finding and extracting hits from clusters module  
digitizer - digitization of Monte Carlo data for detector simulation purposes module   
geometry - various geometry implementations module   
pid - working out the particle ID module  

---

NOTES:  
**TPC API is work in progress and NOT complete yet.**  
- fairTpc directory contains files that will not be ported to API. 
They all inherit directly from FairRoot (except for MpdTpcEDepParams) and will be removed or completely replaced in the future.
- legacy directory are unused old files not included in the build. 
- clusterHitFinder directory (tpcClusterHitFinder library) apart from implementations of hit finding algorithms contains also
hitProducer folder with MpdTpcHitProducer class. This algorithm directly produces hits instead of finding them from digitized data and will not be subject to API port.
