/*
 * NicaMpdInterfaceLinkDef.h
 *
 *  Created on: 28 mar 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class NicaMpdEvent + ;
#pragma link C++ class NicaMpdEventInterface + ;
#pragma link C++ class NicaMpdEventTpcPads + ;
#pragma link C++ class NicaMpdTrack + ;
#pragma link C++ class NicaMpdTrackInterface + ;
#pragma link C++ class NicaMpdTrackTpcPads + ;
#pragma link C++ class NicaTrackTpcPads + ;

#pragma link C++ class NicaMpdDstMCEvent + ;
#pragma link C++ class NicaMpdDstMCEventTpcPads + ;

#pragma link C++ class NicaMpdMcEvent + ;
#pragma link C++ class NicaMpdMcEventInterface + ;
#pragma link C++ class NicaMpdMcTrack + ;
#pragma link C++ class NicaMpdMcTrackInterface + ;

#pragma link C++ class NicaGeneratorReadUnigen + ;
#pragma link C++ class NicaGeneratorWriteUnigen + ;
#pragma link C++ class NicaUnigenEvent + ;
#pragma link C++ class NicaUnigenEventInterface + ;
#pragma link C++ class NicaUnigenSource + ;
#pragma link C++ class NicaUnigenTrack + ;
#pragma link C++ class NicaUnigenTrackInterface + ;

#pragma link C++ class NicaMpdHbtComplexEvent + ;
#pragma link C++ class NicaMpdHbtEvent + ;
#pragma link C++ class NicaMpdHbtTrack + ;

#pragma link C++ class NicaMpdDstMCEvent + ;
#pragma link C++ class NicaMpdMiniDstEvent + ;
#pragma link C++ class NicaMpdMiniDstEventInterface + ;
#pragma link C++ class NicaMpdMiniDstMcEvent + ;
#pragma link C++ class NicaMpdMiniDstMcEventInterface + ;
#pragma link C++ class NicaMpdMiniDstMcTrack + ;
#pragma link C++ class NicaMpdMiniDstMcTrackInterface + ;
#pragma link C++ class NicaMpdMiniDstTrack + ;
#pragma link C++ class NicaMpdMiniDstTrackInterface + ;
#pragma link C++ class NicaMpdMiniDstFullEvent + ;
#pragma link C++ class NicaMiniDstSource + ;

#pragma link C++ class NicaMiniDstFullV0Event + ;
#pragma link C++ class NicaMpdMiniDstEventV0 + ;
#pragma link C++ class NicaMpdMiniDstEventV0Interface + ;


#endif /* INTERFACES_MPDROOT_NICAMPDINTERFACELINKDEF_H_ */
