/*
 * NicaMpdHelpersLinkDef.h
 *
 *  Created on: 19 lip 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class MpdPIDOnTheFly + ;
#pragma link C++ class NicaPairAnaDebugger + ;
#pragma link C++ class MpdSharedHitMap + ;
#pragma link C++ class MpdHitsOnTheFly + ;
#pragma link C++ class NicaTpcSectorGeo + ;
#pragma link C++ namespace MpdDetectorID;
#pragma link C++ namespace NicaMpdConst;
#endif /* INTERFACES_MPDROOT_NICAMPDINTERFACELINKDEF_H_ */
