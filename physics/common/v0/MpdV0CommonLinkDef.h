/*
 * MpdPhysicsCommonLinkDef.h
 *
 *  Created on: 29 paź 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class MpdV0CandidateCut + ;
#pragma link C++ class MpdV0CandidateCutBasic + ;
#pragma link C++ class MpdV0DaughterCut + ;
#pragma link C++ class MpdV0DaughterCutBasic + ;
#pragma link C++ class MpdV0CandidateCutKalman + ;

#pragma link C++ class MpdV0FinderBasic + ;
#pragma link C++ class MpdV0FinderHelix + ;
#pragma link C++ class MpdV0FinderKFPackage + ;
#pragma link C++ class MpdV0Matcher + ;

#pragma link C++ namespace MpdCommonV0;
#pragma link C++ enum MpdCommonV0::EParticleType;
#pragma link C++ enum MpdCommonV0::ESigmaType;

#pragma link C++ class MpdV0Particle + ;
#pragma link C++ class MpdV0Track + ;
#pragma link C++ class MpdSimpleLinks < Int_t> + ;

#endif