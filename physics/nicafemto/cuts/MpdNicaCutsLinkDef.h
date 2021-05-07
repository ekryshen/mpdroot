/*
 * MpdNicaCutsLinkDef.h
 *
 *  Created on: 10 kwi 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef INTERFACES_MPDROOT_CUTS_MPDNICACUTSLINKDEF_H_
#define INTERFACES_MPDROOT_CUTS_MPDNICACUTSLINKDEF_H_

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace MpdPadsFormat;
#pragma link C++ namespace MpdHbtDst;

#pragma link C++ class MpdTrackIsGlobalCut + ;
#pragma link C++ class MpdTrackChargeCut + ;
#pragma link C++ class MpdTrackFirstPointCut + ;
#pragma link C++ class MpdTcpCalibratedCut + ;
#pragma link C++ class MpdPadsEdgeSectorCut + ;
#pragma link C++ class MpdHitsEdgeSectorCut + ;
#pragma link C++ class MpdNSigmaCut + ;
#pragma link C++ class MpdBasicTrackCut + ;

#pragma link C++ class MpdSplittedPairsCut + ;
#pragma link C++ class MpdPairDeltaDCA + ;
#pragma link C++ class MpdPairDeltaPhiStarDeltaEtaMinCut + ;
#pragma link C++ class MpdPairTpcEntrySimpleCut + ;
#pragma link C++ class MpdPairTpcSeparationSimpleCut + ;
#pragma link C++ class MpdPairTpcSimpleCut + ;
#pragma link C++ class MpdCylinderTpcPairCut + ;
#pragma link C++ class MpdPairTcpSeparationSimpleNegCut + ;
#pragma link C++ class MpdPairDeltaPhiStarDeltaEtaCut + ;

#pragma link C++ class MpdPadsFormat::MpdAlicePairQuality + ;
#pragma link C++ class MpdPadsFormat::MpdDeltaEtaDeltaPhiStarAdvancedCut + ;
#pragma link C++ class MpdPadsFormat::MpdHelixSep + ;
#pragma link C++ class MpdPadsFormat::MpdStarPairQualityCut + ;
#pragma link C++ class MpdPadsFormat::MpdModularTpcPairCut + ;
#pragma link C++ class MpdPadsFormat::MpdPairSharedPadsCut + ;
#pragma link C++ class MpdPadsFormat::MpdPairTpcEntranceCut + ;
#pragma link C++ class MpdPadsFormat::MpdTpcPairSeparationCut + ;
#pragma link C++ class MpdPadsFormat::MpdPairTpcEntranceCut2D + ;
#pragma link C++ class MpdPadsFormat::MpdPairSharedHitsCut + ;
#pragma link C++ class MpdPadsFormat::MpdNominalTpcPairPadsDistanceCut + ;
#pragma link C++ class MpdPadsFormat::MpdTpcPadsPairCut + ;

#pragma link C++ class MpdHbtDst::MpdSharedPadsCut + ;
#pragma link C++ class MpdHbtDst::MpdFemtoPairCut + ;
#pragma link C++ class MpdHbtDst::MpdPairTpcEntranceNomininalCut + ;

#pragma link C++ class MpdTpcMonitor + ;
#pragma link C++ class MpdTofMonitor + ;
#pragma link C++ class MpdKinMonitor + ;
#pragma link C++ class MpdDcaMonitor + ;

#endif

#endif /* INTERFACES_MPDROOT_CUTS_MPDNICACUTSLINKDEF_H_ */
