// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file LFResonanceTables.h
/// \brief Definitions of tables of resonance decay candidates
///
/// Inspired by StrangenessTables.h, FemtoDerived.h
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>

#ifndef PWGLF_DATAMODEL_LFRESONANCETABLES_H_
#define PWGLF_DATAMODEL_LFRESONANCETABLES_H_

#include <cmath>

#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/Multiplicity.h"

namespace o2::aod
{
/// Resonance Collisions
namespace resocollision
{
enum {
  kECbegin = 0,
  kINEL = 1,
  kINEL10,
  kINELg0,
  kINELg010,
  kTrig,
  kTrig10,
  kTrigINELg0,
  kTrigINELg010,
  kSel8,
  kSel810,
  kSel8INELg0,
  kSel8INELg010,
  kAllCuts,
  kAllCuts10,
  kAllCutsINELg0,
  kAllCutsINELg010,
  kECend,
};
DECLARE_SOA_COLUMN(Cent, cent, float);             //! Centrality (Multiplicity) percentile (Default: FT0M)
DECLARE_SOA_COLUMN(Spherocity, spherocity, float); //! Spherocity of the event
DECLARE_SOA_COLUMN(EvtPl, evtPl, float);           //! Second harmonic event plane
DECLARE_SOA_COLUMN(EvtPlResAB, evtPlResAB, float); //! Second harmonic event plane resolution of A-B sub events
DECLARE_SOA_COLUMN(EvtPlResAC, evtPlResAC, float); //! Second harmonic event plane resolution of A-C sub events
DECLARE_SOA_COLUMN(EvtPlResBC, evtPlResBC, float); //! Second harmonic event plane resolution of B-C sub events
DECLARE_SOA_COLUMN(BMagField, bMagField, float);   //! Magnetic field
// MC
DECLARE_SOA_COLUMN(IsVtxIn10, isVtxIn10, bool);               //! Vtx10
DECLARE_SOA_COLUMN(IsINELgt0, isINELgt0, bool);               //! INEL>0
DECLARE_SOA_COLUMN(IsTriggerTVX, isTriggerTVX, bool);         //! TriggerTVX
DECLARE_SOA_COLUMN(IsInSel8, isInSel8, bool);                 //! InSel8
DECLARE_SOA_COLUMN(IsInAfterAllCuts, isInAfterAllCuts, bool); //! InAfterAllCuts
DECLARE_SOA_COLUMN(ImpactParameter, impactParameter, float);  //! ImpactParameter

} // namespace resocollision
DECLARE_SOA_TABLE(ResoCollisions, "AOD", "RESOCOLLISION",
                  o2::soa::Index<>,
                  o2::aod::mult::MultNTracksPV,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  resocollision::Cent,
                  resocollision::Spherocity,
                  resocollision::EvtPl,
                  resocollision::EvtPlResAB,
                  resocollision::EvtPlResAC,
                  resocollision::EvtPlResBC,
                  resocollision::BMagField,
                  timestamp::Timestamp,
                  evsel::NumTracksInTimeRange);
using ResoCollision = ResoCollisions::iterator;

DECLARE_SOA_TABLE(ResoMCCollisions, "AOD", "RESOMCCOL",
                  resocollision::IsVtxIn10,
                  resocollision::IsINELgt0,
                  resocollision::IsTriggerTVX,
                  resocollision::IsInSel8,
                  resocollision::IsInAfterAllCuts,
                  resocollision::ImpactParameter);
using ResoMCCollision = ResoMCCollisions::iterator;

// Resonance Daughters
// inspired from PWGCF/DataModel/FemtoDerived.h
namespace resodaughter
{

DECLARE_SOA_INDEX_COLUMN(ResoCollision, resoCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);                                   //! p_T (GeV/c)
DECLARE_SOA_COLUMN(Px, px, float);                                   //! p_x (GeV/c)
DECLARE_SOA_COLUMN(Py, py, float);                                   //! p_y (GeV/c)
DECLARE_SOA_COLUMN(Pz, pz, float);                                   //! p_z (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);                                 //! Eta
DECLARE_SOA_COLUMN(Phi, phi, float);                                 //! Phi
DECLARE_SOA_COLUMN(PartType, partType, uint8_t);                     //! Type of the particle, according to resodaughter::ParticleType
DECLARE_SOA_COLUMN(TempFitVar, tempFitVar, float);                   //! Observable for the template fitting (Track: DCA_xy, V0: CPA)
DECLARE_SOA_COLUMN(Indices, indices, int[2]);                        //! Field for the track indices to remove auto-correlations
DECLARE_SOA_COLUMN(CascadeIndices, cascIndices, int[3]);             //! Field for the track indices to remove auto-correlations (ordered: positive, negative, bachelor)
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                              //! Sign of the track charge
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, uint8_t); //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, uint8_t);             //! Number of TPC clusters found
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, uint8_t);                       //! Number of ITS clusters found
DECLARE_SOA_COLUMN(IsGlobalTrackWoDCA, isGlobalTrackWoDCA, bool);    //! Is global track without DCA
DECLARE_SOA_COLUMN(IsGlobalTrack, isGlobalTrack, bool);              //! Is global track
DECLARE_SOA_COLUMN(IsPrimaryTrack, isPrimaryTrack, bool);            //! Is primary track
DECLARE_SOA_COLUMN(IsPVContributor, isPVContributor, bool);          //! Is primary vertex contributor
DECLARE_SOA_COLUMN(HasITS, hasITS, bool);                            //! Has ITS
DECLARE_SOA_COLUMN(HasTPC, hasTPC, bool);                            //! Has TPC
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);                            //! Has TOF
DECLARE_SOA_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(DaughDCA, daughDCA, float);                               //! DCA between daughters
DECLARE_SOA_COLUMN(CascDaughDCA, cascdaughDCA, float);                       //! DCA between daughters from cascade
DECLARE_SOA_COLUMN(V0CosPA, v0CosPA, float);                                 //! V0 Cosine of Pointing Angle
DECLARE_SOA_COLUMN(CascCosPA, cascCosPA, float);                             //! Cascade Cosine of Pointing Angle
DECLARE_SOA_COLUMN(MLambda, mLambda, float);                                 //! The invariant mass of V0 candidate, assuming lambda
DECLARE_SOA_COLUMN(MAntiLambda, mAntiLambda, float);                         //! The invariant mass of V0 candidate, assuming antilambda
DECLARE_SOA_COLUMN(MK0Short, mK0Short, float);                               //! The invariant mass of V0 candidate, assuming k0s
DECLARE_SOA_COLUMN(MXi, mXi, float);                                         //! The invariant mass of Xi candidate
DECLARE_SOA_COLUMN(TransRadius, transRadius, float);                         //! Transverse radius of the decay vertex
DECLARE_SOA_COLUMN(CascTransRadius, casctransRadius, float);                 //! Transverse radius of the decay vertex from cascade
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);                             //! X position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);                             //! Y position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);                             //! Z position of the decay vertex
DECLARE_SOA_COLUMN(DaughterTPCNSigmaPi1, daughterTPCNSigmaPi1, float);       //! TPC PID of the first daughter as Pion
DECLARE_SOA_COLUMN(DaughterTPCNSigmaKa1, daughterTPCNSigmaKa1, float);       //! TPC PID of the first daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTPCNSigmaPr1, daughterTPCNSigmaPr1, float);       //! TPC PID of the first daughter as Proton
DECLARE_SOA_COLUMN(DaughterTPCNSigmaPi2, daughterTPCNSigmaPi2, float);       //! TPC PID of the second daughter as Pion
DECLARE_SOA_COLUMN(DaughterTPCNSigmaKa2, daughterTPCNSigmaKa2, float);       //! TPC PID of the second daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTPCNSigmaPr2, daughterTPCNSigmaPr2, float);       //! TPC PID of the second daughter as Proton
DECLARE_SOA_COLUMN(DaughterTPCNSigmaPiBach, daughterTPCNSigmaPiBach, float); //! TPC PID of the bachelor daughter as Pion
DECLARE_SOA_COLUMN(DaughterTPCNSigmaKaBach, daughterTPCNSigmaKaBach, float); //! TPC PID of the bachelor daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTPCNSigmaPrBach, daughterTPCNSigmaPrBach, float); //! TPC PID of the bachelor daughter as Proton
DECLARE_SOA_COLUMN(DaughterTOFNSigmaPi1, daughterTOFNSigmaPi1, float);       //! TOF PID of the first daughter as Pion
DECLARE_SOA_COLUMN(DaughterTOFNSigmaKa1, daughterTOFNSigmaKa1, float);       //! TOF PID of the first daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTOFNSigmaPr1, daughterTOFNSigmaPr1, float);       //! TOF PID of the first daughter as Proton
DECLARE_SOA_COLUMN(DaughterTOFNSigmaPi2, daughterTOFNSigmaPi2, float);       //! TOF PID of the second daughter as Pion
DECLARE_SOA_COLUMN(DaughterTOFNSigmaKa2, daughterTOFNSigmaKa2, float);       //! TOF PID of the second daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTOFNSigmaPr2, daughterTOFNSigmaPr2, float);       //! TOF PID of the second daughter as Proton
DECLARE_SOA_COLUMN(DaughterTOFNSigmaPiBach, daughterTOFNSigmaPiBach, float); //! TOF PID of the bachelor daughter as Pion
DECLARE_SOA_COLUMN(DaughterTOFNSigmaKaBach, daughterTOFNSigmaKaBach, float); //! TOF PID of the bachelor daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTOFNSigmaPrBach, daughterTOFNSigmaPrBach, float); //! TOF PID of the bachelor daughter as Proton
// For MC
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! Index of the corresponding MC particle
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);
DECLARE_SOA_COLUMN(ProducedByGenerator, producedByGenerator, bool);
DECLARE_SOA_COLUMN(MothersId, motherId, int);        //! Id of the mother particle
DECLARE_SOA_COLUMN(MotherPDG, motherPDG, int);       //! PDG code of the mother particle
DECLARE_SOA_COLUMN(DaughterPDG1, daughterPDG1, int); //! PDG code of the first Daughter particle
DECLARE_SOA_COLUMN(DaughterPDG2, daughterPDG2, int); //! PDG code of the second Daughter particle
DECLARE_SOA_COLUMN(DaughterID1, daughterId1, int);   //! Id of the first Daughter particle
DECLARE_SOA_COLUMN(DaughterID2, daughterId2, int);   //! Id of the second Daughter particle
DECLARE_SOA_COLUMN(SiblingIds, siblingIds, int[2]);  //! Index of the particles with the same mother
DECLARE_SOA_COLUMN(BachTrkID, bachtrkID, int);       //! Id of the bach track from cascade
DECLARE_SOA_COLUMN(V0ID, v0ID, int);                 //! Id of the V0 from cascade
} // namespace resodaughter
DECLARE_SOA_TABLE(ResoTracks, "AOD", "RESOTRACKS",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::Eta,
                  resodaughter::Phi,
                  resodaughter::Sign,
                  resodaughter::TPCNClsCrossedRows,
                  resodaughter::TPCNClsFound,
                  resodaughter::ITSNCls,
                  o2::aod::track::DcaXY,
                  o2::aod::track::DcaZ,
                  o2::aod::track::X,
                  o2::aod::track::Alpha,
                  resodaughter::HasITS,
                  resodaughter::HasTPC,
                  resodaughter::HasTOF,
                  o2::aod::pidtpc::TPCNSigmaPi,
                  o2::aod::pidtpc::TPCNSigmaKa,
                  o2::aod::pidtpc::TPCNSigmaPr,
                  o2::aod::pidtpc::TPCNSigmaEl,
                  o2::aod::pidtof::TOFNSigmaPi,
                  o2::aod::pidtof::TOFNSigmaKa,
                  o2::aod::pidtof::TOFNSigmaPr,
                  o2::aod::pidtof::TOFNSigmaEl,
                  o2::aod::track::TPCSignal,
                  o2::aod::track::PassedITSRefit,
                  o2::aod::track::PassedTPCRefit,
                  resodaughter::IsGlobalTrackWoDCA,
                  resodaughter::IsGlobalTrack,
                  resodaughter::IsPrimaryTrack,
                  resodaughter::IsPVContributor,
                  resodaughter::TPCCrossedRowsOverFindableCls,
                  o2::aod::track::ITSChi2NCl,
                  o2::aod::track::TPCChi2NCl);
using ResoTrack = ResoTracks::iterator;

DECLARE_SOA_TABLE(ResoV0s, "AOD", "RESOV0S",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::Eta,
                  resodaughter::Phi,
                  resodaughter::Indices,
                  resodaughter::DaughterTPCNSigmaPi1,
                  resodaughter::DaughterTPCNSigmaKa1,
                  resodaughter::DaughterTPCNSigmaPr1,
                  resodaughter::DaughterTPCNSigmaPi2,
                  resodaughter::DaughterTPCNSigmaKa2,
                  resodaughter::DaughterTPCNSigmaPr2,
                  resodaughter::DaughterTOFNSigmaPi1,
                  resodaughter::DaughterTOFNSigmaKa1,
                  resodaughter::DaughterTOFNSigmaPr1,
                  resodaughter::DaughterTOFNSigmaPi2,
                  resodaughter::DaughterTOFNSigmaKa2,
                  resodaughter::DaughterTOFNSigmaPr2,
                  resodaughter::V0CosPA,
                  resodaughter::DaughDCA,
                  v0data::DCAPosToPV,
                  v0data::DCANegToPV,
                  v0data::DCAV0ToPV,
                  resodaughter::MLambda,
                  resodaughter::MAntiLambda,
                  resodaughter::MK0Short,
                  resodaughter::TransRadius,
                  resodaughter::DecayVtxX,
                  resodaughter::DecayVtxY,
                  resodaughter::DecayVtxZ);
using ResoV0 = ResoV0s::iterator;

DECLARE_SOA_TABLE(ResoCascades, "AOD", "RESOCASCADES",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::Eta,
                  resodaughter::Phi,
                  resodaughter::CascadeIndices,
                  resodaughter::DaughterTPCNSigmaPi1,
                  resodaughter::DaughterTPCNSigmaKa1,
                  resodaughter::DaughterTPCNSigmaPr1,
                  resodaughter::DaughterTPCNSigmaPi2,
                  resodaughter::DaughterTPCNSigmaKa2,
                  resodaughter::DaughterTPCNSigmaPr2,
                  resodaughter::DaughterTPCNSigmaPiBach,
                  resodaughter::DaughterTPCNSigmaKaBach,
                  resodaughter::DaughterTPCNSigmaPrBach,
                  resodaughter::DaughterTOFNSigmaPi1,
                  resodaughter::DaughterTOFNSigmaKa1,
                  resodaughter::DaughterTOFNSigmaPr1,
                  resodaughter::DaughterTOFNSigmaPi2,
                  resodaughter::DaughterTOFNSigmaKa2,
                  resodaughter::DaughterTOFNSigmaPr2,
                  resodaughter::DaughterTOFNSigmaPiBach,
                  resodaughter::DaughterTOFNSigmaKaBach,
                  resodaughter::DaughterTOFNSigmaPrBach,
                  resodaughter::V0CosPA,
                  resodaughter::CascCosPA,
                  resodaughter::DaughDCA,
                  resodaughter::CascDaughDCA,
                  cascdata::DCAPosToPV,
                  cascdata::DCANegToPV,
                  cascdata::DCABachToPV,
                  v0data::DCAV0ToPV,
                  cascdata::DCAXYCascToPV,
                  cascdata::DCAZCascToPV,
                  cascdata::Sign,
                  resodaughter::MXi,
                  resodaughter::TransRadius,
                  resodaughter::CascTransRadius,
                  resodaughter::DecayVtxX,
                  resodaughter::DecayVtxY,
                  resodaughter::DecayVtxZ);
using ResoCascade = ResoCascades::iterator;

DECLARE_SOA_TABLE(ResoMCTracks, "AOD", "RESOMCTRACKS",
                  mcparticle::PdgCode,
                  resodaughter::MothersId,
                  resodaughter::MotherPDG,
                  resodaughter::SiblingIds,
                  resodaughter::IsPhysicalPrimary,
                  resodaughter::ProducedByGenerator);
using ResoMCTrack = ResoMCTracks::iterator;

DECLARE_SOA_TABLE(ResoMCV0s, "AOD", "RESOMCV0S",
                  mcparticle::PdgCode,
                  resodaughter::MothersId,
                  resodaughter::MotherPDG,
                  resodaughter::DaughterID1,
                  resodaughter::DaughterID2,
                  resodaughter::DaughterPDG1,
                  resodaughter::DaughterPDG2,
                  resodaughter::IsPhysicalPrimary,
                  resodaughter::ProducedByGenerator);
using ResoMCV0 = ResoMCV0s::iterator;

DECLARE_SOA_TABLE(ResoMCCascades, "AOD", "RESOMCCASCADES",
                  mcparticle::PdgCode,
                  resodaughter::MothersId,
                  resodaughter::MotherPDG,
                  resodaughter::BachTrkID,
                  resodaughter::V0ID,
                  resodaughter::DaughterPDG1,
                  resodaughter::DaughterPDG2,
                  resodaughter::IsPhysicalPrimary,
                  resodaughter::ProducedByGenerator);
using ResoMCCascade = ResoMCCascades::iterator;

DECLARE_SOA_TABLE(ResoMCParents, "AOD", "RESOMCPARENTS",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::McParticleId,
                  mcparticle::PdgCode,
                  resodaughter::DaughterPDG1,
                  resodaughter::DaughterPDG2,
                  resodaughter::IsPhysicalPrimary,
                  resodaughter::ProducedByGenerator,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::Eta,
                  resodaughter::Phi,
                  mcparticle::Y);
using ResoMCParent = ResoMCParents::iterator;

using Reso2TracksExt = soa::Join<aod::FullTracks, aod::TracksDCA>; // without Extra
using Reso2TracksMC = soa::Join<aod::FullTracks, McTrackLabels>;
using Reso2TracksPID = soa::Join<aod::FullTracks, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCEl, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFEl>;
using Reso2TracksPIDExt = soa::Join<Reso2TracksPID, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension>; // Without Extra

} // namespace o2::aod
#endif // PWGLF_DATAMODEL_LFRESONANCETABLES_H_
