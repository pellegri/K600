//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//      ----------------------------------------------------------------
//                      K600 Spectrometer (iThemba Labs)
//      ----------------------------------------------------------------
//
//      Github repository: https://www.github.com/KevinCWLi/K600
//
//      Main Author:    K.C.W. Li
//
//      email: likevincw@gmail.com
//

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

#include "G4UniformMagField.hh"
#include "G4QuadrupoleMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

//#include "TabulatedField3D.hh"
#include "G4PropagatorInField.hh"
#include "G4PropagatorInField.hh"
#include "G4FieldManager.hh"


class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;

///////////////////////////////
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4Isotope;
class G4Element;
class G4LogicalBorderSurface;
class G4LogicalSkinSurface;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class G4ChordFinder;
class G4PropagatorInField;
class G4FieldManager;
class G4UniformMagField;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//////////////////////////////////////////////////////////
//                  DETECTOR ARRAY SETUP                //
//////////////////////////////////////////////////////////

///////////////     VDC DETECTORS     ///////////////////
const G4int     numberOf_VDC = 2;

///////////////     TIARA DETECTORS     ///////////////////
const G4int     numberOf_TIARA = 5;

///////////////     CLOVER DETECTORS     ///////////////////
const G4int     numberOf_CLOVER = 9;
const G4int     numberOf_CLOVER_Shields = 9;


///////////////     TIGRESS DETECTORS     ///////////////////
const G4int     numberOf_TIGRESS = 1;
const G4int     numberOf_TIGRESS_BGO = 1;


///////////////     PADDLE DETECTORS     ///////////////////
const G4int     numberOf_PADDLE = 3;


///////////////     LEPS DETECTORS     ///////////////////
const G4int     numberOf_LEPS = 8;

///////////////     NAIS DETECTORS     ///////////////////
const G4int     numberOf_NAIS = 5;

///////////////     FLATSIDE CHAMBER     ///////////////////
const G4int     numberOf_FLATSIDE = 2;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();
    
public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructField();
    
    // get methods
    //
    //const G4VPhysicalVolume* GetAbsorberPV() const;
    //const G4VPhysicalVolume* GetGapPV() const;
    
private:
    // methods
    //
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
    
    // data members
    //
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;
    // magnetic field messenger
    
    G4VPhysicalVolume*   fAbsorberPV; // the absorber physical volume
    G4VPhysicalVolume*   fGapPV;      // the gap physical volume
    
    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
    
    /////////////////////////////
    //          WORLD
    G4double WorldSize;
    
    
    ////////////////////////////////////
    //      VERTICAL DRIFT CHAMBERS
    ////////////////////////////////////
    
    G4bool              VDC_AllPresent_Override;
    G4bool              VDC_AllAbsent_Override;
    G4bool              VDC_Presence[numberOf_VDC];
    G4RotationMatrix    VDC_rotm[numberOf_VDC];
    G4Transform3D       VDC_transform[numberOf_VDC];
    G4ThreeVector       VDC_CentrePosition[numberOf_VDC];
    G4double            VDC_CentrePositionX[numberOf_VDC];
    G4double            VDC_CentrePositionZ[numberOf_VDC];
    G4double            VDC_RotationY[numberOf_VDC];
    
    //  VDC ASSEMBLY
    G4VPhysicalVolume*  Physical_VDC_ASSEMBLY;
    
    //  VDC SENSE REGION USDS
    G4VPhysicalVolume*  Physical_VDC_SenseRegion_USDS;
    
    //  VDC - STESALIT HV Frames, Upstream, Middle and Downstream
    G4RotationMatrix    VDC_HVFrame_rotm[3];
    G4Transform3D       VDC_HVFrame_transform;
    
    //  VDC - GAS FRAME
    G4RotationMatrix    VDC_GasFrame_rotm;
    G4Transform3D       VDC_GasFrame_transform;
    
    //  VDC - ALUMINIUM OUTER FRAME
    G4RotationMatrix    VDC_Al_Frame_rotm[2];
    G4Transform3D       VDC_Al_Frame_transform;
    
    //      VDC - X WIRES
    G4VPhysicalVolume*  PhysiVDC_X_WIRE;
    G4RotationMatrix    VDC_X_WIRE_rotm;
    G4Transform3D       VDC_X_WIRE_transform;
    
    //      VDC - U WIRES
    G4VPhysicalVolume*  PhysiVDC_U_WIRE;
    G4RotationMatrix    VDC_U_WIRE_rotm;
    G4Transform3D       VDC_U_WIRE_transform;
    
    
    
    /////////////////////////////////////
    //          TIARA DETECTORS
    /////////////////////////////////////
    
    G4bool              TIARA_AllPresent_Override;
    G4bool              TIARA_AllAbsent_Override;
    G4bool              TIARA_Presence[numberOf_TIARA];
    G4RotationMatrix    TIARA_rotm[numberOf_TIARA];
    G4Transform3D       TIARA_transform[numberOf_TIARA];
    G4ThreeVector       TIARA_AA_CentrePosition[numberOf_TIARA];
    G4double            offset_TIARA_BeamAxis;
    
    G4RotationMatrix    TIARA_AA_RS_rotm[8];
    G4Transform3D       TIARA_AA_RS_transform[8];
    
    G4RotationMatrix    TIARA_DL_2M_rotm;
    G4Transform3D       TIARA_DL_transform[2];
    G4Transform3D       TIARA_2M_transform[2];
    
    //  TIARA Active Area - Rings and Sectors
    G4VPhysicalVolume*  PhysiTIARA_AA_RS;
    
    //  TIARA Silicon Wafer - Rings and Sectors
    G4VPhysicalVolume*  PhysiTIARA_SiliconWafer;
    G4RotationMatrix    TIARA_SiliconWafer_rotm;
    G4Transform3D       TIARA_SiliconWafer_transform;
    
    
    ///////////////////////////////
    //      PADDLE DETECTORS
    ///////////////////////////////
    
    G4bool              PADDLE_AllPresent_Override;
    G4bool              PADDLE_AllAbsent_Override;
    G4bool              PADDLE_Presence[numberOf_PADDLE];
    G4RotationMatrix    PADDLE_rotm[numberOf_PADDLE];
    G4Transform3D       PADDLE_transform[numberOf_PADDLE];
    G4ThreeVector       PADDLE_CentrePosition[numberOf_PADDLE];
    G4double            PADDLE_CentrePositionX[numberOf_PADDLE];
    G4double            PADDLE_CentrePositionZ[numberOf_PADDLE];
    G4double            PADDLE_RotationY[numberOf_PADDLE];
    
    G4VPhysicalVolume*  PhysiPADDLE;
    
    /////////////////////////////////////
    //              HAGAR
    /////////////////////////////////////
    
    G4bool              HAGAR_NaICrystal_Presence;
    G4bool              HAGAR_Annulus_Presence;
    G4bool              HAGAR_FrontDisc_Presence;
    
    G4ThreeVector       HAGAR_NaICrystal_CentrePosition;
    G4ThreeVector       HAGAR_Annulus_CentrePosition;
    G4ThreeVector       HAGAR_FrontDisc_CentrePosition;
    G4RotationMatrix    HAGAR_rotm;
    G4Transform3D       HAGAR_transform;
    
    //  HAGAR NaI Crystal
    G4VPhysicalVolume*  PhysiHAGAR_NaICrystal;
    
    //  HAGAR Annulus
    G4VPhysicalVolume*  PhysiHAGAR_Annulus;
    
    //  HAGAR Front Disc
    G4VPhysicalVolume*  PhysiHAGAR_FrontDisc;
    
    
    /////////////////////////////////////
    //          CLOVER DETECTORS
    /////////////////////////////////////
    
    G4bool              CLOVER_AllPresent_Override;
    G4bool              CLOVER_AllAbsent_Override;
    G4bool              CLOVER_Presence[numberOf_CLOVER];
    G4double            CLOVER_Distance[numberOf_CLOVER];
    G4RotationMatrix    CLOVER_rotm[numberOf_CLOVER];
    G4Transform3D       CLOVER_transform[numberOf_CLOVER];
    G4ThreeVector       CLOVER_position[numberOf_CLOVER];
    G4double            CLOVER_phi[numberOf_CLOVER];
    G4double            CLOVER_theta[numberOf_CLOVER];
    
    //  CLOVER HPGe Crystals
    G4VPhysicalVolume*  PhysiCLOVER_HPGeCrystal;
    
    
    /////////////////////////////////////
    //          LEPS DETECTORS
    /////////////////////////////////////
    
    G4bool              LEPS_AllPresent_Override;
    G4bool              LEPS_AllAbsent_Override;
    G4bool              LEPS_Presence[numberOf_LEPS];
    G4double            LEPS_Distance[numberOf_LEPS];
    G4RotationMatrix    LEPS_rotm[numberOf_LEPS];
    G4double            LEPS_phi[numberOf_LEPS];
    G4double            LEPS_theta[numberOf_LEPS];
    
    G4Transform3D       LEPS_transform[numberOf_LEPS];
    G4ThreeVector       LEPS_position[numberOf_LEPS];
    
    G4Transform3D       LEPS_InternalVacuum_transform[numberOf_LEPS];
    G4ThreeVector       LEPS_InternalVacuum_position[numberOf_LEPS];
    
    G4Transform3D       LEPS_Window_transform[numberOf_LEPS];
    G4ThreeVector       LEPS_Window_position[numberOf_LEPS];
    
    //      LEPS HPGe Crystals
    G4VPhysicalVolume*  Physical_LEPS_HPGeCrystal;
    G4Transform3D       LEPS_HPGeCrystal_transform[4];
    G4RotationMatrix    LEPS_HPGeCrystal_rotm[4];
    
    /////////////////////////////////////
    //          NAIS DETECTORS
    /////////////////////////////////////
    
    G4bool              NAIS_AllPresent_Override;
    G4bool              NAIS_AllAbsent_Override;
    G4bool              NAIS_Presence[numberOf_NAIS];
    G4double            NAIS_Distance[numberOf_NAIS];
    G4RotationMatrix    NAIS_rotm[numberOf_NAIS];
    G4double            NAIS_phi[numberOf_NAIS];
    G4double            NAIS_theta[numberOf_NAIS];
    
    G4Transform3D       NAIS_transform[numberOf_NAIS];
    G4ThreeVector       NAIS_position[numberOf_NAIS];
    
    
    //      NAIS NaI Crystals
    G4VPhysicalVolume*  PhysiNAIS_NaICrystal;



    /////////////////////////////////////
    //       FLATSIDE CHAMBER
    /////////////////////////////////////
    
    G4bool              FLATSIDE_AllPresent_Override;
    G4bool              FLATSIDE_AllAbsent_Override;
    G4bool              FLATSIDE_Presence[numberOf_FLATSIDE];
    G4double            FLATSIDE_Distance[numberOf_FLATSIDE];
    G4RotationMatrix    FLATSIDE_rotm[numberOf_FLATSIDE];
    G4double            FLATSIDE_phi[numberOf_FLATSIDE];
    G4double            FLATSIDE_theta[numberOf_FLATSIDE];
    
    G4Transform3D       FLATSIDE_transform[numberOf_FLATSIDE];
    G4ThreeVector       FLATSIDE_position[numberOf_FLATSIDE];
    
    
    /////////////////////////////////////
    //       FLATSIDE CHAMBER
    /////////////////////////////////////
    
    G4bool              PARAFFINTUBE_Presence;
    G4bool              IRONBOX_Presence;
    G4bool              AmBeBOX_Presence;
    
    ///////////////////////////////////////////////////////////////
    //          CLOVER - BGO Shield   (Manufacturer: Cyberstar)
    ///////////////////////////////////////////////////////////////
    
    G4bool              CLOVER_Shield_AllPresent_Override;
    G4bool              CLOVER_Shield_AllAbsent_Override;
    G4bool              CLOVER_Shield_Presence[numberOf_CLOVER_Shields];
    G4ThreeVector       CLOVER_Shield_position[numberOf_CLOVER_Shields];
    G4Transform3D       CLOVER_Shield_transform[numberOf_CLOVER_Shields];
    
    //  Shield BGO Crystal Scintillators
    G4VPhysicalVolume*  PhysiCLOVER_Shield_BGOCrystal;
    
    //  Shield PMT Tubes
    G4VPhysicalVolume*  PhysiCLOVER_Shield_PMT;
    
    
    //////////////////////////////////////
    //          K600 SPECTROMETER
    //////////////////////////////////////
    
    //////////////////////////////////////
    //          K600 - QUADRUPOLE
    G4bool              Ideal_Quadrupole;
    G4bool              Mapped_Quadrupole;
    G4bool              K600_Quadrupole;
    
    G4VPhysicalVolume*  PhysiK600_Quadrupole;
    
    G4ThreeVector       K600_Quadrupole_CentrePosition;
    G4RotationMatrix    K600_Quadrupole_rotm;
    G4Transform3D       K600_Quadrupole_transform;
    
    ////    MAGNETIC FIELD for QUADRUPOLE
    static G4ThreadLocal G4FieldManager* fieldManagerMagneticField_K600_Q;
    static G4ThreadLocal G4QuadrupoleMagField* MagneticField_K600_Q;
    
    G4double                K600_Q_gradient;   // gradient = dB/dr
    G4Mag_UsualEqRhs*       fEquationMagneticField_K600_Q;
    G4MagIntegratorStepper* stepperMagneticField_K600_Q;
    G4ChordFinder*          fChordFinder_K600_Q;
    
    
    //////////////////////////////////////
    //          K600 - DIPOLE 1
    G4bool              K600_Dipole1;
    G4VPhysicalVolume*  PhysiK600_Dipole1;
    
    G4ThreeVector       K600_Dipole1_CentrePosition;
    G4RotationMatrix    K600_Dipole1_rotm;
    G4Transform3D       K600_Dipole1_transform;
    
    ////    MAGNETIC FIELD for DIPOLE 1
    static G4ThreadLocal G4FieldManager* fieldManagerMagneticField_K600_D1;
    static G4ThreadLocal G4UniformMagField* MagneticField_K600_D1;
    
    G4double                K600_Dipole1_BZ;
    G4Mag_UsualEqRhs*       fEquationMagneticField_K600_D1;
    G4double                minStepMagneticField;
    G4MagIntegratorStepper* stepperMagneticField_K600_D1;
    G4ChordFinder*          fChordFinder_K600_D1;
    
    //////////////////////////////////////
    //          K600 - DIPOLE 2
    G4bool              K600_Dipole2;
    G4VPhysicalVolume*  PhysiK600_Dipole2;
    
    G4ThreeVector       K600_Dipole2_CentrePosition;
    G4RotationMatrix    K600_Dipole2_rotm;
    G4Transform3D       K600_Dipole2_transform;
    
    ////    MAGNETIC FIELD for DIPOLE 2
    static G4ThreadLocal G4FieldManager* fieldManagerMagneticField_K600_D2;
    static G4ThreadLocal G4UniformMagField* MagneticField_K600_D2;
    
    G4double                K600_Dipole2_BZ;
    G4Mag_UsualEqRhs*       fEquationMagneticField_K600_D2;
    G4MagIntegratorStepper* stepperMagneticField_K600_D2;
    G4ChordFinder*          fChordFinder_K600_D2;
    
    ////////////////////////////////
    ////        STRUCTURES      ////
    ////////////////////////////////
    
    ////////////////////////////////////////////////
    ////    New K600 Target Chamber - New scattering chamber, both sides on
    G4bool      K600_BACTAR_sidesOn_Presence;
    
    ////////////////////////////////////////////////
    ////    New K600 Target Chamber - New scattering chamber, both sides off
    G4bool      K600_BACTAR_sidesOff_Presence;

    ////////////////////////////////////////////////
    ////    New K600 Target Chamber - New scattering chamber, beam left side side off
    G4bool      K600_BACTAR_beamRightSideOff_Presence;

    ////////////////////////////////////////////////
    ////    New K600 Target Chamber - New scattering chamber, beam right side side off
    G4bool      K600_BACTAR_beamLeftSideOff_Presence;

    
    
    /////////////////////////////////////
    //  K600 Target
    G4bool      K600_Target_Presence;
    
    /////////////////////////////////////
    //  K600 Target Backing
    G4bool      K600_TargetBacking_Presence;
    
};

// inline functions
/*
 inline const G4VPhysicalVolume* DetectorConstruction::GetAbsorberPV() const {
 return fAbsorberPV;
 }
 
 inline const G4VPhysicalVolume* DetectorConstruction::GetGapPV() const  {
 return fGapPV;
 }
 */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

