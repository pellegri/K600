
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

#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Trd.hh"

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4Isotope.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4PhysicalConstants.hh"

#include "G4UniformMagField.hh"
#include "G4QuadrupoleMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"

#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4AutoDelete.hh"

#include "CADMesh.hh"
#include "MagneticFieldMapping.hh"
//#include "G4BlineTracer.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0;

G4ThreadLocal G4QuadrupoleMagField* DetectorConstruction::MagneticField_K600_Q = 0;
G4ThreadLocal G4UniformMagField* DetectorConstruction::MagneticField_K600_D1 = 0;
G4ThreadLocal G4UniformMagField* DetectorConstruction::MagneticField_K600_D2 = 0;

G4ThreadLocal G4FieldManager* DetectorConstruction::fieldManagerMagneticField_K600_Q = 0;
G4ThreadLocal G4FieldManager* DetectorConstruction::fieldManagerMagneticField_K600_D1 = 0;
G4ThreadLocal G4FieldManager* DetectorConstruction::fieldManagerMagneticField_K600_D2 = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
fAbsorberPV(0), fGapPV(0), fCheckOverlaps(false), PhysiCLOVER_HPGeCrystal(0), PhysiCLOVER_Shield_BGOCrystal(0), PhysiCLOVER_Shield_PMT(0), PhysiTIARA_AA_RS(0), PhysiPADDLE(0), PhysiK600_Quadrupole(0), PhysiK600_Dipole1(0), PhysiK600_Dipole2(0), PhysiHAGAR_NaICrystal(0), PhysiHAGAR_Annulus(0), PhysiHAGAR_FrontDisc(0), Physical_LEPS_HPGeCrystal(0)
{
    WorldSize = 15.*m;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    
    //////////////////////////////////////////////////////////
    ////                                                  ////
    ////            DETECTOR ARRAY POSITIONING            ////
    ////                                                  ////
    //////////////////////////////////////////////////////////
    
    
    ////////////////////////////
    ////    VDC SETUP
    VDC_AllPresent_Override = false;
    VDC_AllAbsent_Override = true;
    
    //  VDC 1
    VDC_Presence[0] = true;
    VDC_CentrePositionX[0] = 481.93; // cm
    VDC_CentrePositionZ[0] = 352.050; // cm
    VDC_CentrePosition[0] = G4ThreeVector(VDC_CentrePositionX[0]*cm - 200.*cm, 0., VDC_CentrePositionZ[0]*cm - 100.*cm);
    VDC_RotationY[0] = -14.03; // deg;
    VDC_rotm[0].rotateY(VDC_RotationY[0]*deg);
    
    //  VDC 2
    VDC_Presence[1] = false;
    VDC_CentrePositionX[1] = 531.23; // cm
    VDC_CentrePositionZ[1] = 337.40; // cm
    VDC_CentrePosition[1] = G4ThreeVector(VDC_CentrePositionX[1]*cm - 200.*cm, 0., VDC_CentrePositionZ[1]*cm - 100.*cm);
    VDC_RotationY[1] = -14.03; // deg;
    VDC_rotm[1].rotateY(VDC_RotationY[1]*deg);
    
    
    for (G4int i=0; i<numberOf_VDC; i++)
    {
        if(VDC_AllPresent_Override) VDC_Presence[i] = true;
        if(VDC_AllAbsent_Override) VDC_Presence[i] = false;
        if(VDC_AllPresent_Override && VDC_AllAbsent_Override) VDC_Presence[i] = false;
    }
    
    
    /////////////////////////////
    ////    PADDLE SETUP
    
    PADDLE_AllPresent_Override = false;
    PADDLE_AllAbsent_Override = true;
    
    //  PADDLE 1
    PADDLE_Presence[0] = true;
    PADDLE_CentrePositionX[0] = 513.726; // cm
    PADDLE_CentrePositionZ[0] = 324.573; // cm
    PADDLE_CentrePosition[0] = G4ThreeVector(PADDLE_CentrePositionX[0]*cm - 200.*cm, 0., PADDLE_CentrePositionZ[0]*cm - 100.*cm);
    PADDLE_RotationY[0] = -14.03; // deg
    PADDLE_rotm[0].rotateY(PADDLE_RotationY[0]*deg);
    
    //  PADDLE 2
    PADDLE_Presence[1] = false;
    PADDLE_CentrePositionX[1] = 529.38; // cm
    PADDLE_CentrePositionZ[1] = 320.34; // cm
    PADDLE_CentrePosition[1] = G4ThreeVector(PADDLE_CentrePositionX[1]*cm - 200.*cm, 0., PADDLE_CentrePositionZ[1]*cm - 100.*cm);
    PADDLE_RotationY[1] = -14.03; // deg
    PADDLE_rotm[1].rotateY(PADDLE_RotationY[1]*deg);
    
    //  PADDLE 3
    PADDLE_Presence[2] = true;
    PADDLE_CentrePositionX[2] = 529.38; // cm
    PADDLE_CentrePositionZ[2] = 320.34; // cm
    PADDLE_CentrePosition[2] = G4ThreeVector(PADDLE_CentrePositionX[2]*cm - 200.*cm, 0., PADDLE_CentrePositionZ[2]*cm - 100.*cm);
    PADDLE_RotationY[2] = -14.03; // deg
    PADDLE_rotm[2].rotateY(PADDLE_RotationY[2]*deg);
    
    
    for(G4int i=0; i<3; i++)
    {
        if(PADDLE_AllPresent_Override) PADDLE_Presence[i] = true;
        if(PADDLE_AllAbsent_Override) PADDLE_Presence[i] = false;
        if(PADDLE_AllPresent_Override && PADDLE_AllAbsent_Override) PADDLE_Presence[i] = false;
    }
    
    
    ////////////////////////////
    ////    TIARA SETUP
    
    TIARA_AllPresent_Override = false;
    TIARA_AllAbsent_Override = true;
    
    //offset_TIARA_BeamAxis = -131.0000; // mm
    offset_TIARA_BeamAxis = -131.3217600; // mm
    
    
    //  TIARA 1
    TIARA_Presence[0] = true;
    TIARA_AA_CentrePosition[0] = G4ThreeVector(1.9682273*mm, 0.1740583*mm, offset_TIARA_BeamAxis*mm);
    TIARA_rotm[0].rotateX(22.7683986*deg);
    TIARA_rotm[0].rotateY(-35.0661884*deg);
    TIARA_rotm[0].rotateZ((-54.1496802 + (0*72.))*deg);
    
    ////    This positioning if for viewing the TIARA
    //TIARA_rotm[0].rotateZ(59.63*deg);
    
    //  TIARA 2
    TIARA_Presence[1] = true;
    TIARA_AA_CentrePosition[1] = G4ThreeVector(0.4535687*mm, 1.9147014*mm, offset_TIARA_BeamAxis*mm);
    TIARA_rotm[1].rotateX(22.7683986*deg);
    TIARA_rotm[1].rotateY(-35.0661884*deg);
    TIARA_rotm[1].rotateZ((-54.1496802 + (1*72.))*deg);
    
    
    //  TIARA 3
    TIARA_Presence[2] = true;
    TIARA_AA_CentrePosition[2] = G4ThreeVector(-1.6699365*mm, 1.0120638*mm, offset_TIARA_BeamAxis*mm);
    TIARA_rotm[2].rotateX(22.7683986*deg);
    TIARA_rotm[2].rotateY(-35.0661884*deg);
    TIARA_rotm[2].rotateZ((-54.1496802 + (2*72.))*deg);
    
    
    //  TIARA 4
    TIARA_Presence[3] = true;
    TIARA_AA_CentrePosition[3] = G4ThreeVector(-1.4676763*mm, -1.2864401*mm, offset_TIARA_BeamAxis*mm);
    TIARA_rotm[3].rotateX(22.7683986*deg);
    TIARA_rotm[3].rotateY(-35.0661884*deg);
    TIARA_rotm[3].rotateZ((-54.1496802 + (3*72.))*deg);
    
    
    //  TIARA 5
    //  This TIARA was not present for PR226 and therefore the TIARA_AA_CentrePosition[4] was not yet measured from SolidEdge. Still to do for COMPLETENESS.
    //  In principle, one only has to rotate the TIARA_AA_CentrePosition[] for a particular MMM about the Z-axis for all 5 MMM's in the CAKE array configuration
    TIARA_Presence[4] = false;
    TIARA_AA_CentrePosition[4] = G4ThreeVector(0.90*mm, -1.50*mm, offset_TIARA_BeamAxis*mm);
    TIARA_rotm[4].rotateX(-37.891*deg);
    TIARA_rotm[4].rotateY(17.002*deg);
    TIARA_rotm[4].rotateZ((-38.571 + 4*72)*deg);
    
    for (G4int i=0; i<numberOf_TIARA; i++)
    {
        if(TIARA_AllPresent_Override) TIARA_Presence[i] = true;
        if(TIARA_AllAbsent_Override) TIARA_Presence[i] = false;
        if(TIARA_AllPresent_Override && TIARA_AllAbsent_Override) TIARA_Presence[i] = false;
    }
    
    ////////////////////////////
    ////    HAGAR SETUP
    
    //  HAGAR NaI Crystal
    HAGAR_NaICrystal_Presence = false;
    HAGAR_NaICrystal_CentrePosition = G4ThreeVector(-100.*cm, 0.*cm, 0.*cm);
    HAGAR_rotm.rotateY(-90*deg);
    
    //  HAGAR Annulus
    HAGAR_Annulus_Presence = false;
    HAGAR_Annulus_CentrePosition = G4ThreeVector(-105.*cm, 0.*cm, 0.*cm);
    
    //  HAGAR Front
    HAGAR_FrontDisc_Presence = false;
    HAGAR_FrontDisc_CentrePosition = HAGAR_Annulus_CentrePosition + G4ThreeVector((61./2 +8/2)*cm, 0.*cm, 0.*cm);
    
    
    ////////////////////////////
    ////    CLOVER SETUP
    
    CLOVER_AllPresent_Override = false;
    CLOVER_AllAbsent_Override = false;
    
    CLOVER_Shield_AllPresent_Override = false;
    CLOVER_Shield_AllAbsent_Override = true;
    
    
    //  CLOVER 1
    CLOVER_Presence[0] = true;
    CLOVER_Shield_Presence[0] = false;
    CLOVER_Distance[0] = 10*cm;
    CLOVER_phi[0] = 90*deg;
    CLOVER_theta[0] = 135*deg;
    CLOVER_rotm[0].rotateX(45. *deg);
    
    //  CLOVER 2
    CLOVER_Presence[1] = true;
    CLOVER_Shield_Presence[1] = false;
    CLOVER_Distance[1] = 10*cm;
    CLOVER_phi[1] = 0*deg;
    CLOVER_theta[1] = 135*deg;
    CLOVER_rotm[1].rotateY(-45.0 *deg);
    
    //  CLOVER 3
    CLOVER_Presence[2] = true;
    CLOVER_Shield_Presence[2] = false;
    CLOVER_Distance[2] = 10*cm;
    CLOVER_phi[2] = 270*deg;
    CLOVER_theta[2] = 135*deg;
    CLOVER_rotm[2].rotateX(-45.0 *deg);
    
    //  CLOVER 4
    CLOVER_Presence[3] = true;
    CLOVER_Shield_Presence[3] = false;
    CLOVER_Distance[3] = 10*cm;
    CLOVER_phi[3] = 180*deg;
    CLOVER_theta[3] = 135*deg;
    CLOVER_rotm[3].rotateY(45.0 *deg);
    
    //  CLOVER 5
    CLOVER_Presence[4] = true;
    CLOVER_Shield_Presence[4] = false;
    CLOVER_Distance[4] = 100*cm;
    CLOVER_phi[4] = 45*deg;
    CLOVER_theta[4] = 90*deg;
    CLOVER_rotm[4].rotateY(90.0 *deg);
    CLOVER_rotm[4].rotateZ(-135.0 *deg);
    
    //  CLOVER 6
    CLOVER_Presence[5] = true;
    CLOVER_Shield_Presence[5] = false;
    CLOVER_Distance[5] = 100*cm;
    CLOVER_phi[5] = 0*deg;
    CLOVER_theta[5] = 90*deg;
    CLOVER_rotm[5].rotateY(-90.0 *deg);
    
    //  CLOVER 7
    CLOVER_Presence[6] = true;
    CLOVER_Shield_Presence[6] = false;
    CLOVER_Distance[6] = 100*cm;
    CLOVER_phi[6] = 180*deg;
    CLOVER_theta[6] = 90*deg;
    CLOVER_rotm[6].rotateY(90.0 *deg);
    
    //  CLOVER 8
    CLOVER_Presence[7] = true;
    CLOVER_Shield_Presence[7] = false;
    CLOVER_Distance[7] = 100*cm;
    CLOVER_phi[7] = 135*deg;
    CLOVER_theta[7] = 90*deg;
    CLOVER_rotm[7].rotateY(90.0 *deg);
    CLOVER_rotm[7].rotateZ(-45.0 *deg);
    
    
    //  CLOVER 9
    CLOVER_Presence[8] = false;
    CLOVER_Shield_Presence[8] = false;
    CLOVER_Distance[8] = 10*cm;
    CLOVER_phi[8] = -45*deg;
    CLOVER_theta[8] = 90*deg;
    CLOVER_rotm[8].rotateY(-90.0 *deg);
    CLOVER_rotm[8].rotateZ(-45.0 *deg);
    
    
    for (G4int i=0; i<numberOf_CLOVER; i++)
    {
        if(CLOVER_AllPresent_Override) CLOVER_Presence[i] = true;
        if(CLOVER_AllAbsent_Override) CLOVER_Presence[i] = false;
        if(CLOVER_AllPresent_Override && CLOVER_AllAbsent_Override) CLOVER_Presence[i] = false;
        
        if(CLOVER_Shield_AllPresent_Override) CLOVER_Shield_Presence[i] = true;
        if(CLOVER_Shield_AllAbsent_Override) CLOVER_Shield_Presence[i] = false;
        if(CLOVER_Shield_AllPresent_Override && CLOVER_Shield_AllAbsent_Override) CLOVER_Shield_Presence[i] = false;
    }
    
    
    
    ////////////////////////////////
    ////        LEPS SETUP
    
    LEPS_AllPresent_Override = false;
    LEPS_AllAbsent_Override = true;
    
    
    //  LEPS 1
    LEPS_Presence[0] = false;
    LEPS_Distance[0] = 20*cm;
    LEPS_phi[0] = 90.*deg;
    LEPS_theta[0] = 135.*deg;
    LEPS_rotm[0].rotateX(45.*deg);
    
    //  LEPS 2
    LEPS_Presence[1] = false;
    LEPS_Distance[1] = 14.5*cm;
    LEPS_phi[1] = 0.*deg;
    LEPS_theta[1] = 135.*deg;
    LEPS_rotm[1].rotateY(-45.*deg);
    
    //  LEPS 3
    LEPS_Presence[2] = false;
    LEPS_Distance[2] = 20*cm;
    LEPS_phi[2] = 270.*deg;
    LEPS_theta[2] = 135.*deg;
    LEPS_rotm[2].rotateX(-45.*deg);
    
    //  LEPS 4
    LEPS_Presence[3] = false;
    LEPS_Distance[3] = 14.5*cm;
    LEPS_phi[3] = 180.*deg;
    LEPS_theta[3] = 135*deg;
    LEPS_rotm[3].rotateY(45.*deg);
    
    //  LEPS 5
    LEPS_Presence[4] = false;
    LEPS_Distance[4] = 14.5*cm;
    LEPS_phi[4] = 45*deg;
    LEPS_theta[4] = 90*deg;
    LEPS_rotm[4].rotateY(90.0 *deg);
    LEPS_rotm[4].rotateZ(-135.0 *deg);
    
    //  LEPS 6
    LEPS_Presence[5] = false;
    LEPS_Distance[5] = 14.5*cm;
    LEPS_phi[5] = 0*deg;
    LEPS_theta[5] = 90*deg;
    LEPS_rotm[5].rotateY(-90.*deg);
    
    //  LEPS 7
    LEPS_Presence[6] = true;
    LEPS_Distance[6] = 17*cm;
    LEPS_phi[6] = 180*deg;
    LEPS_theta[6] = 90*deg;
    LEPS_rotm[6].rotateY(90.0 *deg);
    
    //  LEPS 8
    LEPS_Presence[7] = true;
    LEPS_Distance[7] = 17*cm;
    LEPS_phi[7] = 135*deg;
    LEPS_theta[7] = 90*deg;
    LEPS_rotm[7].rotateY(90.0 *deg);
    LEPS_rotm[7].rotateZ(-45.0 *deg);
    
    for (G4int i=0; i<numberOf_LEPS; i++)
    {
        if( LEPS_AllPresent_Override == true ) LEPS_Presence[i] = true;
        if( LEPS_AllAbsent_Override == true ) LEPS_Presence[i] = false;
        if( LEPS_AllPresent_Override == true && LEPS_AllAbsent_Override == true ) LEPS_Presence[i] = false;
    }

    ////////////////////////////////
    ////        FLATSIDES CHAMBER SETUP
    
    
    FLATSIDE_AllPresent_Override = false;
    FLATSIDE_AllAbsent_Override = true;
    
    //  FLATSIDE 1
    FLATSIDE_Presence[0] = true;
    FLATSIDE_Distance[0] = 7*cm;
    FLATSIDE_phi[0] = 180*deg;
    FLATSIDE_theta[0] = 90*deg;
    FLATSIDE_rotm[0].rotateY(85.0 *deg);
    
    //  FLATSIDE 2
    FLATSIDE_Presence[1] = true;
    FLATSIDE_Distance[1] = 7*cm;
    FLATSIDE_phi[1] = 0*deg;
    FLATSIDE_theta[1] = 90*deg;
    FLATSIDE_rotm[1].rotateY(-85.*deg);
    
    
    for (G4int i=0; i<numberOf_FLATSIDE; i++)
    {
        if( FLATSIDE_AllPresent_Override == true ) FLATSIDE_Presence[i] = true;
        if( FLATSIDE_AllAbsent_Override == true ) FLATSIDE_Presence[i] = false;
        if( FLATSIDE_AllPresent_Override == true && FLATSIDE_AllAbsent_Override == true ) FLATSIDE_Presence[i] = false;
    }
    
    ////////////////////////////////
    ////        AmBe Source SETUP
    
    
   /* PARAFFINTUBE_Present_Override = true;
    PARAFFINTUBE_Absent_Override = false;


        if( PARAFFINTUBE_Present_Override == true ) PARAFFINTUBE_Presence = true;
        if( PARAFFINTUBE_Absent_Override == true ) PARAFFINTUBE_Presence = false;
        if( PARAFFINTUBE_Present_Override == true && PARAFFINTUBE_Absent_Override == true ) PARAFFINTUBE_Presence = false;
  */
    
        PARAFFINTUBE_Presence = true;
        IRONBOX_Presence = true;
  
    
    
    
    ////////////////////////////////////////////
    ////                                    ////
    ////        K600 SPECTROMETER SETUP     ////
    ////                                    ////
    ////////////////////////////////////////////
    
    //  K600 Quadrapole
    K600_Quadrupole = false;
    Ideal_Quadrupole = true;
    Mapped_Quadrupole = false;
    //K600_Q_gradient = 0.030*tesla/cm;  // gradient = dB/dr for Ideal Quadrupole
    K600_Q_gradient = 0.001*tesla/cm;  // gradient = dB/dr for Ideal Quadrupole
    K600_Quadrupole_CentrePosition = G4ThreeVector(0.*cm, 0.*cm,100.*cm);
    //K600_Quadrupole_rotm.rotateZ(90.*deg);
    //K600_Quadrupole_rotm.rotateY(90*deg);
    //K600_Quadrupole_rotm2.rotateX(90*deg);
    
    //  K600 Dipole 1
    K600_Dipole1 = false;
    K600_Dipole1_BZ = -2.30*tesla;
    //K600_Dipole1_BZ = -3.30*tesla;
    K600_Dipole1_CentrePosition = G4ThreeVector(75*cm, 0.*cm,250*cm);
    K600_Dipole1_rotm.rotateX(-90.*deg);
    K600_Dipole1_rotm.rotateY(180.*deg);
    
    //  K600 Dipole 2
    K600_Dipole2 = false;
    K600_Dipole2_BZ = -2.40*tesla;
    //K600_Dipole2_BZ = -3.40*tesla;
    K600_Dipole2_CentrePosition = G4ThreeVector(75*cm, 0.*cm,250*cm);
    K600_Dipole2_rotm.rotateX(-90*deg);
    K600_Dipole2_rotm.rotateY(180.*deg);
    
    if(Ideal_Quadrupole && Mapped_Quadrupole || (!Ideal_Quadrupole && !Mapped_Quadrupole))
    {
        Ideal_Quadrupole = false;
        Mapped_Quadrupole = false;
    }
    
    ////////////////////////////////////////////////////
    ////                                            ////
    ////        K600 VAULT - STRUCTURE SETUP        ////
    ////                                            ////
    ////////////////////////////////////////////////////
    
    
    ////////////////////////////////////////////////
    ////    New K600 Target Chamber - New scattering chamber, both sides on
    K600_BACTAR_sidesOn_Presence = false;
    
    ////////////////////////////////////////////////
    ////    New K600 Target Chamber - New scattering chamber, both sides off
    K600_BACTAR_sidesOff_Presence = true;
    
    ////////////////////////////////////////////////
    ////    New K600 Target Chamber - New scattering chamber, beam left side side off
    K600_BACTAR_beamRightSideOff_Presence = false;
    
    ////////////////////////////////////////////////
    ////    New K600 Target Chamber - New scattering chamber, beam right side side off
    K600_BACTAR_beamLeftSideOff_Presence = false;

    
    
    /////////////////////////////////////
    //  K600 Target
    K600_Target_Presence = false;
    
    /////////////////////////////////////
    //  K600 Target Backing
    K600_TargetBacking_Presence = false;
    
    
    // Define materials
    DefineMaterials();
    
    // Define volumes
    return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
    G4NistManager* nistManager = G4NistManager::Instance();
    
    //  NIST Material Database - Materials
    nistManager->FindOrBuildMaterial("G4_Galactic");
    nistManager->FindOrBuildMaterial("G4_AIR");
    nistManager->FindOrBuildMaterial("G4_BGO");
    nistManager->FindOrBuildMaterial("G4_Ge");
    nistManager->FindOrBuildMaterial("G4_Al");
    nistManager->FindOrBuildMaterial("G4_Si");
    nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    nistManager->FindOrBuildMaterial("G4_MYLAR");
    nistManager->FindOrBuildMaterial("G4_W");
    nistManager->FindOrBuildMaterial("G4_Ar");
    nistManager->FindOrBuildMaterial("G4_C");
    nistManager->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
    nistManager->FindOrBuildMaterial("G4_SODIUM_IODIDE");
    nistManager->FindOrBuildMaterial("G4_LITHIUM_CARBONATE");
    nistManager->FindOrBuildMaterial("G4_Be");
    nistManager->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
    nistManager->FindOrBuildMaterial("G4_Fe");
    nistManager->FindOrBuildMaterial("G4_PARAFFIN");

    //  NIST Elementary Material Database - ELEMENTS
    nistManager->FindOrBuildElement("H");
    nistManager->FindOrBuildElement("C");
    nistManager->FindOrBuildElement("N");
    nistManager->FindOrBuildElement("O");
    nistManager->FindOrBuildElement("Fe");
    nistManager->FindOrBuildElement("Co");
    nistManager->FindOrBuildElement("Ni");
    nistManager->FindOrBuildElement("Cu");
    nistManager->FindOrBuildElement("Pb");
    nistManager->FindOrBuildElement("W");
    nistManager->FindOrBuildElement("Li");
    nistManager->FindOrBuildElement("Ar");
    
    
    
    // Print materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
    //////////////////////////////////////
    //          Get Elements            //
    //////////////////////////////////////
    
    G4Element* Hydrogen = G4Element::GetElement("H");
    G4Element* Carbon = G4Element::GetElement("C");
    G4Element* Nitrogen = G4Element::GetElement("N");
    G4Element* Oxygen = G4Element::GetElement("O");
    G4Element* Iron = G4Element::GetElement("Fe");
    G4Element* Cobalt = G4Element::GetElement("Co");
    G4Element* Nickel = G4Element::GetElement("Ni");
    G4Element* Copper = G4Element::GetElement("Cu");
    G4Element* Lead = G4Element::GetElement("Pb");
    G4Element* Tungsten = G4Element::GetElement("W");
    G4Element* Lithium = G4Element::GetElement("Li");
    G4Element* Argon = G4Element::GetElement("Ar");
    G4Element* Americium = G4Element::GetElement("Be");
    G4Element* Beryllium = G4Element::GetElement("Be");
    
    //////////////////////////////////////
    //          Get Materials           //
    //////////////////////////////////////
    
    ////    NIST Defined Elemental Material
    G4Material* G4_Ge_Material = G4Material::GetMaterial("G4_Ge");
    G4Material* G4_Al_Material  = G4Material::GetMaterial("G4_Al");
    G4Material* G4_Si_Material = G4Material::GetMaterial("G4_Si");
    G4Material* G4_W_Material = G4Material::GetMaterial("G4_W");
    G4Material* G4_Ar_Material = G4Material::GetMaterial("G4_Ar");
    G4Material* G4_C_Material = G4Material::GetMaterial("G4_C");
    G4Material* G4_Fe_Material = G4Material::GetMaterial("G4_Fe");
    
    ////    NIST Defined Materials and Compounds
    G4Material* G4_Galactic_Material = G4Material::GetMaterial("G4_Galactic");
    G4Material* G4_AIR_Material = G4Material::GetMaterial("G4_AIR");
    G4Material* G4_BGO_Material = G4Material::GetMaterial("G4_BGO");
    G4Material* G4_PLASTIC_SC_VINYLTOLUENE_Material = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    G4Material* G4_Mylar_Material = G4Material::GetMaterial("G4_MYLAR");
    G4Material* G4_CARBON_DIOXIDE_Material = G4Material::GetMaterial("G4_CARBON_DIOXIDE");
    G4Material* G4_SODIUM_IODIDE_Material = G4Material::GetMaterial("G4_SODIUM_IODIDE");
    G4Material* G4_LITHIUM_CARBONATE_Material = G4Material::GetMaterial("G4_LITHIUM_CARBONATE");
    G4Material* G4_Be_Material = G4Material::GetMaterial("G4_Be");
    G4Material* G4_POLYVINYL_CHLORIDE_Material = G4Material::GetMaterial("G4_POLYVINYL_CHLORIDE");
    G4Material* G4_PARAFFIN_Material = G4Material::GetMaterial("G4_PARAFFIN");

    ////    CLOVER Detector Shield, HEAVIMET Material
    G4Material* Heavimet_Material = new G4Material("Heavimet_Material",19.25*g/cm3, 5);
    Heavimet_Material->AddElement(Tungsten, 94.20*perCent);
    Heavimet_Material->AddElement(Nickel, 4.35*perCent);
    Heavimet_Material->AddElement(Iron, 0.85*perCent);
    Heavimet_Material->AddElement(Cobalt, 0.50*perCent);
    Heavimet_Material->AddElement(Copper, 0.10*perCent);
    
    ////    TIARA Detector Materials
    //G4Material*         TIARA_PCBMaterial;
    
    ////    PADDLE Detector Materials
    G4Material* BC408_Material = new G4Material("BC408_Material",1.032*g/cm3, 2);
    BC408_Material->AddElement(Hydrogen, 8.4748*perCent);
    BC408_Material->AddElement(Carbon, 91.5252*perCent);
    
    ////    VDC Detector Materials
    G4Material* VDC_SR_Gas_Material = new G4Material("VDC_SR_Gas_Material", 0.00184212*g/cm3, 2);
    VDC_SR_Gas_Material->AddMaterial(G4_Ar_Material, 90.*perCent);
    VDC_SR_Gas_Material->AddMaterial(G4_CARBON_DIOXIDE_Material, 10.*perCent);
    
    
    ////    Vacuum - Target Chamber
    G4Material* Vacuum_TC_Material = new G4Material("Vacuum_TC", 1.2e-14*g/cm3, 4);
    Vacuum_TC_Material->AddElement(Carbon, 0.0124*perCent);
    Vacuum_TC_Material->AddElement(Nitrogen, 75.5268*perCent);
    Vacuum_TC_Material->AddElement(Oxygen, 23.1781*perCent);
    Vacuum_TC_Material->AddElement(Argon, 1.2827*perCent);

    ////    AmBeSource_Material
    G4Material* AmBeSource_Material = new G4Material("AmBeSource", 1.0*g/cm3, 2);
    AmBeSource_Material->AddElement(Beryllium, (9.0/250.0)*100*perCent);
    AmBeSource_Material->AddElement(Americium, (241.0/250.0)*100*perCent);

    
    
    
    /*
     4     G4_AIR                     0.00120479        85.7
     6     0.000124
     7     0.755268
     8     0.231781
     18    0.012827
     */
    
    
    ////    Variables for CADMesh
    char    meshPath[512];
    char    meshType[] = "PLY";

    
    //////////////////////////////////////////////////////////
    //                      WORLD                           //
    //////////////////////////////////////////////////////////
    
    G4Box* SolidWorld = new G4Box("World", WorldSize/2, WorldSize/2, WorldSize/2);
    
    G4LogicalVolume*
    LogicWorld = new G4LogicalVolume(SolidWorld,                //its solid
                                     G4_Galactic_Material,      //its material
                                     "World");                  //its name
    G4VPhysicalVolume*
    PhysiWorld = new G4PVPlacement(0,                        //no rotation
                                   G4ThreeVector(),          //at (0,0,0)
                                   LogicWorld,               //its logical volume
                                   "World",                  //its name
                                   0,                        //its mother  volume
                                   false,                    //no boolean operation
                                   0);                       //copy number
    
    
    
    //////////////////////////////////////////////////////////
    //                      VACUUM CHAMBER
    //////////////////////////////////////////////////////////
    
    G4ThreeVector positionVacuumChamber = G4ThreeVector(0,0,0);
    
    G4Box* SolidVacuumChamber = new G4Box("VacuumChamber", (400./2)*cm, (400./2)*cm, (400./2)*cm);
    
    G4LogicalVolume* LogicVacuumChamber = new G4LogicalVolume(SolidVacuumChamber, G4_Galactic_Material,"VacuumChamber",0,0,0);
    
    
    new G4PVPlacement(0,               // no rotation
                      positionVacuumChamber, // at (x,y,z)
                      LogicVacuumChamber,       // its logical volume
                      "VacuumChamber",       // its name
                      LogicWorld,         // its mother  volume
                      false,           // no boolean operations
                      0,               // copy number
                      fCheckOverlaps); // checking overlaps
    
    
   
    
    ////////////////////////////////////////
    ////        FLATSIDE DEFINITION         ////
    ////////////////////////////////////////
    
   // G4double FLATSIDE_thickness = 20.*mm;
    G4double FLATSIDE_thickness = 3.*mm;
    
    G4Box* Solid_FLATSIDE = new G4Box("FlatSide",                             //its name
                               15.*cm, 15.*cm, 0.5*FLATSIDE_thickness);     //its size

    G4LogicalVolume* Logic_FLATSIDE[numberOf_FLATSIDE];
    
    for(G4int i=0; i<numberOf_FLATSIDE; i++)
    {
        Logic_FLATSIDE[i] = new G4LogicalVolume(Solid_FLATSIDE,          //its solid
                                       G4_POLYVINYL_CHLORIDE_Material,           //its material
                                                "LogicFlatSide",0, 0, 0);            //its name
    }
    
    
    
    ////////////////////////////////////////////////////
    //               FLATSIDE INITIALIZATION
    ////////////////////////////////////////////////////
    
    
    for(G4int i=0; i<numberOf_FLATSIDE; i++)
    {
        FLATSIDE_position[i] = (FLATSIDE_Distance[i])*G4ThreeVector( std::sin(FLATSIDE_theta[i]) * std::cos(FLATSIDE_phi[i]), std::sin(FLATSIDE_theta[i]) * std::sin(FLATSIDE_phi[i]), std::cos(FLATSIDE_theta[i]));
        
        FLATSIDE_transform[i] = G4Transform3D(FLATSIDE_rotm[i],FLATSIDE_position[i]);
        
    
        
        /////////////////////////////
        //          FLATSIDE
        if(FLATSIDE_Presence[i] == true)
        {
            
            new G4PVPlacement(FLATSIDE_transform[i],   // transformation matrix
                              Logic_FLATSIDE[i],       // its logical volume
                              "FLATSIDE",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            
        }
        
    }

    
    ////////////////////////////////////////
    ////     AmBe SOURCE                ////
    ////////////////////////////////////////
    
    
   ////************ PARAFFIN TUBE

    G4ThreeVector positionParaffinTube = G4ThreeVector(0.,0.,0.);
    
    G4double ParaffinTube_innerRadius = 0.*cm;
    G4double ParaffinTube_outerRadius = 10.*cm;
    G4double ParaffinTube_hight = 10.*cm;
    G4double startAngleOfTheTube = 0.*deg;
    G4double spanningAngleOfTheTube = 360.*deg;
    
    G4Tubs* Solid_PARAFFINTUBE    = new G4Tubs("ParaffinTube",
                 ParaffinTube_innerRadius,
                 ParaffinTube_outerRadius,
                 ParaffinTube_hight,
                 startAngleOfTheTube,
                 spanningAngleOfTheTube);
    
    
    
    G4double HoleBOX_thickness = 7.*cm;
    G4double HoleBOX_width = 7.*cm;
    G4double HoleBOX_hight = 14*cm;
    
    G4Box* Solid_HoleBOX = new G4Box("HoleBOX",                             //its name
                                     0.5*HoleBOX_width, 0.5*HoleBOX_thickness,0.5*HoleBOX_hight);     //its size
   
    G4ThreeVector translationhole(0.,0.,3*cm);

    G4SubtractionSolid* ParaffinBoxwithHole = new G4SubtractionSolid("ParaffinBoxwithHole",
                                                                     Solid_PARAFFINTUBE,Solid_HoleBOX,0,translationhole);
    
    
    G4LogicalVolume* Logic_PARAFFINBOX = new G4LogicalVolume(ParaffinBoxwithHole,          //its solid ParaffinBoxwithHole
                                                             G4_PARAFFIN_Material,           //it material (G4_PARAFFIN)
                                                             "ParaffinTube",0, 0, 0);            //its name
   
    
    G4RotationMatrix rotm  = G4RotationMatrix();
    G4double theta=-90*deg;
    rotm.rotateX(theta);
    G4Transform3D transform_ParaffinBox = G4Transform3D(rotm,positionParaffinTube);

    
    
    
    if(PARAFFINTUBE_Presence == true)
    {
        new G4PVPlacement(transform_ParaffinBox, // at (x,y,z)
                      Logic_PARAFFINBOX,       // its logical volume
                      "ParaffinBox",       // its name
                      LogicVacuumChamber,         // its mother  volume
                      false,           // no boolean operations
                      0,               // copy number
                      fCheckOverlaps); // checking overlaps

    
    }
    
    /*
    
    if(PARAFFINTUBE_Presence == true)
    {
        new G4PVPlacement(0,
                          positionParaffinTube, // at (x,y,z)
                          Logic_PARAFFINBOX,       // its logical volume
                          "ParaffinBox",       // its name
                          LogicWorld,         // its mother  volume
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
        
        
    }
*/
    
    ////************ IRON BOX


    
    
    G4double IronBox_thickness = 7.*cm;
    G4double IronBox_width = 7.*cm;
    G4double IronBox_hight = 14*cm;
    
    G4Box* Solid_IronBox = new G4Box("IronBox",                             //its name
                                     0.5*IronBox_width,0.5*IronBox_hight, 0.5*IronBox_thickness);     //its size
    
    G4double IronHole_thickness = 3.*cm;
    G4double IronHole_width = 3.*cm;
    G4double IronHole_hight = 12*cm;
    
    G4Box* Solid_IronHole = new G4Box("IronHole",                             //its name
                                     0.5*IronHole_width,0.5*IronHole_hight, 0.5*IronHole_thickness);     //its size
    
    G4ThreeVector translationIronhole(0.,1.0*cm,0.);
    
    G4SubtractionSolid* IronBoxwithHole = new G4SubtractionSolid("IronBoxwithHole",
                                                                     Solid_IronBox,Solid_IronHole,0,translationIronhole);
    
    
    G4LogicalVolume* Logic_IRONBOX = new G4LogicalVolume(IronBoxwithHole,          //its solid
                                                             G4_Fe_Material,           //it material
                                                             "Logic_IronBox",0, 0, 0);            //its name

    
    
    G4ThreeVector positionIronBox = G4ThreeVector(0.,3.0*cm,0.);

    
    if(IRONBOX_Presence == true)
    {
        new G4PVPlacement(0,
                          positionIronBox, // at (x,y,z)
                          Logic_IRONBOX,       // its logical volume
                          "IronBox",       // its name
                          LogicVacuumChamber,         // its mother  volume
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
        
        
    }
    
    
    
    ////************ AmBe Source
    
    
    
    G4ThreeVector positionAmBeSource = G4ThreeVector(0.,0.,0.);
    
    G4Box* Solid_AmBeSource = new G4Box("AmBeSource", 0.01*cm, 0.01*cm, 0.01*cm);     //its size

    G4LogicalVolume* Logic_AmBeSource = new G4LogicalVolume(Solid_AmBeSource,          //its solid
                                                         AmBeSource_Material,           //it material
                                                         "Logic_AmBeSource",0, 0, 0);            //its name

    G4bool AmBeSource_Presence = true;
    
    if(AmBeSource_Presence == true)
    {
        new G4PVPlacement(0,
                          positionIronBox, // at (x,y,z)
                          Logic_AmBeSource,       // its logical volume
                          "AmBeSource",       // its name
                          LogicVacuumChamber,         // its mother  volume
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
        
        
    }

    
    
    //////***************//// PVC WALL = FLAT SIDES OF THE CHAMBER
  /*
    G4double Wthickness = 20.*mm;
    G4double Wall_Frontface_Distance = 50.*mm; // distance of the wall frontface from the source
    
    
    //////////////////// ROTATE THE DETECTOR and Endcup
    G4RotationMatrix rotm  = G4RotationMatrix();
    G4double theta=45*deg;
    rotm.rotateY(theta);
    G4double thetaTrans=theta+90*deg;
    G4ThreeVector uz = G4ThreeVector(std::sin(thetaTrans), 0.,  std::cos(thetaTrans));
    G4ThreeVector Wall_position = Wall_Frontface_Distance*uz;
    G4Transform3D transform_Wall = G4Transform3D(rotm,Wall_position);

    
    /////////////////////////////////
    
    
    
    
    G4ThreeVector posWall = G4ThreeVector(-(Wall_Frontface_Distance + 0.5*Wthickness),0., 0.);

    
    G4Box* solidWall =
    new G4Box("Wall",                       //its name
              30.*cm, 30.*cm, 0.5*Wthickness);     //its size
    
    G4LogicalVolume* logicWall =
    new G4LogicalVolume(solidWall,          //its solid
                        G4_POLYVINYL_CHLORIDE_Material,           //its material
                        "Wall");            //its name
    
    G4VPhysicalVolume* physWall =
    new G4PVPlacement(transform_Wall,       //at (0,0,0)
                      "Wall",               //its name
                      logicWall,            //its logical volume
                      PhysiWorld,             //its mother  volume
                      false,                 //no boolean operation
                      0);
    //                      0,                     //copy number
    //                     checkOverlaps);        //overlaps checking
    
    
    
   */
    
    
    /*********************************************************///////////////
    
    
    
    //////////////////////////////////////////////////////////
    //              Scattering Chamber - CADMesh
    //////////////////////////////////////////////////////////
    
    if(K600_BACTAR_sidesOn_Presence)
    {
        sprintf(meshPath, "../K600/Mesh-Models/STRUCTURES/BACTAR/BACTAR_sidesOn.ply");
    }
    if(K600_BACTAR_sidesOff_Presence)
    {
        sprintf(meshPath, "../K600/Mesh-Models/STRUCTURES/BACTAR/BACTAR_sidesOff.ply");
    }
    if(K600_BACTAR_beamRightSideOff_Presence)
    {
        sprintf(meshPath, "../K600/Mesh-Models/STRUCTURES/BACTAR/BACTAR_beamRightSideOff.ply");
    }
    if(K600_BACTAR_beamLeftSideOff_Presence)
    {
        sprintf(meshPath, "../K600/Mesh-Models/STRUCTURES/BACTAR/BACTAR_beamLightSideOff.ply");
    }

    
    if(K600_BACTAR_sidesOn_Presence || K600_BACTAR_sidesOff_Presence || K600_BACTAR_beamRightSideOff_Presence || K600_BACTAR_beamLeftSideOff_Presence)
    {
        G4ThreeVector offset_BACTAR = G4ThreeVector(0*cm, 0*cm, 0*cm);
        
        CADMesh * mesh_BACTAR = new CADMesh(meshPath, meshType, mm, offset_BACTAR, false);
        
        G4VSolid * SolidBACTAR = mesh_BACTAR->TessellatedMesh();
        
        G4LogicalVolume* LogicBACTAR = new G4LogicalVolume(SolidBACTAR, G4_Al_Material, "BACTAR", 0, 0, 0);
        
        new G4PVPlacement(0,               // no rotation
                          G4ThreeVector(), // at (x,y,z)
                          LogicBACTAR,       // its logical volume
                          "BACTAR",       // its name
                          LogicVacuumChamber,         // its mother  volume
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
        
        
        //  Visualisation
        G4VisAttributes* K600_BACTAR_TC_VisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
        K600_BACTAR_TC_VisAtt->SetForceSolid(true);
        LogicBACTAR ->SetVisAttributes(K600_BACTAR_TC_VisAtt);
    }
    
    
    //////////////////////////////////////////////////////////
    //                  TARGET DEFINITION
    //////////////////////////////////////////////////////////
    
    G4double targetThickness = 2.42; // um
    //G4double targetThickness = 2.18; // um
    //G4double targetThickness = 3.5; // um
    //G4double targetThickness = 4.5; // um
    
    G4ThreeVector positionTarget = G4ThreeVector(0,0,0);
    
    G4Box* SolidTarget = new G4Box("Target", (20./2)*mm, (20./2)*mm, (targetThickness/2)*um);
    
    G4LogicalVolume* LogicTarget = new G4LogicalVolume(SolidTarget, G4_LITHIUM_CARBONATE_Material,"Target",0,0,0);
    
    if(K600_Target_Presence)
    {
        new G4PVPlacement(0,               // no rotation
                          positionTarget, // at (x,y,z)
                          LogicTarget,       // its logical volume
                          "Target",       // its name
                          LogicVacuumChamber,         // its mother  volume
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
    }
    
    
    //////////////////////////////////////////////////////
    //                  TARGET BACKING
    //////////////////////////////////////////////////////
    
    G4double targetBackingThickness = 0.25; // um
    
    G4ThreeVector positionTargetBacking = G4ThreeVector(0,0,-((targetThickness/2) + (targetBackingThickness/2))*um);
    
    G4Box* SolidTargetBacking = new G4Box("TargetBacking", (20./2)*mm, (20./2)*mm, (targetBackingThickness/2)*um);
    
    G4LogicalVolume* LogicTargetBacking = new G4LogicalVolume(SolidTargetBacking, G4_C_Material,"TargetBacking",0,0,0);
    
    if(K600_TargetBacking_Presence)
    {
        new G4PVPlacement(0,               // no rotation
                          positionTargetBacking, // at (x,y,z)
                          LogicTargetBacking,       // its logical volume
                          "TargetBacking",       // its name
                          LogicVacuumChamber,         // its mother  volume
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
    }
    
    
    
    
    /////////////////////////////////////////////
    //             TIARA DEFINITION            //
    /////////////////////////////////////////////
    
    
    /////////////////////////////////////////
    //          TIARA - PCB
    
    G4double thickness_PCB = 1.6e+3; // um
    
    G4Tubs* Solid_TIARA_PCB_filled = new G4Tubs("TIARA_PCB", 0*mm, 170.*mm, (thickness_PCB/2)*um, 0*deg, (57.74)*deg);
    
    G4Box* PCB_razor1 = new G4Box("PCB_razor1", (200/2)*mm, (20.)*mm, (100/2)*mm);
    
    G4ThreeVector   position_TIARA_PCB_razor1[4];
    position_TIARA_PCB_razor1[0] = G4ThreeVector(7.0236770*mm, 2.0832170*mm, 0.*mm);
    position_TIARA_PCB_razor1[1] = G4ThreeVector(171.1783094*mm, 25.6019696*mm, 0.*mm);
    position_TIARA_PCB_razor1[2] = G4ThreeVector(154.2410751*mm, 99.9806850*mm, 0*mm);
    position_TIARA_PCB_razor1[3] = G4ThreeVector(107.0949198*mm, 135.0207195*mm, 0*mm);
    
    G4RotationMatrix* rm1_TIARA_PCB = new G4RotationMatrix();
    G4RotationMatrix* rm2_TIARA_PCB = new G4RotationMatrix();
    G4RotationMatrix* rm3_TIARA_PCB = new G4RotationMatrix();
    G4RotationMatrix* rm4_TIARA_PCB = new G4RotationMatrix();
    
    rm1_TIARA_PCB->rotateZ(59.63*deg);
    rm2_TIARA_PCB->rotateZ(101.2030637*deg);
    rm3_TIARA_PCB->rotateZ(59.63*deg);
    rm4_TIARA_PCB->rotateZ(15.8652730*deg);
    
    G4VSolid* Solid_TIARA_PCB_shear1 = new G4SubtractionSolid("TIARA_PCB_shear1", Solid_TIARA_PCB_filled, PCB_razor1, rm1_TIARA_PCB, position_TIARA_PCB_razor1[0]);
    G4VSolid* Solid_TIARA_PCB_shear2 = new G4SubtractionSolid("TIARA_PCB_shear2", Solid_TIARA_PCB_shear1, PCB_razor1, rm2_TIARA_PCB, position_TIARA_PCB_razor1[1]);
    G4VSolid* Solid_TIARA_PCB_shear3 = new G4SubtractionSolid("TIARA_PCB_shear3", Solid_TIARA_PCB_shear2, PCB_razor1, rm3_TIARA_PCB, position_TIARA_PCB_razor1[2]);
    G4VSolid* Solid_TIARA_PCB_shear4 = new G4SubtractionSolid("TIARA_PCB_shear4", Solid_TIARA_PCB_shear3, PCB_razor1, rm4_TIARA_PCB, position_TIARA_PCB_razor1[3]);
    
    /////////////////////////////////////////
    //          TIARA - PCB PUNCH
    
    G4double thickness_PCBpunch = 2.0e+3; // um
    
    G4Tubs* Solid_TIARA_PCBpunch_filled = new G4Tubs("TIARA_PCBpunch_filled", 0*mm, 160.*mm, (thickness_PCBpunch/2)*um, 0.*deg, (60.)*deg);
    
    G4Box* PCBpunch_shear1 = new G4Box("PCBpunch_shear1", (300/2)*mm, (20.)*mm, (100/2)*mm);
    
    G4ThreeVector   position_TIARA_PCBpunch_shear[6];
    position_TIARA_PCBpunch_shear[0] = G4ThreeVector(7.6909778*mm, 5.8892601*mm, 0.*mm);
    position_TIARA_PCBpunch_shear[1] = G4ThreeVector(82.8083631*mm, -15.9622884*mm, 0.*mm);
    position_TIARA_PCBpunch_shear[2] = G4ThreeVector(152.9061303*mm, 34.2318836*mm, 0.*mm);
    position_TIARA_PCBpunch_shear[3] = G4ThreeVector(133.6409271*mm, 81.5333641*mm, 0.*mm);
    position_TIARA_PCBpunch_shear[4] = G4ThreeVector(105.9034003*mm, 114.8979515*mm, 0.*mm);
    position_TIARA_PCBpunch_shear[5] = G4ThreeVector(26.8609291*mm, 78.1519180*mm, 0.*mm);
    
    G4RotationMatrix* rm_PCBpunch_shear1 = new G4RotationMatrix();
    G4RotationMatrix* rm_PCBpunch_shear2 = new G4RotationMatrix();
    G4RotationMatrix* rm_PCBpunch_shear3 = new G4RotationMatrix();
    G4RotationMatrix* rm_PCBpunch_shear4 = new G4RotationMatrix();
    G4RotationMatrix* rm_PCBpunch_shear5 = new G4RotationMatrix();
    G4RotationMatrix* rm_PCBpunch_shear6 = new G4RotationMatrix();
    
    rm_PCBpunch_shear1->rotateZ(59.6297466*deg);
    rm_PCBpunch_shear2->rotateZ(-3.0000000*deg);
    rm_PCBpunch_shear3->rotateZ(76.7362500*deg);
    rm_PCBpunch_shear4->rotateZ(59.6300000*deg);
    rm_PCBpunch_shear5->rotateZ(42.5237500*deg);
    rm_PCBpunch_shear6->rotateZ(-57.7400000*deg);
    
    G4VSolid* Solid_TIARA_PCB_punch1 = new G4SubtractionSolid("TIARA_PCB_punch1", Solid_TIARA_PCBpunch_filled, PCBpunch_shear1, rm_PCBpunch_shear1, position_TIARA_PCBpunch_shear[0]);
    G4VSolid* Solid_TIARA_PCB_punch2 = new G4SubtractionSolid("TIARA_PCB_punch2", Solid_TIARA_PCB_punch1, PCBpunch_shear1, rm_PCBpunch_shear2, position_TIARA_PCBpunch_shear[1]);
    G4VSolid* Solid_TIARA_PCB_punch3 = new G4SubtractionSolid("TIARA_PCB_punch3", Solid_TIARA_PCB_punch2, PCBpunch_shear1, rm_PCBpunch_shear3, position_TIARA_PCBpunch_shear[2]);
    G4VSolid* Solid_TIARA_PCB_punch4 = new G4SubtractionSolid("TIARA_PCB_punch4", Solid_TIARA_PCB_punch3, PCBpunch_shear1, rm_PCBpunch_shear4, position_TIARA_PCBpunch_shear[3]);
    G4VSolid* Solid_TIARA_PCB_punch5 = new G4SubtractionSolid("TIARA_PCB_punch5", Solid_TIARA_PCB_punch4, PCBpunch_shear1, rm_PCBpunch_shear5, position_TIARA_PCBpunch_shear[4]);
    G4VSolid* Solid_TIARA_PCB_punch6 = new G4SubtractionSolid("TIARA_PCB_punch6", Solid_TIARA_PCB_punch5, PCBpunch_shear1, rm_PCBpunch_shear6, position_TIARA_PCBpunch_shear[5]);
    
    ////    This offset is the offset between the origins of the created PCB and the PCBpunch
    //G4ThreeVector position_TIARA_PCB_perf = G4ThreeVector(1.0437151*mm, -0.2199451*mm, 0.*mm);
    G4ThreeVector position_TIARA_PCB_perf = G4ThreeVector(0.*mm, 0.*mm, 0.*mm);
    
    G4RotationMatrix* rm_perf = new G4RotationMatrix();
    //rm_perf->rotateZ(-3.*deg);
    
    G4VSolid* Solid_TIARA_PCB_punched = new G4SubtractionSolid("TIARA_PCB_punched", Solid_TIARA_PCB_shear4, Solid_TIARA_PCB_punch6, rm_perf, position_TIARA_PCB_perf);
    
    
    /////////////////////////////////////////
    //          TIARA - Silicon Wafer
    
    G4double thickness_SiliconWafer = 400.0; // um
    
    G4double thickness_DL_Front = 2.0; // um
    G4double thickness_DL_Back = 2.0; // um
    
    G4double thickness_2M_Front = 0.3; // um
    G4double thickness_2M_Back = 0.3; // um
    
    ////    Active Area (Silicon) thickness
    G4double thickness_AA = thickness_SiliconWafer - thickness_DL_Front - thickness_DL_Back; // um
    
    G4Tubs* Solid_TIARA_SiliconWafer_filled = new G4Tubs("TIARA_SiliconWafer_filled", 0.*mm, 160.*mm, (thickness_AA/2)*um, 0.*deg, (60.)*deg);
    
    G4Box* SiliconWafer_shear1 = new G4Box("SiliconWafer_shear1", (300/2)*mm, (20.)*mm, (100/2)*mm);
    
    G4ThreeVector   position_TIARA_PCB_shear[6];
    position_TIARA_PCB_shear[0] = G4ThreeVector(7.5760469*mm, 4.1074720*mm, 0.*mm);
    position_TIARA_PCB_shear[1] = G4ThreeVector(86.1695456*mm, -16.7875086*mm, 0.*mm);
    position_TIARA_PCB_shear[2] = G4ThreeVector(152.5181096*mm, 40.2365300*mm, 0.*mm);
    position_TIARA_PCB_shear[3] = G4ThreeVector(135.4631769*mm, 80.4016034*mm, 0.*mm);
    position_TIARA_PCB_shear[4] = G4ThreeVector(105.7796837*mm, 116.3682685*mm, 0.*mm);
    position_TIARA_PCB_shear[5] = G4ThreeVector(28.3491241*mm, 82.3831462*mm, 0.*mm);
    
    G4RotationMatrix* rm_SiliconWafer_shear1 = new G4RotationMatrix();
    G4RotationMatrix* rm_SiliconWafer_shear2 = new G4RotationMatrix();
    G4RotationMatrix* rm_SiliconWafer_shear3 = new G4RotationMatrix();
    G4RotationMatrix* rm_SiliconWafer_shear4 = new G4RotationMatrix();
    G4RotationMatrix* rm_SiliconWafer_shear5 = new G4RotationMatrix();
    G4RotationMatrix* rm_SiliconWafer_shear6 = new G4RotationMatrix();
    
    rm_SiliconWafer_shear1->rotateZ(59.6300000*deg);
    rm_SiliconWafer_shear2->rotateZ(-3.0000000*deg);
    rm_SiliconWafer_shear3->rotateZ(76.7362500*deg);
    rm_SiliconWafer_shear4->rotateZ(59.6300000*deg);
    rm_SiliconWafer_shear5->rotateZ(42.5237500*deg);
    rm_SiliconWafer_shear6->rotateZ(-57.7400000*deg);
    
    G4VSolid* Solid_TIARA_SiliconWafer_shear1 = new G4SubtractionSolid("TIARA_SiliconWafer_shear1", Solid_TIARA_SiliconWafer_filled, SiliconWafer_shear1, rm_SiliconWafer_shear1, position_TIARA_PCB_shear[0]);
    G4VSolid* Solid_TIARA_SiliconWafer_shear2 = new G4SubtractionSolid("TIARA_SiliconWafer_shear2", Solid_TIARA_SiliconWafer_shear1, SiliconWafer_shear1, rm_SiliconWafer_shear2, position_TIARA_PCB_shear[1]);
    G4VSolid* Solid_TIARA_SiliconWafer_shear3 = new G4SubtractionSolid("TIARA_SiliconWafer_shear3", Solid_TIARA_SiliconWafer_shear2, SiliconWafer_shear1, rm_SiliconWafer_shear3, position_TIARA_PCB_shear[2]);
    G4VSolid* Solid_TIARA_SiliconWafer_shear4 = new G4SubtractionSolid("TIARA_SiliconWafer_shear3", Solid_TIARA_SiliconWafer_shear3, SiliconWafer_shear1, rm_SiliconWafer_shear4, position_TIARA_PCB_shear[3]);
    G4VSolid* Solid_TIARA_SiliconWafer_shear5 = new G4SubtractionSolid("TIARA_SiliconWafer_shear3", Solid_TIARA_SiliconWafer_shear4, SiliconWafer_shear1, rm_SiliconWafer_shear5, position_TIARA_PCB_shear[4]);
    G4VSolid* Solid_TIARA_SiliconWafer_Final = new G4SubtractionSolid("TIARA_SiliconWafer_Final", Solid_TIARA_SiliconWafer_shear5, SiliconWafer_shear1, rm_SiliconWafer_shear6, position_TIARA_PCB_shear[5]);
    
    ////    This offset is the offset between the origins of the created PCB and the PCBpunch
    G4ThreeVector position_TIARA_SiliconWafer_shear = G4ThreeVector(0.*mm, 0.*mm, 0.*mm);
    
    G4RotationMatrix* rm_SiliconWafer_shear = new G4RotationMatrix();
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////    Final PCB
    ////    This includes the final PCB shear, which makes space for the Silicon Wafer
    
    // The "+10" for the thickness is to prevent logical issues when using boolean subtraction
    //G4Tubs* Solid_TIARA_SiliconWafer_finalshear_filled = new G4Tubs("TIARA_SiliconWafer_finalshear_filled", 0*mm, 160.*mm, ((thickness_AA/2) + 10)*um, 0.*deg, 60.*deg);
    G4Tubs* Solid_TIARA_SiliconWafer_finalshear_filled = new G4Tubs("TIARA_SiliconWafer_finalshear_filled", 0*mm, 160.*mm, ((400./2) + 10.)*um, 0.*deg, 60.*deg);
    
    G4ThreeVector position_TIARA_SiliconWafer_finalshear = G4ThreeVector(0.*mm, 0.*mm, ((thickness_PCB/2.) - (thickness_SiliconWafer/2.) + 10.)*um);
    
    G4VSolid* Solid_TIARA_finalshear1 = new G4SubtractionSolid("TIARA_finalshear1", Solid_TIARA_SiliconWafer_finalshear_filled, SiliconWafer_shear1, rm_SiliconWafer_shear1, position_TIARA_PCB_shear[0]);
    G4VSolid* Solid_TIARA_finalshear2 = new G4SubtractionSolid("TIARA_finalshear2", Solid_TIARA_finalshear1, SiliconWafer_shear1, rm_SiliconWafer_shear2, position_TIARA_PCB_shear[1]);
    G4VSolid* Solid_TIARA_finalshear3 = new G4SubtractionSolid("TIARA_finalshear3", Solid_TIARA_finalshear2, SiliconWafer_shear1, rm_SiliconWafer_shear3, position_TIARA_PCB_shear[2]);
    G4VSolid* Solid_TIARA_finalshear4 = new G4SubtractionSolid("TIARA_finalshear4", Solid_TIARA_finalshear3, SiliconWafer_shear1, rm_SiliconWafer_shear4, position_TIARA_PCB_shear[3]);
    G4VSolid* Solid_TIARA_finalshear5 = new G4SubtractionSolid("TIARA_finalshear5", Solid_TIARA_finalshear4, SiliconWafer_shear1, rm_SiliconWafer_shear5, position_TIARA_PCB_shear[4]);
    G4VSolid* Solid_TIARA_finalshear6 = new G4SubtractionSolid("TIARA_finalshear6", Solid_TIARA_finalshear5, SiliconWafer_shear1, rm_SiliconWafer_shear6, position_TIARA_PCB_shear[5]);
    
    //  purposely use Solid_TIARA_finalshear5 and not Solid_TIARA_finalshear6 to prevent logical error through boolean operation
    G4VSolid* Solid_TIARA_PCB_Final = new G4SubtractionSolid("TIARA_PCB_Final", Solid_TIARA_PCB_punched, Solid_TIARA_finalshear5, rm_SiliconWafer_shear, position_TIARA_SiliconWafer_finalshear);
    
    G4LogicalVolume* Logic_TIARA_PCB = new G4LogicalVolume(Solid_TIARA_PCB_Final, G4_PLASTIC_SC_VINYLTOLUENE_Material, "TIARA_PCB_Final", 0, 0, 0);
    
    
    
    /////////////////////////////////////////
    //          TIARA - Active Area
    
    G4double rotationphi;
    
    G4Tubs*             Solid_TIARA_RS_filled[16][8];
    G4VSolid*           Solid_TIARA_RS_shear[16][8];
    G4VSolid*           Solid_TIARA_RS[16][8];
    G4LogicalVolume*    Logic_TIARA_AA_RS[16][8];
    
    G4ThreeVector position_TIARA_Sector_shear[8][2];
    G4RotationMatrix* rm_TIARA_Sector_shear[8][2];
    
    for(G4int j=0; j<8; j++)
    {
        for(G4int l=0; l<2; l++)
        {
            rm_TIARA_Sector_shear[j][l] = new G4RotationMatrix();
        }
    }
    
    G4ThreeVector position_TIARA_Sector = G4ThreeVector(1.4448593*mm,  -0.0493525*mm, 0.*mm);
    
    G4Box* TIARA_Sector_razor = new G4Box("TIARA_Sector_razor", (150.)*mm, (150.)*mm, (100/2)*mm);
    
    position_TIARA_Sector_shear[0][0] = G4ThreeVector(93.0470539*mm, -145.2243182*mm, 0.*mm);
    position_TIARA_Sector_shear[0][1] = G4ThreeVector(58.2842550*mm,  162.2461746*mm, 0.*mm);
    rm_TIARA_Sector_shear[0][0]->rotateZ(-3.0000732*deg);
    rm_TIARA_Sector_shear[0][1]->rotateZ(-9.8424272*deg);
    
    position_TIARA_Sector_shear[1][0] = G4ThreeVector(106.8796942*mm, -133.5912908*mm, 0.*mm);
    position_TIARA_Sector_shear[1][1] = G4ThreeVector(35.1278912*mm,  167.0122283*mm, 0.*mm);
    rm_TIARA_Sector_shear[1][0]->rotateZ(-9.8425735*deg);
    rm_TIARA_Sector_shear[1][1]->rotateZ(-16.6849271*deg);
    
    position_TIARA_Sector_shear[2][0] = G4ThreeVector(126.6444045*mm, -118.5243029*mm, 0.*mm);
    position_TIARA_Sector_shear[2][1] = G4ThreeVector(19.1283325*mm,  171.8139616*mm, 0.*mm);
    rm_TIARA_Sector_shear[2][0]->rotateZ(-16.6850728*deg);
    rm_TIARA_Sector_shear[2][1]->rotateZ(-23.5274272*deg);
    
    position_TIARA_Sector_shear[3][0] = G4ThreeVector(138.3487001*mm, -103.2511752*mm, 0.*mm);
    position_TIARA_Sector_shear[3][1] = G4ThreeVector(-1.5603892*mm,  172.8208992*mm, 0.*mm);
    rm_TIARA_Sector_shear[3][0]->rotateZ(-23.5275730*deg);
    rm_TIARA_Sector_shear[3][1]->rotateZ(-30.3700000*deg);
    
    position_TIARA_Sector_shear[4][0] = G4ThreeVector(147.5802745*mm, -87.2542869*mm, 0.*mm);
    position_TIARA_Sector_shear[4][1] = G4ThreeVector(-24.6202415*mm,  169.5199469*mm, 0.*mm);
    rm_TIARA_Sector_shear[4][0]->rotateZ(-30.3700732*deg);
    rm_TIARA_Sector_shear[4][1]->rotateZ(-37.2124996*deg);
    
    position_TIARA_Sector_shear[5][0] = G4ThreeVector(164.4752025*mm, -63.3162228*mm, 0.*mm);
    position_TIARA_Sector_shear[5][1] = G4ThreeVector(-38.4388028*mm,  171.3807948*mm, 0.*mm);
    rm_TIARA_Sector_shear[5][0]->rotateZ(-37.2124996*deg);
    rm_TIARA_Sector_shear[5][1]->rotateZ(-44.0549272*deg);
    
    position_TIARA_Sector_shear[6][0] = G4ThreeVector(169.2760955*mm, -44.7896689*mm, 0.*mm);
    position_TIARA_Sector_shear[6][1] = G4ThreeVector(-59.5849311*mm,  164.3482428*mm, 0.*mm);
    rm_TIARA_Sector_shear[6][0]->rotateZ(-44.0550729*deg);
    rm_TIARA_Sector_shear[6][1]->rotateZ(-50.8974265*deg);
    
    position_TIARA_Sector_shear[7][0] = G4ThreeVector(176.8824451*mm, -20.0264597*mm, 0.*mm);
    position_TIARA_Sector_shear[7][1] = G4ThreeVector(-75.2407310*mm,  161.6241136*mm, 0.*mm);
    rm_TIARA_Sector_shear[7][0]->rotateZ(-50.8974999*deg);
    rm_TIARA_Sector_shear[7][1]->rotateZ(-57.7399268*deg);
    
    G4double innerRadius, outerRadius;
    
    for(G4int j=0; j<16; j++)
    {
        innerRadius = 32.60 + ((6.3125 + 0.1)*j); // mm
        outerRadius = innerRadius + 6.3125; // mm
        
        for(G4int l=0; l<8; l++)
        {
            Solid_TIARA_RS_filled[j][l] = new G4Tubs("TIARA_RS_filled", innerRadius*mm, outerRadius*mm, (thickness_AA/2)*um, 0.*deg, 60.*deg);
            
            Solid_TIARA_RS_shear[j][l] = new G4SubtractionSolid("TIARA_RS_shear", Solid_TIARA_RS_filled[j][l], TIARA_Sector_razor, rm_TIARA_Sector_shear[l][0], position_TIARA_Sector_shear[l][0]);
            
            Solid_TIARA_RS[j][l] = new G4SubtractionSolid("TIARA_RS", Solid_TIARA_RS_shear[j][l], TIARA_Sector_razor, rm_TIARA_Sector_shear[l][1], position_TIARA_Sector_shear[l][1]);
            
            Logic_TIARA_AA_RS[j][l] = new G4LogicalVolume(Solid_TIARA_RS[j][l], G4_Si_Material, "Logic_TIARA_AA_RS", 0, 0, 0);
            
            
        }
    }
    
    
    /////////////////////////////////////////////
    //          TIARA - Active Area Punch
    
    G4Tubs*             Solid_TIARA_RS_punch_filled[16][8];
    G4VSolid*           Solid_TIARA_RS_punch_shear[16][8];
    G4VSolid*           Solid_TIARA_RS_punch[16][8];
    G4LogicalVolume*    Logic_TIARA_RS_punch[16][8];
    
    
    for(G4int j=0; j<16; j++)
    {
        innerRadius = 32.60 + ((6.3125 + 0.1)*j); // mm
        outerRadius = innerRadius + 6.3125; // mm
        
        for(G4int l=0; l<8; l++)
        {
            //  The "+ 300" is to avoid a logical error occuring during the boolean subtraction
            //Solid_TIARA_RS_punch_filled[j][l] = new G4Tubs("TIARA_RS_punch_filled", innerRadius*mm, outerRadius*mm, ((thickness_AA/2) + 1000)*um, 0.*deg, 60.*deg);
            Solid_TIARA_RS_punch_filled[j][l] = new G4Tubs("TIARA_RS_punch_filled", innerRadius*mm, outerRadius*mm, 2.*mm, 0.*deg, 60.*deg);
            
            Solid_TIARA_RS_punch_shear[j][l] = new G4SubtractionSolid("TIARA_RS_punch_shear", Solid_TIARA_RS_punch_filled[j][l], TIARA_Sector_razor, rm_TIARA_Sector_shear[l][0], position_TIARA_Sector_shear[l][0]);
            
            Solid_TIARA_RS_punch[j][l] = new G4SubtractionSolid("TIARA_RS_punch", Solid_TIARA_RS_punch_shear[j][l], TIARA_Sector_razor, rm_TIARA_Sector_shear[l][1], position_TIARA_Sector_shear[l][1]);
            
            Logic_TIARA_RS_punch[j][l] = new G4LogicalVolume(Solid_TIARA_RS_punch[j][l], G4_Si_Material, "Logic_TIARA_RS_punch", 0, 0, 0);
            
        }
    }
    
    
    
    /////////////////////////////////////////////////
    //          TIARA - Inter-strip Region
    
    //TIARA_Interstrip_rotm.rotateZ(3.*deg);
    
    G4ThreeVector position_TIARA_Interstrip = position_TIARA_Sector;
    
    G4RotationMatrix* rm_TIARA_Interstrip = new G4RotationMatrix();
    //rm_TIARA_Interstrip->rotateZ(3.*deg);
    
    //G4Box* Solid_TIARA_Interstrip_filled = new G4Box("TIARA_Interstrip_filled", (40/2)*cm, (40/2)*cm, (thickness_AA/2)*um);
    G4Tubs* Solid_TIARA_Interstrip_filled = new G4Tubs("TIARA_Interstrip_filled", 0.*mm, 160.*mm, (thickness_AA/2)*um, 0.*deg, (60.)*deg);
    
    G4VSolid* Solid_TIARA_Interstrip_Final[16][8];
    
    for(G4int j=0; j<16; j++)
    {
        for(G4int l=0; l<8; l++)
        {
            if(j==0 && l==0)
            {
                //Solid_TIARA_Interstrip_Final[j][l] = new G4SubtractionSolid("Solid_TIARA_Interstrip_Final", Solid_TIARA_SiliconWafer_Final, Solid_TIARA_RS_punch[j][l], rm_TIARA_Interstrip, position_TIARA_Interstrip);
                
                //  Semi Working
                //Solid_TIARA_Interstrip_Final[j][l] = new G4SubtractionSolid("Solid_TIARA_Interstrip_Final", Solid_TIARA_SiliconWafer_filled, Solid_TIARA_RS_punch[j][l], rm_TIARA_Interstrip, position_TIARA_Interstrip);
                
                Solid_TIARA_Interstrip_Final[j][l] = new G4SubtractionSolid("Solid_TIARA_Interstrip_Final", Solid_TIARA_Interstrip_filled, Solid_TIARA_RS_punch[j][l], rm_TIARA_Interstrip, position_TIARA_Interstrip);
                
            }
            if(j!=0 && l==0)
            {
                Solid_TIARA_Interstrip_Final[j][l] = new G4SubtractionSolid("Solid_TIARA_Interstrip_Final", Solid_TIARA_Interstrip_Final[j-1][7], Solid_TIARA_RS_punch[j][l], rm_TIARA_Interstrip, position_TIARA_Interstrip);
                
            }
            if(l!=0)
            {
                Solid_TIARA_Interstrip_Final[j][l] = new G4SubtractionSolid("Solid_TIARA_Interstrip_Final", Solid_TIARA_Interstrip_Final[j][l-1], Solid_TIARA_RS_punch[j][l], rm_TIARA_Interstrip, position_TIARA_Interstrip);
                
            }
        }
    }
    
    //G4VSolid* Solid_TIARA_Interstrip = Solid_TIARA_Interstrip_Final[0][0];
    //G4VSolid* Solid_TIARA_Interstrip = Solid_TIARA_Interstrip_Final[15][7];
    //G4VSolid* Solid_TIARA_Interstrip = new G4SubtractionSolid("Solid_TIARA_Interstrip_Final", Solid_TIARA_Interstrip_shear4, Solid_TIARA_RS_punch[15][7], rm_TIARA_Interstrip, position_TIARA_Interstrip);
    //Solid_TIARA_Interstrip = new G4SubtractionSolid("Solid_TIARA_Interstrip_Final", Solid_TIARA_Interstrip, Solid_TIARA_RS_punch[15][6], rm_TIARA_Interstrip, position_TIARA_Interstrip);
    
    
    ////    STD Shearing attempt
    G4VSolid* Solid_TIARA_Interstrip_shear1 = new G4SubtractionSolid("Solid_TIARA_Interstrip_shear1", Solid_TIARA_Interstrip_Final[15][7], SiliconWafer_shear1, rm_SiliconWafer_shear1, position_TIARA_PCB_shear[0]);
    G4VSolid* Solid_TIARA_Interstrip_shear2 = new G4SubtractionSolid("Solid_TIARA_Interstrip_shear2", Solid_TIARA_Interstrip_shear1, SiliconWafer_shear1, rm_SiliconWafer_shear2, position_TIARA_PCB_shear[1]);
    G4VSolid* Solid_TIARA_Interstrip_shear3 = new G4SubtractionSolid("Solid_TIARA_Interstrip_shear3", Solid_TIARA_Interstrip_shear2, SiliconWafer_shear1, rm_SiliconWafer_shear3, position_TIARA_PCB_shear[2]);
    G4VSolid* Solid_TIARA_Interstrip_shear4 = new G4SubtractionSolid("Solid_TIARA_Interstrip_shear4", Solid_TIARA_Interstrip_shear3, SiliconWafer_shear1, rm_SiliconWafer_shear4, position_TIARA_PCB_shear[3]);
    G4VSolid* Solid_TIARA_Interstrip_shear5 = new G4SubtractionSolid("Solid_TIARA_Interstrip_shear5", Solid_TIARA_Interstrip_shear4, SiliconWafer_shear1, rm_SiliconWafer_shear5, position_TIARA_PCB_shear[4]);
    //G4VSolid* Solid_TIARA_Interstrip = new G4SubtractionSolid("Solid_TIARA_Interstrip", Solid_TIARA_Interstrip_shear5, SiliconWafer_shear1, rm_SiliconWafer_shear6, position_TIARA_PCB_shear[5]);
    
    //  Some progress
    //G4VSolid* Solid_TIARA_Interstrip = new G4SubtractionSolid("Solid_TIARA_Interstrip", Solid_TIARA_Interstrip_shear1, SiliconWafer_shear1, rm_SiliconWafer_shear6, position_TIARA_PCB_shear[5]);
    
    //G4VSolid* Solid_TIARA_Interstrip1 = new G4SubtractionSolid("Solid_TIARA_Interstrip1", Solid_TIARA_Interstrip_shear1, SiliconWafer_shear1, rm_SiliconWafer_shear6, position_TIARA_PCB_shear[5]);
    //G4VSolid* Solid_TIARA_Interstrip = new G4SubtractionSolid("Solid_TIARA_Interstrip", Solid_TIARA_Interstrip1, SiliconWafer_shear1, rm_SiliconWafer_shear2, position_TIARA_PCB_shear[1]);
    
    G4VSolid* Solid_TIARA_Interstrip1 = new G4SubtractionSolid("Solid_TIARA_Interstrip_shear3", Solid_TIARA_Interstrip_Final[15][7], SiliconWafer_shear1, rm_SiliconWafer_shear6, position_TIARA_PCB_shear[5]);
    //G4VSolid* Solid_TIARA_Interstrip = new G4SubtractionSolid("Solid_TIARA_Interstrip_shear3", Solid_TIARA_Interstrip1, SiliconWafer_shear1, rm_SiliconWafer_shear1, position_TIARA_PCB_shear[0]);
    
    
    
    //G4VSolid* Solid_TIARA_Interstrip3 = new G4SubtractionSolid("Solid_TIARA_Interstrip_shear3", Solid_TIARA_Interstrip2, SiliconWafer_shear1, rm_SiliconWafer_shear2, position_TIARA_PCB_shear[1]);
    
    
    
    /*
     ////    Specialised Punch
     G4VSolid* Solid_TIARA_Interstrip_punch1 = new G4UnionSolid("Solid_TIARA_Interstrip_punch1", SiliconWafer_shear1, SiliconWafer_shear1, rm_SiliconWafer_shear1, position_TIARA_PCB_shear[0]);
     G4VSolid* Solid_TIARA_Interstrip_punch2 = new G4UnionSolid("Solid_TIARA_Interstrip_punch1", SiliconWafer_shear1, Solid_TIARA_Interstrip_punch1, rm_SiliconWafer_shear1, position_TIARA_PCB_shear[1]);
     G4VSolid* Solid_TIARA_Interstrip_punch3 = new G4UnionSolid("Solid_TIARA_Interstrip_punch1", SiliconWafer_shear1, Solid_TIARA_Interstrip_punch2, rm_SiliconWafer_shear1, position_TIARA_PCB_shear[2]);
     G4VSolid* Solid_TIARA_Interstrip_punch4 = new G4UnionSolid("Solid_TIARA_Interstrip_punch1", SiliconWafer_shear1, Solid_TIARA_Interstrip_punch3, rm_SiliconWafer_shear1, position_TIARA_PCB_shear[3]);
     G4VSolid* Solid_TIARA_Interstrip_punch5 = new G4UnionSolid("Solid_TIARA_Interstrip_punch1", SiliconWafer_shear1, Solid_TIARA_Interstrip_punch4, rm_SiliconWafer_shear1, position_TIARA_PCB_shear[4]);
     */
    
    /*
     G4Box* TIARA_Interstrip_shear1 = new G4Box("TIARA_Interstrip_shear1", (400/2)*mm, (20.)*mm, (200/2)*mm);
     
     G4VSolid* Solid_TIARA_Interstrip_shear1 = new G4SubtractionSolid("Solid_TIARA_Interstrip_shear1", Solid_TIARA_Interstrip_Final[15][7], TIARA_Interstrip_shear1, rm_SiliconWafer_shear1, position_TIARA_PCB_shear[0]);
     G4VSolid* Solid_TIARA_Interstrip_shear2 = new G4SubtractionSolid("Solid_TIARA_Interstrip_shear2", Solid_TIARA_Interstrip_Final[15][7], TIARA_Interstrip_shear1, rm_SiliconWafer_shear2, position_TIARA_PCB_shear[1]);
     G4VSolid* Solid_TIARA_Interstrip_shear3 = new G4SubtractionSolid("Solid_TIARA_Interstrip_shear3", Solid_TIARA_Interstrip_Final[15][7], TIARA_Interstrip_shear1, rm_SiliconWafer_shear3, position_TIARA_PCB_shear[2]);
     G4VSolid* Solid_TIARA_Interstrip_shear4 = new G4SubtractionSolid("Solid_TIARA_Interstrip_shear4", Solid_TIARA_Interstrip_Final[15][7], TIARA_Interstrip_shear1, rm_SiliconWafer_shear4, position_TIARA_PCB_shear[3]);
     G4VSolid* Solid_TIARA_Interstrip_shear5 = new G4SubtractionSolid("Solid_TIARA_Interstrip_shear5", Solid_TIARA_Interstrip_Final[15][7], TIARA_Interstrip_shear1, rm_SiliconWafer_shear5, position_TIARA_PCB_shear[4]);
     G4VSolid* Solid_TIARA_Interstrip = new G4SubtractionSolid("Solid_TIARA_Interstrip", Solid_TIARA_Interstrip_Final[15][7], SiliconWafer_shear1, rm_SiliconWafer_shear6, position_TIARA_PCB_shear[5]);
     */
    
    /*
     G4RotationMatrix* rm_TIARA_Interstrip_finalpunch = new G4RotationMatrix();
     
     G4Tubs* Solid_TIARA_Interstrip_finalpunch_filled1 = new G4Tubs("TIARA_SiliconWafer_finalshear_filled1", 0.*mm, 180.*mm, 2.*cm, -10.*deg, 80.*deg);
     G4Tubs* Solid_TIARA_Interstrip_finalpunch_filled2 = new G4Tubs("TIARA_SiliconWafer_finalshear_filled2", 0.*mm, 180.*mm, 2.*cm, -10.*deg, 80.*deg);
     
     G4VSolid* Solid_TIARA_Interstrip_finalpunch1 = new G4SubtractionSolid("Solid_TIARA_Interstrip_finalpunch1", Solid_TIARA_Interstrip_finalpunch_filled1, SiliconWafer_shear1, rm_SiliconWafer_shear1, position_TIARA_PCB_shear[0]);
     
     G4VSolid* Solid_TIARA_Interstrip_finalpunch2 = new G4SubtractionSolid("Solid_TIARA_Interstrip_finalpunch2", Solid_TIARA_Interstrip_finalpunch1, SiliconWafer_shear1, rm_SiliconWafer_shear2, position_TIARA_PCB_shear[1]);
     
     G4VSolid* Solid_TIARA_Interstrip_finalpunch3 = new G4SubtractionSolid("Solid_TIARA_Interstrip_finalpunch3", Solid_TIARA_Interstrip_finalpunch2, SiliconWafer_shear1, rm_SiliconWafer_shear3, position_TIARA_PCB_shear[2]);
     
     G4VSolid* Solid_TIARA_Interstrip_finalpunch4 = new G4SubtractionSolid("Solid_TIARA_Interstrip_finalpunch4", Solid_TIARA_Interstrip_finalpunch3, SiliconWafer_shear1, rm_SiliconWafer_shear4, position_TIARA_PCB_shear[3]);
     
     G4VSolid* Solid_TIARA_Interstrip_finalpunch5 = new G4SubtractionSolid("Solid_TIARA_Interstrip_finalpunch5", Solid_TIARA_Interstrip_finalpunch4, SiliconWafer_shear1, rm_SiliconWafer_shear5, position_TIARA_PCB_shear[4]);
     
     G4VSolid* Solid_TIARA_Interstrip_finalpunch6 = new G4SubtractionSolid("Solid_TIARA_Interstrip_finalpunch6", Solid_TIARA_Interstrip_finalpunch5, SiliconWafer_shear1, rm_SiliconWafer_shear6, position_TIARA_PCB_shear[5]);
     
     
     G4VSolid* Solid_TIARA_Interstrip_finalpunch = new G4SubtractionSolid("Solid_TIARA_Interstrip_finalpunch6", Solid_TIARA_Interstrip_finalpunch_filled2, Solid_TIARA_Interstrip_finalpunch6, rm_TIARA_Interstrip_finalpunch, G4ThreeVector(0.*mm, 0.*mm, 0.*mm));
     
     
     //G4VSolid* Solid_TIARA_Interstrip = new G4SubtractionSolid("Solid_TIARA_Interstrip", Solid_TIARA_Interstrip_Final[15][7], Solid_TIARA_Interstrip_finalpunch, rm_TIARA_Interstrip_finalpunch, G4ThreeVector(0.*mm, 0.*mm, 0.*mm));
     //G4VSolid* Solid_TIARA_Interstrip = new G4SubtractionSolid("Solid_TIARA_Interstrip", Solid_TIARA_Interstrip_Final[15][7], Solid_TIARA_Interstrip_finalpunch, rm_TIARA_Interstrip_finalpunch, G4ThreeVector(0.*mm, 0.*mm, 0.*mm));
     
     G4Box* Solid_TIARA_Interstrip_template = new G4Box("TIARA_Interstrip_template", (50/2)*cm, (50/2)*cm, 1.*cm);
     
     G4VSolid* Solid_TIARA_Interstrip_intermediary = new G4SubtractionSolid("Solid_TIARA_Interstrip", Solid_TIARA_Interstrip_template, Solid_TIARA_finalshear6, rm_TIARA_Interstrip_finalpunch, G4ThreeVector(0.*mm, 0.*mm, 0.*mm));
     
     //G4VSolid* Solid_TIARA_Interstrip = new G4IntersectionSolid("TIARA_Interstrip", Solid_TIARA_Interstrip_Final[15][7], Solid_TIARA_finalshear6, rm_TIARA_Interstrip_finalpunch, G4ThreeVector(0.*mm, 0.*mm, 0.*mm));
     G4VSolid* Solid_TIARA_Interstrip = new G4SubtractionSolid("TIARA_Interstrip", Solid_TIARA_Interstrip_Final[15][7], Solid_TIARA_Interstrip_intermediary, rm_TIARA_Interstrip_finalpunch, G4ThreeVector(0.*mm, 0.*mm, 0.*mm));
     */
    
    
    
    
    /*
     G4Tubs* Solid_TIARA_SiliconWafer_finalshear_filled = new G4Tubs("TIARA_SiliconWafer_finalshear_filled", 0*mm, 160.*mm, ((thickness_AA/2) + 10)*um, 0.*deg, 60.*deg);
     
     G4ThreeVector position_TIARA_SiliconWafer_finalshear = G4ThreeVector(0.*mm, 0.*mm, ((thickness_PCB/2) - (thickness_SiliconWafer/2) + 10)*um);
     
     G4VSolid* Solid_TIARA_finalshear1 = new G4SubtractionSolid("TIARA_finalshear1", Solid_TIARA_SiliconWafer_finalshear_filled, SiliconWafer_shear1, rm_SiliconWafer_shear1, position_TIARA_PCB_shear[0]);
     G4VSolid* Solid_TIARA_finalshear2 = new G4SubtractionSolid("TIARA_finalshear2", Solid_TIARA_finalshear1, SiliconWafer_shear1, rm_SiliconWafer_shear2, position_TIARA_PCB_shear[1]);
     G4VSolid* Solid_TIARA_finalshear3 = new G4SubtractionSolid("TIARA_finalshear3", Solid_TIARA_finalshear2, SiliconWafer_shear1, rm_SiliconWafer_shear3, position_TIARA_PCB_shear[2]);
     G4VSolid* Solid_TIARA_finalshear4 = new G4SubtractionSolid("TIARA_finalshear4", Solid_TIARA_finalshear3, SiliconWafer_shear1, rm_SiliconWafer_shear4, position_TIARA_PCB_shear[3]);
     G4VSolid* Solid_TIARA_finalshear5 = new G4SubtractionSolid("TIARA_finalshear5", Solid_TIARA_finalshear4, SiliconWafer_shear1, rm_SiliconWafer_shear5, position_TIARA_PCB_shear[4]);
     */
    
    /*
     G4VSolid* Solid_TIARA_SiliconWafer_shear1 = new G4SubtractionSolid("TIARA_SiliconWafer_shear1", Solid_TIARA_SiliconWafer_filled, SiliconWafer_shear1, rm_SiliconWafer_shear1, position_TIARA_PCB_shear[0]);
     G4VSolid* Solid_TIARA_SiliconWafer_shear2 = new G4SubtractionSolid("TIARA_SiliconWafer_shear2", Solid_TIARA_SiliconWafer_shear1, SiliconWafer_shear1, rm_SiliconWafer_shear2, position_TIARA_PCB_shear[1]);
     G4VSolid* Solid_TIARA_SiliconWafer_shear3 = new G4SubtractionSolid("TIARA_SiliconWafer_shear3", Solid_TIARA_SiliconWafer_shear2, SiliconWafer_shear1, rm_SiliconWafer_shear3, position_TIARA_PCB_shear[2]);
     G4VSolid* Solid_TIARA_SiliconWafer_shear4 = new G4SubtractionSolid("TIARA_SiliconWafer_shear3", Solid_TIARA_SiliconWafer_shear3, SiliconWafer_shear1, rm_SiliconWafer_shear4, position_TIARA_PCB_shear[3]);
     G4VSolid* Solid_TIARA_SiliconWafer_shear5 = new G4SubtractionSolid("TIARA_SiliconWafer_shear3", Solid_TIARA_SiliconWafer_shear4, SiliconWafer_shear1, rm_SiliconWafer_shear5, position_TIARA_PCB_shear[4]);
     G4VSolid* Solid_TIARA_SiliconWafer_Final = new G4SubtractionSolid("TIARA_SiliconWafer_Final", Solid_TIARA_SiliconWafer_shear5, SiliconWafer_shear1, rm_SiliconWafer_shear6, position_TIARA_PCB_shear[5]);
     */
    
    
    //G4VSolid* Solid_TIARA_Interstrip = new G4SubtractionSolid("Solid_TIARA_Interstrip_shear3", Solid_TIARA_SiliconWafer_Final, SiliconWafer_shear1, rm_SiliconWafer_shear1, position_TIARA_PCB_shear[0]);
    
    //G4LogicalVolume* Logic_TIARA_Interstrip = new G4LogicalVolume(Solid_TIARA_Interstrip, G4_Si_Material, "Logic_TIARA_Interstrip", 0, 0, 0);
    
    G4LogicalVolume * Logic_TIARA_SiliconWafer[numberOf_TIARA];
    
    for(G4int i=0; i<numberOf_TIARA; i++)
    {
        Logic_TIARA_SiliconWafer[i] = new G4LogicalVolume(Solid_TIARA_SiliconWafer_Final, G4_Si_Material, "Logic_TIARA_Interstrip", 0, 0, 0);
    }
    
    //G4LogicalVolume* Logic_TIARA_Interstrip = new G4LogicalVolume(Solid_TIARA_SiliconWafer_Final, G4_Si_Material, "Logic_TIARA_Interstrip", 0, 0, 0);
    
    
    
    
    /////////////////////////////////////////////////////
    //          TIARA 2M Windows, Front and Back
    
    G4Tubs* Solid_TIARA_2M_Front = new G4Tubs("TIARA_2M_Front", 0*mm, 140.*mm, (thickness_2M_Front/2)*um, 0.*deg, (54.74)*deg);
    G4Tubs* Solid_TIARA_2M_Back = new G4Tubs("TIARA_2M_Back", 0*mm, 140.*mm, (thickness_2M_Back/2)*um, 0.*deg, (54.74)*deg);
    
    G4LogicalVolume* Logic_TIARA_2M_Front = new G4LogicalVolume(Solid_TIARA_2M_Front, G4_Al_Material, "Logic_TIARA_DL_Front", 0, 0, 0);
    G4LogicalVolume* Logic_TIARA_2M_Back = new G4LogicalVolume(Solid_TIARA_2M_Back, G4_Al_Material, "Logic_TIARA_DL_Back", 0, 0, 0);
    
    /////////////////////////////////////////////////////
    //          TIARA - Dead Layers, Front and Back
    /*
     G4Tubs* Solid_TIARA_DL_Front = new G4Tubs("TIARA_DL_Front", 0*mm, 140.*mm, (thickness_DL_Front/2)*um, 0.*deg, (54.74)*deg);
     G4Tubs* Solid_TIARA_DL_Back = new G4Tubs("TIARA_DL_Back", 0*mm, 140.*mm, (thickness_DL_Back/2)*um, 0.*deg, (54.74)*deg);
     
     G4LogicalVolume* Logic_TIARA_DL_Front = new G4LogicalVolume(Solid_TIARA_DL_Front, G4_Si_Material, "Logic_TIARA_DL_Front", 0, 0, 0);
     G4LogicalVolume* Logic_TIARA_DL_Back = new G4LogicalVolume(Solid_TIARA_DL_Back, G4_Si_Material, "Logic_TIARA_DL_Back", 0, 0, 0);
     */
    
    
    /////////////////////////////////////////
    //          TIARA ASSEMBLY
    
    G4Tubs* Solid_TIARA_Assembly_filled = new G4Tubs("Solid_TIARA_Assembly_filled", 0*mm, 170.*mm, (thickness_PCB/2)*um, 0*deg, (57.74)*deg);
    
    G4VSolid* Solid_TIARA_Assembly_shear1 = new G4SubtractionSolid("TIARA_Assembly_shear1", Solid_TIARA_Assembly_filled, PCB_razor1, rm1_TIARA_PCB, position_TIARA_PCB_razor1[0]);
    G4VSolid* Solid_TIARA_Assembly_shear2 = new G4SubtractionSolid("TIARA_Assembly_shear2", Solid_TIARA_Assembly_shear1, PCB_razor1, rm2_TIARA_PCB, position_TIARA_PCB_razor1[1]);
    G4VSolid* Solid_TIARA_Assembly_shear3 = new G4SubtractionSolid("TIARA_Assembly_shear3", Solid_TIARA_Assembly_shear2, PCB_razor1, rm3_TIARA_PCB, position_TIARA_PCB_razor1[2]);
    G4VSolid* Solid_TIARA_Assembly = new G4SubtractionSolid("Solid_TIARA_Assembly", Solid_TIARA_Assembly_shear3, PCB_razor1, rm4_TIARA_PCB, position_TIARA_PCB_razor1[3]);
    
    G4LogicalVolume*    Logic_TIARA_Assembly[numberOf_TIARA];
    
    for(G4int i=0; i<numberOf_TIARA; i++)
    {
        Logic_TIARA_Assembly[i] = new G4LogicalVolume(Solid_TIARA_Assembly, G4_Galactic_Material, "Logic_TIARA_Assembly", 0, 0, 0);
    }
    
    
    
    
    
    ////////////////////////////////////////////////////
    //               TIARA INITIALIZATION             //
    ////////////////////////////////////////////////////
    
    //G4ThreeVector offset_TIARA_AA = G4ThreeVector(0*cm, 0*cm, (thickness_PCB/2 - thickness_DL_Front - thickness_AA/2)*um);
    //G4ThreeVector offset_TIARA_AA = position_TIARA_Sector + G4ThreeVector(0*cm, 0*cm, (thickness_PCB/2 - thickness_DL_Front - (thickness_AA/2))*um);
    
    //G4ThreeVector offset_TIARA_SiliconWafer = G4ThreeVector(0*cm, 0*cm, (thickness_PCB/2 - thickness_DL_Front - (thickness_AA/2))*um);
    G4ThreeVector offset_TIARA_SiliconWafer = G4ThreeVector(0*cm, 0*cm, (thickness_PCB/2. - (400./2.))*um);
    
    G4ThreeVector offset_TIARA_AA = position_TIARA_Sector + G4ThreeVector(0*cm, 0*cm, ((400./2) - (thickness_AA/2) - thickness_DL_Front)*um);
    
    G4ThreeVector offset_TIARA_DL_Front = offset_TIARA_AA + G4ThreeVector(0*mm, 0*mm, ((thickness_AA/2 + thickness_DL_Front/2))*um);
    G4ThreeVector offset_TIARA_DL_Back = offset_TIARA_AA + G4ThreeVector(0*mm, 0*mm, -((thickness_AA/2 + thickness_DL_Front/2))*um);
    
    G4ThreeVector offset_TIARA_2M_Front = offset_TIARA_DL_Front + G4ThreeVector(0*mm, 0*mm, ((thickness_DL_Front/2 + thickness_2M_Front/2))*um);
    G4ThreeVector offset_TIARA_2M_Back = offset_TIARA_DL_Back + G4ThreeVector(0*mm, 0*mm, -((thickness_DL_Back/2 + thickness_2M_Back/2))*um);
    
    //G4ThreeVector offset_TIARA_Interstrip = position_TIARA_Sector + G4ThreeVector(0*mm, 0*mm, ((thickness_AA/2 + thickness_DL_Front/2))*um);
    //G4ThreeVector offset_TIARA_Interstrip = position_TIARA_Sector + G4ThreeVector(0*cm, 0*cm, (thickness_PCB/2 - thickness_DL_Front - (thickness_AA/2))*um);
    
    //    G4ThreeVector offset_TIARA_Interstrip = G4ThreeVector(2.4748830*mm, 1.3857888*mm, (thickness_PCB/2 - thickness_DL_Front - (thickness_AA/2))*um);
    
    TIARA_SiliconWafer_transform = G4Transform3D(TIARA_SiliconWafer_rotm, offset_TIARA_SiliconWafer);
    
    
    for(G4int i=0; i<numberOf_TIARA; i++)
    {
        TIARA_transform[i] = G4Transform3D(TIARA_rotm[i],TIARA_AA_CentrePosition[i]);
        
        
        if(TIARA_Presence[i])
        {
            ///////////////////////////////////////
            //      TIARA AA - Rings and Sectors
            
            for(G4int j=0; j<16; j++)
            {
                for(G4int l=0; l<8; l++)
                {
                    TIARA_AA_RS_transform[l] = G4Transform3D(TIARA_AA_RS_rotm[l], offset_TIARA_AA);
                    
                    /*
                     PhysiTIARA_AA_RS = new G4PVPlacement(TIARA_AA_RS_transform[l],
                     Logic_TIARA_AA_RS[j][l],
                     "TIARA_AA_RS", // its name
                     Logic_TIARA_Assembly[i],
                     false, // no boolean operations
                     i*128 + l + j*8,  // copy number
                     fCheckOverlaps); // checking overlaps
                     */
                    
                    ////    STD
                    /*
                     PhysiTIARA_AA_RS = new G4PVPlacement(TIARA_AA_RS_transform[l],
                     Logic_TIARA_AA_RS[j][l],
                     "TIARA_AA_RS", // its name
                     Logic_TIARA_SiliconWafer[i],
                     false, // no boolean operations
                     i*128 + l + j*8,  // copy number
                     fCheckOverlaps); // checking overlaps
                     */
                    
                    /*
                     PhysiTIARA_AA_RS = new G4PVPlacement(TIARA_AA_RS_transform[l],
                     Logic_TIARA_AA_RS[j][l],
                     "TIARA_AA_RS", // its name
                     Logic_TIARA_Assembly[i],
                     false, // no boolean operations
                     i*128 + l + j*8,  // copy number
                     fCheckOverlaps); // checking overlaps
                     */
                    
                    PhysiTIARA_AA_RS = new G4PVPlacement(TIARA_AA_RS_transform[l],
                                                         Logic_TIARA_AA_RS[j][l],
                                                         "TIARA_AA_RS", // its name
                                                         Logic_TIARA_SiliconWafer[i],
                                                         false, // no boolean operations
                                                         i*128 + l + j*8,  // copy number
                                                         fCheckOverlaps); // checking overlaps
                    
                }
            }
            
            
            ////////////////////////////////
            //      TIARA INTERSTRIP
            
            PhysiTIARA_SiliconWafer = new G4PVPlacement(TIARA_SiliconWafer_transform,
                                                        Logic_TIARA_SiliconWafer[i],
                                                        "TIARA_SiliconWafer", // its name
                                                        Logic_TIARA_Assembly[i],
                                                        false, // no boolean operations
                                                        i,  // copy number
                                                        fCheckOverlaps); // checking overlaps
            
            
            ////////////////////////////////
            //      TIARA PCB
            
            new G4PVPlacement(0,               // no rotation
                              G4ThreeVector(), // at (x,y,z)
                              Logic_TIARA_PCB,       // its logical volume
                              "TIARA_PCB",       // its name
                              Logic_TIARA_Assembly[i],         // its mother  volume
                              false,           // no boolean operations
                              0,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            /*
             ////////////////////////////////////////////
             //      TIARA 2M Windows, Front and Back
             
             TIARA_2M_transform[0] = G4Transform3D(TIARA_DL_2M_rotm, offset_TIARA_2M_Front);
             TIARA_2M_transform[1] = G4Transform3D(TIARA_DL_2M_rotm, offset_TIARA_2M_Back);
             
             new G4PVPlacement(TIARA_2M_transform[0],
             Logic_TIARA_2M_Front,       // its logical volume
             "Physical_TIARA_2M_Front",       // its name
             Logic_TIARA_Assembly[i],         // its mother  volume
             false,           // no boolean operations
             0,               // copy number
             fCheckOverlaps); // checking overlaps
             
             new G4PVPlacement(TIARA_2M_transform[1],
             Logic_TIARA_2M_Back,       // its logical volume
             "Physical_TIARA_2M_Back",       // its name
             Logic_TIARA_Assembly[i],         // its mother  volume
             false,           // no boolean operations
             0,               // copy number
             fCheckOverlaps); // checking overlaps
             
             
             ////////////////////////////////////////////
             //      TIARA - Dead Layers, Front and Back
             
             TIARA_DL_transform[0] = G4Transform3D(TIARA_DL_2M_rotm, offset_TIARA_DL_Front);
             TIARA_DL_transform[1] = G4Transform3D(TIARA_DL_2M_rotm, offset_TIARA_DL_Back);
             
             new G4PVPlacement(TIARA_DL_transform[0],
             Logic_TIARA_DL_Front,       // its logical volume
             "Physical_TIARA_DL_Front",       // its name
             Logic_TIARA_Assembly[i],         // its mother  volume
             false,           // no boolean operations
             0,               // copy number
             fCheckOverlaps); // checking overlaps
             
             new G4PVPlacement(TIARA_DL_transform[1],
             Logic_TIARA_DL_Front,       // its logical volume
             "Physical_TIARA_DL_Back",       // its name
             Logic_TIARA_Assembly[i],         // its mother  volume
             false,           // no boolean operations
             0,               // copy number
             fCheckOverlaps); // checking overlaps
             */
            ////////////////////////////////
            //      TIARA ASSEMBLY
            
            new G4PVPlacement(TIARA_transform[i], // transformation
                              Logic_TIARA_Assembly[i],
                              "Physical_TIARA_Assembly",
                              LogicVacuumChamber,
                              false,           // no boolean
                              0,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            
        }
    }
    
    
    //////////////////////////////////////
    //          VDC DEFINITION          //
    //////////////////////////////////////
    
    ///////////////////////////////
    //      VDC - ASSEMBLY
    G4Box* Solid_VDC_ASSEMBLY_USDS = new G4Box("Solid_Stesalit_Frame", (936./2)*mm, (240./2)*mm, (24./2)*mm);
    G4Box* Solid_VDC_ASSEMBLY = new G4Box("Solid_Stesalit_Frame", (936./2)*mm, (240./2)*mm, 78./2*mm);
    
    G4LogicalVolume* Logic_VDC_SenseRegion_USDS[numberOf_VDC][2];
    G4LogicalVolume* Logic_VDC_ASSEMBLY[numberOf_VDC];
    
    for(G4int i=0; i<2; i++)
    {
        Logic_VDC_ASSEMBLY[i] = new G4LogicalVolume(Solid_VDC_ASSEMBLY, G4_AIR_Material, "Logic_VDC_ASSEMBLY",0,0,0);
        
        for(G4int j=0; j<2; j++)
        {
            Logic_VDC_SenseRegion_USDS[i][j] = new G4LogicalVolume(Solid_VDC_ASSEMBLY_USDS, VDC_SR_Gas_Material, "Logic_VDC_SenseRegion_USDS",0,0,0);
        }
    }
    
    ///////////////////////////////////////////////
    //      VDC - Stesalit Frame Templates
    G4double thickness_Stesalit_Frame = 8.0; // mm
    G4double thickness_StesalitPCB_Frame = 5.5; // mm
    
    ///////////////////////////////////////////////
    //      VDC - STESALIT STANDARD FRAME
    G4Box* Solid_Stesalit_StdFrame = new G4Box("Solid_Stesalit_StdFrame", (936./2)*mm, (240./2)*mm, (8./2)*mm);
    
    ///////////////////////////////////////////////
    //      VDC - STESALIT STANDARD PCB FRAME
    G4Box* Solid_StesalitPCB_StdFrame = new G4Box("Solid_StesalitPCB_StdFrame", (936./2)*mm, (240./2)*mm, (5.5/2)*mm);
    G4Box* Solid_Wire_subtraction = new G4Box("Solid_Wire_subtraction", (2000./2)*mm, (1000./2)*mm, (10./2)*mm);
    
    ///////////////////////////////////////////////
    //      VDC - PCB FRAME
    G4Box* Solid_PCB_StdFrame = new G4Box("Solid_PCB_StdFrame", (936./2)*mm, (240./2)*mm, (2.5/2)*mm);
    
    ///////////////////////////////////////////////
    //      VDC - ALUMINIUM OUTER FRAME
    G4Box* Solid_Al_StdFrame = new G4Box("Solid_Al_StdFrame", (936./2)*mm, (240./2)*mm, (15./2)*mm);
    
    ///////////////////////////////////////////////
    //      VDC - FRAME FINISHING
    G4ThreeVector   position_punch[4];
    
    G4Box* punch1 = new G4Box("punch1", (825./2)*mm, (100./2)*mm, (20./2)*mm);
    G4Box* punch2 = new G4Box("punch2", (809./2)*mm, (131./2)*mm, (20./2)*mm);
    G4Box* punch3 = new G4Box("punch3", (800./2)*mm, (100./2)*mm, (20./2)*mm);
    G4Box* punch4 = new G4Box("punch4", (852./2)*mm, (100./2)*mm, (20./2)*mm);
    
    ///////////////////////////////////////////////
    //      VDC - STESALIT HV Frames, Upstream and Downstream
    position_punch[0] = G4ThreeVector( (-9.5)*mm, 0.*mm, 0.*mm);
    
    G4VSolid* Solid_Stesalit_HV_USDS_Frame = new G4SubtractionSolid("Solid_Stesalit_HV_USDS_Frame", Solid_Stesalit_StdFrame, punch1, 0, position_punch[0]);
    
    ///////////////////////////////////////////////
    position_punch[0] = G4ThreeVector( -(825./2 - 8.7)*mm, 0.*mm, (0.)*mm);
    
    G4Box* extrusion1 = new G4Box("extrusion1", (830./2)*mm, (113./2)*mm, (20./2)*mm);
    
    G4VSolid* Solid_VDC_HVFrame_USDS2 = new G4SubtractionSolid("Solid_VDC_HVFrame_USDS2", Solid_Stesalit_HV_USDS_Frame, extrusion1, 0, G4ThreeVector( (-9.5-2.5)*mm, 0.*mm, (4+10-2)*mm));
    
    ///////////////////////////////////////////////
    
    G4RotationMatrix* rm_punch = new G4RotationMatrix();
    rm_punch->rotateY(-35.*deg);
    
    G4VSolid* Solid_VDC_HVFrame_USDS = new G4SubtractionSolid("Solid_VDC_HVFrame_USDS", Solid_VDC_HVFrame_USDS2, punch1, rm_punch, position_punch[0]);
    
    G4LogicalVolume*    Logic_VDC_HV_Frame[3];
    
    Logic_VDC_HV_Frame[0] =  new G4LogicalVolume(Solid_VDC_HVFrame_USDS, G4_Al_Material, "Logic_VDC_HV_Frame",0,0,0);
    Logic_VDC_HV_Frame[2] = new G4LogicalVolume(Solid_VDC_HVFrame_USDS, G4_Al_Material, "Logic_VDC_HV_Frame",0,0,0);
    
    
    ///////////////////////////////////////////////
    //      VDC - STESALIT GAS CELL
    position_punch[1] = G4ThreeVector( -1.5*mm, -1.5*mm, 0.*mm);
    
    G4VSolid* Solid_Stesalit_GasFrame = new G4SubtractionSolid("Solid_Stesalit_GasFrame", Solid_Stesalit_StdFrame, punch2, 0, position_punch[1]);
    
    position_punch[1] = G4ThreeVector( -(809./2 + 0.1 - 8.192)*mm, -1.5*mm, 1.3*mm);
    
    G4VSolid* Solid_VDC_GasFrame = new G4SubtractionSolid("Solid_VDC_GasFrame", Solid_Stesalit_GasFrame, punch2, rm_punch, position_punch[1]);
    
    G4LogicalVolume* Logic_VDC_GasFrame = new G4LogicalVolume(Solid_VDC_GasFrame, G4_Al_Material, "Logic_VDC_GasFrame",0,0,0);
    
    
    ///////////////////////////////////////////////
    //      VDC - STESALIT X-WIRE/U-WIRE FRAME
    position_punch[2] = G4ThreeVector( 0.*mm, 0.*mm, 0.*mm);
    
    G4VSolid* Solid_VDC_Stesalit_XU_Frame = new G4SubtractionSolid("Solid_VDC_Stesalit_XU_Frame", Solid_StesalitPCB_StdFrame, punch3, 0, position_punch[2]);
    
    G4VSolid* Solid_Wire_subtraction2 = new G4SubtractionSolid("Solid_VDC_Stesalit_XU_Frame", Solid_Wire_subtraction, punch3, 0, position_punch[2]);
    
    G4LogicalVolume* Logic_VDC_XU_Frame = new G4LogicalVolume(Solid_VDC_Stesalit_XU_Frame, G4_Al_Material, "Logic_VDC_XU_Frame",0,0,0);
    
    
    ///////////////////////////////////////////////
    //      VDC - PCB X-WIRE/U-WIRE FRAME
    position_punch[2] = G4ThreeVector( 0.*mm, 0.*mm, 0.*mm);
    
    G4VSolid* Solid_VDC_PCB_XU_Frame = new G4SubtractionSolid("Solid_VDC_PCB_XU_Frame", Solid_PCB_StdFrame, punch3, 0, position_punch[2]);
    
    G4LogicalVolume* Logic_VDC_XU_PCBFrame = new G4LogicalVolume(Solid_VDC_PCB_XU_Frame, G4_Al_Material, "Logic_VDC_XU_PCBFrame",0,0,0);
    
    ///////////////////////////////////////////////
    //      VDC - STESALIT HV Frame, MIDDLE
    position_punch[3] = G4ThreeVector( 0.*mm, 0.*mm, 0.*mm);
    
    G4VSolid* Solid_Stesalit_HV_MID_Frame = new G4SubtractionSolid("Solid_Stesalit_HV_MID_Frame", Solid_Stesalit_StdFrame, punch4, 0, position_punch[3]);
    
    Logic_VDC_HV_Frame[1] = new G4LogicalVolume(Solid_Stesalit_HV_MID_Frame, G4_Al_Material, "Logic_VDC_HV_Frame",0,0,0);
    
    
    ///////////////////////////////////////////////
    //      VDC - ALUMINIUM OUTER FRAME
    G4Box* AlFrame_punch = new G4Box("AlFrame_punch", (834./2)*mm, (100./2)*mm, (20./2)*mm);
    
    G4Box* Solid_Aluminium_Frame = new G4Box("Solid_Stesalit_Frame", (936./2)*mm, (240./2)*mm, (15./2)*mm);
    
    G4VSolid* Solid_VDC_Al_Frame_shear2perf = new G4SubtractionSolid("Solid_VDC_Al_Frame_shear2perf", Solid_Aluminium_Frame, AlFrame_punch, 0, G4ThreeVector( 0.*mm, 0.*mm, 0.*mm));
    
    G4VSolid* Solid_VDC_Al_Frame = new G4SubtractionSolid("Solid_VDC_Al_Frame", Solid_VDC_Al_Frame_shear2perf, AlFrame_punch, rm_punch, G4ThreeVector( -(834./2 + 14.0 - 8.192)*mm, 0.*mm, 13.0*mm));
    
    G4LogicalVolume* Logic_VDC_Al_Frame = new G4LogicalVolume(Solid_VDC_Al_Frame, G4_Al_Material, "Logic_VDC_Al_Frame",0,0,0);
    
    
    ///////////////////////////////////////////////
    //      VDC - HIGH VOLTAGE PLANE
    G4Box* Solid_VDC_HV_Plane_US = new G4Box("Solid_VDC_HV_Plane_US", (809./2)*mm, (100./2)*mm, (20./2)*um);
    G4Box* Solid_VDC_HV_Plane_MID = new G4Box("Solid_VDC_HV_Plane_MID", (800./2)*mm, (100./2)*mm, (20./2)*um);
    G4Box* Solid_VDC_HV_Plane_DS = new G4Box("Solid_VDC_HV_Plane_DS", (800./2)*mm, (100./2)*mm, (20./2)*um);
    
    G4LogicalVolume*    Logic_VDC_HV_Plane[3];
    
    Logic_VDC_HV_Plane[0] = new G4LogicalVolume(Solid_VDC_HV_Plane_US, G4_Al_Material, "Logic_VDC_HV_Plane_US",0,0,0);
    Logic_VDC_HV_Plane[1] = new G4LogicalVolume(Solid_VDC_HV_Plane_MID, G4_Al_Material, "Logic_VDC_HV_Plane_MID",0,0,0);
    Logic_VDC_HV_Plane[2] = new G4LogicalVolume(Solid_VDC_HV_Plane_DS, G4_Al_Material, "Logic_VDC_HV_Plane_DS",0,0,0);
    
    ///////////////////////////////////////////////
    //      VDC - MYLAR WINDOW
    G4Box* Solid_VDC_MYLAR_Plane = new G4Box("Solid_VDC_MYLAR_Plane", (834./2)*mm, (100./2)*mm, (25./2)*um);
    
    G4LogicalVolume* Logic_VDC_MYLAR_Plane = new G4LogicalVolume(Solid_VDC_MYLAR_Plane, G4_Mylar_Material, "Logic_VDC_MYLAR_Plane",0,0,0);
    
    
    ///////////////////////////////////////////////
    //      VDC - X WIRES
    G4Tubs* Solid_VDC_X_WIRE = new G4Tubs("Solid_VDC_X_WIRE", 0.*um, 20.*um, 100./2*mm, 0.*deg, 360.*deg);
    
    G4LogicalVolume* Logic_VDC_X_WIRE = new G4LogicalVolume(Solid_VDC_X_WIRE, G4_W_Material, "Logic_VDC_X_WIRE",0,0,0);
    
    ///////////////////////////////////////////////
    //      VDC - X GUARD WIRES
    G4Tubs* Solid_VDC_X_GUARDWIRE = new G4Tubs("Solid_VDC_X_GUARDWIRE", 0.*um, 50.*um, 100./2*mm, 0.*deg, 360.*deg);
    
    G4LogicalVolume* Logic_VDC_X_GUARDWIRE = new G4LogicalVolume(Solid_VDC_X_GUARDWIRE, G4_W_Material, "Logic_VDC_X_GUARDWIRE",0,0,0);
    
    
    ///////////////////////////////////////////////
    //      VDC - X GUARD WIRES, THICK
    G4Tubs* Solid_VDC_X_GUARDWIRE_THICK = new G4Tubs("Solid_VDC_X_GUARDWIRE_THICK", 0.*um, 100.*um, 100./2*mm, 0.*deg, 360.*deg);
    
    G4LogicalVolume* Logic_VDC_X_GUARDWIRE_Thick = new G4LogicalVolume(Solid_VDC_X_GUARDWIRE_THICK, G4_W_Material, "Logic_VDC_X_GUARDWIRE_Thick",0,0,0);
    
    ///////////////////////////////////////////////
    //      VDC - U WIRES
    G4Tubs* Solid_VDC_U_WIRE_Full = new G4Tubs("Solid_VDC_U_WIRE", 0.*um, 20.*um, 150./2*mm, 0.*deg, 360.*deg);
    
    G4RotationMatrix* rm_U_WIRE = new G4RotationMatrix();
    rm_U_WIRE->rotateY(90.*deg);
    rm_U_WIRE->rotateZ((90.-40.)*deg);
    
    G4VSolid* Solid_VDC_U_WIRE;
    G4ThreeVector   offset_VDC_U_WIRE_subtraction;
    
    G4LogicalVolume* Logic_VDC_U_WIRE[143];
    
    for(G4int k=0; k<143; k++)
    {
        offset_VDC_U_WIRE_subtraction = G4ThreeVector( 0*mm, (-k + 71)*4*mm, (k - 71)*(4/tan(50*deg))*mm);
        
        Solid_VDC_U_WIRE = new G4SubtractionSolid("Solid_VDC_U_WIRE", Solid_VDC_U_WIRE_Full, Solid_Wire_subtraction2, rm_U_WIRE, offset_VDC_U_WIRE_subtraction);
        
        Logic_VDC_U_WIRE[k] = new G4LogicalVolume(Solid_VDC_U_WIRE, G4_W_Material, "Logic_VDC_U_WIRE",0,0,0);
    }
    
    ///////////////////////////////////////////////
    //      VDC - U GUARD WIRES
    G4Tubs* Solid_VDC_U_GUARDWIRE_Full = new G4Tubs("Solid_VDC_U_GUARDWIRE", 0.*um, 50.*um, 150./2*mm, 0.*deg, 360.*deg);
    
    G4VSolid* Solid_VDC_U_GUARDWIRE;
    G4LogicalVolume* Logic_VDC_U_GUARDWIRE[144];
    
    for(G4int k=0; k<144; k++)
    {
        offset_VDC_U_WIRE_subtraction = G4ThreeVector( 0*mm, (-k + 71.5)*4*mm, (k - 71.5)*(4/tan(50*deg))*mm);
        
        Solid_VDC_U_GUARDWIRE = new G4SubtractionSolid("Solid_VDC_U_WIRE", Solid_VDC_U_GUARDWIRE_Full, Solid_Wire_subtraction2, rm_U_WIRE, offset_VDC_U_WIRE_subtraction);
        
        Logic_VDC_U_GUARDWIRE[k] = new G4LogicalVolume(Solid_VDC_U_GUARDWIRE, G4_W_Material, "Logic_VDC_U_GUARDWIRE",0,0,0);
        
    }
    
    ////////////////////////////////////////////////
    //      VDC - U GUARD WIRES, Thick
    G4Tubs* Solid_VDC_U_GUARDWIRE_Thick_Full = new G4Tubs("Solid_VDC_U_GUARDWIRE_Thick", 0.*um, 100.*um, 150./2*mm, 0.*deg, 360.*deg);
    
    G4VSolid* Solid_VDC_U_GUARDWIRE_Thick;
    
    G4LogicalVolume*    Logic_VDC_U_GUARDWIRE_Thick[2];
    
    for(G4int k=0; k<2; k++)
    {
        offset_VDC_U_WIRE_subtraction = G4ThreeVector( 0*mm, (-k*144. + 72.)*4*mm, (k*144. - 72.)*(4/tan(50*deg))*mm);
        
        Solid_VDC_U_GUARDWIRE_Thick = new G4SubtractionSolid("Solid_VDC_U_WIRE_Thick", Solid_VDC_U_GUARDWIRE_Thick_Full, Solid_Wire_subtraction2, rm_U_WIRE, offset_VDC_U_WIRE_subtraction);
        
        Logic_VDC_U_GUARDWIRE_Thick[k] = new G4LogicalVolume(Solid_VDC_U_GUARDWIRE_Thick, G4_W_Material, "Logic_VDC_U_GUARDWIRE_Thick",0,0,0);
    }
    
    
    ////////////////////////////////////////////////
    //             VDC INITIALISATION             //
    ////////////////////////////////////////////////
    
    G4ThreeVector   offset_VDC_SenseRegion_USDS[3];
    offset_VDC_SenseRegion_USDS[0] = G4ThreeVector( 0*mm, 0*mm, 1.5*thickness_Stesalit_Frame*mm);
    offset_VDC_SenseRegion_USDS[1] = G4ThreeVector( 0*mm, 0*mm, -1.5*thickness_Stesalit_Frame*mm);
    
    G4ThreeVector   offset_VDC_HVFrame[3];
    offset_VDC_HVFrame[0] = G4ThreeVector( 0*mm, 0*mm, thickness_Stesalit_Frame*mm);
    offset_VDC_HVFrame[1] = G4ThreeVector( 0*mm, 0*mm, thickness_Stesalit_Frame*mm);
    offset_VDC_HVFrame[2] = G4ThreeVector( 0*mm, 0*mm, -thickness_Stesalit_Frame*mm);
    VDC_HVFrame_rotm[2].rotateY(180.*deg);
    
    G4ThreeVector   offset_VDC_HV_Plane[3];
    offset_VDC_HV_Plane[0] = G4ThreeVector( -1.5*mm, 0.*mm, (4000.)*um);
    offset_VDC_HV_Plane[1] = G4ThreeVector( 0.*mm, 0.*mm, (12000. - 10.)*um);
    offset_VDC_HV_Plane[2] = G4ThreeVector( 0.*mm, 0.*mm, -(3990.)*um);
    
    G4ThreeVector   offset_VDC_GasFrame;
    offset_VDC_GasFrame = G4ThreeVector( 0*mm, 0*mm, 0.*mm);
    
    G4ThreeVector   offset_VDC_XU_Frame[2];
    offset_VDC_XU_Frame[0] = G4ThreeVector( 0*mm, 0.*mm, (-1.5*thickness_Stesalit_Frame + 5.5/2)*mm);
    offset_VDC_XU_Frame[1] = G4ThreeVector( 0*mm, 0.*mm, -(2.5/2)*mm);
    
    G4ThreeVector   offset_VDC_XU_PCBFrame[2];
    offset_VDC_XU_PCBFrame[0] = G4ThreeVector( 0*mm, 0.*mm, (-thickness_Stesalit_Frame/2 - 2.5/2)*mm);
    offset_VDC_XU_PCBFrame[1] = G4ThreeVector( 0*mm, 0.*mm, (thickness_Stesalit_Frame/2 - 2.5/2)*mm);
    
    G4ThreeVector   offset_VDC_Al_Frame[2];
    offset_VDC_Al_Frame[0] = G4ThreeVector( 0*mm, 0.*mm, (3.*thickness_Stesalit_Frame  + 7.5)*mm);
    offset_VDC_Al_Frame[1] = G4ThreeVector( 0*mm, 0.*mm, -(3.*thickness_Stesalit_Frame + 7.5)*mm);
    VDC_Al_Frame_rotm[1].rotateY(180.*deg);
    
    G4ThreeVector   offset_VDC_MYLAR_Plane[2];
    offset_VDC_MYLAR_Plane[0] = G4ThreeVector( 0*mm, 0.*mm, (24000. + 25./2)*um);
    offset_VDC_MYLAR_Plane[1] = G4ThreeVector( 0*mm, 0.*mm, -(24000. + 25./2)*um);
    
    G4RotationMatrix* rm_GAS_VDC_HV = new G4RotationMatrix();
    rm_GAS_VDC_HV->rotateY(180.*deg);
    
    G4int usds;
    
    G4ThreeVector   offset_VDC_X_WIRE;
    VDC_X_WIRE_rotm.rotateX(90.*deg);
    
    G4ThreeVector   offset_VDC_U_WIRE;
    VDC_U_WIRE_rotm.rotateX(-90.*deg);
    VDC_U_WIRE_rotm.rotateY(-90.*deg);
    VDC_U_WIRE_rotm.rotateZ(40.*deg);
    
    
    for(G4int i=0; i<2; i++)
    {
        if(VDC_Presence[i])
        {
            ///////////////////////////////////////////
            //      VDC - STESALIT STANDARD FRAME
            for(G4int j=0; j<3; j++)
            {
                VDC_HVFrame_transform = G4Transform3D(VDC_HVFrame_rotm[j], offset_VDC_HVFrame[j]);
                
                if(j==0) usds = 0;
                if(j==1 || j ==2) usds = 1;
                
                
                new G4PVPlacement(VDC_HVFrame_transform,
                                  Logic_VDC_HV_Frame[j],
                                  "VDC_HV_Frame",
                                  Logic_VDC_SenseRegion_USDS[i][usds],
                                  false,    // no boolean operations
                                  i*2 + j,    // copy number
                                  fCheckOverlaps); // checking overlaps
                
                new G4PVPlacement(0,
                                  offset_VDC_HV_Plane[j],
                                  Logic_VDC_HV_Plane[j],
                                  "VDC_HV_Plane",
                                  Logic_VDC_SenseRegion_USDS[i][usds],
                                  false,    // no boolean operations
                                  i*2 + j,    // copy number
                                  fCheckOverlaps); // checking overlaps
                
            }
            
            ///////////////////////////////////////
            //      VDC - GAS CELL
            VDC_GasFrame_transform = G4Transform3D(VDC_GasFrame_rotm, offset_VDC_GasFrame);
            
            new G4PVPlacement(VDC_GasFrame_transform,
                              Logic_VDC_GasFrame,
                              "VDC_GasFrame",
                              Logic_VDC_SenseRegion_USDS[i][0],
                              false,    // no boolean operations
                              0,    // copy number
                              fCheckOverlaps); // checking overlaps
            
            //////////////////////////////////////////////////////
            //     PLACING OBJECTS within each SENSE REGION
            for(G4int j=0; j<2; j++)
            {
                if(j==0) usds = 0;
                if(j==1) usds = 1;
                
                if(j==1)
                {
                    new G4PVPlacement(0,    // no rotation
                                      offset_VDC_XU_PCBFrame[j],
                                      Logic_VDC_XU_PCBFrame,
                                      "VDC_XU_PCBFrame",
                                      Logic_VDC_SenseRegion_USDS[i][usds],
                                      false,    // no boolean operations
                                      i*2 + j,    // copy number
                                      fCheckOverlaps); // checking overlaps
                    
                }
                
                /////////////////////////////////////////////////
                //      VDC - STESALIT X-WIRE/U-WIRE FRAME
                new G4PVPlacement(0,    // no rotation
                                  offset_VDC_XU_Frame[j],
                                  Logic_VDC_XU_Frame,
                                  "VDC_XU_Frame",
                                  Logic_VDC_SenseRegion_USDS[i][usds],
                                  false,    // no boolean operations
                                  i*2 + j,    // copy number
                                  fCheckOverlaps); // checking overlaps
                
                ///////////////////////////////////////////
                //      VDC - PCB X-WIRE/U-WIRE FRAME
                new G4PVPlacement(0,    // no rotation
                                  offset_VDC_XU_PCBFrame[j],
                                  Logic_VDC_XU_PCBFrame,
                                  "VDC_XU_PCBFrame",
                                  Logic_VDC_SenseRegion_USDS[i][usds],
                                  false,    // no boolean operations
                                  i*2 + j,    // copy number
                                  fCheckOverlaps); // checking overlaps
                
                if(j==0)
                {
                    //////////////////////////////////////////////
                    //      VDC - U WIRES
                    for(G4int k=0; k<143; k++)
                    {
                        //offset_VDC_U_WIRE = G4ThreeVector( (k - 71)*4*tan(40*deg)*mm, 0.*mm, (4000. - 20./2)*um);
                        offset_VDC_U_WIRE = G4ThreeVector( (k - 71)*(4/sin(50.*deg))*mm, 0.*mm, (-4000. - 20./2)*um);
                        
                        VDC_U_WIRE_transform = G4Transform3D(VDC_U_WIRE_rotm, offset_VDC_U_WIRE);
                        
                        PhysiVDC_U_WIRE = new G4PVPlacement(VDC_U_WIRE_transform,
                                                            Logic_VDC_U_WIRE[k],
                                                            "VDC_U_WIRE",
                                                            Logic_VDC_SenseRegion_USDS[i][j],
                                                            false,    // no boolean operations
                                                            i*143 + k,    // copy number
                                                            fCheckOverlaps); // checking overlaps
                        
                    }
                    
                    //////////////////////////////////////////////
                    //      VDC - U GUARD WIRES
                    for(G4int k=0; k<144; k++)
                    {
                        offset_VDC_U_WIRE = G4ThreeVector( (k - 71.5)*(4/sin(50.*deg))*mm, 0.*mm, (-4000. - 20./2)*um);
                        VDC_U_WIRE_transform = G4Transform3D(VDC_U_WIRE_rotm, offset_VDC_U_WIRE);
                        
                        new G4PVPlacement(VDC_U_WIRE_transform,
                                          Logic_VDC_U_GUARDWIRE[k],
                                          "VDC_U_GUARDWIRE",
                                          Logic_VDC_SenseRegion_USDS[i][j],
                                          false,    // no boolean operations
                                          i*144 + k,    // copy number
                                          fCheckOverlaps); // checking overlaps
                        
                    }
                    
                    //////////////////////////////////////////////
                    //      VDC - U GUARD WIRES, Thick
                    for(G4int k=0; k<2; k++)
                    {
                        offset_VDC_U_WIRE = G4ThreeVector( (k*144. - 72.)*(4/sin(50.*deg))*mm, 0.*mm, (-4000. - 20./2)*um);
                        VDC_U_WIRE_transform = G4Transform3D(VDC_U_WIRE_rotm, offset_VDC_U_WIRE);
                        
                        new G4PVPlacement(VDC_U_WIRE_transform,
                                          Logic_VDC_U_GUARDWIRE_Thick[k],
                                          "VDC_U_GUARDWIRE_Thick",
                                          Logic_VDC_SenseRegion_USDS[i][j],
                                          false,    // no boolean operations
                                          i*2 + k,    // copy number
                                          fCheckOverlaps);
                        
                    }
                }
                
                if(j==1)
                {
                    //////////////////////////////////////////////
                    //      VDC - X WIRES
                    for(G4int k=0; k<198; k++)
                    {
                        offset_VDC_X_WIRE = G4ThreeVector( ( (k*4.) - 394.)*mm, 0.*mm, (4000. - 20./2)*um);
                        VDC_X_WIRE_transform = G4Transform3D(VDC_X_WIRE_rotm, offset_VDC_X_WIRE);
                        
                        PhysiVDC_X_WIRE = new G4PVPlacement(VDC_X_WIRE_transform,
                                                            Logic_VDC_X_WIRE,
                                                            "VDC_X_WIRE",
                                                            Logic_VDC_SenseRegion_USDS[i][j],
                                                            false,    // no boolean operations
                                                            i*198 + k,    // copy number
                                                            fCheckOverlaps); // checking overlaps
                        
                    }
                    
                    //////////////////////////////////////////////
                    //      VDC - X GUARD WIRES
                    for(G4int k=0; k<199; k++)
                    {
                        offset_VDC_X_WIRE = G4ThreeVector( ( (k*4.) - 396.)*mm, 0.*mm, (4000. - 20./2)*um);
                        VDC_X_WIRE_transform = G4Transform3D(VDC_X_WIRE_rotm, offset_VDC_X_WIRE);
                        
                        new G4PVPlacement(VDC_X_WIRE_transform,
                                          Logic_VDC_X_GUARDWIRE,
                                          "VDC_X_GUARDWIRE",
                                          Logic_VDC_SenseRegion_USDS[i][j],
                                          false,    // no boolean operations
                                          i*199 + k,    // copy number
                                          fCheckOverlaps); // checking overlaps
                        
                    }
                    
                    //////////////////////////////////////////////
                    //      VDC - X GUARD WIRES, Thick
                    for(G4int k=0; k<2; k++)
                    {
                        offset_VDC_X_WIRE = G4ThreeVector( ( (k*2.*398.) - 398.)*mm, 0.*mm, (4000. - 20./2)*um);
                        VDC_X_WIRE_transform = G4Transform3D(VDC_X_WIRE_rotm, offset_VDC_X_WIRE);
                        
                        new G4PVPlacement(VDC_X_WIRE_transform,
                                          Logic_VDC_X_GUARDWIRE_Thick,
                                          "VDC_X_GUARDWIRE_Thick",
                                          Logic_VDC_SenseRegion_USDS[i][j],
                                          false,    // no boolean operations
                                          i*2 + k,    // copy number
                                          fCheckOverlaps);
                        
                    }
                }
            }
        }
    }
    
    
    for(G4int i=0; i<2; i++)
    {
        if(VDC_Presence[i])
        {
            VDC_transform[i] = G4Transform3D(VDC_rotm[i],VDC_CentrePosition[i]);
            
            for(G4int j=0; j<2; j++)
            {
                
                ////////////////////////////////////////////////////
                //      VDC - US and DS components of ASSEMBLY
                
                Physical_VDC_SenseRegion_USDS = new G4PVPlacement(0,
                                                                  offset_VDC_SenseRegion_USDS[j],
                                                                  Logic_VDC_SenseRegion_USDS[i][j],
                                                                  "VDC_SenseRegion_USDS",
                                                                  Logic_VDC_ASSEMBLY[i],
                                                                  false,           // no boolean
                                                                  i*2 + j,               // copy number
                                                                  fCheckOverlaps); // checking overlaps
                
                ////////////////////////////////////////
                //      VDC - ALUMINIUM FRAME
                VDC_Al_Frame_transform = G4Transform3D(VDC_Al_Frame_rotm[j], offset_VDC_Al_Frame[j]);
                
                new G4PVPlacement(VDC_Al_Frame_transform,
                                  Logic_VDC_Al_Frame,
                                  "VDC_Al_Frame",
                                  Logic_VDC_ASSEMBLY[i],
                                  false,    // no boolean operations
                                  i*2 + j,    // copy number
                                  fCheckOverlaps); // checking overlaps
                
                ////////////////////////////////////////
                //      VDC - MYLAR WINDOW
                
                new G4PVPlacement(0,
                                  offset_VDC_MYLAR_Plane[j],
                                  Logic_VDC_MYLAR_Plane,
                                  "VDC_MYLAR_Plane",
                                  Logic_VDC_ASSEMBLY[i],
                                  false,    // no boolean operations
                                  i*2 + j,    // copy number
                                  fCheckOverlaps); // checking overlaps
                
            }
            
            ////////////////////////////////
            //      VDC ASSEMBLY
            
            Physical_VDC_ASSEMBLY = new G4PVPlacement(VDC_transform[i], // transformation
                                                      Logic_VDC_ASSEMBLY[i],
                                                      "VDC_ASSEMBLY",
                                                      LogicWorld,
                                                      false,           // no boolean
                                                      i,               // copy number
                                                      fCheckOverlaps); // checking overlaps
            
        }
    }
    
    
    
    
    
    
    
    //////////////////////////////////
    //      PADDLE DEFINITION       //
    //////////////////////////////////
    
    G4Box* Solid_PADDLE_1 = new G4Box("Solid_PADDLE_1", (122/2)*cm, (10.2/2)*cm, (0.3175/2)*cm);
    G4Box* Solid_PADDLE_2 = new G4Box("Solid_PADDLE_2", (122/2)*cm, (10.2/2)*cm, (0.6350/2)*cm);
    G4Box* Solid_PADDLE_3 = new G4Box("Solid_PADDLE_3", (122/2)*cm, (10.2/2)*cm, (1.2700/2)*cm);
    
    G4LogicalVolume * Logic_PADDLE[numberOf_PADDLE];
    
    Logic_PADDLE[0] = new G4LogicalVolume(Solid_PADDLE_1, BC408_Material,"PADDLE_1",0,0,0);
    Logic_PADDLE[1] = new G4LogicalVolume(Solid_PADDLE_2, BC408_Material,"PADDLE_2",0,0,0);
    Logic_PADDLE[2] = new G4LogicalVolume(Solid_PADDLE_3, BC408_Material,"PADDLE_3",0,0,0);
    
    
    //////////////////////////////////////
    //      PADDLE INITIALISATION       //
    //////////////////////////////////////
    
    
    for(G4int i=0; i<3; i++)
    {
        if(PADDLE_Presence[i])
        {
            PADDLE_transform[i] = G4Transform3D(PADDLE_rotm[i],PADDLE_CentrePosition[i]);
            
            PhysiPADDLE = new G4PVPlacement(PADDLE_transform[i],
                                            Logic_PADDLE[i],    // its logical volume
                                            "PADDLE",           // its name
                                            LogicWorld,         // its mother  volume
                                            false,              // no boolean operations
                                            i,                  // copy number
                                            fCheckOverlaps);    // checking overlaps
            
        }
    }
    
    
    
    //////////////////////////////////
    //      HAGAR DEFINITION        //
    //////////////////////////////////
    
    ////////////////////////////////
    //      HAGAR - NaI Crystal
    G4Tubs* Solid_HAGAR_NaICrystal = new G4Tubs("HAGAR_NaICrystal", 0.*cm, (23.8/2)*cm, (35.6/2)*cm, 0.*deg, 360.*deg);
    
    G4LogicalVolume* Logic_HAGAR_NaICrystal = new G4LogicalVolume(Solid_HAGAR_NaICrystal, G4_SODIUM_IODIDE_Material, "Logic_HAGAR_NaICrystal", 0, 0, 0);
    
    ////////////////////////////////
    //      HAGAR - Annulus
    G4Tubs* Solid_HAGAR_Anulus = new G4Tubs("HAGAR_Anulus", (28.86/2)*cm, ((28.86/2) + 9.84)*cm, (61/2)*cm, 0.*deg, 360.*deg);
    
    G4LogicalVolume* Logic_HAGAR_Annulus = new G4LogicalVolume(Solid_HAGAR_Anulus, BC408_Material, "Logic_HAGAR_Annulus", 0, 0, 0);
    
    ////////////////////////////////
    //      HAGAR - Front Disc
    G4Tubs* Solid_HAGAR_FrontDisc = new G4Tubs("HAGAR_FrontDisc", 0.*cm, (48.58/2)*cm, (8/2)*cm, 0.*deg, 360.*deg);
    
    G4LogicalVolume* Logic_HAGAR_FrontDisc = new G4LogicalVolume(Solid_HAGAR_FrontDisc, BC408_Material, "Logic_HAGAR_FrontDisc", 0, 0, 0);
    
    
    
    ////////////////////////////////////////////////////
    //               HAGAR INITIALIZATION             //
    ////////////////////////////////////////////////////
    
    
    ///////////////////////////////
    //      HAGAR - NaI Crystal
    if(HAGAR_NaICrystal_Presence)
    {
        HAGAR_transform = G4Transform3D(HAGAR_rotm, HAGAR_NaICrystal_CentrePosition);
        
        PhysiHAGAR_NaICrystal = new G4PVPlacement(HAGAR_transform,   // transformation matrix
                                                  Logic_HAGAR_NaICrystal,       // its logical volume
                                                  "HAGAR_NaICrystal",       // its name
                                                  LogicWorld,         // its mother  volume
                                                  false,           // no boolean operations
                                                  0,               // copy number
                                                  fCheckOverlaps); // checking overlaps
    }
    
    /////////////////////////////
    //      HAGAR - Annulus
    if(HAGAR_Annulus_Presence)
    {
        HAGAR_transform = G4Transform3D(HAGAR_rotm, HAGAR_Annulus_CentrePosition);
        
        PhysiHAGAR_Annulus = new G4PVPlacement(HAGAR_transform,   // transformation matrix
                                               Logic_HAGAR_Annulus,       // its logical volume
                                               "HAGAR_Annulus",       // its name
                                               LogicWorld,         // its mother  volume
                                               false,           // no boolean operations
                                               0,               // copy number
                                               fCheckOverlaps); // checking overlaps
    }
    
    ///////////////////////////////
    //      HAGAR - Front Disc
    if(HAGAR_FrontDisc_Presence)
    {
        HAGAR_transform = G4Transform3D(HAGAR_rotm, HAGAR_FrontDisc_CentrePosition);
        
        PhysiHAGAR_FrontDisc = new G4PVPlacement(HAGAR_transform,   // transformation matrix
                                                 Logic_HAGAR_FrontDisc,       // its logical volume
                                                 "HAGAR_FrontDisc",       // its name
                                                 LogicWorld,         // its mother  volume
                                                 false,           // no boolean operations
                                                 0,               // copy number
                                                 fCheckOverlaps); // checking overlaps
    }
    
    
    
    
    
    
    
    /////////////////////////////////////////////
    //             CLOVER DEFINITION           //
    /////////////////////////////////////////////
    
    G4double CLOVERtoShield_displacement = 10;  // cm
    
    G4ThreeVector offset_CLOVERInternalVacuum = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    G4ThreeVector offset_CLOVEREncasement = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    G4ThreeVector offset_CLOVERHPGeCrystal1 = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    G4ThreeVector offset_CLOVERHPGeCrystal2 = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    G4ThreeVector offset_CLOVERHPGeCrystal3 = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    G4ThreeVector offset_CLOVERHPGeCrystal4 = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    
    G4LogicalVolume * Logic_CLOVER_InternalVacuum[numberOf_CLOVER];
    G4LogicalVolume * Logic_CLOVER_Encasement;
    G4LogicalVolume * Logic_CLOVER_HPGeCrystal[4];
    
    if( CLOVER_Presence[0] || CLOVER_Presence[1] || CLOVER_Presence[2] || CLOVER_Presence[3] || CLOVER_Presence[4] || CLOVER_Presence[5] || CLOVER_Presence[6] || CLOVER_Presence[7] )
    {
        //////////////////////////////////////////////////////////
        //              CLOVER Internal Vacuum - CADMesh
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/CLOVER-InternalVacuum/CloverInternalVacuum.ply");

        CADMesh * mesh_CLOVERInternalVacuum = new CADMesh(meshPath, meshType, mm, offset_CLOVERInternalVacuum, false);
        
        G4VSolid * Solid_CLOVERInternalVacuum = mesh_CLOVERInternalVacuum->TessellatedMesh();
        
        for(G4int i=0; i<numberOf_CLOVER; i++)
        {
            Logic_CLOVER_InternalVacuum[i] = new G4LogicalVolume(Solid_CLOVERInternalVacuum, G4_Galactic_Material, "LogicCLOVERInternalVacuum", 0, 0, 0);
        }
        
        
        ///////////////////////////////////////////////////////
        //              CLOVER Encasement - CADMesh
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Clover-Encasement/CloverEncasement.ply");

        CADMesh * mesh_CLOVEREncasement = new CADMesh(meshPath, meshType, mm, offset_CLOVEREncasement, false);
        
        G4VSolid * Solid_CLOVEREncasement = mesh_CLOVEREncasement->TessellatedMesh();
        
        Logic_CLOVER_Encasement = new G4LogicalVolume(Solid_CLOVEREncasement, G4_Al_Material, "LogicCLOVERCloverEncasement", 0, 0, 0);
        
        
        //////////////////////////////////////////////////////////
        //              CLOVER HPGeCrystals - CADMesh
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe-Crystal1.ply");
        CADMesh * mesh_CLOVERHPGeCrystal1 = new CADMesh(meshPath, meshType, mm, offset_CLOVERHPGeCrystal1, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe-Crystal2.ply");
        CADMesh * mesh_CLOVERHPGeCrystal2 = new CADMesh(meshPath, meshType, mm, offset_CLOVERHPGeCrystal2, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe-Crystal3.ply");
        CADMesh * mesh_CLOVERHPGeCrystal3 = new CADMesh(meshPath, meshType, mm, offset_CLOVERHPGeCrystal3, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe-Crystal4.ply");
        CADMesh * mesh_CLOVERHPGeCrystal4 = new CADMesh(meshPath, meshType, mm, offset_CLOVERHPGeCrystal4, false);
        
        G4VSolid * Solid_HPGeCrystal1 = mesh_CLOVERHPGeCrystal1->TessellatedMesh();
        G4VSolid * Solid_HPGeCrystal2 = mesh_CLOVERHPGeCrystal2->TessellatedMesh();
        G4VSolid * Solid_HPGeCrystal3 = mesh_CLOVERHPGeCrystal3->TessellatedMesh();
        G4VSolid * Solid_HPGeCrystal4 = mesh_CLOVERHPGeCrystal4->TessellatedMesh();
        
        Logic_CLOVER_HPGeCrystal[0] = new G4LogicalVolume(Solid_HPGeCrystal1, G4_Ge_Material,"LogicCLOVERHPGeCrystal",0,0,0);
        Logic_CLOVER_HPGeCrystal[1] = new G4LogicalVolume(Solid_HPGeCrystal2, G4_Ge_Material,"LogicCLOVERHPGeCrystal",0,0,0);
        Logic_CLOVER_HPGeCrystal[2] = new G4LogicalVolume(Solid_HPGeCrystal3, G4_Ge_Material,"LogicCLOVERHPGeCrystal",0,0,0);
        Logic_CLOVER_HPGeCrystal[3] = new G4LogicalVolume(Solid_HPGeCrystal4, G4_Ge_Material,"LogicCLOVERHPGeCrystal",0,0,0);
        
    }
    
    
    ////////////////////////////////////////////
    //              CLOVER SHIELD
    ////////////////////////////////////////////
    
    G4ThreeVector offset_CLOVER_Shield_Body = G4ThreeVector(0*cm, 0*cm, 0*cm);
    G4ThreeVector offset_CLOVER_Shield_Heavimet = G4ThreeVector(0*cm, 0*cm, 0*cm);
    G4ThreeVector offset_CLOVER_Shield_BGOCrystals = G4ThreeVector(0*cm, 0*cm, 0*cm);
    G4ThreeVector offset_CLOVER_Shield_PMT = G4ThreeVector(0*cm, 0*cm, 0*cm);
    G4ThreeVector offset_CLOVER_Shield_PMTConArray = G4ThreeVector(0*cm, 0*cm, 0*cm);
    
    G4LogicalVolume* Logic_CLOVER_Shield_Body;
    G4LogicalVolume* Logic_CLOVER_Shield_Heavimet;
    G4LogicalVolume* Logic_CLOVER_Shield_PMTConArray;
    G4LogicalVolume* Logic_CLOVER_Shield_BGOCrystal[16];
    G4LogicalVolume* Logic_CLOVER_Shield_PMT[16];
    
    if( CLOVER_Shield_Presence[0] || CLOVER_Shield_Presence[1] || CLOVER_Shield_Presence[2] || CLOVER_Shield_Presence[3] || CLOVER_Shield_Presence[4] || CLOVER_Shield_Presence[5] || CLOVER_Shield_Presence[6] || CLOVER_Shield_Presence[7] || CLOVER_Shield_Presence[8])
    {
        ///////////////////////////////////////////////////////
        //              CLOVER Shield Body - CADMesh
        ///////////////////////////////////////////////////////
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/Body/Body.ply");
        CADMesh * mesh_CLOVER_Shield_Body = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_Body, false);
        
        G4VSolid * Solid_CLOVER_Shield_Body = mesh_CLOVER_Shield_Body->TessellatedMesh();
        
        Logic_CLOVER_Shield_Body = new G4LogicalVolume(Solid_CLOVER_Shield_Body, G4_Al_Material, "LogicCLOVERShieldBody", 0, 0, 0);
        
        ///////////////////////////////////////////////////////
        //              CLOVER Shield Heavimet - CADMesh
        ///////////////////////////////////////////////////////
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/Heavimet-Shield/HeavimetShield.ply");
        CADMesh * mesh_CLOVER_Shield_Heavimet = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_Heavimet, false);
        
        G4VSolid * Solid_CLOVER_Shield_Heavimet = mesh_CLOVER_Shield_Heavimet->TessellatedMesh();
        
        Logic_CLOVER_Shield_Heavimet = new G4LogicalVolume(Solid_CLOVER_Shield_Heavimet, Heavimet_Material, "LogicCLOVERShieldHeavimet", 0, 0, 0);
        
        ///////////////////////////////////////////////////////
        //      CLOVER Shield PMT Connecter Array - CADMesh
        ///////////////////////////////////////////////////////
        
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/PMT-Connectors/PMT-ConnecterArray.ply");
        CADMesh * mesh_CLOVER_Shield_PMTConArray = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMTConArray, false);
        
        G4VSolid * Solid_CLOVER_Shield_PMTConArray = mesh_CLOVER_Shield_PMTConArray->TessellatedMesh();
        
        Logic_CLOVER_Shield_PMTConArray = new G4LogicalVolume(Solid_CLOVER_Shield_PMTConArray, G4_Al_Material, "LogicCLOVERShieldHeavimet", 0, 0, 0);
        
        
        ///////////////////////////////////////////////////////
        //      CLOVER Shield BGO Crystals - CADMesh
        ///////////////////////////////////////////////////////
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal1.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal1 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal2.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal2 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal3.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal3 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal4.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal4 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal5.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal5 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal6.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal6 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal7.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal7 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal8.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal8 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal9.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal9 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal10.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal10 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal11.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal11 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal12.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal12 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal13.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal13 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal14.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal14 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal15.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal15 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal16.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal16 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        G4VSolid * Solid_CLOVER_Shield_BGOCrystal[16];
        
        Solid_CLOVER_Shield_BGOCrystal[0] = mesh_CLOVER_Shield_BGOCrystal1->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[1] = mesh_CLOVER_Shield_BGOCrystal2->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[2] = mesh_CLOVER_Shield_BGOCrystal3->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[3] = mesh_CLOVER_Shield_BGOCrystal4->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[4] = mesh_CLOVER_Shield_BGOCrystal5->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[5] = mesh_CLOVER_Shield_BGOCrystal6->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[6] = mesh_CLOVER_Shield_BGOCrystal7->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[7] = mesh_CLOVER_Shield_BGOCrystal8->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[8] = mesh_CLOVER_Shield_BGOCrystal9->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[9] = mesh_CLOVER_Shield_BGOCrystal10->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[10] = mesh_CLOVER_Shield_BGOCrystal11->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[11] = mesh_CLOVER_Shield_BGOCrystal12->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[12] = mesh_CLOVER_Shield_BGOCrystal13->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[13] = mesh_CLOVER_Shield_BGOCrystal14->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[14] = mesh_CLOVER_Shield_BGOCrystal15->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[15] = mesh_CLOVER_Shield_BGOCrystal16->TessellatedMesh();
        
        
        for (G4int k=0; k<16; k++)
        {
            Logic_CLOVER_Shield_BGOCrystal[k] = new G4LogicalVolume(Solid_CLOVER_Shield_BGOCrystal[k], G4_BGO_Material,"LogicCLOVERShieldBGOCrystal",0,0,0);
        }
        
        ////////////////////////////////////
        //      CLOVER Shield PMT's
        ////////////////////////////////////
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT1.ply");
        CADMesh * mesh_CLOVER_Shield_PMT1 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT2.ply");
        CADMesh * mesh_CLOVER_Shield_PMT2 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT3.ply");
        CADMesh * mesh_CLOVER_Shield_PMT3 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT4.ply");
        CADMesh * mesh_CLOVER_Shield_PMT4 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT5.ply");
        CADMesh * mesh_CLOVER_Shield_PMT5 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT6.ply");
        CADMesh * mesh_CLOVER_Shield_PMT6 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT7.ply");
        CADMesh * mesh_CLOVER_Shield_PMT7 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT8.ply");
        CADMesh * mesh_CLOVER_Shield_PMT8 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT9.ply");
        CADMesh * mesh_CLOVER_Shield_PMT9 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT10.ply");
        CADMesh * mesh_CLOVER_Shield_PMT10 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT11.ply");
        CADMesh * mesh_CLOVER_Shield_PMT11 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT12.ply");
        CADMesh * mesh_CLOVER_Shield_PMT12 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT13.ply");
        CADMesh * mesh_CLOVER_Shield_PMT13 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT14.ply");
        CADMesh * mesh_CLOVER_Shield_PMT14 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT15.ply");
        CADMesh * mesh_CLOVER_Shield_PMT15 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT16.ply");
        CADMesh * mesh_CLOVER_Shield_PMT16 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        
        G4VSolid * Solid_CLOVER_Shield_PMT[16];
        
        Solid_CLOVER_Shield_PMT[0] = mesh_CLOVER_Shield_PMT1->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[1] = mesh_CLOVER_Shield_PMT2->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[2] = mesh_CLOVER_Shield_PMT3->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[3] = mesh_CLOVER_Shield_PMT4->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[4] = mesh_CLOVER_Shield_PMT5->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[5] = mesh_CLOVER_Shield_PMT6->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[6] = mesh_CLOVER_Shield_PMT7->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[7] = mesh_CLOVER_Shield_PMT8->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[8] = mesh_CLOVER_Shield_PMT9->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[9] = mesh_CLOVER_Shield_PMT10->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[10] = mesh_CLOVER_Shield_PMT11->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[11] = mesh_CLOVER_Shield_PMT12->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[12] = mesh_CLOVER_Shield_PMT13->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[13] = mesh_CLOVER_Shield_PMT14->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[14] = mesh_CLOVER_Shield_PMT15->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[15] = mesh_CLOVER_Shield_PMT16->TessellatedMesh();
        
        for (G4int k=0; k<16; k++)
        {
            Logic_CLOVER_Shield_PMT[k] = new G4LogicalVolume(Solid_CLOVER_Shield_PMT[k], G4_Al_Material,"LogicCLOVERShieldPMT",0,0,0);
        }
    }
    
    
    ////////////////////////////////////////////////////
    //               CLOVER INITIALIZATION            //
    ////////////////////////////////////////////////////
    
    
    for(G4int i=0; i<numberOf_CLOVER; i++)
    {
        CLOVER_position[i] = (CLOVER_Distance[i])*G4ThreeVector( sin(CLOVER_theta[i]) * cos(CLOVER_phi[i]), sin(CLOVER_theta[i]) * sin(CLOVER_phi[i]), cos(CLOVER_theta[i]));
        
        CLOVER_transform[i] = G4Transform3D(CLOVER_rotm[i],CLOVER_position[i]);
        
        CLOVER_Shield_position[i] = (CLOVER_Distance[i])*G4ThreeVector( sin(CLOVER_theta[i]) * cos(CLOVER_phi[i]), sin(CLOVER_theta[i]) * sin(CLOVER_phi[i]), cos(CLOVER_theta[i]));
        
        CLOVER_Shield_transform[i] = CLOVER_transform[i];
        
        /////////////////////////////
        //          CLOVER
        if(CLOVER_Presence[i])
        {
            new G4PVPlacement(CLOVER_transform[i],   // transformation matrix
                              Logic_CLOVER_Encasement,       // its logical volume
                              "CLOVER_Encasement",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            
            
            
            
            for (int j=0; j<4; j++)
            {
                PhysiCLOVER_HPGeCrystal = new G4PVPlacement(0,               // no rotation
                                                            G4ThreeVector(0,0,0), // at (x,y,z)
                                                            Logic_CLOVER_HPGeCrystal[j],
                                                            "CLOVER_HPGeCrystal", // its name
                                                            Logic_CLOVER_InternalVacuum[i],
                                                            false,           // no boolean operations
                                                            i*4 + j,               // copy number
                                                            fCheckOverlaps); // checking overlaps
                
            }
            
            
            new G4PVPlacement(CLOVER_transform[i],
                              Logic_CLOVER_InternalVacuum[i],
                              "CLOVER_InternalVacuum",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
        }
        
        /////////////////////////////
        //      CLOVER SHIELD
        if(CLOVER_Shield_Presence[i])
        {
            new G4PVPlacement(CLOVER_Shield_transform[i],
                              Logic_CLOVER_Shield_Body,
                              "CLOVER_Shield_Body",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            for (int j=0; j<16; j++)
            {
                PhysiCLOVER_Shield_BGOCrystal = new G4PVPlacement(CLOVER_Shield_transform[i],
                                                                  Logic_CLOVER_Shield_BGOCrystal[j],
                                                                  "CLOVER_Shield_BGOCrystal",
                                                                  LogicVacuumChamber,
                                                                  false,
                                                                  i*16 + j,
                                                                  fCheckOverlaps);
                
                PhysiCLOVER_Shield_PMT = new G4PVPlacement(CLOVER_Shield_transform[i],
                                                           Logic_CLOVER_Shield_PMT[j],
                                                           "CLOVER_Shield_PMT", // its name
                                                           LogicVacuumChamber,
                                                           false, // no boolean operations
                                                           i*16 + j,  // copy number
                                                           fCheckOverlaps); // checking overlaps
                
            }
            
            new G4PVPlacement(CLOVER_Shield_transform[i],
                              Logic_CLOVER_Shield_PMTConArray,
                              "CLOVER_Shield_PMTConArray",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            new G4PVPlacement(CLOVER_Shield_transform[i],
                              Logic_CLOVER_Shield_Heavimet,
                              "CLOVER_Shield_Heavimet",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
        }
    }
    
    
    ////////////////////////////////////////
    ////        LEPS DEFINITION         ////
    ////////////////////////////////////////
    
    //////////////////////////////////////////////////////////
    //              LEPS Internal Vacuum - CADMesh
    //////////////////////////////////////////////////////////
    
    G4Tubs* Solid_LEPS_InternalVacuum = new G4Tubs("Solid_LEPSInternalVacuum", 0.*mm, 38.4*mm, 45.0*mm, 0.*deg, 360*deg);
    G4LogicalVolume* Logic_LEPS_InternalVacuum[numberOf_LEPS];
    
    for(G4int i=0; i<numberOf_LEPS; i++)
    {
        Logic_LEPS_InternalVacuum[i] = new G4LogicalVolume(Solid_LEPS_InternalVacuum, G4_Galactic_Material, "LogicLEPSInternalVacuum", 0, 0, 0);
    }
    
    ///////////////////////////////////////////////////////
    //              LEPS Encasement - CADMesh
    ///////////////////////////////////////////////////////
    
    G4Tubs* Solid_LEPS_Encasement = new G4Tubs("Solid_LEPSEncasement", 38.5*mm, 40.0*mm, (90./2)*mm, 0.*deg, 360*deg);
    
    G4LogicalVolume* Logic_LEPS_Encasement = new G4LogicalVolume(Solid_LEPS_Encasement, G4_Al_Material, "LogicLEPSLEPSEncasement", 0, 0, 0);
    
    
    ///////////////////////////////////////////////////////
    //              LEPS Beryllium Window - CADMesh
    ///////////////////////////////////////////////////////
    
    G4Tubs* Solid_LEPS_Window = new G4Tubs("Solid_LEPSWindow", 0.*mm, 38.5*mm, (0.3/2)*mm, 0.*deg, 360*deg);
    
    G4LogicalVolume* Logic_LEPS_Window = new G4LogicalVolume(Solid_LEPS_Window, G4_Be_Material, "Logic_LEPS_Window", 0, 0, 0);
    
    
    //////////////////////////////////////////////////////////
    //              LEPS HPGeCrystals - CADMesh
    //////////////////////////////////////////////////////////
    
    G4Tubs* Solid_HPGeCrystal = new G4Tubs("Solid_HPGeCrystal1", 0.*mm, 33.0*mm, 5.5*mm, 0.*deg, 90.*deg);
    
    G4LogicalVolume* Logic_LEPS_HPGeCrystal = new G4LogicalVolume(Solid_HPGeCrystal, G4_Ge_Material,"LogicLEPSHPGeCrystal",0,0,0);
    
    LEPS_HPGeCrystal_rotm[0].rotateZ(0.*deg);
    LEPS_HPGeCrystal_rotm[1].rotateZ(90.*deg);
    LEPS_HPGeCrystal_rotm[2].rotateZ(180.*deg);
    LEPS_HPGeCrystal_rotm[3].rotateZ(270.*deg);
    
    for(G4int i=0; i<4; i++)
    {
        LEPS_HPGeCrystal_transform[i] = G4Transform3D(LEPS_HPGeCrystal_rotm[i], G4ThreeVector(0,0,(29.0-0.5)*mm));
    }
    
    ////////////////////////////////////////////////////
    //               LEPS INITIALIZATION
    ////////////////////////////////////////////////////
    
    
    for(G4int i=0; i<numberOf_LEPS; i++)
    {
        LEPS_position[i] = (LEPS_Distance[i] + 4.5*cm)*G4ThreeVector( std::sin(LEPS_theta[i]) * std::cos(LEPS_phi[i]), std::sin(LEPS_theta[i]) * std::sin(LEPS_phi[i]), std::cos(LEPS_theta[i]));
        
        LEPS_transform[i] = G4Transform3D(LEPS_rotm[i],LEPS_position[i]);
        
        LEPS_InternalVacuum_position[i] = (LEPS_Distance[i]+ 4.5*cm + 0.5*mm)*G4ThreeVector( std::sin(LEPS_theta[i]) * std::cos(LEPS_phi[i]), std::sin(LEPS_theta[i]) * std::sin(LEPS_phi[i]), std::cos(LEPS_theta[i]));
        LEPS_InternalVacuum_transform[i] = G4Transform3D(LEPS_rotm[i],LEPS_InternalVacuum_position[i]);
        
        LEPS_Window_position[i] = (LEPS_Distance[i] + 4.5*cm +(-45.0+0.15)*mm)*G4ThreeVector( std::sin(LEPS_theta[i]) * std::cos(LEPS_phi[i]), std::sin(LEPS_theta[i]) * std::sin(LEPS_phi[i]), std::cos(LEPS_theta[i]));
        LEPS_Window_transform[i] = G4Transform3D(LEPS_rotm[i],LEPS_Window_position[i]);
        
        /////////////////////////////
        //          LEPS
        if(LEPS_Presence[i] == true)
        {
            
            new G4PVPlacement(LEPS_transform[i],   // transformation matrix
                              Logic_LEPS_Encasement,       // its logical volume
                              "LEPSEncasement",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            new G4PVPlacement(LEPS_Window_transform[i],   // transformation matrix
                              Logic_LEPS_Window,       // its logical volume
                              "LEPSWindow",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            
            
            for(int j=0; j<4; j++)
            {
                Physical_LEPS_HPGeCrystal = new G4PVPlacement(LEPS_HPGeCrystal_transform[j],
                                                              Logic_LEPS_HPGeCrystal,       // its logical volume
                                                              "LEPSHPGeCrystal",       // its name
                                                              Logic_LEPS_InternalVacuum[i],    // its mother  volume
                                                              false,           // no boolean operations
                                                              j + (i*4),               // copy number
                                                              fCheckOverlaps); // checking overlaps
                
            }
            
            new G4PVPlacement(LEPS_InternalVacuum_transform[i],
                              Logic_LEPS_InternalVacuum[i],
                              "LEPSInternalVacuum",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            
        }
        
    }
    
    
    
    
    //////////////////////////////////////////////////
    //      K600 SPECTROMETER INITIIALIZATION       //
    //////////////////////////////////////////////////
    
    
    // A field object is held by a field manager
    // Find the global Field Manager
    G4TransportationManager* tmanagerMagneticField = G4TransportationManager::GetTransportationManager();
    tmanagerMagneticField->GetPropagatorInField()->SetLargestAcceptableStep(1*mm);
    G4double minStepMagneticField = 0.0025*mm ;
    
    
    //////////////////////////////////////////////////////
    //              K600 - QUADRUPOLE
    //////////////////////////////////////////////////////
    
    if(K600_Quadrupole)
    {
        
        G4RotationMatrix* K600_Q_MagField_rotm = new G4RotationMatrix;
        //K600_Q_MagField_rotm->rotateX(90.*deg);
        
        K600_Quadrupole_transform = G4Transform3D(K600_Quadrupole_rotm, K600_Quadrupole_CentrePosition);
        
        G4Box* Solid_K600_Quadrupole = new G4Box("Solid_K600_Quadrupole", (50./2)*cm, (50./2)*cm, (30./2)*cm);
        
        G4LogicalVolume* Logic_K600_Quadrupole = new G4LogicalVolume(Solid_K600_Quadrupole, G4_Galactic_Material,"Logic_K600_Quadrupole",0,0,0);
        
        ////    IDEAL MAGNETIC FIELD for QUADRUPOLE
        if(Ideal_Quadrupole)
        {
            MagneticField_K600_Q = new G4QuadrupoleMagField(K600_Q_gradient, K600_Quadrupole_CentrePosition, K600_Q_MagField_rotm);
            fEquationMagneticField_K600_Q = new G4Mag_UsualEqRhs(MagneticField_K600_Q);
            
            fieldManagerMagneticField_K600_Q = new G4FieldManager(MagneticField_K600_Q);
            
            stepperMagneticField_K600_Q = new G4ClassicalRK4( fEquationMagneticField_K600_Q );
            fieldManagerMagneticField_K600_Q -> SetDetectorField(MagneticField_K600_Q);
            
            fChordFinder_K600_Q = new G4ChordFinder( MagneticField_K600_Q, minStepMagneticField, stepperMagneticField_K600_Q);
            
            Logic_K600_Quadrupole -> SetFieldManager(fieldManagerMagneticField_K600_Q, true) ;
            
            
            //G4BlineTracer* theBlineTool = new G4BlineTracer();
            //theBlineTool->ComputeBlines();
        }
        
        ////    MAPPED MAGNETIC FIELD for QUADRUPOLE
        if(Mapped_Quadrupole)
        {
            
            G4double z_Q_Offset = 4.4*mm+ 100*cm;
            
            G4MagneticField* PurgMagField = new MagneticFieldMapping("../K600/MagneticFieldMaps/Quadrupole_MagneticFieldMap.TABLE", z_Q_Offset);
            fEquationMagneticField_K600_Q = new G4Mag_UsualEqRhs(PurgMagField);
            
            fieldManagerMagneticField_K600_Q = new G4FieldManager(PurgMagField);
            
            stepperMagneticField_K600_Q = new G4ClassicalRK4( fEquationMagneticField_K600_Q );
            fieldManagerMagneticField_K600_Q -> SetDetectorField(PurgMagField);
            
            fChordFinder_K600_Q = new G4ChordFinder( PurgMagField, minStepMagneticField, stepperMagneticField_K600_Q);
            
            Logic_K600_Quadrupole -> SetFieldManager(fieldManagerMagneticField_K600_Q, true) ;
            
            
            //G4BlineTracer* theBlineTool = new G4BlineTracer();
            
        }
        
        PhysiK600_Quadrupole = new G4PVPlacement(K600_Quadrupole_transform,
                                                 Logic_K600_Quadrupole,       // its logical volume
                                                 "K600_Quadrupole",       // its name
                                                 LogicWorld,         // its mother  volume
                                                 false,           // no boolean operations
                                                 0,               // copy number
                                                 fCheckOverlaps); // checking overlaps
        
    }
    
    //////////////////////////////////////////////////////
    //              K600 - DIPOLE 1
    //////////////////////////////////////////////////////
    
    /*
     // magnetic field ----------------------------------------------------------
     MagneticField_K600_D1 = new G4UniformMagField(G4ThreeVector(0., K600_Dipole1_BZ, 0.));
     fieldManagerMagneticField_K600_D1 = new G4FieldManager();
     fieldManagerMagneticField_K600_D1->SetDetectorField(MagneticField_K600_D1);
     fieldManagerMagneticField_K600_D1->CreateChordFinder(MagneticField_K600_D1);
     G4bool forceToAllDaughters = true;
     */
    
    if(K600_Dipole1)
    {
        ////    MAGNETIC FIELD for DIPOLE 1
        MagneticField_K600_D1 = new G4UniformMagField(G4ThreeVector(0., K600_Dipole1_BZ, 0.));
        fEquationMagneticField_K600_D1 = new G4Mag_UsualEqRhs(MagneticField_K600_D1);
        
        fieldManagerMagneticField_K600_D1 = new G4FieldManager(MagneticField_K600_D1);
        
        stepperMagneticField_K600_D1 = new G4ClassicalRK4( fEquationMagneticField_K600_D1 );
        fieldManagerMagneticField_K600_D1 -> SetDetectorField(MagneticField_K600_D1);
        
        fChordFinder_K600_D1 = new G4ChordFinder( MagneticField_K600_D1, minStepMagneticField, stepperMagneticField_K600_D1);
        
        
        /////////////////////////////////////////////
        
        K600_Dipole1_transform = G4Transform3D(K600_Dipole1_rotm, K600_Dipole1_CentrePosition);
        
        //G4Box* Solid_K600_Dipole1 = new G4Box("Solid_K600_Dipole1", (50./2)*cm, (50./2)*cm, (30./2)*cm);
        //G4Tubs* Solid_K600_Dipole1 = new G4Tubs("Solid_K600_Dipole1", 50.*cm, 100.0*cm, 30.*cm, 0.*deg, 40.*deg);
        G4Tubs* Solid_K600_Dipole1 = new G4Tubs("Solid_K600_Dipole1", 30.*cm, 150.0*cm, 30.*cm, 0.*deg, 40.*deg);
        
        G4LogicalVolume* Logic_K600_Dipole1 = new G4LogicalVolume(Solid_K600_Dipole1, G4_Galactic_Material,"Logic_K600_Dipole1",0,0,0);
        Logic_K600_Dipole1 -> SetFieldManager(fieldManagerMagneticField_K600_D1, true) ;
        
        // Register the field and its manager for deleting
        //G4AutoDelete::Register(MagneticField_K600_D1);
        //G4AutoDelete::Register(fieldManagerMagneticField_K600_D1);
        
        PhysiK600_Dipole1 = new G4PVPlacement(K600_Dipole1_transform,
                                              Logic_K600_Dipole1,       // its logical volume
                                              "K600_Dipole1",       // its name
                                              LogicWorld,         // its mother  volume
                                              false,           // no boolean operations
                                              0,               // copy number
                                              fCheckOverlaps); // checking overlaps
        
    }
    
    
    
    //////////////////////////////////////////////////////
    //              K600 - DIPOLE 2
    //////////////////////////////////////////////////////
    
    if(K600_Dipole2)
    {
        
        ////    MAGNETIC FIELD for DIPOLE 2
        MagneticField_K600_D2 = new G4UniformMagField(G4ThreeVector(0., K600_Dipole2_BZ, 0.));
        fEquationMagneticField_K600_D2 = new G4Mag_UsualEqRhs(MagneticField_K600_D2);
        
        fieldManagerMagneticField_K600_D2 = new G4FieldManager(MagneticField_K600_D2);
        
        stepperMagneticField_K600_D2 = new G4ClassicalRK4( fEquationMagneticField_K600_D2 );
        fieldManagerMagneticField_K600_D2 -> SetDetectorField(MagneticField_K600_D2);
        
        fChordFinder_K600_D2 = new G4ChordFinder( MagneticField_K600_D2, minStepMagneticField, stepperMagneticField_K600_D2);
        
        
        /////////////////////////////////////////////
        
        K600_Dipole2_transform = G4Transform3D(K600_Dipole2_rotm, K600_Dipole2_CentrePosition);
        
        //G4Tubs* Solid_K600_Dipole2 = new G4Tubs("Solid_K600_Dipole2", 50.*cm, 100.0*cm, 30.*cm, 50.*deg, 70.*deg);
        G4Tubs* Solid_K600_Dipole2 = new G4Tubs("Solid_K600_Dipole2", 30.*cm, 150.0*cm, 30.*cm, 50.*deg, 70.*deg);
        
        
        G4LogicalVolume* Logic_K600_Dipole2 = new G4LogicalVolume(Solid_K600_Dipole2, G4_Galactic_Material,"Logic_K600_Dipole2",0,0,0);
        Logic_K600_Dipole2 -> SetFieldManager(fieldManagerMagneticField_K600_D2, true) ;
        
        PhysiK600_Dipole2 = new G4PVPlacement(K600_Dipole2_transform,
                                              Logic_K600_Dipole2,       // its logical volume
                                              "K600_Dipole2",       // its name
                                              LogicWorld,         // its mother  volume
                                              false,           // no boolean operations
                                              0,               // copy number
                                              fCheckOverlaps); // checking overlaps
        
        
    }
    
    
    
    
    //////////////////////////////////////////////////////////
    //                      VISUALISATION
    //////////////////////////////////////////////////////////
    
    G4VisAttributes* World_VisAtt= new G4VisAttributes(G4Colour(0., 0., 0.));
    World_VisAtt->SetVisibility(true);
    LogicWorld->SetVisAttributes(World_VisAtt);
    
    //G4VisAttributes* Target_VisAtt= new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    //LogicTarget->SetVisAttributes(Target_VisAtt);
    
    
    ////////////////////////////
    //  Vacuum Chamber
    
    G4VisAttributes* VacuumChamber_VisAtt= new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    //VacuumChamber_VisAtt->SetVisibility(false);
    LogicVacuumChamber->SetVisAttributes(VacuumChamber_VisAtt);
    
    
    ////////////////////////////////
    //      VDC VISUALIZATION
    ////////////////////////////////
    
    //  VDC - ASSEMBLY
    G4VisAttributes* VDC_ASSEMBLY_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
    VDC_ASSEMBLY_VisAtt->SetVisibility(false);
    
    //  VDC ASSEMBLY USDS
    G4VisAttributes* VDC_ASSEMBLY_USDS_VisAtt = new G4VisAttributes(G4Colour(0.0, 0.7, 0.0));
    VDC_ASSEMBLY_USDS_VisAtt->SetVisibility(false);
    
    //  VDC - STESALIT HV Frames, Upstream, Middle and Downstream
    G4VisAttributes* VDC_HV_Frame_VisAtt = new G4VisAttributes(G4Colour(0.8, 0.5, 0.35));
    VDC_HV_Frame_VisAtt->SetForceSolid(true);
    //VDC_HV_Frame_VisAtt->SetVisibility(false);
    
    //  VDC - HIGH VOLTAGE PLANE
    G4VisAttributes* VDC_HV_Plane_VisAtt = new G4VisAttributes(G4Colour(0.8, 0.5, 0.35));
    //VDC_HV_Plane_VisAtt->SetForceSolid(true);
    //VDC_HV_Plane_VisAtt->SetVisibility(false);
    
    //  VDC - GAS FRAME
    G4VisAttributes* VDC_GasFrame_VisAtt = new G4VisAttributes(G4Colour(0.8, 0.5, 0.35));
    VDC_GasFrame_VisAtt->SetForceSolid(true);
    //VDC_GasFrame_VisAtt->SetVisibility(false);
    
    //  VDC - STESALIT X-WIRE/U-WIRE FRAME
    G4VisAttributes* VDC_XU_Frame_VisAtt = new G4VisAttributes(G4Colour(0.8, 0.5, 0.35));
    VDC_XU_Frame_VisAtt->SetForceSolid(true);
    //VDC_XU_Frame_VisAtt->SetVisibility(false);
    
    //  VDC - PCB X-WIRE/U-WIRE FRAME
    G4VisAttributes* VDC_XU_PCBFrame_VisAtt = new G4VisAttributes(G4Colour(0.8, 0.5, 0.35));
    VDC_XU_PCBFrame_VisAtt->SetForceSolid(true);
    //VDC_XU_PCBFrame_VisAtt->SetVisibility(false);
    
    //  VDC - ALUMINIUM OUTER FRAME
    G4VisAttributes* VDC_Al_Frame_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    VDC_Al_Frame_VisAtt->SetForceSolid(true);
    
    //  VDC - MYLAR PLANE
    G4VisAttributes* VDC_MYLAR_Plane_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    //VDC_MYLAR_Plane_VisAtt->SetForceSolid(true);
    
    //  VDC - X WIRES
    G4VisAttributes* VDC_X_WIRE_VisAtt = new G4VisAttributes(G4Colour(0., 0.7, 0.7));
    VDC_X_WIRE_VisAtt->SetForceSolid(true);
    
    //  VDC - X GUARD WIRES
    G4VisAttributes* VDC_X_GUARDWIRE_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    VDC_X_GUARDWIRE_VisAtt->SetVisibility(false);
    VDC_X_GUARDWIRE_VisAtt->SetForceSolid(true);
    
    //  VDC - U WIRES
    G4VisAttributes* VDC_U_WIRE_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.));
    VDC_U_WIRE_VisAtt->SetForceSolid(true);
    
    //  VDC - U GUARD WIRES
    G4VisAttributes* VDC_U_GUARDWIRE_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    VDC_U_GUARDWIRE_VisAtt->SetVisibility(false);
    VDC_U_GUARDWIRE_VisAtt->SetForceSolid(true);
    
    Logic_VDC_GasFrame->SetVisAttributes(VDC_GasFrame_VisAtt);
    Logic_VDC_XU_Frame->SetVisAttributes(VDC_XU_Frame_VisAtt);
    Logic_VDC_XU_PCBFrame->SetVisAttributes(VDC_XU_PCBFrame_VisAtt);
    Logic_VDC_Al_Frame->SetVisAttributes(VDC_Al_Frame_VisAtt);
    Logic_VDC_MYLAR_Plane->SetVisAttributes(VDC_MYLAR_Plane_VisAtt);
    Logic_VDC_X_WIRE->SetVisAttributes(VDC_X_WIRE_VisAtt);
    Logic_VDC_X_GUARDWIRE->SetVisAttributes(VDC_X_GUARDWIRE_VisAtt);
    Logic_VDC_X_GUARDWIRE_Thick->SetVisAttributes(VDC_X_GUARDWIRE_VisAtt);
    
    for(G4int k=0; k<143; k++)
    {
        Logic_VDC_U_WIRE[k]->SetVisAttributes(VDC_U_WIRE_VisAtt);
    }
    
    for(G4int k=0; k<144; k++)
    {
        Logic_VDC_U_GUARDWIRE[k]->SetVisAttributes(VDC_U_GUARDWIRE_VisAtt);
    }
    
    for(G4int k=0; k<2; k++)
    {
        Logic_VDC_U_GUARDWIRE_Thick[k]->SetVisAttributes(VDC_U_GUARDWIRE_VisAtt);
    }
    
    for(G4int k=0; k<3; k++)
    {
        Logic_VDC_HV_Frame[k]->SetVisAttributes(VDC_HV_Frame_VisAtt);
        Logic_VDC_HV_Plane[k]->SetVisAttributes(VDC_HV_Plane_VisAtt);
    }
    
    for(G4int i=0; i<2; i++)
    {
        for(G4int j=0; j<2; j++)
        {
            Logic_VDC_SenseRegion_USDS[i][j]->SetVisAttributes(VDC_ASSEMBLY_USDS_VisAtt);
        }
    }
    
    for(G4int i=0; i<numberOf_VDC; i++)
    {
        Logic_VDC_ASSEMBLY[i]->SetVisAttributes(VDC_ASSEMBLY_VisAtt);
    }
    
    
    
    
    //////////////////////////////////
    //      PADDLE VISUALIZATION
    //////////////////////////////////
    
    G4VisAttributes* scintillatorVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
    for(G4int i=0; i<3; i++)
    {Logic_PADDLE[i] -> SetVisAttributes(scintillatorVisAtt);}
    
    
    
    //////////////////////////////////
    //      HAGAR VISUALIZATION
    //////////////////////////////////
    
    //  HAGAR - NaI Crystal
    G4VisAttributes* HAGAR_NaICrystal_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
    HAGAR_NaICrystal_VisAtt->SetForceSolid(true);
    Logic_HAGAR_NaICrystal->SetVisAttributes(HAGAR_NaICrystal_VisAtt);
    
    //  HAGAR - Annulus
    G4VisAttributes* HAGAR_Annulus_VisAtt = new G4VisAttributes(G4Colour(0.0, 0.7, 0.0));
    Logic_HAGAR_Annulus->SetVisAttributes(HAGAR_Annulus_VisAtt);
    
    //  HAGAR - Front Disc
    G4VisAttributes* HAGAR_FrontDisc_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
    HAGAR_FrontDisc_VisAtt->SetForceSolid(true);
    Logic_HAGAR_FrontDisc->SetVisAttributes(HAGAR_FrontDisc_VisAtt);
    
    
    /////////////////////////////////
    //      TIARA VISUALIZATION
    /////////////////////////////////
    
    //  ASSEMBLY
    G4VisAttributes* TIARA_Assembly_VisAtt = new G4VisAttributes(G4Colour(0.0, 0.7, 0.0));
    TIARA_Assembly_VisAtt->SetVisibility(false);
    //TIARA_Assembly_VisAtt->SetForceSolid(false);
    //TIARA_Assembly_VisAtt->SetForceWireframe(false);
    
    //  PCB
    G4VisAttributes* TIARA_PCB_VisAtt = new G4VisAttributes(G4Colour(0.0, 0.7, 0.0));
    //TIARA_PCB_VisAtt->SetForceSolid(true);
    
    //  AA - RS
    //G4VisAttributes* TIARA_AA_RS_VisAtt = new G4VisAttributes(G4Colour(0., 0.7, 0.7));
    G4VisAttributes* TIARA_AA_RS_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    //TIARA_AA_RS_VisAtt->SetForceSolid(true);
    TIARA_AA_RS_VisAtt->SetForceLineSegmentsPerCircle(100);
    
    //  Silicon Wafer
    G4VisAttributes* TIARA_SiliconWafer_VisAtt = new G4VisAttributes(G4Colour(0., 0.2, 0.2));
    //TIARA_SiliconWafer_VisAtt->SetForceSolid(true);
    TIARA_SiliconWafer_VisAtt->SetForceLineSegmentsPerCircle(100);
    
    
    //  DL
    G4VisAttributes* TIARA_DL_VisAtt = new G4VisAttributes(G4Colour(1., 0., 0.));
    //TIARA_DL_VisAtt->SetForceSolid(true);
    
    //  2M
    G4VisAttributes* TIARA_2M_VisAtt = new G4VisAttributes(G4Colour(1., 1., 1.));
    TIARA_2M_VisAtt->SetVisibility(false);
    //TIARA_2M_VisAtt->SetForceSolid(false);
    //TIARA_2M_VisAtt->SetForceWireframe(false);
    
    
    Logic_TIARA_PCB->SetVisAttributes(TIARA_PCB_VisAtt);
    //Logic_TIARA_DL_Front->SetVisAttributes(TIARA_DL_VisAtt);
    //Logic_TIARA_DL_Back->SetVisAttributes(TIARA_DL_VisAtt);
    //Logic_TIARA_2M_Front->SetVisAttributes(TIARA_2M_VisAtt);
    //Logic_TIARA_2M_Back->SetVisAttributes(TIARA_2M_VisAtt);
    
    for(G4int i=0; i<numberOf_TIARA; i++)
    {
        Logic_TIARA_Assembly[i]->SetVisAttributes(TIARA_Assembly_VisAtt);
        Logic_TIARA_SiliconWafer[i]->SetVisAttributes(TIARA_SiliconWafer_VisAtt);
        
    }
    
    for(G4int j=0; j<16; j++)
    {
        for(G4int l=0; l<8; l++)
        {
            Logic_TIARA_AA_RS[j][l]->SetVisAttributes(TIARA_AA_RS_VisAtt);
            Logic_TIARA_RS_punch[j][l]->SetVisAttributes(TIARA_AA_RS_VisAtt);
            
        }
    }
    
    //
    //always return the physical World
    //
    return PhysiWorld;
    
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructField()
{
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue = G4ThreeVector();
    //G4ThreeVector fieldValue = G4ThreeVector(0., 5*tesla, 0.);
    
    fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    fMagFieldMessenger->SetVerboseLevel(1);
    
    // Register the field messenger for deleting
    G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
