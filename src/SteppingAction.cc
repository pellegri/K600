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

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(const DetectorConstruction* detectorConstruction, EventAction* eventAction)
: G4UserSteppingAction(),
fDetConstruction(detectorConstruction),
fEventAction(eventAction)
{
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
    G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
    
    // get particle name/definition
    G4String particleName = aStep->GetTrack()->GetDefinition()->GetParticleName();
    //G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition();
    
    // get particle lifetime
    //G4double lifetime = aStep->GetTrack()->GetDefinition()->GetIonLifeTime()/ns;
    
    // get interaction time of the current step
    interactiontime = preStepPoint->GetGlobalTime()/ns;
    
    // get volume of the current step
    G4VPhysicalVolume* volume = theTouchable->GetVolume();
    
    // get volume name of the current step
    volumeName = volume->GetName();
    
    
    //G4cout << "Here is the particleName    "<< particleName << G4endl;
    //G4cout << "Here is the lifetime    "<< lifetime << G4endl;
    //G4cout << "Here is the interactiontime    "<< interactiontime << G4endl;
    //G4cout << "Here is the TEST    "<< TEST << G4endl;
    //G4cout << "Here is the interactiontime    "<< interactiontime << G4endl;
    //G4cout << "                                "<< G4endl;
    
    
    //G4ParticleDefinition* particleOI = G4Gamma::Gamma();
    
    
    ////////////////////////////////////////////
    //              TIARA ARRAY
    ////////////////////////////////////////////
    
    if(interactiontime < TIARA_TotalSampledTime && volumeName == "TIARA_AA_RS")
    {
        edepTIARA_AA = aStep->GetTotalEnergyDeposit()/MeV;
        
        if(edepTIARA_AA != 0.)
        {
            //G4cout << "Here we are in the Stepping Action" << G4endl;
            
            channelID = volume->GetCopyNo();
            
            TIARANo = channelID/128;
            TIARA_RowNo = (channelID - (TIARANo*128))/8;
            TIARA_SectorNo = (channelID - (TIARANo*128))%8;
            
            iTS = interactiontime/TIARA_SamplingTime;
            edepTIARA_AA = aStep->GetTotalEnergyDeposit()/MeV;
            
            if(fEventAction->GetVar_TIARA_AA(TIARANo, TIARA_RowNo, TIARA_SectorNo, 0, iTS)==0)
            {
                worldPosition = preStepPoint->GetPosition();
                
                xPosW = worldPosition.x()/m;
                yPosW = worldPosition.y()/m;
                zPosW = worldPosition.z()/m;
                
                normVector = pow(pow(xPosW,2) + pow(yPosW,2) + pow(zPosW,2) , 0.5);
                theta = acos(zPosW/normVector)/deg;
                
                if(xPosW==0)
                {
                    if(yPosW==0) phi = 0;
                    if(yPosW>0) phi = 90;
                    if(yPosW<0) phi = 270;
                }
                else
                {
                    phi = atan(yPosW/xPosW)/deg;
                    
                    if(xPosW>0 && yPosW>0) phi = phi; // deg
                    if(xPosW<0 && yPosW>0) phi = phi + 180.; // deg
                    if(xPosW<0 && yPosW<0) phi = phi + 180.; // deg
                    if(xPosW>0 && yPosW<0) phi = phi + 360.; // deg
                }
                
                fEventAction->SetVar_TIARA_AA(TIARANo, TIARA_RowNo, TIARA_SectorNo, 1, iTS, theta);
                fEventAction->SetVar_TIARA_AA(TIARANo, TIARA_RowNo, TIARA_SectorNo, 2, iTS, phi);
                
            }
            
            fEventAction->FillVar_TIARA_AA(TIARANo, TIARA_RowNo, TIARA_SectorNo, 0, iTS, edepTIARA_AA);
        }
    }
    
    ////////////////////////////////////////////////
    //              VDC DETECTORS
    ////////////////////////////////////////////////
    
    if(interactiontime < VDC_TotalSampledTime)
    {
        if(volumeName == "VDC_SenseRegion_USDS")
        {
            WireChamberNo = volume->GetCopyNo();
            
            iTS = interactiontime/PADDLE_SamplingTime;
            edepVDC = aStep->GetTotalEnergyDeposit()/keV;
            
            worldPosition = preStepPoint->GetPosition();
            localPosition = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPosition);
            
            G4int cellNo = 0;
            G4int bufferNo = 0;
            G4bool CompletedVDCFilling = false;
            
            
            //  X WireChamber
            if( (WireChamberNo==0) || (WireChamberNo==2) )
            {
                xPosL = localPosition.x()/mm;
                yPosL = localPosition.y()/mm;
                zPosL = localPosition.z()/mm + 4.0;
                
                if(abs(zPosL)>8) CompletedVDCFilling = true;
                
                while(cellNo<198 && !CompletedVDCFilling)
                {
                    if( (xPosL > (-99+cellNo)*4) && (xPosL <= (-98+cellNo)*4) )
                    {
                        if(WireChamberNo==0) channelID = cellNo;
                        if(WireChamberNo==2) channelID = cellNo + 341;
                        //G4cout << "Here is the WireChamberNo     -->     "<< WireChamberNo << G4endl;
                        //G4cout << "Here is the X WireChamber Triggered Cell     -->     "<< i << G4endl;
                        //fEventAction->FillVDC_Observables(channelID, edepVDC, edepVDC*zPosL, edepVDC*interactiontime);
                        
                        while(bufferNo<hit_buffersize && !CompletedVDCFilling)
                        {
                            hit_StoredChannelNo = fEventAction->GetVDC_ObservablesChannelID(bufferNo);
                            
                            if( (hit_StoredChannelNo < 0) || (hit_StoredChannelNo == channelID) )
                            {
                                fEventAction->FillVDC_Observables(bufferNo, channelID, edepVDC, edepVDC*zPosL, edepVDC*interactiontime);
                                
                                CompletedVDCFilling = true;
                            }
                            
                            bufferNo++;
                        }
                    }
                    
                    cellNo++;
                }
            }
            
            //  U WireChamber
            if( (WireChamberNo==1) || (WireChamberNo==3) )
            {
                xPosL = localPosition.x()/mm;
                yPosL = localPosition.y()/mm;
                zPosL = localPosition.z()/mm - 4.0;
                
                if(abs(zPosL)>8) CompletedVDCFilling = true;
                
                xOffset = -(1/tan(50))*yPosL;
                
                while(cellNo<143 && !CompletedVDCFilling)
                {
                    if( (xPosL > (-71.5+cellNo)*abs(xShift) + xOffset) && (xPosL <= (-70.5+cellNo)*abs(xShift) + xOffset) )
                    {
                        if(WireChamberNo==1) channelID = cellNo + 198;
                        if(WireChamberNo==3) channelID = cellNo + 539;
                        
                        while(bufferNo<hit_buffersize && !CompletedVDCFilling)
                        {
                            hit_StoredChannelNo = fEventAction->GetVDC_ObservablesChannelID(bufferNo);
                            
                            if( (hit_StoredChannelNo < 0) || (hit_StoredChannelNo == channelID) )
                            {
                                fEventAction->FillVDC_Observables(bufferNo, channelID, edepVDC, edepVDC*zPosL, edepVDC*interactiontime);
                                
                                CompletedVDCFilling = true;
                            }
                            
                            bufferNo++;
                        }
                    }
                    
                    cellNo++;
                }
                
            }
            
            ////    The PRE-point
            if(zPosL<0. && aStep->GetTrack()->GetParentID()==0)
            {
                fEventAction->SetVDC_WireplaneTraversePos(WireChamberNo, 0, 0, xPosL);
                fEventAction->SetVDC_WireplaneTraversePos(WireChamberNo, 0, 1, yPosL);
                fEventAction->SetVDC_WireplaneTraversePos(WireChamberNo, 0, 2, zPosL);
            }
            
            ////    The POST-point
            if(zPosL>0. && aStep->GetTrack()->GetParentID()==0 && fEventAction->GetVDC_WireplaneTraversePOST(WireChamberNo)==false)
            {
                fEventAction->SetVDC_WireplaneTraversePOST(WireChamberNo, true);
                fEventAction->SetVDC_WireplaneTraversePos(WireChamberNo, 1, 0, xPosL);
                fEventAction->SetVDC_WireplaneTraversePos(WireChamberNo, 1, 1, yPosL);
                fEventAction->SetVDC_WireplaneTraversePos(WireChamberNo, 1, 2, zPosL);
 
            }
            
        }
    }
    
    
    
    ////////////////////////////////////////////////
    //              PADDLE DETECTORS
    ////////////////////////////////////////////////
    
    if (interactiontime < PADDLE_TotalSampledTime)
    {
        if (volumeName == "PADDLE")
        {
            channelID = volume->GetCopyNo();
            
            PADDLENo = channelID;
            
            iTS = interactiontime/PADDLE_SamplingTime;
            edepPADDLE = aStep->GetTotalEnergyDeposit()/MeV;
            
            worldPosition = preStepPoint->GetPosition();
            localPosition = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPosition);
            
            fEventAction->AddEnergy_PADDLE( PADDLENo, iTS, edepPADDLE);
            fEventAction->TagTOF_PADDLE(PADDLENo, iTS, interactiontime);
            fEventAction->AddEWpositionX_PADDLE( PADDLENo, iTS, edepPADDLE*localPosition.x());
            fEventAction->AddEWpositionY_PADDLE( PADDLENo, iTS, edepPADDLE*localPosition.y());
            
            //if(fEventAction->Get_PADDLE_Trig(i) == false) fEventAction->Set_PADDLE_Trig(i, true);
        }
    }
    
    
    
    ////////////////////////////////////////////////
    //                  CLOVERS
    ////////////////////////////////////////////////
    
    if(interactiontime < CLOVER_TotalSampledTime)
    {
       // if(volumeName == "CLOVER_HPGeCrystal" && particleName=="neutron")
       // if(volumeName == "CLOVER_HPGeCrystal" && particleName=="gamma")
        if(volumeName == "CLOVER_HPGeCrystal" && particleName == "gamma")
        {
            //G4cout << "particleName:    " << particleName <<  G4endl;

            channelID = volume->GetCopyNo();
            
            CLOVERNo = channelID/4;
            CLOVER_HPGeCrystalNo = channelID%4;
            
            /*
             G4cout << "Here is the copyNo    "<< copyNo << G4endl;
             G4cout << "Here is the CLOVERNo    "<< CLOVERNo << G4endl;
             G4cout << "Here is the CLOVER_HPGeCrystalNo    "<< CLOVER_HPGeCrystalNo << G4endl;
             G4cout << " "<< G4endl;
             */
            
            iTS = interactiontime/CLOVER_SamplingTime;
            edepCLOVER_HPGeCrystal = aStep->GetTotalEnergyDeposit()/keV;
            
            fEventAction->AddEnergyCLOVER_HPGeCrystal(CLOVERNo, CLOVER_HPGeCrystalNo, iTS, edepCLOVER_HPGeCrystal);
            
            
            

         //   G4cout << "line 330  ******************* "<< fEventAction->GetCLOVER_iEDep(CLOVERNo) <<G4endl;
            
            if(fEventAction->GetCLOVER_iEDep(CLOVERNo)==0)
            {
                
             //   G4cout << "line 335  &&&&&&&&&&&&&&&&&&& "<< fEventAction->GetCLOVER_iEDep(CLOVERNo) <<G4endl;
                
                G4double initialE = aStep->GetPreStepPoint()->GetKineticEnergy()/keV;
                fEventAction->SetCLOVER_iEDep(CLOVERNo, initialE);
            
              //  G4cout << " CLOVER incident E     "<< initialE <<  "		particleName     "<<particleName << G4endl;
		
            //    if(initialE>0.0) G4cout << "HELOOOOOOO:    " << initialE << G4endl;
            //    G4cout << "line 340 -------------- "<< initialE <<G4endl;
            }

        }
    }
    
    /*
     if (interactiontime < CLOVER_Shield_BGO_TotalSampledTime)
     {
     for(G4int i=0; i<8; i++)
     {
     for(G4int l=0; l<16; l++)
     {
     if ( volume == fDetConstruction->GetVolume_CLOVER_Shield_BGOCrystal(i, l) && interactiontime <     CLOVER_Shield_BGO_TotalSampledTime )
     {
     iTS = interactiontime/CLOVER_Shield_BGO_SamplingTime;
     edepCLOVER_BGOCrystal = aStep->GetTotalEnergyDeposit()/keV;
     
     fEventAction->AddEnergyBGODetectors(i, l, iTS, edepCLOVER_BGOCrystal);
     //G4cout << "Here is the edepCLOVER_BGOCrystal    "<< edepBGO << G4endl;
     }
     }
     }
     }
     */

    
    ////////////////////////////////////////////////
    //              PARAFFIN BOX
    ////////////////////////////////////////////////
    
 G4double ParaffinBoxInitialE;
    
 const G4String& processName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

    // if(volumeName == "ParaffinBox" && particleName == "gamma")
    if(volumeName == "ParaffinBox" && particleName == "gamma")
    {
        //if(particleName=="neutron") G4cout << "particleName:    " << particleName <<  G4endl;
        //G4cout << "particleName:    " << particleName <<  G4endl;

        
        edepParaffinBox = aStep->GetTotalEnergyDeposit()/keV;
        //G4cout << "line 412 ------------ E depos in Paraffin   "<< edepParaffinBox << G4endl;
        //G4cout << "line 412 ------------ E depos in Paraffin   "<< edepParaffinBox << G4endl;
        
        fEventAction->AddEnergyParaffinBox(edepParaffinBox);
	/*
           if(ParaffinBoxFinalE>5000){
		 G4cout << " ParaffinBoxFinalE     "<< ParaffinBoxFinalE <<  "			   particleName     "<<particleName << G4endl;
		}
       */ 
       		/*if(particleName!="neutron"){ G4cout << "particleName     "<< particleName << "  Process name " << aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
            }*/
            //G4cout << "line 377  ******************* "<< ParaffinBoxInitialE <<G4endl;
            //G4cout << "line 378  §§§§§§§§§§§ Parrafin box energy initial BEFORE loop "<< fEventAction->GetPARAFFINBOX_iEDep() <<G4endl;
            if(fEventAction->GetPARAFFINBOX_iEDep() == 0.0){
                //G4cout << "line 379  +++++++++++++ Parrafin box energy initial AFTER loop "<< fEventAction->GetPARAFFINBOX_iEDep() << G4endl;
                 ParaffinBoxInitialE = aStep->GetPreStepPoint()->GetKineticEnergy()/keV;
		/*if(ParaffinBoxInitialE>5000){
		 G4cout << " PARAFFIN incident E     "<< ParaffinBoxInitialE <<  "		particleName     "<< particleName << "Process name " << aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
		}*/
                //G4cout << "&&&&&&&&&&&&&&&&&&&&&   Parrafin box energy initial SET  "<< ParaffinBoxInitialE << G4endl;
                fEventAction->SetPARAFFINBOX_iEDep(ParaffinBoxInitialE);
            }
        
        
     /*
        G4double flagParticle=-1.;
        G4double flagProcess=-1.;
        G4double x,y,z,xp,yp,zp;
        const G4String& particleName = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
        const G4String& processName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
        // Particle
        if (particleName == "neutron")       flagParticle = 1;
        else if (particleName == "e-")       flagParticle = 2;
        else if (particleName == "proton")   flagParticle = 3;
        else if (particleName == "deuteron") flagParticle = 4;
        else if (particleName == "e+")       flagParticle = 5;
        else if (particleName == "helium")   flagParticle = 6;
        else if (particleName == "12C")      flagParticle = 7;
     
        // Processes
        if (processName=="hadElastic")    flagProcess =11;
        else if (processName=="eIoni")    flagProcess =12;
        else if (processName=="hIoni")    flagProcess =13;
        else if (processName=="ionIoni")    flagProcess =14;
        else if (processName=="RadioactiveDecay")  flagProcess =15;
        

            x=aStep->GetPreStepPoint()->GetPosition().x()/nanometer;
            y=aStep->GetPreStepPoint()->GetPosition().y()/nanometer;
            z=aStep->GetPreStepPoint()->GetPosition().z()/nanometer;
            xp=aStep->GetPostStepPoint()->GetPosition().x()/nanometer;
            yp=aStep->GetPostStepPoint()->GetPosition().y()/nanometer;
            zp=aStep->GetPostStepPoint()->GetPosition().z()/nanometer;
      
        
            // get analysis manager
         //   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
            
            // fill ntuple
            G4AnalysisManager::Instance()->FillH1(flagParticle);
            //analysisManager->FillNtupleDColumn(1, flagProcess);
            //analysisManager->FillNtupleDColumn(2, x);
          //  analysisManager->FillNtupleDColumn(3, y);
         //   analysisManager->FillNtupleDColumn(4, z);
          //  analysisManager->FillNtupleDColumn(2, aStep->GetTotalEnergyDeposit()/eV);
         //   analysisManager->FillNtupleDColumn(6,std::sqrt((x-xp)*(x-xp)+y-yp)*(y-yp)+(z-zp)*(z-zp))/nm);
         //   analysisManager->FillNtupleDColumn(3,(aStep->GetPreStepPoint()->GetKineticEnergy()-aStep->GetPostStepPoint()->GetKineticEnergy())/eV );
           // analysisManager->FillNtupleIColumn(4, G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID());
            //analysisManager->AddNtupleRow();

       */
        
    }
    
    
    ////////////////////////////////////////////////
    //              IRON BOX
    ////////////////////////////////////////////////
    
    G4double IronBoxInitialE;
    
    // if(volumeName == "IronBox" && particleName == "neutron")
   // if(volumeName == "IronBox" && particleName == "gamma")
    if(volumeName == "IronBox" && particleName == "gamma")
    {
        
        
        edepIronBox = aStep->GetTotalEnergyDeposit()/keV;
        //G4cout << "line 412 ------------ E depos in Iron   "<< edepIronBox << G4endl;
        
        
        fEventAction->AddEnergyIronBox(edepIronBox);

        
        //G4cout << "line 418  §§§§§§§§§§ Iron box energy initial BEFORE loop "<< fEventAction->GetIRONBOX_iEDep() <<G4endl;
        
        if(fEventAction->GetIRONBOX_iEDep() == 0.0){
            //G4cout << "line 421  +++++++++++++ Iron box energy initial AFTER loop "<< fEventAction->GetIRONBOX_iEDep() << G4endl;
            IronBoxInitialE = aStep->GetPreStepPoint()->GetKineticEnergy()/keV;
	/*	if(IronBoxInitialE>5000){
		 G4cout << "IRON incident E       "<< IronBoxInitialE <<  "		particleName     "<<particleName << G4endl;
		}*/
            fEventAction->SetIRONBOX_iEDep(IronBoxInitialE);
        }
        
        
    }
    
    
    edepParaffinBox = aStep->GetTotalEnergyDeposit()/keV;
    if(edepParaffinBox>0.0)
    {
    //    G4cout << "Energy deposited:   "<< edepParaffinBox << G4endl;
        fEventAction->SetTotalEnergyDeposition(edepParaffinBox);
    }
    
    
    ////////////////////////////////////////////////
    //              LEPS DETECTOR ARRAY
    ////////////////////////////////////////////////
    
    
    if((interactiontime < LEPS_TotalSampledTime) && (volumeName == "LEPSHPGeCrystal"))
    {
        channelID = volume->GetCopyNo();
        
        LEPSNo = channelID/4;
        LEPS_HPGeCrystalNo = channelID%4;
        
        iTS = interactiontime/LEPS_SamplingTime;
        edepLEPS_HPGeCrystal = aStep->GetTotalEnergyDeposit()/keV;
        
        fEventAction->AddEnergyLEPS_HPGeCrystals(LEPSNo, LEPS_HPGeCrystalNo, iTS, edepLEPS_HPGeCrystal);
        
    }
    
    
    ////////////////////////////////////////////////
    //              NAIS DETECTOR ARRAY
    ////////////////////////////////////////////////
    
    
    if((interactiontime < NAIS_TotalSampledTime) && (volumeName == "NAISNaICrystal"))
    {
        channelID = volume->GetCopyNo();
        
        NAISNo = channelID;
        
        iTS = interactiontime/NAIS_SamplingTime;
        edepNAIS_NaICrystal = aStep->GetTotalEnergyDeposit()/keV;
        
        fEventAction->AddEnergyNAIS_NaICrystals(NAISNo, iTS, edepNAIS_NaICrystal);
        
    }
    
    
    //if (volumeName=="TIARA_Assembly" && !fEventAction->GA_GetLineOfSight() )   G4cout << "Here is the TIARA_Assembly Hit!" << G4endl;
    
    ////////////////////////////////////////////
    //              TIARA ARRAY
    ////////////////////////////////////////////
        
    if(GA_MODE)
    {
        
        if(((volumeName=="TIARA_AA_RS" || volumeName=="TIARA_SiliconWafer") && ((GA_LineOfSightMODE && fEventAction->GA_GetLineOfSight()==true) || !GA_LineOfSightMODE)) || (volumeName == "World" && GA_GenInputVar))
        //if((((volumeName=="TIARA_AA_RS" || volumeName=="TIARA_SiliconWafer") && ((GA_LineOfSightMODE && fEventAction->GA_GetLineOfSight()==true) || !GA_LineOfSightMODE)) || (volumeName == "World" && GA_GenInputVar)) && particleName == "gamma")
        {
            channelID = volume->GetCopyNo();
            worldPosition = preStepPoint->GetPosition();
            //worldPosition = worldPosition.unit();
            
            xPosW = worldPosition.x()/m;
            yPosW = worldPosition.y()/m;
            zPosW = worldPosition.z()/m;
            
            /*
             if(volumeName == "TIARA_SiliconWafer")
             {
             G4cout << " " << G4endl;
             G4cout << "Here is the TIARA_SiliconWafer HIT" << G4endl;
             }
             
             if(volumeName == "TIARA_AA_RS")
             {
             G4cout << " " << G4endl;
             G4cout << "Here is the TIARA_AA_RS HIT" << G4endl;
             }
             */
            
            if(volumeName == "TIARA_AA_RS")
            {
                fEventAction->FillGA_TIARAstor(channelID, 0, xPosW);
                fEventAction->FillGA_TIARAstor(channelID, 1, yPosW);
                fEventAction->FillGA_TIARAstor(channelID, 2, zPosW);
                fEventAction->FillGA_TIARAstor(channelID, 3, 1.);
            }
            
            if(GA_GenAngDist && fEventAction->GetGA_TIARA(channelID, 0)==0)
            {
                normVector = pow(pow(xPosW,2) + pow(yPosW,2) + pow(zPosW,2) , 0.5);
                theta = acos(zPosW/normVector)/deg;
                
                if(xPosW==0)
                {
                    if(yPosW==0) phi = 0;
                    if(yPosW>0) phi = 90;
                    if(yPosW<0) phi = 270;
                }
                else
                {
                    phi = atan(yPosW/xPosW)/deg;
                    
                    if(xPosW>0 && yPosW>0) phi = phi; // deg
                    if(xPosW<0 && yPosW>0) phi = phi + 180.; // deg
                    if(xPosW<0 && yPosW<0) phi = phi + 180.; // deg
                    if(xPosW>0 && yPosW<0) phi = phi + 360.; // deg
                }
                
                if(volumeName == "TIARA_AA_RS")
                {
                    fEventAction->SetGA_TIARA(channelID, 0, 1);
                    fEventAction->SetGA_TIARA(channelID, 1, theta);
                    fEventAction->SetGA_TIARA(channelID, 2, phi);
                }
            }
            
            if(GA_GenInputVar && volumeName == "World")
            {
                fEventAction->SetInputDist(0, theta);
                fEventAction->SetInputDist(1, phi);
            }
            //G4cout << "Here is the geantino Hit!     -->     " << G4endl;
        }
    }
    
    
    ////    Here, one declares the volumes that one considers will block the particles of interest and effectively mask the relevant volume of interest.
    if (GA_LineOfSightMODE && (volumeName == "TIARA_AA_RS" || volumeName=="TIARA_PCB" || volumeName=="TIARA_SiliconWafer"))
    {
        //G4cout << "Here is the volumeName    "<< volumeName << G4endl;
        fEventAction->GA_SetLineOfSight(false);
    }
    
    
    
    /*
     // Collect energy and track length step by step
     
     // get volume of the current step
     G4VPhysicalVolume* volume
     = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
     
     // energy deposit
     G4double edep = step->GetTotalEnergyDeposit();
     
     
     // step length
     G4double stepLength = 0.;
     if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
     stepLength = step->GetStepLength();
     }
     
     if ( volume == fDetConstruction->GetAbsorberPV() ) {
     fEventAction->AddAbs(edep,stepLength);
     }
     
     if ( volume == fDetConstruction->GetGapPV() ) {
     fEventAction->AddGap(edep,stepLength);
     }
     */
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
