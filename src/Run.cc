/// \brief Implementation of the Run class
//

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4EmCalculator.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

#include "N6SetUp.hh"

#include <iomanip>


Run::Run(DetectorConstruction* det)
: G4Run(),
  fDetector(det), 
  fParticle(0), fEkin(0.), 
  fHistoManager(0), fTreeCutX(1.0*mm), fTreeCutY(1.0*mm),
  fTreeParticle(0), fSkipEvents(0)
{
  fEnergyDeposit  = fEnergyDeposit2  = 0.;
  fTrakLenCharged = fTrakLenCharged2 = 0.;
  fTrakLenNeutral = fTrakLenNeutral2 = 0.;
  fNbStepsCharged = fNbStepsCharged2 = 0.;
  fNbStepsNeutral = fNbStepsNeutral2 = 0.;
  fMscProjecTheta = fMscProjecTheta2 = 0.;

  fNbGamma = fNbElect = fNbPosit = 0;

  fTransmit[0] = fTransmit[1] = fReflect[0] = fReflect[1] = 0;
  
  fMscEntryCentral = 0;
  
  fEnergyLeak[0] = fEnergyLeak[1] = fEnergyLeak2[0] = fEnergyLeak2[1] = 0.;
}



Run::~Run()
{ }



void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin = energy;

}
 


void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  // pass information about primary particle
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;

  // accumulate sums
  //
  fEnergyDeposit   += localRun->fEnergyDeposit;  
  fEnergyDeposit2  += localRun->fEnergyDeposit2;  
  fTrakLenCharged  += localRun->fTrakLenCharged;
  fTrakLenCharged2 += localRun->fTrakLenCharged2;   
  fTrakLenNeutral  += localRun->fTrakLenNeutral;  
  fTrakLenNeutral2 += localRun->fTrakLenNeutral2;
  fNbStepsCharged  += localRun->fNbStepsCharged;
  fNbStepsCharged2 += localRun->fNbStepsCharged2;
  fNbStepsNeutral  += localRun->fNbStepsNeutral;
  fNbStepsNeutral2 += localRun->fNbStepsNeutral2;
  fMscProjecTheta  += localRun->fMscProjecTheta;
  fMscProjecTheta2 += localRun->fMscProjecTheta2;

    
  fNbGamma += localRun->fNbGamma;
  fNbElect += localRun->fNbElect;      
  fNbPosit += localRun->fNbPosit;
  
  fTransmit[0] += localRun->fTransmit[0];  
  fTransmit[1] += localRun->fTransmit[1];
  fReflect[0]  += localRun->fReflect[0];
  fReflect[1]  += localRun->fReflect[1];
  
  fMscEntryCentral += localRun->fMscEntryCentral;
  
  fEnergyLeak[0]  += localRun->fEnergyLeak[0];
  fEnergyLeak[1]  += localRun->fEnergyLeak[1];
  fEnergyLeak2[0] += localRun->fEnergyLeak2[0];
  fEnergyLeak2[1] += localRun->fEnergyLeak2[1];
                  
  G4Run::Merge(run); 
} 



void Run::EndOfRun()
{
  // compute mean and rms
  //
  G4int TotNbofEvents = numberOfEvent;
  if (TotNbofEvents == 0) return;
  
  G4double EnergyBalance = fEnergyDeposit + fEnergyLeak[0] + fEnergyLeak[1];
  EnergyBalance /= TotNbofEvents;

  fEnergyDeposit /= TotNbofEvents; fEnergyDeposit2 /= TotNbofEvents;
  G4double rmsEdep = fEnergyDeposit2 - fEnergyDeposit*fEnergyDeposit;
  if (rmsEdep>0.) rmsEdep = std::sqrt(rmsEdep/TotNbofEvents);
  else            rmsEdep = 0.;

  fTrakLenCharged /= TotNbofEvents; fTrakLenCharged2 /= TotNbofEvents;
  G4double rmsTLCh = fTrakLenCharged2 - fTrakLenCharged*fTrakLenCharged;
  if (rmsTLCh>0.) rmsTLCh = std::sqrt(rmsTLCh/TotNbofEvents);
  else            rmsTLCh = 0.;

  fTrakLenNeutral /= TotNbofEvents; fTrakLenNeutral2 /= TotNbofEvents;
  G4double rmsTLNe = fTrakLenNeutral2 - fTrakLenNeutral*fTrakLenNeutral;
  if (rmsTLNe>0.) rmsTLNe = std::sqrt(rmsTLNe/TotNbofEvents);
  else            rmsTLNe = 0.;

  fNbStepsCharged /= TotNbofEvents; fNbStepsCharged2 /= TotNbofEvents;
  G4double rmsStCh = fNbStepsCharged2 - fNbStepsCharged*fNbStepsCharged;
  if (rmsStCh>0.) rmsStCh = std::sqrt(rmsTLCh/TotNbofEvents);
  else            rmsStCh = 0.;

  fNbStepsNeutral /= TotNbofEvents; fNbStepsNeutral2 /= TotNbofEvents;
  G4double rmsStNe = fNbStepsNeutral2 - fNbStepsNeutral*fNbStepsNeutral;
  if (rmsStNe>0.) rmsStNe = std::sqrt(rmsTLCh/TotNbofEvents);
  else            rmsStNe = 0.;

  G4double transmit[2];
  transmit[0] = 100.*fTransmit[0]/TotNbofEvents;
  transmit[1] = 100.*fTransmit[1]/TotNbofEvents;

  G4double reflect[2];
  reflect[0] = 100.*fReflect[0]/TotNbofEvents;
  reflect[1] = 100.*fReflect[1]/TotNbofEvents;

  G4double rmsMsc = 0., tailMsc = 0.;
  if (fMscEntryCentral > 0) {
    fMscProjecTheta /= fMscEntryCentral; fMscProjecTheta2 /= fMscEntryCentral;
    rmsMsc = fMscProjecTheta2 - fMscProjecTheta*fMscProjecTheta;
    if (rmsMsc > 0.) { rmsMsc = std::sqrt(rmsMsc); }
    if(fTransmit[1] > 0.0) {
      tailMsc = 100.- (100.*fMscEntryCentral)/(2*fTransmit[1]);
    }    
  }
  
  fEnergyLeak[0] /= TotNbofEvents; fEnergyLeak2[0] /= TotNbofEvents;
  G4double rmsEl0 = fEnergyLeak2[0] - fEnergyLeak[0]*fEnergyLeak[0];
  if (rmsEl0>0.) rmsEl0 = std::sqrt(rmsEl0/TotNbofEvents);
  else           rmsEl0 = 0.;
  
  fEnergyLeak[1] /= TotNbofEvents; fEnergyLeak2[1] /= TotNbofEvents;
  G4double rmsEl1 = fEnergyLeak2[1] - fEnergyLeak[1]*fEnergyLeak[1];
  if (rmsEl1>0.) rmsEl1 = std::sqrt(rmsEl1/TotNbofEvents);
  else           rmsEl1 = 0.;    
  
   
  // normalize histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  G4int ih = 1;
  G4double binWidth = analysisManager->GetH1Width(ih);
  G4double unit     = analysisManager->GetH1Unit(ih);  
  G4double fac = unit/(TotNbofEvents*binWidth);
  analysisManager->ScaleH1(ih,fac);

  ih = 10;
  binWidth = analysisManager->GetH1Width(ih);
  unit     = analysisManager->GetH1Unit(ih);  
  fac = unit/(TotNbofEvents*binWidth);
  analysisManager->ScaleH1(ih,fac);

//< ob
  ih = 20;
  binWidth = analysisManager->GetH1Width(ih);
  unit     = analysisManager->GetH1Unit(ih);  
  fac = unit/(TotNbofEvents*binWidth);
  analysisManager->ScaleH1(ih,fac);
//> ob
  
  ih = 12;
  analysisManager->ScaleH1(ih,1./TotNbofEvents);
                    
//< ob
  ih = 22;
  analysisManager->ScaleH1(ih,1./TotNbofEvents);
//> ob

  // reset default precision
  //G4cout.precision(prec);
}   



void Run::SetHistoManager (HistoManager *hmng) 
{ 
  fHistoManager = hmng; 
  if (fHistoManager) {
    fTreeCutX = fHistoManager->GetTreeCutX();
    fTreeCutY = fHistoManager->GetTreeCutY();

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    fTreeParticle = particleTable->FindParticle(fHistoManager->GetTreeParticle());
    
    std::cout << "TTree X cut: " << G4BestUnit(fTreeCutX, "Length") << std::endl; 
    std::cout << "TTree Y cut: " << G4BestUnit(fTreeCutY, "Length") << std::endl;
    std::cout << "TTree particles: " << fHistoManager->GetTreeParticle() << std::endl;
    if (fTreeParticle) std::cout << "   it was found in G4ParticleTable with PDG: " << fTreeParticle->GetPDGEncoding() << std::endl;
    else std::cout << "   it was not found in G4ParticleTable, so all particles will be saved to the TTree\n";
  }
}



void Run::ScanGeometry(const G4VPhysicalVolume *topvol, int verbose)
{
  std::vector<const G4VPhysicalVolume*> pvol;
  std::vector<int> nd;
  std::vector<int> id;

  std::cout << "\n***************** Scanning geometry: ***************** \n";

  pvol.push_back(topvol);
  nd.push_back(topvol->GetLogicalVolume()->GetNoDaughters());
  id.push_back(0);
  int l = 0;

  while (l >= 0 ) {

    if (id[l] == nd[l]) {
      pvol.pop_back();
      nd.pop_back();
      id.pop_back();
      --l;
      continue;
    }

    if (verbose > 1) std::cout << std::string(2*l, ' ') << pvol[l]->GetName() << std::endl;

    const G4LogicalVolume *lv = pvol[l]->GetLogicalVolume();
    const G4VPhysicalVolume *pv = lv->GetDaughter(id[l]);
    if (verbose > 1) {
      G4LogicalVolume  *lvsens = pv->GetLogicalVolume();
      const G4Material  *lvmat = lvsens->GetMaterial();
      std::cout << std::string(2*(l+1), ' ') << pv->GetName() << "  instance:  " << pv->GetCopyNo()
                << "  mass: " << lvsens->GetMass()/kg
                << "  material: " << lvmat->GetName() << "  density: " << lvmat->GetDensity()/(kg/m3)
                << std::endl;
      }

    const auto &itr = fSensitiveVol.find(pv->GetName());
    if (itr != fSensitiveVol.end()) {
      if (verbose)  std::cout << std::string(2*(l+1), ' ') << "Detector: " << pv->GetName()
                              << "  instance:  " << pv->GetCopyNo() <<   std::endl;
      SaveDetectorTransformation(pv, pvol, l);
//       ++id[l];
//       continue;
    }

    const auto &vitr = fInterceptVol.find(pv->GetName());
    if (vitr != fInterceptVol.end()) {
      if (verbose)  std::cout << std::string(2*(l+1), ' ') << "Intercept volume: " << pv->GetName()
                              << "  instance:  " << pv->GetCopyNo() <<   std::endl;
      SaveVolumeTransformation(pv, pvol, l);
    }

    ++id[l];

    int nn = pv->GetLogicalVolume()->GetNoDaughters();
    if (nn > 0) {
      pvol.push_back(pv);
      nd.push_back(nn);
      id.push_back(0);
      ++l;
    }
  }
}


// Produces transformation from local coordinates of pv volume to global
void Run::SaveDetectorTransformation(const G4VPhysicalVolume *pv, std::vector<const G4VPhysicalVolume*> pvol, const int l)
{
  G4AffineTransform transform(pv->GetRotation(), pv->GetTranslation());
  for (int jj = l; jj > 0; --jj) {
    transform *= G4AffineTransform(pvol[jj]->GetRotation(), pvol[jj]->GetTranslation());
  }

  G4RotationMatrix mtrx = transform.NetRotation();
  G4ThreeVector trns = transform.NetTranslation();

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4int dettype(-1), subdepth(1), ldepth(0);
  G4int ncellx, ncelly;
  G4double size_x, size_y, size_z;
  G4String detname = pv->GetName();

  try {
    dettype = std::get<0>(fSensitiveVol.at(detname));
    subdepth = std::get<1>(fSensitiveVol.at(detname));
    ldepth = std::get<2>(fSensitiveVol.at(detname));
    const auto &sseg = fDetector->GetSensorSegmentation(detname);
    size_x = std::get<0>(sseg) /mm;
    size_y = std::get<1>(sseg) /mm;
    ncellx = std::get<2>(sseg);
    ncelly = std::get<3>(sseg);
    size_z = 0.0 /mm;
  } catch (const std::out_of_range& oor) {
    G4String msgstr("Something is wrong! Segmentation for the sensitive detector should be defined at this stage!");
    G4Exception("Run::", "SaveDetectorTransformation()", FatalException, msgstr.c_str());
  }

  G4int detid;
  if (subdepth <= 0) detid = pv->GetCopyNo();
  else detid = pvol[l-subdepth+1]->GetCopyNo();
  detid += dettype;
  G4int layerid;
  if (ldepth <= 0) layerid = pv->GetCopyNo();
  else layerid = pvol[l-ldepth+1]->GetCopyNo();

  analysisManager->FillNtupleIColumn(4, 0, detid);
  analysisManager->FillNtupleSColumn(4, 1, pv->GetName());
  analysisManager->FillNtupleIColumn(4, 2, layerid);
  analysisManager->FillNtupleIColumn(4, 3, ncellx);
  analysisManager->FillNtupleIColumn(4, 4, ncelly);
  analysisManager->FillNtupleDColumn(4, 5, size_x);
  analysisManager->FillNtupleDColumn(4, 6, size_y);
  analysisManager->FillNtupleDColumn(4, 7, size_z);
  analysisManager->FillNtupleDColumn(4, 8, trns.x() );
  analysisManager->FillNtupleDColumn(4, 9, trns.y() );
  analysisManager->FillNtupleDColumn(4, 10, trns.z() );
  analysisManager->FillNtupleDColumn(4, 11, mtrx.getPhi() );
  analysisManager->FillNtupleDColumn(4, 12, mtrx.getTheta() );
  analysisManager->FillNtupleDColumn(4, 13, mtrx.getPsi() );

  G4LogicalVolume  *lvsens = pv->GetLogicalVolume();
  const G4Material  *lvmat = lvsens->GetMaterial();
  analysisManager->FillNtupleDColumn(4, 14, lvsens->GetMass()/kg);
  analysisManager->FillNtupleSColumn(4, 15, lvmat->GetName());
  analysisManager->FillNtupleDColumn(4, 16, lvmat->GetDensity()/(kg/m3));

  analysisManager->AddNtupleRow(4);

}


////////////////////////////////////////////////////////////////////////////
// TRotation drt;
// drt.SetXEulerAngles(phi, theta, psi);
// 
// TVector3 trns(translation.x, translation.y, translation.z);
// xl = (cellx+0.5) * sizex / static_cast<double>(ncellx) - 0.5*sizex;
// yl = (celly+0.5) * sizey / static_cast<double>(ncelly) - 0.5*sizey;
// zl = 0.5*sizez;
// TVector3 hitl(xl, yl, zl);
// 
// hitg = drt * hitl + trns;
////////////////////////////////////////////////////////////////////////////


void Run::SaveVolumeTransformation(const G4VPhysicalVolume *pv, std::vector<const G4VPhysicalVolume*> pvol, const int l)
{
  G4AffineTransform transform(pv->GetRotation(), pv->GetTranslation());
  for (int jj = l; jj > 0; --jj) {
    transform *= G4AffineTransform(pvol[jj]->GetRotation(), pvol[jj]->GetTranslation());
  }

  G4RotationMatrix mtrx = transform.NetRotation();
  G4ThreeVector trns = transform.NetTranslation();

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4int ncellx, ncelly;
  G4double size_x, size_y, size_z;
  G4String detname = pv->GetName();

  const auto &volids = fInterceptVol[detname];

  G4int detid = volids[0] + pv->GetCopyNo();
  G4int bdec = 1000;
  for (auto iv = volids.cbegin()+1; iv != volids.cend(); ++iv) {
    detid += bdec * pvol[l-(*iv)+1]->GetCopyNo();
    bdec *= 10;
  }

  G4int layerid = -1;
  ncellx = ncelly = -1;
  size_x = size_y = size_z = -1.0;

  G4LogicalVolume  *lvsens = pv->GetLogicalVolume();
  const G4VSolid *solsens = lvsens->GetSolid();
  G4String soltype(solsens->GetEntityType());

  if (soltype == "G4Box") {
    size_x = 2.0 * (dynamic_cast<const G4Box*>(solsens))->GetXHalfLength();
    size_y = 2.0 * (dynamic_cast<const G4Box*>(solsens))->GetYHalfLength();
    size_z = 2.0 * (dynamic_cast<const G4Box*>(solsens))->GetZHalfLength();
  } else if (soltype == "G4Tubs") {
    size_x = dynamic_cast<const G4Tubs*>(solsens)->GetInnerRadius();
    size_y = dynamic_cast<const G4Tubs*>(solsens)->GetOuterRadius();
    size_z = 2.0 * dynamic_cast<const G4Tubs*>(solsens)->GetZHalfLength();
  }

  analysisManager->FillNtupleIColumn(4, 0, detid);
  analysisManager->FillNtupleSColumn(4, 1, pv->GetName());
  analysisManager->FillNtupleIColumn(4, 2, layerid);
  analysisManager->FillNtupleIColumn(4, 3, ncellx);
  analysisManager->FillNtupleIColumn(4, 4, ncelly);
  analysisManager->FillNtupleDColumn(4, 5, size_x);
  analysisManager->FillNtupleDColumn(4, 6, size_y);
  analysisManager->FillNtupleDColumn(4, 7, size_z);
  analysisManager->FillNtupleDColumn(4, 8, trns.x() );
  analysisManager->FillNtupleDColumn(4, 9, trns.y() );
  analysisManager->FillNtupleDColumn(4, 10, trns.z() );
  analysisManager->FillNtupleDColumn(4, 11, mtrx.getPhi() );
  analysisManager->FillNtupleDColumn(4, 12, mtrx.getTheta() );
  analysisManager->FillNtupleDColumn(4, 13, mtrx.getPsi() );

  const G4Material  *lvmat = lvsens->GetMaterial();
  analysisManager->FillNtupleDColumn(4, 14, lvsens->GetMass()/kg);
  analysisManager->FillNtupleSColumn(4, 15, lvmat->GetName());
  analysisManager->FillNtupleDColumn(4, 16, lvmat->GetDensity()/(kg/m3));

  analysisManager->AddNtupleRow(4);

}

