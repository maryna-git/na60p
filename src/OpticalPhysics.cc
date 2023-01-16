#include "OpticalPhysics.hh"

// #include "G4OpticalParameters.hh" // for newer installs / versions of G4 at least G4.10.07.p01

#include "G4PhysicsListHelper.hh"
#include "G4ProcessManager.hh"
#include "CLHEP/Units/SystemOfUnits.h"

using namespace CLHEP;

OpticalPhysics::OpticalPhysics(const G4String& name):
G4VPhysicsConstructor(name)
{
	G4cout << "Optical physics" << G4endl;
}

OpticalPhysics::~OpticalPhysics()
{
	
}

void OpticalPhysics::ConstructParticle()
{
	G4OpticalPhoton::OpticalPhotonDefinition();
}

void OpticalPhysics::ConstructProcess()
{
  // 	 G4OpticalParameters* params = G4OpticalParameters::Instance();  // use this construct for newer installs of G4, at least G4.10.07.p01
 

  G4Cerenkov* fCerenkovProcess = new G4Cerenkov("Cerenkov");
	fCerenkovProcess->SetMaxNumPhotonsPerStep(1000000);
	fCerenkovProcess->SetMaxBetaChangePerStep(10.0);
	fCerenkovProcess->SetTrackSecondariesFirst(true);
	G4Scintillation* fScintillationProcess = new G4Scintillation("Scintillation");
	fScintillationProcess->SetScintillationByParticleType(false);
	fScintillationProcess->SetScintillationYieldFactor(0.1);
	fScintillationProcess->SetTrackSecondariesFirst(true);
	G4cout << "Optical physics!: scint Yield factor: " << fScintillationProcess->GetScintillationYieldFactor() << G4endl;

	fCerenkovProcess->SetVerboseLevel(0);
	fScintillationProcess->SetVerboseLevel(0);


	// Use Birks Correction in the Scintillation process
	if(G4Threading::IsMasterThread())
	 {
	 G4EmSaturation* emSaturation =
		G4LossTableManager::Instance()->EmSaturation();
		fScintillationProcess->AddSaturation(emSaturation);
	 }

	auto particleIterator=GetParticleIterator();
	particleIterator->reset();
	while( (*particleIterator)() )
	{
		G4ParticleDefinition* particle = particleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();
		if (fCerenkovProcess->IsApplicable(*particle)) 
		{
			pmanager->AddProcess(fCerenkovProcess);
			pmanager->SetProcessOrdering(fCerenkovProcess,idxPostStep);
		}
		if (fScintillationProcess->IsApplicable(*particle)) 
		{
			pmanager->AddProcess(fScintillationProcess);
			pmanager->SetProcessOrderingToLast(fScintillationProcess, idxAtRest);
			pmanager->SetProcessOrderingToLast(fScintillationProcess, idxPostStep);
		}
		if (particleName == "opticalphoton") 
		{

		}
	}
}
