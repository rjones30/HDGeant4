//
// GlueXRunAction class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXRunAction.hh"

#include "G4Run.hh"

GlueXRunAction::GlueXRunAction()
{}

GlueXRunAction::~GlueXRunAction()
{}

void GlueXRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
}

void GlueXRunAction::EndOfRunAction(const G4Run*)
{ }
