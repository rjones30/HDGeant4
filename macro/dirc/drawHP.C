#define glx__sim
#include "../../../../sim-recon/master/src/plugins/Analysis/pid_dirc/DrcEvent.h"
#include "glxtools.C"

void drawHP(TString infile="drc.root"){
  if(!glx_init(infile,1,"data/drawHP")) return;

  DrcHit hit;
  for (Int_t ievent=0; ievent<glx_ch->GetEntries(); ievent++){
    glx_nextEvent(ievent,10);
    
    for(Int_t h=0; h<glx_event->GetHitSize(); h++){
      hit = glx_event->GetHit(h);
      Int_t pmt = hit.GetPmtId();
      Int_t pix = hit.GetPixelId();
      TVector3 gpos = hit.GetPosition();
      Double_t time = hit.GetLeadTime();
      if(pmt<102) glx_hdigi[pmt]->Fill(pix%8, 7-pix/8);
    }
  }
  glx_drawDigi("m,p,v\n",0);
  glx_canvasAdd(glx_cdigi);
  // glx_canvasSave(1,0);
}
