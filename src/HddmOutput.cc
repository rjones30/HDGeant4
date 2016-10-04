//
// class implementation for HddmOutput
//
// author: richard.t.jones at uconn.edu
// version: september 28, 2016

#include "HddmOutput.hh"

int HddmOutput::fRunNo = 0;
int HddmOutput::fEventNo = 0;
std::ofstream *HddmOutput::fHDDMoutfile = 0;
hddm_s::ostream *HddmOutput::fHDDMostream = 0;

HddmOutput::HddmOutput(const std::string &filename)
{
   if (fHDDMostream != 0) {
      delete fHDDMostream;
      delete fHDDMoutfile;
   }
   fHDDMoutfile = new std::ofstream(filename);
   fHDDMostream = new hddm_s::ostream(*fHDDMoutfile);
}

HddmOutput::~HddmOutput()
{
   if (fHDDMostream != 0) {
      delete fHDDMostream;
      delete fHDDMoutfile;
      fHDDMostream = 0;
      fHDDMoutfile = 0;
   }
}

HddmOutput::HddmOutput(HddmOutput &src)
{}

HddmOutput& HddmOutput::operator=(HddmOutput &src)
{
   return *this;
}

void HddmOutput::WriteOutputHDDM(hddm_s::HDDM &record)
{
   if (fHDDMostream != 0) {
      *fHDDMostream << record;
   }
}
