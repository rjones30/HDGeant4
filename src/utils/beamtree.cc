//
// beamtree.cc : main program for GlueXBremsstrahlung beam generation
//
// author: richard.t.jones at uconn.edu
// version: february 4, 2017
//

#include <GlueXBremsstrahlungGenerator.hh>
#include <unistd.h>
#include <iostream>

void usage()
{
   std::cout << std::endl
          << "Usage: beamtree [options]" << std::endl
          << " where options are:" << std::endl
          << "   -n # : generate # events, default is 1000" << std::endl;
   exit(9);
}

int main(int argc,char** argv)
{
   int nevents = 1000;

   char c;
   while ((c = getopt(argc, argv, "n:")) != -1) {
      if (c == 'n') {
         nevents = atoi(optarg);
      }
      else {
         usage();
      }
   }

   GlueXBremsstrahlungGenerator gen;
   gen.GenerateBeamPhotons(nevents);
   return 0;
}
