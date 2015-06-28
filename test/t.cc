//
// test_GlueXUserOptions.cc
//
// purpose: Test program for GlueXUserOptions class
// author: richard.t.jones at uconn.edu
// version: june 27, 2015
//

#include "../src/GlueXUserOptions.hh"
#include <iostream>

int main()
{
   GlueXUserOptions options;
   options.ReadControl_in();

   std::map<int, std::string> svals;
   int retvals = options.Find("beam", svals);
   std::cout << "string retval was " << retvals << std::endl;
   std::cout << "svals was";
   for (int i=1; i < 99; ++i) {
      if (svals.find(i) != svals.end())
         std::cout << " " << svals[i];
   }
   std::cout << std::endl;

   std::map<int, double> dvals;
   int retvald = options.Find("beam", dvals);
   std::cout << "double retval was " << retvald << std::endl;
   std::cout << "dvals was";
   for (int i=1; i < 99; ++i) {
      if (dvals.find(i) != dvals.end())
         std::cout << " " << dvals[i];
   }
   std::cout << std::endl;

   std::map<int, int> ivals;
   int retvali = options.Find("beam", ivals);
   std::cout << "int retval was " << retvali << std::endl;
   std::cout << "ivals was";
   for (int i=1; i < 99; ++i) {
      if (ivals.find(i) != ivals.end())
         std::cout << " " << ivals[i];
   }
   std::cout << std::endl;

   int retvalx = options.Find("beam2", ivals);
   std::cout << "wrong retval was " << retvalx << std::endl;
}
