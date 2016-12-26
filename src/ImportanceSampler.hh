//
// ImportanceSampler class header
//
// author: richard.t.jones at uconn.edu
// version: december 24, 2016
//
// Utility class for importance sampling of random variables.

#ifndef ImportanceSampler_H
#define ImportanceSampler_H

#include <vector>

class ImportanceSampler {
 public:
   ImportanceSampler()
    : Psum(0), Pcut(1), Pmax(0), Ntested(0), Npassed(0) {}

   std::vector<double> randvar;
   std::vector<double> density;
   std::vector<double> integral;
   double Psum;
   double Pcut;
   double Pmax;
   long int Ntested;
   long int Npassed;

   static unsigned int search(double u, const std::vector<double> &list)
   {
      // Perform a fast search through non-decreasing list
      // for the index of the first element not less than u.

      int imin = -1;
      int imax = list.size() - 1;
      while (imax > imin + 1) {
         int i = (imin + imax) / 2;
         if (list[i] >= u)
            imax = i;
         else
            imin = i;
      }
      return imax;
   }
   unsigned int search(double u) const
   {
      return search(u, integral);
   }

};


#endif
