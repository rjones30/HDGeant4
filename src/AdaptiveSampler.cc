//
// AdaptiveSampler - provide an adaptive algorithm for improved 
//                   importance sampling in Monte Carlo integration
//
// file: AdaptiveSampler.cc (see AdaptiveSampler.hh for header, usage)
// author: richard.t.jones at uconn.edu
// original release: february 23, 2016
//
// new feature release: may 30,2 2020
//    - see new notes 7,8 below
//    - richard.t.jones at uconn.edu
//
// implementation notes:
// 1) See the description at the top of AdaptiveSampler.hh for an overview
//    of the problem that this class is designed to help solve, and how to
//    use it.
//
// 2) If an AdaptiveSampler object is being shared by more than one thread,
//    only one of the threads can be calling the feedback method in the
//    sample loop at a time. Simultaneous entry from more than one thread
//    to any of the mutating methods (feedback, adapt, reset_stats, etc.)
//    produces undefined behavior. Either you can share a single instance
//    of AdaptiveSampler among many threads with the calls to mutating
//    methods protected by a mutex, or you can construct an AdapativeSample
//    instance for each thread and then pause at regular intervals to 
//    combine their statistics using save/merge/restore.
//
// 3) Users must provide their uniform random number generator function to
//    the AdaptiveSampler constructor. If the sample method is to be called
//    from multiple threads, this user-written function must be reentrant.
//    The user is responsible for initializing the state of the generator
//    before beginning the Monte Carlo loop.
//
// 4) Successful use of AdaptiveSampler requires that the user understand
//    some basic things about opportunistic importance sampling. The basic
//    strategy it uses is to look for hot spots in the integration domain
//    and zoom in to oversample in those regions. By implication, there
//    are regions in the integration domain that are undersampled, and 
//    this can lead to systematic underbias in the result. Because of this
//    bias, importance sampling is really a variance-reduction technique:
//    only by repeated application under a variety of conditions does one
//    reach any level of confidence that the result is correct. Tuning
//    parameters are provided to allow the user to vary the conditions
//    for adaptation and check the robustness of their result.
//      *) Adaptation_sampling_threshold - count the effective number 
//         of samples that need to have been generated within a given
//         cell before the statistical information within that cell is
//         sufficiently precise for the adaptation algorithm to decide
//         if and how it should be split. Before adaptaion is permitted
//         within any given cell, the product of the cell hit count and
//         its average hit efficiency (<wI>**2 / (N <wI**2>) must be
//         above this threshold.
//      *) Adaptation_efficiency_target - this is the goal that is set
//         for the adaptation algorithm to shoot for. Adaptations are
//         allowed if and only if the improvements to the efficiency
//         that they entail have a reasonable chance of achieving the
//         target value before Adaptation_maximum_cells (see below) is
//         reached.
//      *) Adaptation_maximum_depth - this is the maximum depth of the
//         adaptation tree. The finest subdivision of the domain of
//         integration will be a cell of volume 3**maximum_depth. This
//         prevents the algorithm from tuning on noise if the inherent
//         resolution in the integrand function is has a finite limit.
//         Because of the intrinsic limitations of double precision
//         floats, this parameter should normally not exceed 35.
//      *) Adaptation_maximum_cells - maximum number of cells that
//         should be allowed by the adaptation algorithm in building
//         out the sampling tree. The sampling tree is a hierarchical
//         subdivision of the integration domain that selects certain
//         regions in the ndim-dimensional hypercube for oversampling.
//         Each cell in the tree requires memory resources to hold the
//         statistics that it collects for the hits that are generated
//         within that cell. All cells in the tree have the same size,
//         approximately 40 + (8 * 3^ndim) for ndim dimensions. The
//         number of cells that are needed to achieve a certain target
//         efficiency depends on the problem. Generally one should
//         plan to accommodate as many cells as memory resources can
//         support. The algorithm attempts to achieve the efficiency
//         target using as few cells as it can.
//
// 5) Generally, one can expect to achieve a factor 3 improvement in
//    sampling efficiency each time the maximum depth of the tree
//    advances by ndim. Generically it is the maximum depth of the
//    tree and not the number of cells it contains that indicates how
//    much improvement the sampler has achieved over unbiased sampling.
//    The number of cells required to advance the depth another ndim
//    levels may increase or decrease with depth depending on the
//    fractal dimension of the integrand. It really depends on the
//    problem; there is no way in advance to predict how this will
//    behave. Obviously integrands of low fractal dimension will be
//    more susceptible to the generic subdivision strategy employed
//    by the AdaptiveSampling algorithm than those with high dimension.
//
// 6) Here are rules of thumb that will keep you from common blunders.
//     a) Do not be too aggressive in adapting, especially at first
//        when the AdaptiveSampler is new and undifferentiated. Keep
//        the sampling_threshold value as high as you can (several
//        hundred is usually safe), only lower it if you cannot get
//        the adaptation to progress. Keep in mind that the lower
//        this threshold is set, the more susceptible you are to
//        tuning on noise and introducing uncontrolled underbias.
//     b) Bias is generically on the low side of the true result.
//        This means that if you repeat a measurement 25 times and
//        one repetition gives an answer that is many sigma above
//        the rest, there is a low-visibility region in integration
//        space that this repetition has found. Based on this, you
//        should modify your adaptation strategy to increase the
//        visiblity of that region early on during development. See
//        point (c) for how that can be done.
//     c) If you know there are low-visiblity regions of the domain
//        of integration that tend to get excluded by the adaptation
//        algorithm, you can "train" the AdaptiveSampler by running it
//        initially for several thousand events on a integrand that is
//        biased in favor of these regions. Mixing in these events with
//        a larger sample of the true integrand over several rounds of
//        adapt/reset_stats/repeat will build sensitivity to these
//        regions into the tree that should persist throughout the
//        remainder of its evolution.
//     d) Run large batches of MC jobs in parallel during adaptation
//        and use saveState/mergeState, then adapt on their combined
//        results. Do not adapt separately in each individual thread
//        unless your goal is to generate as diverse a set of trees
//        as possible. The adaptation algorithm in AdaptiveSampler is
//        designed to deterministically produce the optimum sampling
//        tree in the large-N limit, so the best trees are those that
//        are based on the largest statistical sample. Do not try to
//        generate a super-high performance tree by producing many of
//        them and chosing the one with the highest efficiency -- it
//        may have high efficiency, but by doing this you have almost
//        guaranteed that your result will have a large bias.
//
// 7) The user can now set up the problem with higher dimension than
//    the number of random numbers needed per event, with the extra
//    dimensions allocated to user-generated random variables. These
//    must be transformed by the user to fall within the domain [0,1]
//    and passed to sample() and feedback() in the first nfixed
//    elements of the u vector (first argument), with nfixed supplied
//    as the last argument.
//
// 8) If sample() is called with nfixed > 0 then a pass through the
//    entire cell tree needs to be made to eliminate all cells from
//    the subset sum that do not intersect with the hyperplane that
//    is defined by the first nfixed components of u. Those cells
//    that do intersect with the fixed-u hyperplane must have their
//    subsets converted into a density projected onto the hyperplane
//    and the sum computed to normalize the total probability to 1.
//    These updated subsets are computed for a given fixed u-vector
//    and saved in the ss elements of the cells in the tree upon
//    each event. In this way, the existing algorithms for sample()
//    and feedback() still work without modification, and findCell()
//    is the only thing that needs substantial modification.

#include <assert.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <exception>
#include "AdaptiveSampler.hh"

int AdaptiveSampler::verbosity = 3;

AdaptiveSampler::AdaptiveSampler(int dim, Uniform01 user_generator)
 : fNdim(dim),
   fMaximum_depth(30),
   fMaximum_cells(1000000),
   fSampling_threshold(25),
   fEfficiency_target(0.9),
   fMinimum_sum_wI2_delta(0),
   fNfixed(0)
{
   fFixed_u0 = new double[dim];
   fFixed_u1 = new double[dim];
   fRandom = user_generator;
   fTopCell = new Cell(dim);
   fTopCell->subset = 1;
   reset_stats();
}

AdaptiveSampler::AdaptiveSampler(const AdaptiveSampler &src)
{
   fNfixed = 0;
   fNdim = src.fNdim;
   fRandom = src.fRandom;
   fFixed_u0 = new double[fNdim];
   fFixed_u1 = new double[fNdim];
   fTopCell = new Cell(*src.fTopCell);
   fSampling_threshold = src.fSampling_threshold;
   fEfficiency_target = src.fEfficiency_target;
   fMaximum_cells = src.fMaximum_cells;
   fMaximum_depth = src.fMaximum_depth;
   fMinimum_sum_wI2_delta = src.fMinimum_sum_wI2_delta;
}

AdaptiveSampler::~AdaptiveSampler()
{
   delete [] fFixed_u0;
   delete [] fFixed_u1;
   delete fTopCell;
}

AdaptiveSampler AdaptiveSampler::operator=(const AdaptiveSampler &src)
{
   return AdaptiveSampler(src);
}

double AdaptiveSampler::sample(double *u, int nfixed)
{
   for (int i=0; i < nfixed; ++i) {
      if (u[i] < 0 || u[i] >= 1) {
         std::cerr << "AdaptiveSampler::sample error - "
                   << "fixed parameter " << i  << " is "
                   << "outside the allowed interval [0,1), "
                   << "value is " << u[i] << std::endl;
         return 0;
      }
   }
   int depth;
   double *u0 = new double[fNdim];
   double *u1 = new double[fNdim];
   double *uu = new double[fNdim + 1];
   (*fRandom)(fNdim - nfixed + 1, uu + nfixed);
   Cell *cell = findCell(uu[fNdim], depth, u0, u1, u, nfixed);
   for (int i=nfixed; i < fNdim; ++i) {
      u[i] = uu[i] * u0[i] + (1 - uu[i]) * u1[i];
   }
   delete [] u0;
   delete [] u1;
   delete [] uu;
   // WARNING: strong condition is assumed here!
   // Provided the rule is respected that all splits
   // along axes 0..nfixed are subsetted equally into
   // thirds, the following weight is correct even in
   // the case where nfixed > 0.
   double dNu = (depth > 0)? pow(1/3., depth) : 1;
   return dNu / cell->subset;
}

double AdaptiveSampler::sum_subsets(const double *u, int nfixed)
{
   double *u0 = new double[fNdim];
   double *u1 = new double[fNdim];
   std::fill(u0, u0 + fNdim, 0);
   std::fill(u1, u1 + fNdim, 1);
   std::fill(fFixed_u0, fFixed_u0 + fNdim, 0);
   std::fill(fFixed_u1, fFixed_u1 + fNdim, 1);

   class Branch {
    public:
      Cell *cell;
      int subdiv;
      int ndim;
      double *cell_u0;
      double *cell_u1;
      Branch(Cell *c, int s, double *u0, double *u1, int dim)
       : cell(c), subdiv(s), ndim(dim)
      {
         cell_u0 = new double[ndim];
         cell_u1 = new double[ndim];
         for (int i=0; i < ndim; ++i) {
            cell_u0[i] = u0[i];
            cell_u1[i] = u1[i];
         }
      }
      Branch(const Branch &src)
       : cell(src.cell), subdiv(src.subdiv), ndim(src.ndim)
      {
         cell_u0 = new double[ndim];
         cell_u1 = new double[ndim];
         for (int i=0; i < ndim; ++i) {
            cell_u0[i] = src.cell_u0[i];
            cell_u1[i] = src.cell_u1[i];
         }
      }
      ~Branch() {
         delete [] cell_u0;
         delete [] cell_u1;
      }
   };

   std::vector<Branch> tree;
   tree.push_back(Branch(fTopCell,0,u0,u1,fNdim));
   double ss = 0;
   while (tree.size() > 0) {
      int diva = tree.back().cell->divAxis;
      if (diva >= 0 && diva < nfixed) {
         double du = (u1[diva] - u0[diva]) / 3;
         int s = int((u[diva] - u0[diva]) / du);
         assert (s >= 0 && s <= 2);
         u0[diva] += s * du;
         u1[diva] = u0[diva] + du;
         if (u0[diva] > fFixed_u0[diva])
            fFixed_u0[diva] = u0[diva];
         if (u1[diva] < fFixed_u1[diva])
            fFixed_u1[diva] = u1[diva];
         tree.back().subdiv = 3;
         Cell *c = tree.back().cell->subcell[s];
         tree.push_back(Branch(c,0,u0,u1,fNdim));
      }
      else if (diva >= 0) {
         double du = (u1[diva] - u0[diva]) / 3;
         int s = tree.back().subdiv++;
         u0[diva] += s * du;
         u1[diva] = u0[diva] + du;
         Cell *c = tree.back().cell->subcell[s];
         tree.push_back(Branch(c,0,u0,u1,fNdim));
      }
      else {
         ss += tree.back().cell->subset;
         tree.back().cell->ss = ss;
         tree.pop_back();
         while (tree.size() > 0 && tree.back().subdiv == 3) {
            tree.back().cell->ss = ss;
            tree.pop_back();
         }
         if (tree.size() > 0) {
            for (int i=0; i < fNdim; ++i) {
               u0[i] = tree.back().cell_u0[i];
               u1[i] = tree.back().cell_u1[i];
            }
         }
      }
   }
   delete [] u0;
   delete [] u1;
   return ss;
}

AdaptiveSampler::Cell *AdaptiveSampler::findCell(double ucell,
                                                 int &depth, 
                                                 double *u0,
                                                 double *u1,
                                                 const double *u,
                                                 int nfixed)
{
   if (fTopCell->ss == 0 || fNfixed != nfixed)
      sum_subsets(u, nfixed);
   for (int i=0; i < nfixed; ++i) {
      if (u[i] < fFixed_u0[i] || u[i] > fFixed_u1[i])
         sum_subsets(u, nfixed);
   }
   
   std::fill(u0, u0 + fNdim, 0);
   std::fill(u1, u1 + fNdim, 1);
   depth = 0;
   Cell *sel = fTopCell;
   double usearch = ucell * sel->ss;
   while (sel->divAxis > -1) {
      int i;
      for (i=0; i < 2; i++)
         if (usearch < sel->subcell[i]->ss)
            break;
      assert (usearch < sel->subcell[i]->ss);
      int j = sel->divAxis;
      double du = (u1[j] - u0[j]) / 3;
      u0[j] = u0[j] + du * i;
      u1[j] = u0[j] + du;
      sel = sel->subcell[i];
      ++depth;
   }
   return sel;
}

AdaptiveSampler::Cell *AdaptiveSampler::findCell(const double *u,
                                                 int &depth,
                                                 double *u0,
                                                 double *u1) const
{
   std::fill(u0, u0 + fNdim, 0);
   std::fill(u1, u1 + fNdim, 1);
   depth = 0;
   double ucell = 0;
   Cell *sel = fTopCell;
   while (sel->divAxis > -1) {
      int j = sel->divAxis;
      double du = (u1[j] - u0[j]) / 3;
      int i = floor((u[j] - u0[j]) / du);
      i = (i < 0)? 0 : (i < 3)? i : 2;
      u0[j] = u0[j] + du * i;
      u1[j] = u0[j] + du;
      ucell += (i > 0)? sel->subcell[0]->subset : 0;
      ucell += (i > 1)? sel->subcell[1]->subset : 0;
      if (sel->subcell[i] == 0) {
         std::cerr << "bad romans!" << std::endl;
      }
      sel = sel->subcell[i];
      ++depth;
   }
   return sel;
}

void AdaptiveSampler::feedback(const double *u, double wI)
{
   int depth;
   double *u0 = new double[fNdim];
   double *u1 = new double[fNdim];
   Cell *cell = findCell(u, depth, u0, u1);
   cell->nhit += 1;
   cell->sum_wI += wI;
   double wI2 = wI*wI;
   cell->sum_wI2 += wI2;
   cell->sum_wI4 += wI2*wI2;
   int jbin = 0;
   for (int n = fNdim - 1; n >= 0; --n) {
      double i = floor(3 * (u[n] - u0[n]) / (u1[n] - u0[n]));
      i = (i < 0)? 0 : (i < 3)? i : 2;
      jbin = jbin * 3 + i;
   }
   cell->sum_wI2u[jbin] += wI2;
   delete [] u0;
   delete [] u1;
}

int AdaptiveSampler::adapt()
{
   int ncells = fTopCell->sum_stats(fNfixed);
   double sum_wI2_target = pow(fTopCell->sum_wI, 2) / 
                              (fTopCell->nhit * fEfficiency_target);
   double sum_wI2_error = sqrt(fTopCell->sum_wI4) + 1e-99;
   double gain_sig = (fTopCell->sum_wI2 - sum_wI2_target) / sum_wI2_error;
   if (gain_sig > 3) {
      if (verbosity > 0)
         std::cout << "adapt - sample statistics indicate that significant"
                   << " gains (" << gain_sig << " sigma)" << std::endl
                   << "in sampling efficiency can be achieved through"
                   << " further optimization, beginning optimization pass."
                   << std::endl;
   }
   else {
      if (verbosity > 0)
         std::cout << "adapt - sample statistics indicate that significant"
                   << " gains in sampling efficiency cannot be achieved through"
                   << " further optimization,"
                   << "(" << gain_sig << " sigma)"
                   << std::endl;
      return ncells;
   }

   fMinimum_sum_wI2_delta = (fTopCell->sum_wI2 - sum_wI2_target) / 
                            (fMaximum_cells - ncells);
   std::vector<int> cellIndex;
   int newcells = recursively_update(cellIndex);
   if (verbosity > 0) {
      if (newcells) {
         if (verbosity > 0)
            std::cout << "adapt - " << newcells << " added this pass,"
                      << " new count is " << getNcells()
                      << std::endl;
      }
      else {
         if (verbosity > 0)
            std::cout << "adapt - no changes this pass,"
                      << " count remains unchanged at " << getNcells()
                      << std::endl;
      }
   }
   if (verbosity > 1) {
      std::cout << "present efficiency " << getEfficiency()
                << ", new value predicted to be " << getEfficiency(true)
                << std::endl;
   }
   return newcells;
}

int AdaptiveSampler::recursively_update(std::vector<int> index)
{
   Cell *cell = fTopCell;
   unsigned int depth;
   std::stringstream cellpath;
   for (depth=0; depth < index.size(); ++depth) {
      cellpath << "/" << cell->divAxis << ":" << index[depth];
      cell = cell->subcell[index[depth]];
   }

   int count = 0;
   double efficiency = 0;
   double sum_wI2_sig = 0;
   if (cell->nhit > 0 && cell->sum_wI2 > 0) {
      sum_wI2_sig = cell->sum_wI2 / sqrt(cell->sum_wI4);
      efficiency = pow(cell->sum_wI, 2) / (cell->nhit * cell->sum_wI2);
   }
   if (cell->divAxis < 0) {

      // This is the adaptation algorithm: all of the power
      // of AdaptiveSampler lies in this little bit of code.

      if (verbosity > 3) {
         std::cout << "checking if we should partition " << cellpath.str()
                   << " with sum_wI2=" << cell->sum_wI2
                   << " +/- " << sqrt(cell->sum_wI4)
                   << std::endl;
      }
      if (sum_wI2_sig > fSampling_threshold &&
          efficiency < fEfficiency_target &&
          depth < fMaximum_depth)
      {
         double new_sum_wI2s = 0;
         class Slice {
          public:
            double sum_wI2[3];
            double sum_wI4[3];
            Slice() {
               for (int i=0; i < 3; ++i) {
                  sum_wI2[i] = 0;
                  sum_wI4[i] = 0;
               }
            }
         };
         Slice *slice = new Slice[fNdim];
         int *jbase3 = new int[fNdim];
         std::fill(jbase3, jbase3 + fNdim, 0);
         for (int j=0; j < cell->n3dim; ++j) {
            new_sum_wI2s += sqrt(cell->sum_wI2u[j]);
            for (int n=0; n < fNdim; ++n) {
               slice[n].sum_wI2[jbase3[n]] += cell->sum_wI2u[j];
               slice[n].sum_wI4[jbase3[n]] += pow(cell->sum_wI2u[j], 2) *
                                              cell->n3dim / cell->nhit;
            }
            for (int n=0; n < fNdim; ++n) {
               jbase3[n] += 1;
               if (jbase3[n] > 2)
                  jbase3[n] = 0;
               else
                  break;
            }
         }
         double new_sum_wI2 = new_sum_wI2s * new_sum_wI2s / cell->n3dim;
         double sum_wI2_sigma = sqrt(cell->sum_wI4);
         if (new_sum_wI2 > cell->sum_wI2 + sum_wI2_sigma * 1e-3) {
            std::cerr << "AdaptiveSampler::recursive_update error - "
                      << "new_sum_wI2 > old sum_wI2, this cannot be!"
                      << std::endl
                      << "  old sum_wI2 = " << cell->sum_wI2 << std::endl
                      << "  new sum_wI2 = " << new_sum_wI2 << std::endl
                      << std::endl;
            exit(8);
         }
         if (verbosity > 3) {
            std::cout << "checking if to split at level " << depth
                      << " with sum_wI2 = " << cell->sum_wI2
                      << " +/- " << sqrt(cell->sum_wI4)
                      << " compared to " << new_sum_wI2
                      << " with fMinimum_sum_wI2_delta = "
                      << fMinimum_sum_wI2_delta
                      << std::endl;
         }
         double sum_wI2_delta = cell->sum_wI2 - new_sum_wI2 -
                                fSampling_threshold * sqrt(cell->sum_wI4);
         if (sum_wI2_delta > fMinimum_sum_wI2_delta) {
            int best_axis = -1;
            double best_sum_wI2 = 1e99;
            for (int n=0; n < fNdim; ++n) {
               double sum_wI2 = pow(sqrt(slice[n].sum_wI2[0]) +
                                     sqrt(slice[n].sum_wI2[1]) +
                                     sqrt(slice[n].sum_wI2[2]), 2) / 3;
               if (sum_wI2 < best_sum_wI2) {
                  best_sum_wI2 = sum_wI2;
                  best_axis = n;
               }
            }
            if (best_sum_wI2 > cell->sum_wI2) {
               std::cerr << "AdaptiveSampler::recursive_update error - "
                         << "best_sum_wI2 > sum_wI2, this cannot be!"
                         << std::endl;
            }
            else if (best_sum_wI2 > 0) {
               cell->divAxis = best_axis;
               int ibase = 1;
               for (int i=0; i < best_axis; ++i)
                  ibase *= 3;
               for (int n=0; n < 3; ++n) {
                  cell->subcell[n] = new Cell(fNdim);
                  cell->subcell[n]->nhit = cell->nhit / 3;
                  cell->subcell[n]->sum_wI = cell->sum_wI / 3;
                  cell->subcell[n]->sum_wI2 = slice[best_axis].sum_wI2[n];
                  cell->subcell[n]->sum_wI4 = slice[best_axis].sum_wI4[n];
                  cell->subcell[n]->subset = cell->subset / 3;
                  std::fill(jbase3, jbase3 + fNdim, 0);
                  int jstride = 1;
                  for (int j=0; j < cell->n3dim; ++j) {
                     int jsrc = j + jstride * (n - jbase3[best_axis]);
                     cell->subcell[n]->sum_wI2u[j] = cell->sum_wI2u[jsrc] / 3;
                     for (int i=0; i < fNdim; ++i) {
                        jbase3[i] += 1;
                        if (jbase3[i] > 2)
                           jbase3[i] = 0;
                        else
                           break;
                     }
                  }
               }
               cell->nhit = 0;
               cell->sum_wI = 0;
               cell->sum_wI2 = 0;
               cell->sum_wI4 = 0;

               if (verbosity > 2) {
                  std::cout << "splitting this cell[";
                  for (unsigned int i=0; i < depth; i++) {
                     std::cout << ((i > 0)? "," : "") << index[i];
                  }
                  std::cout << "] along axis " << best_axis << std::endl;
               }
               count += 1;
            }
            else {
               if (verbosity > 3) {
                  std::cout << "nope, missing statistics to find best axis!"
                            << std::endl;
               }
            }
         }
         else {
            if (verbosity > 3) {
               std::cout << "nope, fails to make the threshold!" << std::endl;
            }
         }
         delete [] jbase3;
         delete [] slice;
      }
      else if (sum_wI2_sig < fSampling_threshold) {
         if (verbosity > 3) {
            std::cout << "nope, fails the statistical significance test!"
                      << std::endl;
         }
      }
      else if (efficiency >= fEfficiency_target) {
         if (verbosity > 3) {
            std::cout << "nope, efficiency is " << efficiency
                      << ", already meets the target value!"
                      << std::endl;
         }
      }
      else {
         if (verbosity > 3) {
            std::cout << "nope, this cell cannot be split any further!"
                      << std::endl;
         }
      }
   }
   else {
      index.push_back(0);
      for (int i=0; i < 3; ++i) {
         index[depth] = i;
         count += recursively_update(index);
      }
      index.pop_back();
   }
   return count;
}

void AdaptiveSampler::reset_stats()
{
   fTopCell->reset_stats();
}

long int AdaptiveSampler::getNsample() const
{
   fTopCell->sum_stats(fNfixed);
   return fTopCell->nhit;
}

double AdaptiveSampler::getWItotal() const
{
   fTopCell->sum_stats(fNfixed);
   return fTopCell->sum_wI;
}

double AdaptiveSampler::getWI2total(bool optimized)
{
   fTopCell->sum_stats(fNfixed);
   if (optimized) {
      optimize_tree();
      return fTopCell->opt_wI2;
   }
   return fTopCell->sum_wI2;
}

double AdaptiveSampler::getEfficiency(bool optimized)
{
   fTopCell->sum_stats(fNfixed);
   if (fTopCell->nhit == 0)
      return 0;
   else if (optimized) {
      optimize_tree();
      return pow(fTopCell->sum_wI, 2) / 
             (fTopCell->opt_wI2 * fTopCell->opt_nhit + 1e-99);
   }
   return pow(fTopCell->sum_wI, 2) / 
          (fTopCell->sum_wI2 * fTopCell->nhit);
}

double AdaptiveSampler::getResult(double *error,
                                  double *error_uncertainty)
{
   fTopCell->sum_stats(fNfixed);
   double eff = pow(fTopCell->sum_wI, 2) /
                   (fTopCell->sum_wI2 * fTopCell->nhit);
   double result = (eff > 0)? fTopCell->sum_wI / fTopCell->nhit : 0;
   if (error) {
      if (eff > 0)
         *error = sqrt((1 - eff) * fTopCell->sum_wI2) / fTopCell->nhit;
      else
         *error = 0;
   }
    
   if (error_uncertainty) {
     if (eff > 0)
         *error_uncertainty = sqrt((1 - eff) / fTopCell->sum_wI2) / 2 *
                              sqrt(fTopCell->sum_wI4) / fTopCell->nhit;
     else
         *error_uncertainty = 0;
   }
   return result;
}

double AdaptiveSampler::getReweighted(double *error,
                                      double *error_uncertainty)
{
   fTopCell->sum_stats(fNfixed);
   optimize_tree();
   double eff = pow(fTopCell->sum_wI, 2) / 
                (fTopCell->opt_wI2 * fTopCell->opt_nhit + 1e-99);
   double result = (eff > 0)? fTopCell->sum_wI / fTopCell->opt_nhit : 0;
   if (error) {
      if (eff > 0)
         *error = sqrt((1 - eff) * fTopCell->opt_wI2) / 
                  (fTopCell->opt_nhit + 1e-99);
      else
         *error = 0;
   }
   if (error_uncertainty) {
      if (eff > 0)
         *error_uncertainty = sqrt((1 - eff) / fTopCell->opt_wI2) / 2 *
                              sqrt(fTopCell->opt_wI4) /
                              (fTopCell->opt_nhit + 1e-99);
      else
         *error_uncertainty = 0;
   }
   return result;
}

int AdaptiveSampler::getNcells() const
{
   return fTopCell->sum_stats(fNfixed);
}

void AdaptiveSampler::setAdaptation_sampling_threshold(double threshold)
{
   fSampling_threshold = threshold;
}

double AdaptiveSampler::getAdaptation_sampling_threshold() const
{
   return fSampling_threshold;
}

void AdaptiveSampler::setAdaptation_efficiency_target(double target)
{
   fEfficiency_target = target;
}

double AdaptiveSampler::getAdaptation_efficiency_target() const
{
   return fEfficiency_target;
}

void AdaptiveSampler::setAdaptation_maximum_depth(int depth)
{
   fMaximum_depth = depth;
}

int AdaptiveSampler::getAdaptation_maximum_depth() const
{
   return fMaximum_depth;
}

void AdaptiveSampler::setAdaptation_maximum_cells(int ncells)
{
   fMaximum_cells = ncells;
}

int AdaptiveSampler::getAdaptation_maximum_cells() const
{
   return fMaximum_cells;
}

void AdaptiveSampler::optimize_tree()
{
   fTopCell->sum_stats(fNfixed);
   if (fTopCell->sum_wI > 0) {
      fTopCell->opt_nhit = fTopCell->nhit;
      fTopCell->opt_subset = 1;
      fTopCell->optimize(fNfixed);
   }
}

void AdaptiveSampler::display_tree(bool optimized)
{
   double *u0 = new double[fNdim];
   double *u1 = new double[fNdim];
   std::fill(u0, u0 + fNdim, 0);
   std::fill(u1, u1 + fNdim, 1);
   display_tree(fTopCell, 1, 0, u0, u1, optimized);
   delete [] u0;
   delete [] u1;
}

double AdaptiveSampler::display_tree(Cell *cell, double subset, int level,
                                     double *u0, double *u1, bool optimized)
{
   for (int i=0; i < level; ++i)
      std::cout << " ";

   char numeric[80];
   std::cout << level << ": ";
   double cell_subset = (optimized)? cell->opt_subset : cell->subset;
   snprintf(numeric, 80, "%9.5e", cell_subset);
   std::cout << numeric;
   if (cell_subset > subset * 99.95) {
      std::cout << "(100%) ";
   }
   else {
      snprintf(numeric, 30, "(%4.1f) ", 100 * cell_subset / subset);
      std::cout << numeric;
   }

   if (cell->divAxis > -1) {
      double ssum = 0;
      double u0m = u0[cell->divAxis];
      double u1m = u1[cell->divAxis];
      double du = (u1m - u0m) / 3.;
      snprintf(numeric, 80, " [%15.13f,%15.13f] ", u0m, u1m);
      std::cout << " split along axis " << cell->divAxis << numeric;
      double dlev = -log(du) / log(3.);
      snprintf(numeric, 80, " (1/3**%5.3f) ", dlev);
      std::cout << numeric;
      std::cout << ((optimized)? cell->opt_nhit : cell->nhit) << std::endl;
      for (int n=0; n < 3; ++n) {
         u0[cell->divAxis] = u0m + du * n;
         u1[cell->divAxis] = u0m + du * (n + 1);
         ssum += display_tree(cell->subcell[n], cell_subset, 
                              level+1, u0, u1, optimized);
      }
      if (fabs(ssum - cell_subset) > 1e-15 * cell_subset) {
         std::cerr << "Error in AdaptiveSampler::display_tree - "
                      "subcell subsets fail to obey the sum rule, "
                      "tree is invalid !!!" << std::endl;
      }
      u0[cell->divAxis] = u0m;
      u1[cell->divAxis] = u1m;
   }
   else {
      std::cout << ((optimized)? cell->opt_nhit : cell->nhit )
                << " " << cell->sum_wI
                << " " << ((optimized)? cell->opt_wI2 : cell->sum_wI2)
                << " " << ((optimized)? cell->opt_wI4 : cell->sum_wI4)
                << std::endl;
   }
   return cell_subset;
}

int AdaptiveSampler::saveState(const std::string filename, bool optimized) const
{
   std::ofstream fout(filename);
   fout << "fNdim=" << fNdim << std::endl;
   fout << "fNfixed=" << fNfixed << std::endl;
   fout << "fSampling_threshold=" << fSampling_threshold << std::endl;
   fout << "fMaximum_depth=" << fMaximum_depth << std::endl;
   fout << "fMaximum_cells=" << fMaximum_cells << std::endl;
   fout << "fEfficiency_target=" << fEfficiency_target << std::endl;
   fout << "=" << std::endl;
   int ncells = fTopCell->serialize(fout, optimized);
   return (ncells > 0);
}

int AdaptiveSampler::mergeState(const std::string filename)
{
   std::ifstream fin(filename);
   if (!fin.is_open()) {
      std::cerr << "AdaptiveSampler::mergeState error - "
                << "unable to open " << filename << " for input."
                << std::endl;
      return 0;
   }
   std::map<std::string,double> keyval;
   while (true) {
      std::string key;
      std::getline(fin, key, '=');
      if (key.size() == 0) {
         std::getline(fin, key);
         break;
      }
      fin >> keyval[key];
      std::getline(fin, key);
   }
   if (fNdim != keyval["fNdim"]) {
      std::cerr << "AdaptiveSampler::mergeState error - "
                << "cannot merge with state saved in " << filename
                << " because they have different dimensions."
                << std::endl;
      return 0;
   }
   try {
      fNfixed = keyval.at("fNfixed");
      fSampling_threshold = keyval.at("fSampling_threshold");
      fMaximum_depth = keyval.at("fMaximum_depth");
      fMaximum_cells = keyval.at("fMaximum_cells");
      fEfficiency_target = keyval.at("fEfficiency_target");
   }
   catch (std::exception e) {
      std::cerr << "AdaptiveSampler::mergeState warning - "
                << "required keyword missing in " << filename
                << std::endl;
   }
   int ncells = fTopCell->deserialize(fin);
   return (ncells > 0);
}

int AdaptiveSampler::restoreState(const std::string filename)
{
   if (fTopCell->divAxis > -1) {
      for (int i=0; i < 3; ++i) {
         delete fTopCell->subcell[i];
         fTopCell->subcell[i] = 0;
      }
   }
   fTopCell->divAxis = -1;
   reset_stats();
   return mergeState(filename);
}

int AdaptiveSampler::getVerbosity()
{
   return verbosity;
}

void AdaptiveSampler::setVerbosity(int verbose)
{
   verbosity = verbose;
}
