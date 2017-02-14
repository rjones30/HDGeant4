//
// GlueXTimer class header
//
// author: richard.t.jones at uconn.edu
// version: january 30, 2017
//
// This class is an alternative to G4Timer that is included in the
// G4 library. I created this one so I could get access to the total
// cpu time used on a per-thread basis. I also wanted to be able to
// suspend and resume the clock and have it accummulate the total
// time used while it was running. It only keeps track of one kind
// of time, which is the cpu time (user + system) consumed.
//
// In the context of the Geant4 event-level multithreading model,
// this class is "shared", ie. is common to all threads.
// This code is written to be thread-safe.
 
#ifndef GlueXTimer_h
#define GlueXTimer_h 1

#define MAX_THREADS 128

#include <time.h>
#include <unistd.h>
#include <assert.h>
#include <string>
#include <map>

#include "G4Threading.hh"
#include "G4AutoLock.hh"

// GlueXTimer is a general tool for finding efficiency bottlenecks
// in a multi-threaded application. All you need to do is to add
// lines like the following to your existing c++ code.
//
//  GLUEXTIMER_START("unique user-defined timer description")
//     ... code to be time goes here ...
//  GLUEXTIMER_STOP("unique user-defined timer description")
//
// Repeated sections of this kind with the same name in user code
// will all be accummulated together in a single aggregated timer.
// The user should make sure than any execution thread that passes
// through a GLUEXTIMER_START also hits the matching GLUEXTIMER_STOP
// or else the timer will go into a runaway state. Repeated execution
// of GLUEXTIME_START without an intervening GLUEXTIMER_STOP, or of
// GLUEXTIMER_STOP without a prior GLUEXTIMER_START will result in
// a warning message being printed to stderr, but it will not 
// otherwise affect the flow of the code.
//
// Of course one may also instantiate GlueXTimer objects explicitly
// and control them by calling Start/Stop and other methods at the
// appropriate time in the code. The above macros are defined to make
// it simple to manage these objects. To print out a summary report
// at the end of execution, simply call GlueXTimer::PrintAll().

#define GLUEXTIMER_START(NAME) {\
   GlueXTimer *timer = GlueXTimer::GetInstance(NAME);\
   if (timer == 0) {\
      timer = new GlueXTimer(NAME);\
   }\
   if (timer->IsStarted())\
      timer->Resume();\
   else\
      timer->Start();\
}
#define GLUEXTIMER_STOP(NAME) {\
   std::string timername(NAME);\
   GlueXTimer *timer = GlueXTimer::GetInstance(NAME);\
   assert(timer != 0);\
   timer->Suspend();\
}
   

class GlueXTimer {
 public:
   GlueXTimer(std::string name);
   ~GlueXTimer();

   void Reset();   // clears the timer, can be invoked at ANY time
   void Start();   // resets the time, starts a new measurement
   void Stop();    // same as Suspend, but prevents further Resume
   void Suspend(); // ends a measurement interval, adds the time to total
   void Resume();  // starts a new measurement interval

   static GlueXTimer *GetInstance(std::string name) {
      if (fInstance.find(name) != fInstance.end())
         return fInstance[name];
      return 0;
   }
   std::string GetName() const {
      return fName;
   }
   int IsStarted() const {
      int thread = G4Threading::G4GetThreadId() + 1;
      return fIsStarted[thread];
   }
   int IsRunning() const {
      int thread = G4Threading::G4GetThreadId() + 1;
      return fIsRunning[thread];
   }
   G4double GetCPUseconds() const {
      int thread = G4Threading::G4GetThreadId() + 1;
      return fRunningTotal[thread];
   }
   G4double GetCPUseconds(int threadId) const {
      return fRunningTotal[threadId+1];
   }
   G4double GetProcessTotal() const;

   static void PrintAll();
   
 protected:
   std::string fName;
   int fIsStarted[MAX_THREADS];
   int fIsRunning[MAX_THREADS];
   double fTimeLastGo[MAX_THREADS];
   double fRunningTotal[MAX_THREADS];
   struct timespec fClock[MAX_THREADS];

   static std::map<std::string, GlueXTimer*> fInstance;
   static G4Mutex fMutex;

 private:
   GlueXTimer(const GlueXTimer &src);
   GlueXTimer &operator=(const GlueXTimer &src);
};

#endif
