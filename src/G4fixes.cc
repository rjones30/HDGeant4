#include <iostream>
#include <boost/python.hpp>

void G4fixes_init(int verbose=0)
{
   if (verbose)
      std::cout << "G4fixes_init - successfully loaded" << std::endl;
}

BOOST_PYTHON_MODULE(libG4fixes)
{
   boost::python::def("G4fixes_init", &G4fixes_init);
}
