#ifndef MYTIME_H
#define MYTIME_H

//#include <ctime>
#include <chrono>
#include <val_basics.h>

namespace val
{



class DLL_PUBLIC ChronoClass
{
private:
    std::chrono::time_point<std::chrono::system_clock> a;
public:
    ChronoClass() : a( std::chrono::system_clock::now()) {}
    double operator() (); // Time difference between last call of () and object creation.
};

} // end namespace val


#endif
