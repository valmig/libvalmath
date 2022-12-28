#include <MyTime.h>

namespace val
{


double ChronoClass::operator() ()
{
        std::chrono::time_point<std::chrono::system_clock> b(std::chrono::system_clock::now());
        unsigned s,ms;

        s=std::chrono::duration_cast<std::chrono::seconds>(b-a).count();
        if (s>=60) {
            a=std::chrono::system_clock::now();
            return double(s);
        }
        else {
            ms=std::chrono::duration_cast<std::chrono::milliseconds>(b-a).count();
            a=std::chrono::system_clock::now();
            return (double (ms))/1000.0;
        }
}

} // end namespace val
