#ifndef RAND_H
#define RAND_H

// class to generate random numbers

#include <cstdlib>
#include <val_basics.h>


namespace val
{

class integer;
DLL_PUBLIC void initialize_random();

DLL_PUBLIC unsigned random();
DLL_PUBLIC int random(int a);        // number between 0 and a;
DLL_PUBLIC int random(int a,int b);  // number between a and b

class DLL_PUBLIC RandomClass
{
public:
    RandomClass();
    unsigned random();
    int random(int a);        // number between 0 and a;
    int random(int a,int b);  // number between a and b

protected:
    static int randominitialized;
};


class DLL_PUBLIC integerRandClass : public val::RandomClass
{
public:
    val::integer integerrandom(int n=1);     // random number of type integer of length n
};


} // end namespace val




#endif
