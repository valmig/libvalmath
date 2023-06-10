#include <rand.h>
#include <ctime>
#include <integer.h>

namespace val
{

namespace hilfrand
{
 unsigned numberofbits=8*sizeof(unsigned int);
 int max(int a, int b) {if (a>b) return a; else return b;}
 int min(int a,int b) {if (a<b) return a; else return b;}
}


void initialize_random()
{
    std::srand( unsigned(std::time( NULL )) );
}

unsigned random()
{
 unsigned r,bit=0,z=0;

 while (bit<hilfrand::numberofbits) {
	 r = rand()% 2;  // => 0 <= r < 2  Pseudo random number;
	 z=z<<1;
	 z+=r;
	 bit++;
 }
 return z;
}


int random(int a)
{
    if (!a) return 0;
    int b=a,z;

    if (a<0) b*=-1;
    z=random();

    if (z<0) z*=-1;
    z%=b;
    if (a<0) z*=-1;
    return z;
}

int random(int a,int b)
{
    int c=hilfrand::min(a,b),d=hilfrand::max(a,b),z=0,k;

    if (c==d) return c;

    k=d-c;
    z=random(k);
    z+=c;
    return z;
}



int RandomClass::randominitialized=0;



RandomClass::RandomClass()
{
    if (randominitialized==0) {
        initialize_random();
        randominitialized=1;
    }
}


unsigned RandomClass::random()
{
    return val::random();
}

int RandomClass::random(int a)
{
    return val::random(a);
}


int RandomClass::random(int a,int b)
{
    return val::random(a,b);
}



#ifdef INTEGER_H

val::integer integerRandClass::integerrandom(int n)
{
    if (n<=0) return val::integer(0);

    int i;
    unsigned z=random();

    val::integer a(z);

    for (i=1;i<n;i++) {
        a=val::shift(a,1);
        z=random();
        a+=val::integer(z);
    }
    i=random(0,2);
    if (i) a.changesign();

    return a;
}

#endif // INTEGER_H



}// end namespace val
