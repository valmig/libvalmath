#ifndef ANALYSIS_H_INCLUDED
#define ANALYSIS_H_INCLUDED

#include <pol.h>

#define __int64 long long

namespace val
{

template <class T> class vector;
template <class T> class complex_type;
typedef complex_type<double> complex;


DLL_PUBLIC void sign_exp_mantissa(const double&,int&,int&,__int64&);
DLL_PUBLIC int isNaN(const double&);
DLL_PUBLIC void normalize(const double& a,double& m,int& exp); // a = (2^exp) * m , mit 0.5 <= |m| < 1.


typedef double (d_function)(const double &);




DLL_PUBLIC double integral(const pol<double> &p,const double& a,const double& b);


// Real zeros of a double-polynomial in [x1,x2]. If x1 = x2 = 0 all zeros will be computed.
// zeros will be stored in vector Roots. A polynomial will be returned that could eventually be used in function roots.
DLL_PUBLIC pol<double> realRoots(const pol<double> &g,val::vector<double> &Roots,const double& eps=1e-9,const double &x1=0,const double &x2=0);


// Returns number of computed roots of f. Function fails sometimes in finding all complex zeros of f.
// Real roots are stored in vector Real, complex roots (without conjugates) in vector Comp, so that
// number of computed roots = Real.dimension + 2* Comp.dimension.  
DLL_PUBLIC int roots(const pol<double> &f,val::vector<double> &Real,val::vector<val::complex> &Comp,const double &eps=1e-9);


// Computes local maxima and minima of f in [x1,x2]. Returns number of extrema.
DLL_PUBLIC int extreme(const pol<double>& g,val::vector<double> &Maxima,val::vector<double> &Minima,const double &eps=1e-9,const double &x1=0,const double& x2=0);

extern DLL_PUBLIC const double PI;

//double sqrt(const double& a);
//int abs(int);
//double abs(const double&);
DLL_PUBLIC double power(const double&,int);
DLL_PUBLIC double exp(const double&); // Exponential-function
DLL_PUBLIC double log(const double&); // logarithmus naturalis
DLL_PUBLIC double exp(const double& a,const double &x); //  = a^x
DLL_PUBLIC double log(const double& a,const double &x); // = log_a(x)
DLL_PUBLIC double sin(const double&);
DLL_PUBLIC double cos(const double&);
DLL_PUBLIC double tan(const double&);
DLL_PUBLIC double arctan(const double&);
DLL_PUBLIC double arcsin(const double&);
DLL_PUBLIC double arccos(const double&);
DLL_PUBLIC double sinh(const double&);
DLL_PUBLIC double cosh(const double&);
DLL_PUBLIC double tanh(const double&);
DLL_PUBLIC double arsinh(const double&);
DLL_PUBLIC double arcosh(const double&);
DLL_PUBLIC double artanh(const double&);




// Iterative secant method for the computation of a root of the double functor f.
// Return values: -2 number of iteration n not reached, method failed.
//                -1 number of iteration is reached but |f(x1) > eps|
//                 i = number of iterations and |f(x1) <= eps|
template <class T>
int SecantMethod(const T& f,double& x0,double& x1,const double& eps=1e-9,int n=15);





template <class T>
int SecantMethod(const T& f,double& x0,double& x1,const double& eps,int n)
{
    int i;
    double y0,y1,x2;

    y0=f(x0);
    if (abs(y0)<eps) {x1=x0;return 0;}

    for (i=0;i<n;i++) {
        //std::cout<<" i = "<<i<<" , x1  = "<<x1<<std::endl;
        //std::cout<<" f(x1) = "<<f(x1)<<std::endl;
        if (abs(y1=f(x1))<eps) {
            //std::cout<<"\nAnzahl Iterationen in secantmethod: "<<i+1<<std::endl;
            return i;
        }
        //std::cout<<" evaluated!"<<std::endl;
        if (y0==y1) return -2;
        x2= x1-y1 * (x1-x0)/(y1-y0);
        y0=y1;
        x0=x1;
        x1=x2;
    }
    return -1;
}


}


#endif // ANALYSIS_H_INCLUDED
