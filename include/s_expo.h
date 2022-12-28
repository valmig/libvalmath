#ifndef S_EXPO_H_INCLUDED
#define S_EXPO_H_INCLUDED

// class for multivariate exponents (terms) of polynomials with static dimension.
// Order relation defined by a quadratic int - matrix, or lexicographic (ordtype = -1), or DegRevLex ( ordtype = -2 ),
// or DegLex ( ordtype = 0 ).

#include <matrix.h>


namespace val
{

class s_expo;
template <class T> class s_polynom;

DLL_PUBLIC s_expo gcd(const s_expo&,const s_expo&);
DLL_PUBLIC s_expo lcm(const s_expo&,const s_expo&);
DLL_PUBLIC int expocompare(const s_expo&,const s_expo&);


class DLL_PUBLIC s_expo
{
private:
    int *coeff;
    static int dim;
    static int ordtype;
    static matrix<int> ordmatrix;
    static const matrix<int> &c_ordmatrix;
public:
    //
    s_expo() {coeff = new int[dim];}         // without initialization
    explicit s_expo(int c);                  // *this = (c,...,c)
    s_expo(const s_expo&);
    s_expo(s_expo&&);
    //
    ~s_expo() {delete[] coeff;coeff=NULL;}
    //
    const s_expo& operator =(const s_expo&);
    const s_expo& operator =(s_expo&&);
    //
    static int getdim() {return dim;}
    static int dimension() {return dim;}
    static int getordtype() {return ordtype;}
    static void setordtype(int ord) {ordtype=ord;}
    //static const matrix<int>& getordmatrix() {return ordmatrix;}
    static matrix<int>& getordmatrix() {return ordmatrix;}
    static void setdim(int n) {if (n<=0) return; else dim=n;}
    static void setordmatrix(const matrix<int> &M) {ordmatrix=M;}
    static void setordmatrix(matrix<int> &&M) {ordmatrix=std::move(M);}
    //
    int operator ==(const s_expo&) const;
    int operator !=(const s_expo& Y) const {return !operator==(Y);}
    int operator <(const s_expo&) const;
    int operator >(const s_expo& Y) const {return Y<*this;}
    //
    int& operator[](int);
    int operator[](int) const;
    //
    s_expo operator *(const s_expo&) const;
    const s_expo& operator *=(const s_expo&);
    int operator |(const s_expo&) const;
    s_expo operator /(const s_expo&) const;
    const s_expo& operator /=(const s_expo&);
    friend s_expo gcd(const s_expo&,const s_expo&);
    friend s_expo lcm(const s_expo&,const s_expo&);
    friend int expocompare(const s_expo&,const s_expo&);
    //
    int totdeg() const;
    //
    friend DLL_PUBLIC std::istream& operator >>(std::istream&,s_expo&);
    friend DLL_PUBLIC std::ostream& operator <<(std::ostream&,const s_expo&);
    //
    template <class T> friend class s_polynom;
};



} // end namespace



#endif // S_EXPO_H_INCLUDED
