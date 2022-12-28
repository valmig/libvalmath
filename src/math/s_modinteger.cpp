
#include <s_modinteger.h>
#include <error.h>

namespace val
{

integer s_modinteger::q=integer(1);

s_modinteger::s_modinteger(const integer& a)
{
    x = a%q;
    if (x.signum()==-1) x+=q;
}

s_modinteger s_modinteger:: operator+(const s_modinteger& a) const
{
 s_modinteger b;
 b.x=x+a.x;
 if (b.x>=q) b.x-=q;
 return b;
}

s_modinteger s_modinteger:: operator -() const
{
 s_modinteger b;
 b.x = q-x;
 return b;
}

s_modinteger s_modinteger:: operator -(const s_modinteger& a) const
{
 s_modinteger b;
 b.x = x -a.x;
 if (b.x.signum()==-1) b.x+=q;
 return b;
}

s_modinteger s_modinteger::operator *(const s_modinteger& a) const
{
    s_modinteger b;
    b.x=x*a.x;
    if (b.x>q) b.x%=q;
    return b;
}


s_modinteger inv(const s_modinteger &a)
{
 integer x0,x1,one(1),zero;
 s_modinteger b;

 x1=euclid(s_modinteger::q,a.x,x0,b.x);
 if ((x1 != one) || (a.x==zero)){
     Error::error("\nERROR s_modinteger: inv(a): a is not a unit!!");
 }
 if (b.x.signum()==-1) b.x+=s_modinteger::q;
 return b;
}


s_modinteger s_modinteger:: operator /(const s_modinteger& a) const
{
 integer x1,y1,d,zero;
 d=euclid(a.x,q,x1,y1);
 if (x%d!=zero) {
    Error::error("\nERROR: s_modinteger::operator/ : division by zero!");
 }

 return s_modinteger(x1)*s_modinteger(EDIV(x,d));
}


std::istream& operator >>(std::istream& is, s_modinteger& a)
{
 is>>a.x;
 a.x%=s_modinteger::q;
 if (a.x.signum()==-1) a.x+=s_modinteger::q;
 return is;
}

std::ostream& operator <<(std::ostream& os,const s_modinteger &a)
{
 os<<a.x;
 return os;
}


} // end namespace val


