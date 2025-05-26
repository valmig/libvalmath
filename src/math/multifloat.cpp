

#include <multifloat.h>

#define __int64 long long


namespace val
{

int hilfmultifloat::highestbit(unsigned x)
{
    if (x==0) return 0;
    int i= 8*sizeof(unsigned);
    unsigned bit = 1<< (i-1);
    for (;i>0;i--) {
        if (bit&x) break;
        bit>>=1;
    }
    return i;
}

int hilfmultifloat::abs(int x)
{
    if (x<0) return -x;
    else return x;
}


multifloat multifloat::NaN(double(0.0/0.0));
multifloat multifloat::InfPos(1.0/0.0);
multifloat multifloat::InfNeg(-1.0/0.0);

int multifloat::prec=1;

const int multifloat::numberofbits = 8*sizeof(unsigned);


multifloat::multifloat(const multifloat& x) : exp(x.exp),laenge(x.laenge)
{
    if (x.mantissa!=nullptr) {
        if (laenge==0) {   // NaN
            mantissa = new unsigned[1];
            mantissa[0] = 0;
        }
        else {
            int i,l=x.abslength();
            mantissa= new unsigned[l];
            for (i=0;i<l;i++) mantissa[i]=x.mantissa[i];
        }

    }
    else mantissa=nullptr;
}


multifloat::multifloat(multifloat&& x)
{
    mantissa=x.mantissa;
    exp=x.exp;laenge=x.laenge;
    x.mantissa=nullptr;x.exp=x.laenge=0;
}



multifloat::multifloat(int x)
{
    if (x==0) {
        exp=laenge=0;
        mantissa=nullptr;
    }
    else {
        exp=0;laenge=1;
        mantissa = new unsigned[1];
        mantissa[0] = hilfmultifloat::abs(x);
        // Normalisiere:
        while (mantissa[0]%2==0) {mantissa[0]>>=1;exp++;}
        if (x<0) laenge*=-1;
    }
}


multifloat::multifloat(const double& a)
{
    //double c=a;
    int sign=1,i;
    //unsigned bit=1;
    __int64 eins(1),bit, b;// b =  *((__int64*) &a);

    val::MemCopy(b, a);

    exp=0;

    // Zeige Mantissa:


    //sign:
    if ((eins<<63) & b) sign=-1;
    //std::cout<<"\nsign= "<<sign;
    // Exponent:
    for (i=62;i>51;i--) {
        if ((eins<<i) &b) exp+=1;
        exp<<=1;
    }
    exp>>=1;
    //std::cout<<"\n exp = "<<exp;
    //Sonderfälle:
    if (exp==0) {laenge=0;mantissa=nullptr;return;} // subnormale werden als 0 behandelt.

    if (exp==2047) {   //+ oder - inf oder nan
        mantissa=new unsigned[1];
        mantissa[0]=0;
        int nan=0;
        for (i=51;i>=0;i--) {
            if ((eins<<i) & b) {nan=1;break;}
        }
        if (nan) {laenge=0;return;}
        else {laenge=sign;}
        //kuerze();
        return;
    }
    exp-=1023;
    //std::cout<<"\nexp in double : "<<exp;
    exp-=52;
    __int64 man(1);
    man<<=52;   // man = 2^52;
    //Bestimme Mantissa von a

    bit=1;
    for (i=0;i<52;i++,bit=bit<<1) {
        //std::cout<<(h&b)<<"   ";
        if ((bit&b) !=0) {man|=bit;}
    }

    //std::cout<<"\nman = "<<man;
    //b<<=12;
    //b>>=12;
    //man|=b;
    //std::cout<<"\nVor Normalisierung man , exp , b "<<man<<"  "<<exp<<"   "<<b;
    //Normalisiere man:
    bit=1;
    i=0;

    while (!(bit&man)) {bit<<=1;i++;exp++;}
    if (i) man>>=i;

    //while ((man%2)==0) {man>>=1;i++;exp++;}

    //std::cout<<"\nNach Normalisierung man , exp "<<man<<"  "<<exp;
    if (52-i<numberofbits) {
        laenge=1;
        mantissa=new unsigned[1];
        mantissa[0]= unsigned(man);
        laenge*=sign;
        return;
    }
    else if (prec>=2) {
        int k=52-i;
        laenge=2;
        mantissa= new unsigned[2];
        mantissa[0]=mantissa[1]=0;
        bit=1;
        for (i=0;i<numberofbits;i++,bit<<=1) mantissa[0] |= unsigned(bit&man);
        man=(__int64) (man>>numberofbits);
        bit=1;
        for (;i<=k;i++,bit<<=1) mantissa[1]|=unsigned(bit&man);
        laenge*=sign;
    }
    else {
        int k=52-i;
        laenge=1;
        mantissa= new unsigned[1];
        // man muss um k-31 bits nach rechts verschoben werden, dann runden und normalisieren.
        exp+=k-(numberofbits-1);
        bit=1<<(k+1-numberofbits);   // (k-31)
        bit&=man;
        man>>=(k+1-numberofbits);    //(k-31)
        mantissa[0]=unsigned(man);
        if (bit) { //nach oben runden;
            mantissa[0]+=1;
            if (mantissa[0]==0) {exp+=32;mantissa[0]=1;}
            else {   // Normalisiere:
                bit=1;i=0;
                while (!(bit&mantissa[0])) {
                    bit<<=1;i++;exp++;
                }
                if (i) mantissa[0]>>=i;
            }
        }
        laenge*=sign;
    }
}



const multifloat& multifloat::operator =(const multifloat& x)
{
    if (mantissa == x.mantissa) return *this;
    if (mantissa!= nullptr) delete[] mantissa;

    exp=x.exp;laenge=x.laenge;
    if (x.mantissa==nullptr) {
        mantissa=nullptr;
        return *this;
    }
    if (laenge==0) { // Fall Nan
        mantissa= new unsigned[1];
        mantissa[0]=0;
        return *this;
    }

    int i,l=hilfmultifloat::abs(laenge);

    mantissa = new unsigned[l];
    for (i=0;i<l;i++) mantissa[i]=x.mantissa[i];
    return *this;
}


const multifloat& multifloat::operator=(multifloat&& x)
{
    if (mantissa == x.mantissa) return *this;
    if (mantissa!=nullptr) delete[] mantissa;
    mantissa=x.mantissa;
    exp=x.exp;laenge=x.laenge;
    x.mantissa=nullptr;x.exp=x.laenge=0;
    return *this;
}



int multifloat::operator==(const multifloat &x) const
{
    if (mantissa==nullptr) return (x.mantissa==nullptr);
    if (isNaN() || isInfNeg() || isInfPos()) return 0;

    if (laenge!=x.laenge) return 0;
    if (exp!=x.exp) return 0;

    int i,l=abslength();
    for (i=0;i<l;i++) if (mantissa[i]!=x.mantissa[i]) return 0;
    return 1;
}


int multifloat::operator <(const multifloat &x) const
{
    if (isNaN() || x.isNaN()) return 0;
    if (x.isInfPos()) {
        if (isInfPos()) return 0;
        else return 1;
    }
    if (isInfNeg()) {
        if (x.isInfNeg()) return 0;
        else return 1;
    }
    if (isInfPos() || x.isInfNeg()) return 0;

    int sign=signum(),signx=x.signum(),hexpo=highestexp(),hexpox=x.highestexp(),l=abslength(),lx=x.abslength();
    int ind,indx,first,last,firstx,lastx;
    unsigned m,mx;

    if (sign!=signx) return 1+signx;
    // Nun sind Vorzeichen gleich.
    if (mantissa==nullptr) return 0;    // denn dann ist auch x.mantissa=nullptr

    if (hexpo<hexpox) return 1 + sign;
    if (hexpox<hexpo) return 1 - sign;
    // Nun sind höchste Exponente gleich
    first=hilfmultifloat::highestbit(mantissa[l-1]);
    firstx=hilfmultifloat::highestbit(x.mantissa[lx-1]);
    last=numberofbits-first; lastx=numberofbits-firstx;
    ind=l-1;indx=lx-1;
    while (ind!=0 && indx!=0) {
        m=mx=0;
        m=(mantissa[ind]<<last) | (mantissa[ind-1]>>first);
        mx=(x.mantissa[indx]<<lastx) | (x.mantissa[indx-1]>>firstx);
        if (m<mx) return 1+sign;
        else if (mx<m) return 1-sign;
        ind--;indx--;
    }
    if (ind==indx) {
        m=mantissa[0];mx=x.mantissa[0];
    }
    else if (ind==0) {  // indx>0
        m = mantissa[ind]<<last;
        mx=(x.mantissa[indx]<<lastx) | (x.mantissa[indx-1]>>firstx);
    }
    else {
        m=(mantissa[ind]<<last) | (mantissa[ind-1]>>first);
        mx=x.mantissa[indx]<<lastx;
    }
    if (m<mx) return 1+sign;
    else if (mx==m) return 0;
    else return 1-sign;
}




int multifloat::operator <=(const multifloat &x) const
{
    if (isNaN() || x.isNaN()) return 0;
    if (x.isInfPos()) {
        if (isInfPos()) return 0;
        else return 1;
    }
    if (isInfNeg()) {
        if (x.isInfNeg()) return 0;
        else return 1;
    }
    if (isInfPos() || x.isInfNeg()) return 0;

    int sign=signum(),signx=x.signum(),hexpo=highestexp(),hexpox=x.highestexp(),l=abslength(),lx=x.abslength();
    int ind,indx,first,last,firstx,lastx;
    unsigned m,mx;

    if (sign!=signx) return 1+signx;
    // Nun sind Vorzeichen gleich.
    if (mantissa==nullptr) return 1;    // denn dann ist auch x.mantissa=nullptr

    if (hexpo<hexpox) return 1 + sign;
    if (hexpox<hexpo) return 1 - sign;
    // Nun sind höchste Exponente gleich
    first=hilfmultifloat::highestbit(mantissa[l-1]);
    firstx=hilfmultifloat::highestbit(x.mantissa[lx-1]);
    last=numberofbits-first; lastx=numberofbits-firstx;
    ind=l-1;indx=lx-1;
    while (ind!=0 && indx!=0) {
        m=mx=0;
        m=(mantissa[ind]<<last) | (mantissa[ind-1]>>first);
        mx=(x.mantissa[indx]<<lastx) | (x.mantissa[indx-1]>>firstx);
        if (m<mx) return 1+sign;
        else if (mx<m) return 1-sign;
        ind--;indx--;
    }
    if (ind==indx) {
        m=mantissa[0];mx=x.mantissa[0];
    }
    else if (ind==0) {  // indx>0
        m = mantissa[ind]<<last;
        mx=(x.mantissa[indx]<<lastx) | (x.mantissa[indx-1]>>firstx);
    }
    else {
        m=(mantissa[ind]<<last) | (mantissa[ind-1]>>first);
        mx=x.mantissa[indx]<<lastx;
    }
    if (m<mx) return 1+sign;
    else if (mx==m) return 1;
    else return 1-sign;
}

int multifloat::operator >=(const multifloat &x) const
{
    if (isNaN() || x.isNaN()) return 0;
    return !(*this<x);
}


int multifloat::operator >(const multifloat &x) const
{
    if (isNaN() || x.isNaN()) return 0;
    return !(*this<=x);
}



void multifloat::setprecision(int p)
{
    if (p<=0) prec=1;
    else prec=p;
}

int multifloat::isNaN() const
{
    if (mantissa==nullptr) return 0;
    else return (laenge==0 && mantissa[0]==0);
}


int multifloat::isInfPos() const
{
    if (mantissa==nullptr) return 0;
    else return (laenge==1 && mantissa[0]==0);
}

int multifloat::isInfNeg() const
{
    if (mantissa==nullptr) return 0;
    else return (laenge==-1 && mantissa[0]==0);
}



int multifloat::signum() const
{
    if (laenge==0) return 0;
    else if (laenge<0) return -1;
    else return 1;
}


int multifloat::highestexp() const
{
    if (mantissa==nullptr || isNaN() || isInfPos() || isInfNeg()) return 0;
    int exponent=exp,l=abslength();

    exponent+=(l-1)*numberofbits + hilfmultifloat::highestbit(mantissa[l-1]);
    return exponent;
}


multifloat multifloat::round(int precisionbits) const
{
    if (mantissa==nullptr || isNaN() || isInfPos() || isInfNeg()) return *this;
    precisionbits = hilfmultifloat::abs(precisionbits);
    int l=abslength(),k,i,wert,maxlength,highest,rest,thisbits,pos,s;
    unsigned bit(1),*m=nullptr;

    multifloat x;
    if (precisionbits==0) return x;

    if (precisionbits>numberofbits*prec) precisionbits=numberofbits*prec;

    maxlength = precisionbits/numberofbits;
    if ((precisionbits%numberofbits)) maxlength++;

    // x.mantissa ist höchstens maxlength lang.

    //if (maxlength>l) return *this;
    // Nun maxlength<=l:
    highest = hilfmultifloat::highestbit(mantissa[l-1]);
    thisbits = (l-1)*numberofbits + highest;   // Anzahl der signifikanten bits in mantissa
    /*
    if (maxlength==l) {
        if (highest <= rest) return *this;
    }
    */
    if (thisbits<=precisionbits) return *this;
    // Fall precisionbits < Anzahl der bits in mantissa :

    x.laenge=maxlength;
    x.exp=exp;
    m = new unsigned[maxlength];
    for (i=0;i<maxlength;i++) m[i]=0;
    pos=l-1;
    k=precisionbits-highest;
    pos-= (k / numberofbits);
    if ((rest=(k%numberofbits))) pos--;
    s=thisbits-precisionbits;
    x.exp+=s;
    //
    if (rest==0) {
        bit=1<<(numberofbits-1);
        wert= mantissa[pos-1] & bit;
        for (i=0;i<maxlength;i++) m[i]=mantissa[pos+i];
    }
    else {
        bit=1<<(numberofbits-rest-1);
        wert=mantissa[pos] &bit;
        for (i=0;i<maxlength-1;i++) m[i]=(mantissa[pos+i]>>(numberofbits-rest)) | (mantissa[pos+i+1]<<rest);
        m[i] = (mantissa[pos+i]>>(numberofbits-rest));
        if (pos+i<l-1) m[i] |= (mantissa[pos+i+1]<<rest);
    }
    // runden:
    for (i=0;i<maxlength;i++) {
        m[i]+=wert;
        if (m[i]) break;
    }
    x.exp+=numberofbits*i;
    if (i==maxlength) {
        x.laenge=1;
        x.mantissa= new unsigned[1];
        x.mantissa[0]=1;
        if (laenge<0) x.laenge*=-1;
        delete[] m;
        return x;
    }
    //Normalisiere:     m[i]!=0;
    bit=1;
    k=0;
    while (!(bit&m[i])) {bit<<=1;x.exp++;k++;}
    x.laenge-=i;
    s=i;
    if (k) {
        for (;i<maxlength-1;i++) {
            m[i]>>=k;
            m[i] |= m[i+1] << (numberofbits-k);
        }
        m[i]>>=k;
        if (m[i]==0) {x.laenge--;}
    }
    if (x.laenge==maxlength) {
        x.mantissa=m;
        m=nullptr;
    }
    else {
        x.mantissa=new unsigned[x.laenge];
        for (i=0;i<x.laenge;i++) x.mantissa[i]=m[i+s];
        delete[] m;
    }
    if (laenge<0) x.laenge*=-1;

    return x;

}


multifloat::operator double() const
{
    if (mantissa==nullptr) return 0.0;
    if (isInfPos()) return (1.0/0.0);
    if (isInfNeg()) return (-1.0/0.0);
    if (isNaN()) return (0.0/0.0);

    int l,exponent=highestexp()-1,h;
    if (exponent>1024) {
        if (laenge<0) return (-1.0/0.0);
        else return (1.0/0.0);
    }
    else if (exponent+1023<=0) return 0.0;

    multifloat x(round(53));

    double b=0.0;
    __int64 Exp(0),*bint= (__int64 *) &b;

    l=x.abslength();
    exponent+=1023;

    if (laenge<0) Exp+=1;
    Exp<<=11;
    Exp|=(__int64) exponent;
    Exp<<=52;

    h=hilfmultifloat::highestbit(x.mantissa[l-1]);
    *bint = (__int64)(x.mantissa[l-1]<<(numberofbits+1-h)) >> (numberofbits+1-h);
    //std::cout<<"\n *bint = "<<*bint<<"    ";
    if (l==2) {
        *bint<<=numberofbits;
        *bint|=(__int64) x.mantissa[0];
         h+=32;
    }
    *bint<<=(53-h);
    //std::cout<<"\n *bint = "<<*bint<<"    ";
    /*
    //bint= (x.mantissa[l-1]<<(numberofbits+1-h)) >> (numberofbits+1-h);
    *bint= (x.mantissa[l-1]<<(numberofbits+1-h));
    if (l==2) {
        *bint<<=numberofbits;
        *bint|=x.mantissa[0];
    }
    */
    *bint|=Exp;

    return b;
}

// z.mantissa = nullptr and z.exp inititalized, l is length of c. Get z.mantissa from c
void multifloat::normalizemantissa(unsigned *c, int l, multifloat &z)
{
    unsigned wert;
    int i, start = l - prec;
    unsigned highon = 1 << (numberofbits -1); //highest bit is 1, rest 0.
                                              //
    if (z.exp > (z.exp + start * numberofbits)) {
        z = InfPos;
        return;
    }
    z.exp += start * numberofbits;
    // round:
    wert = 0;
    if (c[start-1] & highon) wert = 1;


    for (i = start; i < l; ++i) {
        c[i] += wert;
        if ((!c[i])) {
            start++;
            wert = 1;
            if (z.exp > 0 && (z.exp > z.exp + numberofbits)) {
                z = InfPos;
                return;
            }
            z.exp += numberofbits;
        }
        else {
            wert = 0;
            break;
        }
    }
    if (wert) {
        z.mantissa = new unsigned[1];
        z.mantissa[0] = 1;
        z.laenge = 1;
        return;
    }

    // eliminate 0 at end:
    //
    unsigned last;
    int r, s;

    r = 0;
    wert = c[start];
    while (wert % 2 == 0) {
        ++r;
        ++z.exp;
        wert /= 2;
        if (z.exp == 0) {
            z = InfPos;
            return;
        }
    }

    last = c[l-1] >> r;
    if (last) z.laenge = l - start;
    else z.laenge = l - 1 - start;
    s = numberofbits - r;

    z.mantissa = new unsigned[z.laenge];

    for (i = start; i < l-1; ++i) {
        wert = c[i+1] << s;
        z.mantissa[i-start] = wert | (c[i] >> r);
    }
    if (last) z.mantissa[z.laenge - 1] = last;

    return;
}

// Additive functions:

// ahexp = a.highestexp(), bhexp = b.highestexp()   |ahexp| >= |bhexp|.
multifloat multifloat::add(const multifloat& a,const multifloat& b, int ahexp, int bhexp)  // |a| + |b|
{
    int la = a.abslength(), lb = b.abslength(), aexp = a.exp, bexp = b.exp, mexp = Min(aexp,bexp),
        digits = ahexp - mexp, maxdigits = numberofbits * prec,l;
    multifloat z;

    std::cout << "\n digits = " << digits << "; maxdigits = " << maxdigits;
    digits = Min(maxdigits,digits);
    l = digits / numberofbits;
    if (digits % numberofbits) ++l;
    mexp = ahexp - digits;
    std::cout << "\n New length: " << l << "; Minimal exp = " << mexp;
    if ((bhexp + 1) < mexp) return a;
    return z;
}


multifloat multifloat::operator +(const multifloat &y) const
{
    if (mantissa == nullptr) return y;
    if (y.mantissa == nullptr) return *this;
    if (laenge == 0 || y.laenge == 0) return NaN;

    int lx = abslength(), ly = y.abslength(), signx = signum(), signy = y.signum();

    if (lx == 1 && mantissa[0] == 0) {
        if (ly == 1 && y.mantissa[0] == 0) {
            if (signx * signy < 0) return NaN;
            else return (*this);
        }
        else return *this;
    }
    if (ly == 1 && y.mantissa[0] == 0) return y;

    multifloat z;

    if (signx == signy) {
        int ahexp = highestexp(), bhexp = y.highestexp();
        if (ahexp >= bhexp) {
            z = add(*this, y, ahexp, bhexp);
        }
        else z = add(y,*this, bhexp, ahexp);
        z.laenge *= signx;
    }

    return z;
}



// Multiplicative functions:


void multifloat::multunsigned(unsigned a, unsigned b, unsigned &high, unsigned &low)
{
    __int64 c;
    c=(__int64) (a) * (__int64)(b);
    high=unsigned(c>>numberofbits);
    low=unsigned(c);
}


int multifloat::multbyunsigned(unsigned *dat,int l,unsigned x,unsigned *c)
{
    unsigned high,low,wert=0;
    int i;

    for (i = 0; i < l; ++i, ++c) {
        multunsigned(dat[i], x, high, low);
        *c = wert + low;
        if (*c < wert) wert = 1;
        else wert = 0;
        wert += high;
    }
    if (wert) {
        *c = wert;
        l++;
    }
    else *c = 0;

    return l;
}


multifloat multifloat::operator *(const multifloat &y) const
{
    multifloat z;
    int expz = exp + y.exp, signz = signum() * y.signum();

    if (isNaN() || y.isNaN()) return NaN;
    if (iszero() && (y.isInfNeg() || y.isInfPos())) return NaN;
    if (y.iszero() && (isInfNeg() || isInfPos())) return NaN;
    if (isInfPos()) {
        if (y.signum() == 1) return InfPos;
        else return InfNeg;
    }
    if (isInfNeg()) {
        if (y.signum() == 1) return  InfNeg;
        else return InfPos;
    }
    if (iszero() || y.iszero()) return z;
    if (exp > 0 && y.exp > 0) {
        if (expz < exp || expz < y.exp) {
            if (signz > 0) return  InfPos;
            else return  InfNeg;
        }
    }
    if (exp < 0 && y.exp < 0)  {
        if (expz > exp || expz > y.exp) return z;
    }

    int i, j, lx = hilfmultifloat::abs(laenge), ly = hilfmultifloat::abs(y.laenge), la, lb, l = lx + ly;
    unsigned wert, *c, *h = nullptr, *a, *b, *ph;

    if (lx<ly) {
        la = ly; lb = lx;
        a = y.mantissa ; b=mantissa;
    }
    else {
        la = lx; lb = ly;
        a = mantissa; b = y.mantissa;
    }
    c = new unsigned int[l];
    c[l-1] = 0;

    if (lb > 1) {
        h = new unsigned int[la+1];
    }
    // 1.:
    multbyunsigned(a, la, b[0], c);

    for (i = 1; i < lb; ++i) {
        multbyunsigned(a, la, b[i], h);
        // add now:
        wert = 0;
        ph = h;
        for (j = 0; j < la; ++j, ++ph) {
            *ph += wert;
            c[j+i] += *ph;
            if (wert && !(*ph)) wert = 1;
            else if (c[j+i] < *ph) wert = 1;
            else wert = 0;
        }
        c[i+la] = wert + *ph;
    }
    if (lb>1) delete[] h;

    if (!c[l-1]) l--;

    // since both mantissas are odd, so their product.

    if (l <= prec) {
        z.mantissa = c;
        c = nullptr;
        z.laenge = signz * l;
        z.exp = expz;
        return z;
    }

    // now l > prec:

    normalizemantissa(c, l, z);

    z.laenge *= signz;
    delete [] c;

    return z;
}


std::istream& operator >>(std::istream& is ,multifloat& x)
{
    //Vorläufig:
    double b;
    is>>b;
    x = multifloat(b);
    return is;
}


std::ostream& operator <<(std::ostream& os,const multifloat& x)
{
    //Vorläufig:
    if (x.isNaN()) {
        os << "nan";
        return os;
    }
    if (x.isInfNeg()) {
        os << "-inf";
        return os;
    }
    if (x.isInfPos()) {
        os << "inf";
        return os;
    }
    if (x.iszero()) {
        os << "0";
        return os;
    }
    double b(x);
    os<<b;
    return os;
}


void multifloat::write() const
{
    int i,l=abslength();
    if (mantissa==nullptr) std::cout<<"0 0   0";
    for (i=l-1;i>=0;i--) std::cout<<"  "<<mantissa[i];
    //std::cout<<"    "<<exp;
}





}
