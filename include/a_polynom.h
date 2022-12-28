#ifndef A_POLYNOM_H_INCLUDED
#define A_POLYNOM_H_INCLUDED

#include <n_polynom.h>

namespace val
{

template <class T> class a_polynom;
template <class T> class a_polynomIterator;
class integer;

template <class T> std::ostream& operator <<(std::ostream&,const a_polynom<T>&);


template <class T>
class a_polynom
{
private:
   //
   /*
   struct monomial {
        T coeff;
        n_expo X;
        //
        monomial() = default;
        monomial(const T &wert): coeff(wert) {}
        monomial(T&& wert) : coeff(std::move(wert)) {}
        monomial(const T& wert,const n_expo& Z) : coeff(wert), X(Z) {}
        monomial(T&& wert,const n_expo& Z) : coeff(std::move(wert)), X(Z) {}
        monomial(T&& wert,n_expo&& Z) : coeff(std::move(wert)), X(std::move(Z)) {}
        monomial& operator = (monomial &m) {coeff = m.coeff; X = m.X; return *this;}
        const T& actualcoef() const {return coeff;}
        const n_expo& actualterm() const {return X;}
    };
    */
    //
    //monomial *head = nullptr;
    T *coeff = nullptr;
    n_expo *X = nullptr;
    int len = 0;
    //
    //void reserve();
    //void resize();
    //
    static const T zero;
    static const n_expo neutral;
    //
public:
    //
    a_polynom() = default;
    a_polynom(const T& c) : coeff(new(c) T[1]), X(new n_expo[1]), len(1) {}
    a_polynom(T&& c) : coeff(new T[1]{std::move(c)}), X(new n_expo[1]), len(1) {}
    a_polynom(const T& c, const n_expo &X) : coeff(new T[1]{c}),X(new n_expo[1](X)), len(1) {}
    a_polynom(T&& c, const n_expo &X) : coeff(new T[1]{std::move(c)}), X(new n_expo[1](X)), len(1) {}
    a_polynom(T&& c,n_expo &&X) : coeff(new T[1]{std::move(c)}), X(new n_expo[1]{std::move(X)}), len(1) {}
    a_polynom(const a_polynom<T> &f);
    a_polynom(a_polynom<T> &&f);
    a_polynom(const n_polynom<T> &f);
    //
    ~a_polynom() {delete[] coeff; delete[] X;}
    void del() {delete[] coeff; delete[] X; coeff = X = nullptr; len = 0;}
    //
    a_polynom<T>& operator =(const a_polynom<T> &f);
    a_polynom<T>& operator =(a_polynom<T> &&f);
    //
    int operator ==(const a_polynom<T> &f) const;
    int operator !=(const a_polynom<T> &f) const {return !(*this == f);}
    int operator <(const a_polynom<T>&) const;
    //
    void malmonom(const T&,const n_expo&);
    void operator *=(const n_expo&);
    void operator *=(const T&);
    //
    void minusmalmonom(const a_polynom<T>& g,const T& a,const n_expo& Y,int index = 0);     // *this -= a*Y*g
    //
    int iszero() const {return (coeff == nullptr);}
    const n_expo& LT() const {return (coeff == nullptr) ? neutral : X[0];}
    const T& LC() const {return (coeff == nullptr) ? zero : coeff[0];}
    int length() const {return len;}
    //
    a_polynomIterator<T> begin() const {return a_polynomIterator<T>(*this);}
    a_polynomIterator<T> end() const {a_polynomIterator<T> it(*this); if (coeff != nullptr) it.index = len; return it;}
    //
    static int expocompare(const n_expo&,const n_expo&);
    static void changesign(T &a);    //
    friend std::ostream& operator <<<T>(std::ostream&,const a_polynom&);
    friend class a_polynomIterator<T>;
};


template <class T>
class a_polynomIterator
{
private:
    T *coeff = nullptr;
    n_expo *X = nullptr;
    int len = 0, index = 0;
public:
    a_polynomIterator() = default;
    a_polynomIterator(const a_polynom<T>& f) : coeff(f.coeff), X(f.X), len(f.len), index(0) {}
    const a_polynomIterator& operator = (const n_polynom<T>& f) {coeff = f.coeff; X = f.X; len = f.len; index = 0; return *this;}
    operator int() const {return (index < len);}
    void operator++(int) {index++;}
    void operator++() {index++;}
    const a_polynomIterator<T>& operator* () const {return *this;}
    const T& actualcoef() const {if (index >= len) return a_polynom<T>::zero;else return coeff[index];}
    const n_expo& actualterm() const {if (index >= len) return a_polynom<T>::neutral;else return X[index];}
    int actualvalid() const {return (index < len);}
    int moveactualtolast() {if (coeff != nullptr) index = len-1; return index;}    // point to last monomial
    friend class a_polynom<T>;
};



// Def.

template <class T>
const T a_polynom<T>::zero(val::zero_element<T>());

template <class T>
const n_expo a_polynom<T>::neutral;


//

template <class T>
int a_polynom<T>::expocompare (const n_expo& X,const n_expo& Y)
{
    if (X==Y) return 0;
    else if (X<Y) return -1;
    else return 1;
}


//
template <class T>
a_polynom<T>::a_polynom(const a_polynom<T> &f)
{
    if (f.coeff == nullptr) return;
    coeff = new T[len = f.len]; X = new n_expo[len];
    for (int i = 0; i < len; ++i) {
        coeff[i] = f.coeff[i]; X[i] = f.X[i];
    }
}

template <class T>
a_polynom<T>::a_polynom(a_polynom<T> &&f)
{
    if (f.coeff == nullptr) return;
    coeff = f.coeff; X = f.X; len = f.len;
    f.coeff = f.X = nullptr; f.len = 0;
}

template <class T>
a_polynom<T>::a_polynom(const n_polynom<T> &f)
{
    if (f.iszero()) return;
    int i = 0;

    len = f.length();
    coeff = new T[len]; X = new n_expo[len];

    for (const auto& m : f) {
        coeff[i] = m.actualcoef();
        X[i] = m.actualterm();
        ++i;
    }
}

//
template <class T>
a_polynom<T>& a_polynom<T>::operator =(const a_polynom<T> &f)
{
    if (coeff == f.coeff) return *this;
    del();
    if (f.coeff == nullptr) {
        return *this;
    }

    coeff = new T[len = f.len]; X = new n_expo[len];


    for(int i = 0; i < len; ++i) {
        coeff[i] = f.coeff[i];
        X[i] = f.X[i];
    }

    return *this;
}


template <class T>
a_polynom<T>& a_polynom<T>::operator =(a_polynom<T> &&f)
{
    if (f.coeff == coeff) return *this;
    del();
    len = f.len; coeff = f.coeff; X = f.X;
    f.len = 0; f.coeff = f.X = nullptr;
    return *this;
}

//

template <class T>
int a_polynom<T>::operator ==(const a_polynom<T> &f) const
{
	if (coeff == f.coeff) return 1;
	if (len != f.len) return 0;
	for (int i = 0; i < len; ++i) {
		if (coeff[i] != f.coeff[i]) return 0;
		if (X[i] != f.X[i]) return 0;
	}
	return 1;
}

template <class T>
int a_polynom<T>::operator <(const a_polynom<T> &f) const
{
    if (f.coeff == nullptr) return 0;

    for (int i = 0; i < f.len; ++i) {
       if (i == len) return 1;
       if (X[i] < f.X[i]) return 1;
       else if (X[i] == f.X[i]) continue;
       else return 0;
    }
    return 0;
}

// -----------------------------------------------------------------------------------------------------------

template <class T>
void a_polynom<T>::malmonom(const T& a, const n_expo& Y)
{
    if (coeff == nullptr) return;
    if (a == zero) {del(); return;}

    for (int i = 0; i < len; ++i) {
        coeff[i] *= a;
        X[i] *= Y;
    }
}


template <class T>
void a_polynom<T>::operator *=(const n_expo& Y)
{
    if (coeff == nullptr) return;
    for (int i = 0; i < len; ++i) {
        X[i] *= Y;
    }
}

template <class T>
void a_polynom<T>::operator *=(const T &a)
{
    if (coeff == nullptr) return;
    if (a == zero) {
        delete[] coeff; delete[] X; len = 0;
        return;
    }
    for (int i = 0; i < len ; ++i) coeff[i] *= a;
}

// -----------------------------------------------------------------------------------------------------------

template <class T>
void a_polynom<T>::minusmalmonom(const a_polynom<T> &g,const T &a,const n_expo &Y,int index)
{
    int c, i=index, i_g = 0, i_n = 0, n_len = len + g.len;
    T b, *n_coeff = nullptr;
    n_expo Z, *n_X = nullptr;

    if (g.coeff == nullptr || a==zero) return;

    if (i == 0) {
        // 1. Set head:
        for (i_g = 0; i_g < g.len; ++i_g, --n_len) {
            if (i == len) break;
            Z = g.X[i_g] * Y;
            c = expocompare(Z,X[i]);
            if (c==-1) {                                // case: Z < head->X
                n_coeff = new T[n_len]; n_X = new n_expo[n_len];
                n_coeff[0] = std::move(coeff[i]); n_X[0] = std::move(X[i]);
                ++i; i_n = 1;
                break;
            }
            else if (c==1) {                            // case: head->X < Z
                b = g.coeff[i_g] * a;
                changesign(b);
                n_coeff = new T[n_len];
                n_X = new n_expo[n_len];
                n_coeff[0] = std::move(b); n_X[0] = std::move(Z);
                ++i_g; ++i_n;
                break;
            }
            else {
                b = g.coeff[i_g] * a;
                changesign(b);
                coeff[i] += b; //h
                if (coeff[i] == zero) {//dat==NULL) {
                    --n_len; ++i;
                }
                else {
                    ++i_g;
                    n_coeff = new T[n_len]; n_X = new n_expo[n_len];
                    n_coeff[0] = std::move(coeff[i]); n_X[0] = std::move(X[i]);
                    ++i; ++i_n;
                    break;
                }
            }
        }
    }
    else {
        n_coeff = new T[n_len]; n_X = new n_expo[n_len];
        for (i = 0; i < index; ++i) {
            n_coeff[i] = std::move(coeff[i]); n_X[i] = std::move(X[i]);
        }
        i = i_n = index;
    }


    if (n_coeff == nullptr) {
        if (i_g < g.len || i < len) {
            n_coeff = new T[n_len]; n_X = new n_expo[n_len];
        }
        for (; i < len; ++i, ++i_n ) {
            n_coeff[i_n] = std::move(coeff[i]); n_X[i_n] = std::move(X[i]);
        }
        for (; i_g < g.len; ++i_g, ++i_n) {
            b = g.coeff[i_g] * a;
            changesign(b);
            Z  = g.X[i_g] * Y;
            n_coeff[i_n] = std::move(b); n_X[i_n] = std::move(Z);
        }
    }

    // Invariant: X[i-1] > g.X[i_g] * Y
    while (i_g < g.len && i < len) {
        Z = g.X[i_g] * Y;
        c = expocompare(Z,X[i]);
        if (c == -1) {                                                      // case: Z < X[i]
            n_coeff[i_n] = std::move(coeff[i]); n_X[i_n] = std::move(X[i]);
            ++i;
        }
        else if (c == 1) {                                                  // case: X[i] < Z
            b = g.coeff[i_g] * a;
            changesign(b);
            n_coeff[i_n] = std::move(b); n_X[i_n] = std::move(Z);
            ++i_g;
        }
        else {
            b = g.coeff[i_g] * a;
            changesign(b);
            coeff[i] += b;
            if (coeff[i] != zero) {
                n_coeff[i_n] = std::move(coeff[i]); n_X[i_n] = std::move(X[i]);
            }
            else {--n_len; --i_n;}
            ++i; ++i_g;
        }
        ++i_n;
    }
    // Append eventually rest.
    for (; i< len; ++i, ++i_n) {
        n_coeff[i_n] = std::move(coeff[i]); n_X[i_n] = std::move(X[i]);
    }
    for (; i_g < g.len; ++i_g, ++i_n) {
        b = g.coeff[i_g] * a;
        changesign(b);
        Z = g.X[i_g] * Y;
        n_coeff[i_n] = std::move(b); n_X[i_n] = std::move(Z);
    }

    delete[] coeff; delete[] X;
    coeff = n_coeff; X = n_X; len = n_len;

    return;
}


// ------------------------------------------------------------------------------------------------------------

template <class T>
std::ostream& operator <<(std::ostream& os,const a_polynom<T>& f)
{
    for (int i = 0; i < f.len; ++i) {
        os<<f.coeff[i]<<std::endl;
        os<<f.X[i];
        os<<std::endl;
    }
    os<<0<<std::endl;
    return os;
}


}// end namespace val

#endif // A_POLYNOM_H_INCLUDED
