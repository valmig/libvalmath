#ifndef S_GROEBNER_H_INCLUDED
#define S_GROEBNER_H_INCLUDED

// Header for groebner-basis computation.

#include <string>
#include <val_basics.h>

namespace val
{

typedef void (strg_function)(const std::string &);
typedef void (void_function) (void);

DLL_PUBLIC void SetWriteText(strg_function &);
DLL_PUBLIC void SetClear(void_function &);
DLL_PUBLIC void SetMyMessage(strg_function&);
DLL_PUBLIC void WriteText(const std::string& s);

DLL_PUBLIC void SetPrimlistPath(const std::string&);

DLL_PUBLIC void SetNumberofThreads(int);
//

template <class T> class Glist;
template <class T> class s_polynom;
template <class T> class pol;
template <class T> class matrix;
class modq; class integer; class rational;

//
DLL_PUBLIC int getFileType(const std::string& name);


template <class T>
DLL_PUBLIC void saveGBasis(const std::string &arg,const val::Glist<val::s_polynom<T> >& P);


template <class T>
DLL_PUBLIC int readPolynomBasis(const std::string &arg,Glist<s_polynom<T> > &G,int onlytotdeghom=0);


DLL_PUBLIC int groebner(const std::string &name,val::Glist< val::s_polynom<val::modq> > &G);

DLL_PUBLIC int groebner(const std::string &name,val::Glist< val::s_polynom<val::integer> > &G);

DLL_PUBLIC int groebner(val::Glist<val::s_polynom<val::integer>> &G,int comment = 0);

DLL_PUBLIC int groebner(val::Glist<val::s_polynom<val::modq>> &G,int comment=0);

DLL_PUBLIC int primitiv_groebner(val::Glist<val::s_polynom<val::integer>> &G);


template <class T>
DLL_PUBLIC int isgroebner(const std::string &name,Glist<s_polynom<T> > &G,int anzeige=0);



template <class T>
DLL_PUBLIC int homgroebner(const std::string &name,Glist<s_polynom<T> > &G);


template <class T>
DLL_PUBLIC int hilbertconversion(const std::string &name,Glist<s_polynom<T>> &G,int order,const matrix<int> &M=matrix<int>());


DLL_PUBLIC pol<rational> Hilbertpolynomial(const std::string &name,int affin=1);


DLL_PUBLIC pol<rational> Hilbertpolynomial(const Glist<s_polynom<integer>>& G,int affin=1);
DLL_PUBLIC pol<rational> Hilbertpolynomial(const Glist<s_polynom<modq>>& G,int affin=1);


DLL_PUBLIC pol<rational> minimalpolynom(const Glist<s_polynom<integer> > &G,int k,int comment=1);

DLL_PUBLIC pol<modq> minimalpolynom(const Glist<s_polynom<modq> > &G,int k,int comment=1);


template <class T>
DLL_PUBLIC void groebnerwalk(val::Glist<val::s_polynom<T> > &G,int neworder,int k=2,int l=2,const val::matrix<int> &M=val::matrix<int>(),int comment=1);

template <class T>
DLL_PUBLIC void groebnerwalk(const std::string& name,val::Glist<val::s_polynom<T> > &G,int neworder,int k=2,int l=2,const val::matrix<int> &M=val::matrix<int>());


DLL_PUBLIC int isinG(const std::string& file,const std::string& gbfile);


DLL_PUBLIC int Iddim(const std::string& file);


template <class T>
DLL_PUBLIC int reduceGroebner(const std::string& name,Glist<s_polynom<T> >& G);


DLL_PUBLIC int homogenize(const std::string& i_name,const std::string& o_name,int comment=0);


DLL_PUBLIC int dehomogenize(const std::string& i_name,const std::string& o_name,int comment=0);


DLL_PUBLIC int preparenormalpos(const std::string& i_name,const std::string& o_name,int comment=0);



//

} //end namespace val


#endif // S_GROEBNER_H_INCLUDED
