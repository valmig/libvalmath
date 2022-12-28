#ifndef BB_INT_H_INCLUDED
#define BB_INT_H_INCLUDED


#include <s_groebner/common_bb.h>

namespace bb_int
{

int readfromfile(const std::string &name,val::Glist<val::s_polynom<val::integer> > &G,int &order,val::matrix<int> &Mold,int onlytotdegcompatible=0);

int homogenize(val::Glist<val::s_polynom<val::integer>> &G,int &order,val::matrix<int> &Mold);

void updateG(val::Glist< val::s_polynom<val::integer> > &G,int &nG,val::Glist< val::s_polynom<val::integer> > &Gd,
             val::Glist<val::s_polynom<val::modq> > &Gp,val::Glist<val::s_polynom<val::modq> > &Gpd,val::Glist<common_bb::spair> &lspair);

int HomGroebner(val::Glist<val::s_polynom<val::integer> > &H,val::Glist<val::s_polynom<val::integer> > &G,int nH);

int bbhommod(val::Glist<val::s_polynom<val::integer> > &H,int order,val::matrix<int> &M,int nH=-1);

int bbhommod(val::Glist<val::s_polynom<val::integer>> &H,int comment = 0);

int bbhommodmain(const std::string &argv,val::Glist<val::s_polynom<val::integer> > &G);

}  // end namespace

#endif // BB_INT_H_INCLUDED
