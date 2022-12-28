#ifndef BB_MODQ_H_INCLUDED
#define BB_MODQ_H_INCLUDED

#include <s_groebner/common_bb.h>

namespace bb_modq
{

extern int nGrest,nGrestold,nGall;


void update(val::s_polynom<val::modq> &f,val::Glist< val::s_polynom<val::modq> > &G,
            val::Glist< val::s_polynom<val::modq> > &Grest,val::Glist< common_bb::spair > &lspair,int &m);

int Groebner(val::Glist< val::s_polynom<val::modq> > &G,int comment=1);

int bbgccmain(const std::string &argv,val::Glist< val::s_polynom<val::modq> > &G);

} // end namespace val


#endif // BB_MODQ_H_INCLUDED
