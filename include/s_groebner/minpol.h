#ifndef MINPOL_H_INCLUDED
#define MINPOL_H_INCLUDED


#include <s_groebner/common_bb.h>
#include <pol.h>

namespace minpol
{


template <class T>
class m_polynom : public val::s_polynom<T>
{
public:
    using val::s_polynom<T>::s_polynom;
    const m_polynom<T>& operator =(const val::s_polynom<T> &f) {val::s_polynom<T>::operator =(f);return *this;}
    const m_polynom<T>& operator =(val::s_polynom<T> &&f) {val::s_polynom<T>::operator =(std::move(f));return *this;}
    void minusmalconst(const m_polynom<T>&,const T&);
    template <class S> int linearreduction(const val::Glist<m_polynom<T> >&,const val::Glist<val::pol<S> >&,val::pol<S> &);
    template <class S> int normalform(const val::Glist<m_polynom<T> >&,S &);
    T normalize_cont();
};


template <class T>
int convertanddelete(val::Glist<val::s_polynom<T> > &G,val::Glist<m_polynom<T> > &H);

template <class T>
int convert(const val::Glist<val::s_polynom<T> > &G,val::Glist<m_polynom<T> > &H);



template <class T>
int readfromfile(char* name,val::Glist<m_polynom<T> > &G);

template<class T>
int iszerodimensional(const val::Glist< m_polynom<T> > &G,val::vector<int>&);



val::pol<val::rational> minimalpolynom(const val::Glist<m_polynom<val::integer> > &G,int k,int comment=1);
val::pol<val::modq> minimalpolynom(const val::Glist<m_polynom<val::modq> > &G,int k,int comment=1);


val::pol<val::rational> minimalpolynom(const val::Glist<val::s_polynom<val::integer> > &G,int k,int comment=1);
val::pol<val::modq> minimalpolynom(const val::Glist<val::s_polynom<val::modq> > &G,int k,int comment=1);




template <class T>
int convertanddelete(val::Glist<val::s_polynom<T> > &G,val::Glist<m_polynom<T> > &H)
{
    int nG=0;
    if (!G.isempty()) G.dellist();
    for (H.resetactual();H.actualvalid();) {
        G.inserttoend(std::move(H.actualvalue()));
        H.skiphead();
        nG++;
    }
    return nG;
}


template <class T>
int convert(const val::Glist<val::s_polynom<T> > &G,val::Glist<m_polynom<T> > &H)
{
    using namespace val;
    int i=0;
    GlistIterator<s_polynom<T> > ItG;
    m_polynom<T> h;

    if (!H.isempty()) H.dellist();
    for (ItG=G;ItG;ItG++,i++) {
         h=ItG();
         H.inserttoend(std::move(h));
    }

    return i;
}


template <class T>
int readfromfile(char* name,val::Glist<m_polynom<T> > &G)
{
    val::Glist<val::s_polynom<T> > hG;
    int n;
    m_polynom<T> h;

    n=common_bb::readfromfile(name,hG);

    for (hG.resetactual();hG.actualvalid();hG.moveactual()) {
         h=std::move(hG.actualvalue());
         G.inserttoend(std::move(h));
    }
    return n;
}


template<class T>
int iszerodimensional(const val::Glist< m_polynom<T> > &G,val::vector<int> &degX)
{
    using namespace val;
    int i,j,n=s_expo::getdim(),is=1;
    GlistIterator<m_polynom<T> > ItG;

    if (G.isempty()) return 0;

    for (i=0;i<n;i++) degX[i]=0;

    for (ItG.settohead(G);ItG.actualvalid();ItG.moveactual()) {
        for (i=0;i<n;i++) {
            if (degX[i]) continue;
            else if (ItG.getelement().LT()[i]!=0) {
                degX[i]=ItG.getelement().LT()[i];
                for (j=0;j<n;j++) {
                    if (j==i) continue;
                    if (ItG.getelement().LT()[j]!=0) {
                        degX[i]=0;
                        break;
                    }
                }
            }
        }
    }
    for (i=0;i<n;i++) {
        if (!degX[i]) {
            is=0;
            break;
        }
    }
    return is;
}



} //end namespace

#endif // MINPOL_H_INCLUDED
