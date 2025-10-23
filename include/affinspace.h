#ifndef AFFINSPACE_H_INCLUDED
#define AFFINSPACE_H_INCLUDED

#include <LA.h>


// --------------------------------------------------------------------------------------------

namespace val
{

template <class T> class affinspace;
template <class T> matrix<T> getprojectionmatrix(const affinspace<T> &A, const matrix<T>& SP= matrix<T>());


template <class T>
class affinspace
{
private:
    val::vector<T> Point;
    val::vector<val::vector<T> > Vspace,OrthBasis,OrthSpace;
    val::matrix<T> LES;
    int Globaldim=0;
    int Dimension=0;
    void complete(const val::matrix<T>& SP=val::matrix<T>());

public:
    affinspace() = default;
    //affinspace(const std::string& s);
    explicit affinspace(val::matrix<T>&& A,const val::matrix<T>& SP=val::matrix<T>()) : LES(std::move(A)) {complete(SP);}   //{LES=std::move(A);complete(SP);}
    explicit affinspace(const val::matrix<T>& A,const val::matrix<T>& SP=val::matrix<T>()) : LES(A) {complete(SP);}   //{LES=std::move(A);complete(SP);}
    explicit affinspace(val::vector<val::vector<T> >&& V,const val::matrix<T> &SP=val::matrix<T>());
    explicit affinspace(const val::vector<val::vector<T> >& V,const val::matrix<T> &SP=val::matrix<T>());

    enum Intersection_Type {DISJOINT,WINDSCHIEF,PARALLEL,SUPSPACE,SUBSPACE,INTERSECTS,EQUAL};

    int globaldim() const {return Globaldim;}
    int dimension() const {return Dimension;}
    const val::vector<T>& getPoint() const {return Point;}
    const val::vector<val::vector<T> >& getVspace() const {return Vspace;}
    const val::vector<val::vector<T> >& getOrthBasis() const {return OrthBasis;}
    const val::vector<val::vector<T> >& getOrthSpace() const {return OrthSpace;}
    const val::matrix<T>& getLES() const {return LES;}
    Intersection_Type intersection(const affinspace<T>& A,val::matrix<T>& X) const;
    //P, Q sind gemeinsame Lote P gehoert zu this, Q zu B
    T squaredistance(const affinspace<T>& B,val::vector<T> &P,val::vector<T> &Q,const val::matrix<T>& SP=val::matrix<T>()) const;
    //matrix<T> getprojectionmatrix(const matrix<T>& SP= matrix<T>()) const;
};





template <class T>
affinspace<T>::affinspace(val::vector<val::vector<T> >&& V,const val::matrix<T>& SP)
{
    if (V.isempty()) {Globaldim=Dimension=0; return;}

    Globaldim=V(0).dimension();
    if (Globaldim==0) {Dimension=0; return;}

    int i,m=V.dimension()-1;
    Point = std::move(V(0));
    if (m==0) { complete();return;}

    Vspace= val::vector<val::vector<T> >(m);

    for (i=1;i<V.dimension();i++) Vspace(i-1) = std::move(V(i));

    complete(SP);
}


template <class T>
affinspace<T>::affinspace(const val::vector<val::vector<T> >& V,const val::matrix<T>& SP)
{
    if (V.isempty()) {Globaldim=Dimension=0; return;}

    Globaldim=V(0).dimension();
    if (Globaldim==0) {Dimension=0; return;}

    int i,m=V.dimension()-1;
    Point = V(0);
    if (m==0) { complete();return;}

    Vspace= val::vector<val::vector<T> >(m);

    for (i=1;i<V.dimension();i++) Vspace(i-1) = V(i);

    complete(SP);
}





template <class T>
void affinspace<T>::complete(const val::matrix<T>& SP)
{
    T zero = zero_element<T>();
    if (Point.isempty() && Vspace.isempty() && LES.isempty()) {
        Globaldim=Dimension=0;
        OrthBasis=OrthSpace=val::vector<val::vector<T> >();
        return;
    }
    else if (!LES.isempty()) {
        val::matrix<T> X;
        T det;
        int r=val::les(LES,X,det);

        if (r==0) {
            LES=val::matrix<T>();
            Point=val::vector<T>();
            Vspace=val::vector<val::vector<T> >();
            Globaldim=Dimension=0;
            OrthBasis=OrthSpace=val::vector<val::vector<T> >();
            return;
        }
        int i,j;

        Globaldim=LES.numberofcolumns()-1;
        Dimension=X.numberofrows()-1;
        Point = val::vector<T>(Globaldim);
        for (i=0;i<Globaldim;i++) Point(i) = std::move(X(0,i));

        if (Dimension>0) {
            Vspace=val::vector<val::vector<T> >(val::vector<T>(Globaldim),Dimension);
            for (i=1;i<=Dimension;i++) {
                for (j=0;j<Globaldim;j++) Vspace(i-1)(j) = std::move(X(i,j));
            }
        }
        else Vspace=val::vector<val::vector<T> >();
    }
    else {
        if (Vspace.isempty()) {
            Globaldim=Point.dimension();
            Dimension=0;
            LES=val::matrix<T>(Globaldim,Globaldim+1);
            LES.make_identity();
            for (int i=0;i<Globaldim;i++) LES(i,Globaldim) = Point(i);
            //return;
        }
        else {
            Globaldim=Vspace(0).dimension();
            if (Globaldim==0) {
                Point=val::vector<T>();
                Vspace=val::vector<val::vector<T> > ();
                Dimension=0;
                OrthBasis=OrthSpace=val::vector<val::vector<T> >();
                return;
            }
            val::spanvspace(Vspace);
            Dimension=Vspace.dimension();
            if (Point.isempty() || Point.dimension()!=Globaldim) Point = val::vector<T>(zero, Globaldim);
            if (Dimension==0) {
                LES=val::matrix<T>(Globaldim,Globaldim+1);
                LES.make_identity();
                for (int i=0;i<Globaldim;i++) LES(i,Globaldim) = Point(i);
                //return;
            }
            else {
                val::matrix<T> A(Dimension,Globaldim+1),X;
                int i,j,k,r;
                T det;
                for (i=0;i<Dimension;i++) {
                    for (j=0;j<Globaldim;j++) A(i,j) = std::move(Vspace(i)(j));
                    A(i,Globaldim) = zero;
                }

                r=val::les(A,X,det);
                r--;
                LES = val::matrix<T>(r,Globaldim+1);

                for (i=0;i<r;i++) {
                    for (j=0;j<Globaldim;j++) LES(i,j) = std::move(X(i+1,j));
                    LES(i,Globaldim) = zero;
                    for (k=0;k<Globaldim;k++) LES(i,Globaldim) += LES(i,k) * Point(k);
                }

                for (i=0;i<Dimension;i++) {
                    for (j=0;j<Globaldim;j++) Vspace(i)(j) = std::move(A(i,j));
                }
            }
        }
    }
    //Konstruiere Orthogonal-Raum und Orthogonal-Basis:

    if (Globaldim==0) {
        OrthBasis=OrthSpace=val::vector<val::vector<T> >();
        return;
    }
    if (Vspace.isempty()) {
        val::matrix<T> B(Globaldim);
        B.make_identity();
        val::movefrommatrixtovector(B,OrthSpace);
        OrthBasis=val::vector<val::vector<T> >();
        return;
    }
    else if (Dimension==Globaldim) {
        val::matrix<T> B(Globaldim);
        B.make_identity();
        val::movefrommatrixtovector(B,OrthBasis);
        OrthSpace=val::vector<val::vector<T> >();
        return;
    }
    else {
        int i;
        OrthBasis=Vspace;
        val::orthogonalize(OrthBasis,SP);
        for (i=0;i<OrthBasis.dimension();i++) val::makeprimitiv(OrthBasis(i));
        OrthSpace=val::orthogonalspace(Vspace,SP);
        for (i=0;i<OrthSpace.dimension();i++) val::makeprimitiv(OrthSpace(i));
        //WriteText(val::ToString(OrthSpace.dimension()),0);
    }

}


 template <class T>
 T affinspace<T>::squaredistance(const affinspace<T>& B,val::vector<T> &P,val::vector<T> &Q,const val::matrix<T>& SP) const
{
    T zero = zero_element<T>();
    P=Q=val::vector<T>();

    if (Globaldim==0 || B.Globaldim==0) return zero;
    if (Globaldim!=B.Globaldim) return zero;
    if (Dimension==0 && B.Dimension==0) {
        P=Point;
        Q=B.Point;
    }
    else {
        int i,j,m=Dimension+B.Dimension,r;
        T det;
        val::matrix<T> LGS(zero,m,m+1),X;
        val::vector<T> diff = Point-B.Point;

        // Stelle LGS auf:
        for (i=0;i<Dimension;i++) {
             LGS(i,i) = val::innerproduct(OrthBasis(i),OrthBasis(i),SP);
             for (j=Dimension;j<m;j++) LGS(i,j) = val::innerproduct(OrthBasis(i),B.OrthBasis(j-Dimension),SP);
             LGS(i,m) = val::innerproduct(OrthBasis(i),diff,SP);
        }
        for (;i<m;i++) {
            LGS(i,i) = val::innerproduct(B.OrthBasis(i-Dimension),B.OrthBasis(i-Dimension),SP);
            for (j=0;j<Dimension;j++) LGS(i,j) = val::innerproduct(B.OrthBasis(i-Dimension),OrthBasis(j),SP);
            LGS(i,m)=val::innerproduct(B.OrthBasis(i-Dimension),diff,SP);
        }
        r=val::les(LGS,X,det);
        if (r==0) return zero;
        P=Point; Q= B.Point;
        for (i=0;i<Dimension;i++)
            for (j=0;j<Globaldim;j++) P(j)-=X(0,i)*OrthBasis(i)(j);
        for (;i<m;i++)
            for (j=0;j<Globaldim;j++) Q(j)+=X(0,i)*B.OrthBasis(i-Dimension)(j);
    }
    val::vector<T> diff = P-Q;
    return val::innerproduct(diff,diff,SP);
}


template<class T>
typename affinspace<T>::Intersection_Type affinspace<T>::intersection(const affinspace<T>& A,val::matrix<T>& X) const
{

    X=val::matrix<T>();
    if (Globaldim==0 || A.Globaldim==0) {return DISJOINT;}

    if (Globaldim!=A.Globaldim) {
        return DISJOINT;
    }


    if (LES.isempty() || A.LES.isempty()) {return DISJOINT;}

    int i,j,m1=LES.numberofrows(),m2=A.LES.numberofrows(),m=m1+m2,n=LES.numberofcolumns(),r;
    val::matrix<T> B(m,n);
    T det;

    for (i=0;i<m1;i++)
        for (j=0;j<n;j++) B(i,j) = LES(i,j);
    for (;i<m;i++)
        for (j=0;j<n;j++) B(i,j) = A.LES(i-m1,j);

    r=val::les(B,X,det);

    if (r==0) {
        if (Dimension>0 && A.Dimension>0) {
            // Bestimmme Rang B:
            T zero = zero_element<T>();
            int k,dimhom;

            if (B.numberofcolumns()<B.numberofrows()) k=B.numberofcolumns();
            else k = B.numberofrows();
            r=0;
            for (i=0;i<k;i++) if (B(i,i)!=zero) r++;
            dimhom=Globaldim-r;
            if (dimhom>=Dimension || dimhom>=A.Dimension) return PARALLEL;
            else return WINDSCHIEF;
        }
        else return DISJOINT;
    }

    if (Dimension==r-1 && A.Dimension==r-1) {
        return EQUAL;
    }
    else if (Dimension==r-1) {
        return SUBSPACE;
    }
    else if (A.Dimension==r-1) {
        return SUPSPACE;
    }
    else {
        return INTERSECTS;
    }
}


template <class T>
matrix<T> getprojectionmatrix(const affinspace<T> &Aff,const matrix<T>& SP)
{
    if (!Aff.globaldim() || !Aff.dimension()) return matrix<T>();

    const vector<vector<T>> &OrthBasis = Aff.getOrthBasis();
    int i,j,m=OrthBasis.dimension(),Globaldim = Aff.globaldim();
    vector<T> v(Globaldim),prod(m),e(Globaldim);
    T one = val::unity_element<T>(), zero = val::zero_element<T>();
    matrix<T> A(Globaldim);

    //std::cout<<"\n Globaldim = "<<Globaldim<<" . m = "<<m;


    for (i=0;i<m;++i) prod(i) = innerproduct(OrthBasis(i),OrthBasis(i),SP);

    for (j=0;j<Globaldim;j++) {
        //std::cout<<"\n j = "<<j<<std::endl;
        v.make_zero();
        if (j>0) e(j-1) = zero;
        e(j) = one;
        for (i=0;i<m;++i) {
            v+=(innerproduct(e,OrthBasis(i),SP)/prod(i))*OrthBasis(i);
        }
        for (i=0;i<Globaldim;++i) {
            //std::cout<<"\n i = "<<i<<std::endl;
            A(i,j) = std::move(v(i));
            //std::cout<<"\n A nun : \n"<<A;
        }
    }
    return A;
}


} // end namespace val

// -------------------------------------------------------------------------------------------


#endif // AFFINSPACE_H_INCLUDED
