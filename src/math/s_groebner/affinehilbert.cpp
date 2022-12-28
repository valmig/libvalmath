
#include <s_groebner/affinehilbert.h>

namespace a_hilbert
{

int readfromfile(char* name,val::Glist<val::s_polynom<val::modq> > &G,int &order,val::matrix<int> &Mold,int onlytotdegcompatible)
{
	int i=0,n,nG=0;
	val::s_polynom<val::modq> f;
	std::string line;
	val::matrix<int> Mnew;

	std::ifstream file(name,std::ios::in); 
	if (!file) {                              
		common_bb::MyMessage("FILE DOES NOT EXIST");
        return 0;
	}
    do {
        getline(file,line);
        if (line!="") i++;
        else break;
	}
	while(file);

	if (i!=3) {
        common_bb::MyMessage("\nWrong type of file!");
        return 0;
	}

	file.clear();
	file.seekg(0,std::ios::beg);

	file>>val::modq::q>>n>>order;
	if (order==-1000) {
		int j;
		Mold=val::matrix<int>(n);
		Mnew = val::matrix<int>(n+1);
		val::s_expo::setordtype(-1000);
		for (i=0;i<=n;i++) Mnew(0,i)=1;
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) {
				file>>Mold(i,j);
				Mnew(i+1,j)=Mold(i,j);
			}
			Mnew(i+1,n)=0;
		}
		if (onlytotdegcompatible) {
            for (i=0;i<n;i++)
                if (Mold(0,1)!=1) {
                   common_bb::MyMessage("\nOrder is not totdeg-compatible!");
                   return 0;
                }
		}
		val::s_expo::setordmatrix(Mnew);
	}

	else if (order==-2) {
		val::s_expo::setordtype(-2);
	}
	else if (order==-1) {
        if (onlytotdegcompatible) {common_bb::MyMessage("\nOrder is not totdeg-compatible!"); return 0;}
		val::s_expo::setordtype(0);
	}
	else if (order==0) {
        Mnew = val::matrix<int>(0,n+1,n+1);
        for (i=0;i<n;i++) Mnew(0,i) = Mnew(1,i) = 1;
        Mnew(0,n) =1;
        for (i=2;i<=n;i++) Mnew(i,i-2) =1;
        val::s_expo::setordmatrix(Mnew);
	}
	else val::s_expo::setordtype(-2);

	val::s_expo::setdim(n+1);
	val::s_expo X(0);
	val::modq coeff,zero(0);

	do {
		file>>coeff;
		if (coeff!=zero) {
			for (i=0;i<n;i++) {
				file>>X[i];
			}
			f.insert(std::move(coeff),X);
		}
		else break;
		do {
			file>>coeff;
			if (coeff!=zero) {
				for (i=0;i<n;i++) {
					file>>X[i];
				}
				f.insert(std::move(coeff),X);
			}
			else break;
		}
		while (1);

		if (!f.iszero()){
			f.homogenize();
			f.normalize();
			f.reord();
			G.sinsert(std::move(f));
			nG++;
		}
		else break;
	}
	while (file);
	file.close();
	common_bb::WriteText("\nRed list G of file.\nG has " + val::ToString(nG) + " elements.");
	return nG;
}



template <>
int affinehilbertconversionmain(char* name,val::Glist<val::s_polynom<val::modq> >& G,int order,const val::matrix<int> &M)
{
    using namespace val;
    int nG,s_order;
    matrix<int> Ms;

    if (!G.isempty()) G.dellist();

    readfromfile(name,G,s_order,Ms,1);
    nG=affinehilbertconversion(G,order,M);
    return nG;
}


template <>
int affinehilbertconversionmain(char* name,val::Glist<val::s_polynom<val::integer> >& G,int order,const val::matrix<int> &M)
{
    using namespace val;
    int nG,s_order;
    matrix<int> Ms;

    if (!G.isempty()) G.dellist();

    bb_int::readfromfile(name,G,s_order,Ms,1);
    nG=affinehilbertconversion(G,order,M);
    common_bb::WriteText("\nMaximal integer-length: " + val::ToString(val::integer::GetMaxlength()));
    common_bb::WriteText("\nMaximal integer-length in G: " + val::ToString(common_bb::Maximalinteger(G)));
    return nG;
}


} //end namespace
