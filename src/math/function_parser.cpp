#include <function_parser.h>
#include <val_basics.h>
#include <val_utils.h>
#include <analysis.h>
#include <pol.h>
#include <rational.h>
#include <fraction.h>
#include <error.h>
#include <n_polynom.h>
#include <s_polynom.h>
#include <pol_arithmetic.h>
//#include <s_modinteger.cpp>
//#include <polfactor.h>
#include <s_groebner.h>


// Computing lcm of multivariate polynomials via groebner-bases.
val::n_polynom<val::rational> lcm(const val::n_polynom<val::rational>& f,const val::n_polynom<val::rational>& g)
{
    using namespace val;
    if (f.iszero()) return f;
    if (g.iszero()) return g;

    n_polynom<rational> h;
    n_polynom<integer> f1 = val::primitivpart(f), g1=val::primitivpart(g);
    s_polynom<integer> f2,g2;
    Glist<s_polynom<integer>> G;

    int i,n = val::Max(f1.getdim(),g1.getdim());

    s_expo::setordtype(-1);
    s_expo::setdim(n+1);
    s_expo X,Y(0);
    // f2 = X0 * f1;
    for (const auto &monom : f1) {
        for (i=0;i<n;++i) X[i+1] = monom.actualterm()[i];
        X[0] = 1;
        f2.insert(monom.actualcoef(),X);
    }
    // g2 = (X0-1)*g1
    for (const auto& monom : g1) {
        for (i=0;i<n;++i) {
            X[i+1] = Y[i+1] = monom.actualterm()[i];
        }
        X[0] = 1;
        g2.insert(monom.actualcoef(),X);
        g2.insert(-monom.actualcoef(),Y);
    }
    G.sinsert(std::move(f2));
    G.sinsert(std::move(g2));
    if (!primitiv_groebner(G)) {
        std::cout<<"\nGroebner-Basis computation failed!! G.length() = "<<G.length();
        exit(-1);
    }
    n_expo Z(n);
    for (auto monom=G[0].begin();monom;++monom) {
        for (i=0;i<n;++i) Z[i] = monom.actualterm()[i+1];
        h.insert(rational(monom.actualcoef()),Z);
    }
    h.normalize();
    return h;
}

val::n_polynom<val::rational> gcd(const val::n_polynom<val::rational> &f,const val::n_polynom<val::rational> &g)
{
    using namespace val;

    if (f.iszero()) return g;
    if (g.iszero()) return f;

    n_polynom<rational> h,r;

    divrem(f*g,lcm(f,g),h,r);
    if (!r.iszero()) {
        std::cout<<"\ngcd computation failed!, r = \n"<<r;
        exit(-1);
    }
    h.normalize();
    return h;
}


namespace fparser
{
val::d_array<std::string> sfunctionlist({"sqrt", "exp", "log", "abs", "sinh", "cosh", "tanh", "arsinh", "arcosh", "artanh", 
										"sin", "cos", "tan", "arcsin", "arccos", "arctan"});

val::d_array<val::GPair<std::string,val::d_function*>> 
functionpairs ({ {"sqrt",val::sqrt}, {"abs",val::abs}, {"exp", val::exp}, {"log", val::log}, {"sin", val::sin}, {"cos", val::cos}, 
					{"tan", val::tan}, {"arcsin", val::arcsin}, {"arccos", val::arccos}, {"arctan", val::arctan}, {"sinh", val::sinh},
					{"cosh", val::cosh}, {"tanh", val::tanh}, {"arsinh", val::arsinh}, {"arcosh", val::arcosh}, {"artanh", val::artanh}	});


std::string getstringfunction(const std::string &s, int j)
{
	int n = s.length() - j, m, found, i;
	
	for (const auto& sf : sfunctionlist) {
		m = sf.length();
		found = 0;
		if (n >= m) {
			found = 1;
			for (i = 0; i < m; ++i) {
				if (sf[i] != s[i+j]) {
					found = 0;
					break;
				}
			}
		}
		if (found) return sf;
	}

	return "";
}

int getindexoffunction(const std::string &sf)
{
	if (sf =="") return -1;
	int i = 0;
	
	for (const auto &f : functionpairs) {
		if (f.x == sf) return i;
		++i;
	} 
	
	return -1;
}



int isinteger(const double &a,const double &epsilon=1e-9)
{
    int h=int (val::integer(val::rational(a)));
    double d_h=double(h);
    if (val::abs(d_h-a)<epsilon) return 1;
    else return 0;
}


int has_point(const std::string &s)
{
    int n=s.length();

    for (int i=0;i<n;++i)
        if (s[i]==',' || s[i]=='.') return 1;
    return 0;
}


int checkbrackets(const std::string& s)
{
    int i=0,n=s.size(), k=0;

    for (i=0;i<n;++i) {
        if (s[i] == '(') k++;
        else if (s[i]==')') k--;
        else continue;
        if (k<0) return 0;
    }

    if (k!=0) return 0;

    return 1;
}


std::string findnumber(const std::string &s, int &i)
{
    int n=s.size(),e_set=0,p_set=0;
    std::string out;

    for (;i<n;i++) {
        if (s[i]>='0'  && s[i]<='9') out+=s[i];
        else if (s[i]=='e' && i<n-1 && s[i+1]!='x') {
            if (e_set) {i++; break;}
            if (i==n-1) {i++;break;}
            if (s[i+1] == '-' || s[i+1]=='+') {
                out+=s[i];i++;out+=s[i];e_set=1;
            }
            else {out+=s[i];e_set=1;}
        }
        else if (s[i]==',' || s[i]=='.') {
            if (p_set) {i++;break;}
            out+=s[i];
            p_set=1;
        }
        else {break;}
    }

    return out;
}

int foundpattern(const std::string &s,const std::string &pattern, int& i)
{
    int is=1,n=s.size(),m=pattern.size(),j;

    if (i+m>n) return 0;

    for(j=0;i<n && j<m;++i,++j) {
        if (s[i]!=pattern[j]) {return 0;}
    }

    return is;
}


int is_sum(const std::string& s)
{
    int n=s.length();
    for (int i=0;i<n;++i) {
        if (s[i]=='+' || s[i]=='-') return 1;
    }
    return 0;
}

int is_arithmeticterm(const std::string &s)
{
    int n=s.length();
    for (int i=0;i<n;++i) {
        if (s[i]=='+' || s[i]=='-' || s[i]=='*' || s[i]=='/' || s[i]=='^') return 1;
    }
    return 0;
}

/*
int is_functionterm(const std::string &s)
{
    int j=0;
    if (foundpattern(s,"sin",j)) return 1;
    j=0;
    if (foundpattern(s,"cos",j)) return 1;
    j=0;
    if (foundpattern(s,"tan",j)) return 1;
    j=0;
    if (foundpattern(s,"arcsin",j)) return 1;
    j=0;
    if (foundpattern(s,"arccos",j)) return 1;
    j=0;
    if (foundpattern(s,"arctan",j)) return 1;
    j=0;
    if (foundpattern(s,"exp",j)) return 1;
    j=0;
    if (foundpattern(s,"log",j)) return 1;
    j=0;
    if (foundpattern(s,"sqrt",j)) return 1;
    j=0;
    if (foundpattern(s,"abs",j)) return 1;
    j=0;
    if (foundpattern(s,"sinh",j)) return 1;
    j=0;
    if (foundpattern(s,"cosh",j)) return 1;
    j=0;
    if (foundpattern(s,"tanh",j)) return 1;
    j=0;
    if (foundpattern(s,"arsinh",j)) return 1;
    j=0;
    if (foundpattern(s,"arcosh",j)) return 1;
    j=0;
    if (foundpattern(s,"artanh",j)) return 1;
    j=0;
    return 0;
}
*/

int is_sum_operator(const std::string& s)
{
    if (s=="") return 0;
    if (s[0]=='+' || s[0]=='-') return 1;
    return 0;
}

int is_arithmetic_operator(const std::string& s)
{
    if (s=="") return 0;
    if (s[0]=='+' || s[0]=='-' || s[0]=='*' || s[0]=='/' || s[0]=='^') return 1;
    return 0;
}


val::rationalfunction derive(const val::rationalfunction &f)
{
    return val::rationalfunction(f.nominator().derive()*f.denominator()-f.nominator()*f.denominator().derive(),f.denominator()*f.denominator());
}

int issimilar(const double &x,const double &y,const double &eps=1-9)
{
    if (val::abs(x-y)<=eps) return 1;
    else return 0;
}

val::Glist<double> unify(const val::Glist<double> &H1, const val::Glist<double> &H2,const double& eps)
{
    if (H1.isempty()) return H2;
    if (H2.isempty()) return H1;
    val::Glist<double> G;
    int i1,i2,n1=H1.length(),n2=H2.length();

    for (i1=0,i2=0; i1<n1 && i2 < n2;) {
        if (H1[i1]<H2[i2]) {
            G.push_back(H1[i1]);
            if (issimilar(H1[i1],H2[i2],eps)) {
                i2++;
            }
            i1++;
        }
        else {
            G.push_back(H2[i2]);
            if (issimilar(H1[i1],H2[i2],eps)) i1++;
            i2++;
        }
    }
    for (;i1<n1;i1++) G.push_back(H1[i1]);
    for (;i2<n2;i2++) G.push_back(H2[i2]);
    return G;
}

val::Glist<val::GPair<double>> unify(const val::Glist<val::GPair<double>> &X,const val::Glist<val::GPair<double>> &Y,const double &epsilon)
{
    val::Glist<val::GPair<double>> Z;
    if (X.isempty()) return Y;
    if (Y.isempty()) return X;
    int nx=X.length(),ny=Y.length(),ready=0,ix=0,iy=0,i,n;
    val::d_array<int> xready(0,nx),yready(0,ny);
    double x1,y1;
    val::GPair<double> const *Cand=nullptr;
    val::GPair<double> P;

    while (ready<nx+ny) {
        x1=y1=val::Inf;
        Cand=nullptr;
        for (i=0;i<nx;++i) {
            if (xready[i]) continue;
            if (X[i].x<x1) {
                ix=i;
                x1 = X[i].x;
            }
        }
        for (i=0;i<ny;++i) {
            if (yready[i]) continue;
            if (Y[i].x<y1) {
                iy=i;
                y1 = Y[i].x;
            }
        }
        if (x1==val::Inf) {
            Cand = &Y[iy];
            yready[iy]=1;
        }
        else if (y1==val::Inf) {
            Cand = &X[ix];
            xready[ix] = 1;
        }
        else if (issimilar(x1,y1,epsilon)) {
            if (X[ix].y>Y[iy].y) Cand = &X[ix];
            else Cand = &Y[iy];
            xready[ix] = yready[iy] = 1;
            ready++;
        }
        else if (x1<y1) {
            if (issimilar(X[ix].y,y1,epsilon) || y1<X[ix].y){
                P.x=x1; P.y = val::Max(X[ix].y,Y[iy].y);
                Cand = &P;
                xready[ix] = yready[iy] = 1;
                ready++;
            }
            else {
                Cand = &X[ix];
                xready[ix] = 1;
            }
        }
        else if (x1>=y1) {
            if (issimilar(Y[iy].y,x1,epsilon) || x1 < Y[iy].y) {
                P.x=y1; P.y = val::Max(X[ix].y,Y[iy].y);
                Cand = &P;
                xready[ix] = yready[iy] = 1;
                ready++;
            }
            else {
                Cand = &Y[iy];
                yready[iy]=1;
            }
        }
        if (Cand==nullptr) {
            val::Error::error("Nullpointer in union of Glist<GPair<double>>!\n");
        }
        n = Z.length();
        if (Z.isempty()) Z.push_back(*Cand);
        else if (Z[n-1].y>=Cand->x) {
            Z[n-1].y = val::Max(Z[n-1].y,Cand->y);
        }
        else if (issimilar(Z[n-1].y,Cand->x,epsilon)) {
            Z[n-1].y = Cand->y;
        }
        else Z.push_back(*Cand);
        ready++;
    }
    return Z;
}

/*
std::string factorize(const val::pol<val::rational> &f)
{
    if (f.degree()==0 || f.iszero()) return val::ToString(f.LC());
    if (f.degree()==1) return val::PolToString(f);
    std::string s;
    val::d_array<val::pol<val::rational>> factors = val::polfactor(f);
    val::pol<val::rational> h;
    val::rational cont = val::content(f);
    int e,brackets,r,signf,sign;

    signf = f.LC().signum();
    sign = cont.signum();

    for (const auto& g : factors) sign*=g.LC().signum();

    h = f;
    if (sign != signf) cont.changesign();
    if (cont == val::rational(-1)) s += "-";
    else if (cont != val::rational(1)) s += val::ToString(cont);
    if (factors.length()>0 && s!="" && s!="-") s+="*";
    r=factors.length();
    for (int i=0;i<r;++i) {
        e=0;

        while ((h%factors[i]).iszero())
        {
            h/=factors[i];
            ++e;
        }
        brackets=0;
        if ((r > 1 || e > 1) && factors[i].length()>1) {
            s+="(";
            brackets=1;
        }
        s+=val::PolToString(factors[i]);
        if (brackets) s+=")";
        if (e>1) {
            s+="^"+val::ToString(e);
        }
        if (i<r-1) s+=" ";
    }
    return s;
}

std::string factorize(const val::rationalfunction &F)
{
	int bracket=0;
	std::string sf;

	if (F.nominator().length()>1) {sf+='(';bracket=1;}
	sf+=factorize(F.nominator());


	if (bracket) sf+=')';
	bracket=0;
	if (F.denominator()!=val::pol<val::rational>(1,0)) {
		sf+='/';
		if (F.denominator().length()>1) {
			bracket=1;
			sf+='(';
		}
		sf+=factorize(F.denominator());
		if (bracket) sf+=')';
	}

	return sf;
}


std::string PolfractionToString(const val::rationalfunction &F)
{
	int bracket=0;
	std::string sf;
	val::pol<val::rational> zero,one(1,0),minusone(-1,0);

	if (F.nominator()==zero) return "0";
	if (F.denominator()==one) return factorize(F.nominator());
	if (F.denominator()==minusone) return factorize(-F.nominator());

	val::pol<val::integer> nom,denom;
	val::rational cont,c1;

	val::primitivpart(F.nominator(),nom,cont);
	val::primitivpart(F.denominator(),denom,c1);
	cont/=c1;
	nom*=val::nominator(cont);
	denom*=val::denominator(cont);

	if (F.nominator().length()>1) {sf+='(';bracket=1;}
	sf+=factorize(val::toRationalPolynom(nom));
	if (bracket) sf+=')';
	bracket=0;
	if (F.denominator()!=one) {
		sf+='/';
		if (F.denominator().length()>1) {
			bracket=1;
			sf+='(';
		}
		sf+=factorize(val::toRationalPolynom(denom));
		if (bracket) sf+=')';
	}

	return sf;
}

*/

} //end namespace fparser


namespace val
{

const std::string valfunction::zero_string="0";

const std::string& valfunction::getinfixnotation() const
{
    if (s_infix=="") return zero_string;
    return s_infix;
}

std::string valfunction::getfirstoperator() const
{
    if (Gdat.isempty()) return "";
    int n=Gdat.length();
    return Gdat[n-1].data;
}


void valfunction::squeeze(d_array<token> &f,const d_array<token> &g,int a,int b)
{
    int n=f.length(),m=g.length(),k;
    if (m == 0 || a<0 || a>=b || b>n) return;

    if (m>b-a) {
        f.resize(n+m-b+a);
        for (k=n-1;k>=b;--k) f[k+m-b+a] = std::move(f[k]);
        for (k=0;k<m;++k) f[k+a]=g[k];
    }
    else {
        for (k=0;k<m;++k,++a) f[a] = g[k];
        m=b-a;
        if (m) for (;a<n-m;++a) f[a] = std::move(f[a+m]);
        f.resize(n-m);
    }
}

void valfunction::to_double(d_array<token> &f)
{
    int n=f.length(),j;
    d_array<token> h(1);
    rational a,b;
    for (j=0;j<n;++j) {
        if (f[j].data=="/" && f[j+1].type==0 && f[j+2].type==0) {
            a = FromString<rational>(f[j+2].data);
            b = FromString<rational>(f[j+1].data);
            h[0] = token(ToString(double(a/b)),0);
            squeeze(f,h,j,j+3);
            n-=2;
        }
    }
}

int valfunction::has_variable(const d_array<token> &f)
{
    for (const auto& value : f) {
        if (value.data[0]=='x') return 1;
    }
    return 0;
}

int valfunction::has_operator(const d_array<token> &f,const std::string &op)
{
    for (const auto& value : f) {
        if (value.data==op) return 1;
    }
    return 0;
}

void valfunction::subst_var_t_pi(d_array<token> &f_t,d_array<d_array<token>> &toklist,int &nx,int nvar)
{
    int found,par=0,i,realnvar=0;
    std::string svar;
    d_array<int> invar(0,nvar);

    for (const auto& value : f_t) {
		if (value.type==0 && value.data=="PI") {
			toklist.push_back(d_array<token>{token("PI",0)});
			par++;
			break;
		}
	}

    for (const auto& value : f_t) {
		if (value.type==0 && value.data=="t") {
			toklist.push_back(d_array<token>{token("t",0)});
			par++;
			break;
		}
	}
    for (i=1;i<=nvar;++i) {
		svar="x" + val::ToString(i);
		for (const auto& value : f_t) {
			if (value.data==svar || (i==1 && value.data=="x")) {
				toklist.push_back(d_array<token>{token(svar,1)});
				invar[i-1]=1;
                realnvar++;
				break;
			}
		}
	}
	//std::cout<<"\n realnvar = "<<realnvar;
	//std::cout<<"\n invar = ";
	for (i=nvar-1;i>=0;--i) {
        if (invar[i]) {
            invar[i]=realnvar;
            realnvar--;
        }
        //std::cout<<invar[i]<<"  ";
	}
	
	d_array<int> replaced(0,f_t.length());
    for (i=nvar;i>=1;--i) {
		if (!invar[i-1]) continue;
		svar="x"+ val::ToString(i);
		for (int j = 0; j < f_t.length(); ++j) {
			if (replaced[j]) continue;
			if (f_t[j].data==svar || (i==1 && f_t[j].data=="x")) {
				//std::cout<<"\n original: "<<value.data;
				f_t[j].data= "x" + val::ToString(invar[i-1]+par);
				replaced[j] = 1;
				//std::cout<<", replaced with: "<<value.data;
			}
		}
	}
	
	//std::cout<<"\n After renaming vars, f_t = ";
	//for (const auto& v : f_t) std::cout<<v.data<<"  ";

    // Find and replace PI:
    nx=toklist.length();
    found=0;
    par=0;
    for (auto& value : f_t) {
        if (value.type==0 && value.data=="PI") {
            value.data="x" + val::ToString(par+1) ;value.type=1;
            found=1;
        }
    }
    if (found) par++;

    // Find and replace parameter t:
    found=0;
    for (auto& value : f_t) {
        if (value.type==0 && value.data=="t") {
            value.data="x" + val::ToString(par+1);value.type=1;
            found=1;
        }
    }
    if (found) par++;
}

void valfunction::back_subst(d_array<token> & f_t,const d_array<d_array<token>> &toklist,int nx)
{
    int n,i,k,d,l;
    std::string ns;

    d=toklist.length();
    n=f_t.length();
    

    /*
    for (k = 0; k < n; ) {
		if (f_t[k].data[0] == 'x') {
			l = f_t[k].data.length();
			i = FromString<int>(tailofstring(f_t[k].data,l-1)) -1;
			if (i >= nx && i < d) {
				squeeze(f_t,toklist[i],k,k+1);
				k+=toklist[i].length();
				n = f_t.length();
				continue;
			}
		}
		++k;
	}
    */
    
    
    for (i=d-1;i>=nx;--i) {
        ns="x" + ToString(i+1);
        l=toklist[i].length();
        for (k=0;k<n;) {
            if (f_t[k].data==ns) {
                squeeze(f_t,toklist[i],k,k+1);
                k+=l;
                n=f_t.length();
            }
            else ++k;
        }
        //std::cout<<"\n Substitution i = "<<i+1<<"  f_t = ";
        //for (auto& value : f_t) std::cout<<value.data<<" ";
    }
    
  
    for (k = 0; k < n; ) {
		if (f_t[k].data[0] == 'x') {
			l = f_t[k].data.length();
			i = FromString<int>(tailofstring(f_t[k].data,l-1)) -1;
			if (i >= 0 && i < nx) {
				squeeze(f_t,toklist[i],k,k+1);
				k+=toklist[i].length();
				n = f_t.length();
				continue;
			}
		}
		++k;
	}
    
    /*
    for (i=0;i<nx;++i) {
        ns="x" + ToString(i+1);
        l=toklist[i].length();
        for (k=0;k<n;) {
            if (f_t[k].data==ns) {
                //std::cout<<"\n data = "<<ns;
                squeeze(f_t,toklist[i],k,k+1);
                //std::cout<<", replaced with: "<<f_t[k].data;
                k+=l;
                n=f_t.length();
            }
            else ++k;
        }
        //std::cout<<"\n Substitution i = "<<i+1<<"  f_t = ";
        //for (auto& value : f_t) std::cout<<value.data<<" ";
    }
    */
    
}



void valfunction::simplify_exp(d_array<token> &f_t,int nvar,int prod)
{
    int i,k,m,m1,m2,l,n=f_t.length(),found,nx=0;
    d_array<token> h_t,tok,h1,tok1,tok2;
    d_array<d_array<token>> toklist;

    // log(exp):
    for (i=0;i<n-1;++i) {
		if (f_t[i].data=="log" && f_t[i+1].data=="exp") {
			k=i+2;
			tok=splitfunction(f_t,k);
			squeeze(f_t,tok,i,k);
			n=f_t.length();
			--i;
		}
	}

    for (i=0;i<n-1;++i) {
        if (f_t[i].data=="^") {  // powers
            k=i+1;
            h_t=splitfunction(f_t,k);
            tok=splitfunction(f_t,k);
            if (tok[0].data!="exp") continue;
            m1=h_t.length(); m2=tok.length(); m=m1+m2;
            tok[0].data="*";
            tok.resize(m);
            for (l=0;l<m1;++l) tok[m2+l] = h_t[l];
            f_t[i].data="exp";
            squeeze(f_t,tok,i+1,k);
            i=k-1;
            n=f_t.length();
        }
        else if (f_t[i].data == "exp" && f_t[i+1].data == "0") {  // epx(0)
            k = i + 2;
            tok.reserve(1); tok.push_back(token("1",0));
            squeeze(f_t,tok,i,k);
            i=k-1;
            n=f_t.length();
        }
    }


    // Divisions:
    for (i=0;i<n-2;++i) {
        if (f_t[i].data=="/" && f_t[i+1].data=="exp") {
            h1.del();
            k=i+2;
            tok=splitfunction(f_t,k);
            m=tok.length()+1;
            h1.reserve(m);
            for (l=1;l<m;++l) h1[l] = tok[l-1];
            h1[0] = token("m",2);
            f_t[i] = token("*",2);
            squeeze(f_t,h1,i+2,k);
            n=f_t.length();
        }
    }
    //

    if (!prod) return;

    // Put all exp-products together:
    subst_var_t_pi(f_t,toklist,nx,nvar);
    // Replace exp:
    for (i=n-1;i>=0;--i) {
        if (f_t[i].data=="exp") {
            k=i;
            splitfunction(f_t,k);
            tok.del();
            tok.reserve(k-i);
            for (int j=i;j<k;++j) tok.push_back(f_t[j]);
            l=toklist.length()+1;
            toklist.push_back(std::move(tok));
            //Replace operators from i to k-1 by xl:
            f_t[i] = token("x"+val::ToString(l),1);
            for (int j=k;j<n;++j) f_t[j-k+i+1]=std::move(f_t[j]);
            n-=k-i-1;
            f_t.resize(n);
            
        }
    }
    //std::cout<<"\n Nach exp-Subst. f_t : ";
    //for (auto& value : f_t) std::cout<<value.data<<" ";

    // Replace other operators:
    n = f_t.length();
    for (i=n-1;i>=0;--i) {
        if (f_t[i].type==2 && f_t[i].data!="+" && f_t[i].data!="-" && f_t[i].data!="*" && f_t[i].data!="/" && f_t[i].data!="m" && f_t[i].data[0]!='x') {
            k=i;

            if (f_t[i].data=="^") {
                if (i>= n-1) continue;
                if (f_t[i+1].type==0 && f_t[i+1].data!="t" && isinteger(f_t[i+1].data)) continue;
            }

            tok = splitfunction(f_t,k);
            //std::cout<<"\n opertor: "<<f_t[i].data<<", replace: ";
            //for (const auto &v : tok) std::cout<<v.data<<" ";
            //tok.del();
            //tok.reserve(k-i);
            //for (int j=i;j<k;++j) tok.push_back(f_t[j]);
            l=toklist.length()+1;
            toklist.push_back(std::move(tok));
            //Replace operators from i to k-1 by xl:
            tok.del(); tok.reserve(1); tok.push_back(token("x"+val::ToString(l),1));
            squeeze(f_t,tok,i,k);
            n = f_t.length();
            //std::cout<<"\n After replace, f_t = ";
            //for (const auto &v : f_t) std::cout<<v.data<<"  ";
            /*
            f_t[i] = token("x"+val::ToString(l),2);
            for (int j=k;j<n;++j) f_t[j-k+i+1]=std::move(f_t[j]);
            n-=k-i-1;
            f_t.resize(n);
            */
        }
    }
    //std::cout<<"\n Nach äußere Subst. f_t : ";
    //for (auto& value : f_t) std::cout<<value.data<<" ";
    simplifypolynomial(f_t);

    back_subst(f_t,toklist,nx);
    n=f_t.length();
    //std::cout<<"\n Nach exp-simplifypolynomial, f_t : ";
    //for (auto& value : f_t) std::cout<<value.data<<" ";
    // All exp-products are together.

    // exp-products:
    do {
        found=0;
        for (i=0;i<n;++i) {
            if (f_t[i].data=="*" || f_t[i].data=="/") {
                k=i+1;
                h_t=splitfunction(f_t,k);
                tok=splitfunction(f_t,k);
                if (h_t[0].data=="exp" && tok[0].data=="exp") {
                    m1=h_t.length(),m2=tok.length(),m=m1+m2-1;
                    h1.del();
                    h1.reserve(m);
                    if (f_t[i].data=="*") h1[0] = token("+",2);
                    else h1[0] = token("-",2);
                    for (l=1;l<m1;++l) h1[l] = h_t[l];
                    for (l=1;l<m2;++l) h1[l+m1-1] = tok[l];
                    f_t[i].data="exp";
                    squeeze(f_t,h1,i+1,k);
                    n=f_t.length();
                    found=1;
                    //std::cout<<"\nexp simplified f_t = ";
                    //for (const auto& value : f_t) std::cout<<value.data + " ";
                    break;
                }
                if (f_t[i].data=="*" &&  h_t[0].data=="exp" && tok[0].data=="*" && tok[1].data=="exp") {
					l=1;
					m1=h_t.length();m2=tok.length();
					h1.del();
					h1.reserve(m=m1+m2-1);
					tok1 = splitfunction(tok,l);
					tok2= splitfunction(tok,l);
					h1.push_back(h_t[0]); h1.push_back(token("+",2));
					for (l=1;l<m1;++l) h1.push_back(h_t[l]);
					for (l=1;l<tok1.length();++l) h1.push_back(tok1[l]);
					for (l=0;l<tok2.length();++l) h1.push_back(tok2[l]);
					squeeze(f_t,h1,i+1,k);
					n=f_t.length();
					found=1;
					break;
				}
            }
        }
    }
    while (found);
}


void valfunction::simplify_log(d_array<token> &f_t,int extended)
{
    int i,k,m,m1,m2,n=f_t.length(),p;
    d_array<token> h_t,tok,tok1;

    // exp(log):
    if (extended) {
        for (i=0;i<n-1;++i) {
            if (f_t[i].data=="exp" && f_t[i+1].data=="log") {
                k=i+2;
                tok=splitfunction(f_t,k);
                squeeze(f_t,tok,i,k);
                n=f_t.length();
                --i;
            }
        }
    }

    for (i = 0; i < n-1; ++i) {
        if (f_t[i].data != "log") continue;
        if (f_t[i+1].data == "1") {
            k = i + 2;
            tok.reserve(1); tok.push_back(token("0",0));
            squeeze(f_t,tok,i,k);
            i = k-1;
            n = f_t.length();
        }
        else if (f_t[i+1].data == "^") {                        // powers
            if (i<n-2 && f_t[i+2].data == "m") {
                f_t[i].data = "m"; f_t[i+1].data = "log"; f_t[i+2].data = "^";
                continue;
            }
            k=i+2;
            h_t=splitfunction(f_t,k);
            tok=splitfunction(f_t,k);
            m1 = h_t.length(); m2 = tok.length();
            if (!extended && m1 == 1 && h_t[0].type == 0) {
                if (val::isinteger(h_t[0].data)) {
                    p = val::FromString<int>(h_t[0].data);
                    //std::cout<<"\n p = "<<p<<std::endl;
                    if (p % 2 == 0) {
                        int q = p / 2;
                        if (q == 1) continue;
                        else {
                            h_t[0].data = val::ToString(2);
                            tok1.reserve(m2+4);
                            f_t[i].data = "*";
                            tok1.push_back(token("log",2));
                            tok1.push_back(token("^",2));
                            tok1.append(std::move(h_t));
                            tok1.append(std::move(tok));
                            tok1.push_back(token(val::ToString(q),0));
                            squeeze(f_t,tok1,i+1,k);
                            i = k -1;
                            n = f_t.length();
                            continue;
                        }
                    }
                }
            }
            m = m1 + m2;
            tok1.reserve(m+1);
            f_t[i].data = "*";
            tok1.push_back(token("log",2));
            tok1.append(std::move(tok));
            tok1.append(std::move(h_t));
            squeeze(f_t,tok1,i+1,k);
            i=k-1;
            n=f_t.length();
        }
        else if (f_t[i+1].data == "/") {                        // log(1/h(x))
            k = i+2;
            h_t=splitfunction(f_t,k);
            tok=splitfunction(f_t,k);
            if (tok[0].data != "1") continue;
            f_t[i].data = "m"; f_t[i+1].data = "log";
            squeeze(f_t,h_t,i+2,k);
            i = k-1;
            n = f_t.length();
        }
    }
}


void valfunction::simplify_sqrt(d_array<token> &f_t, int nvar, int prod)
{
    int i,k,m,n=f_t.length(),p,q, found,m1,m2,l, nx =0;
    d_array<token> h_t,tok,tok1,tok11,tok12,tok2,h1;
    d_array<d_array<token>> toklist;
    
    
    // sqrt(1)
    for (i = 0; i < n; ++i) {
		if (f_t[i].data != "sqrt") continue;
		k = i+1;
		h_t = splitfunction(f_t,k);
		if (h_t[0].type != 0 || h_t[0].data != "1") continue;
		tok.del();
		tok.push_back(token("1",0));
		squeeze(f_t,tok,i,k);
		n = f_t.length();
	}
	
    if (prod) {
	    // Put all sqrt-products together:
	    subst_var_t_pi(f_t,toklist,nx,nvar);
	    // Replace sqrt:
	    for (i=n-1;i>=0;--i) {
	        if (f_t[i].data=="sqrt") {
	            k=i;
	            splitfunction(f_t,k);
	            tok.del();
	            tok.reserve(k-i);
	            for (int j=i;j<k;++j) tok.push_back(f_t[j]);
	            l=toklist.length()+1;
	            toklist.push_back(std::move(tok));
	            //Replace operators from i to k-1 by xl:
	            f_t[i] = token("x"+val::ToString(l),1);
	            for (int j=k;j<n;++j) f_t[j-k+i+1]=std::move(f_t[j]);
	            n-=k-i-1;
	            f_t.resize(n);
	        }
	    }
	    //std::cout<<"\n Nach sqrt-Subst. f_t : ";
	    //for (auto& value : f_t) std::cout<<value.data<<" ";
	
	    // Replace other operators:
	    n = f_t.length();
	    for (i=n-1;i>=0;--i) {
	        if (f_t[i].type==2 && f_t[i].data!="+" && f_t[i].data!="-" && f_t[i].data!="*" && f_t[i].data!="/" && f_t[i].data!="m" && f_t[i].data[0]!='x') {
	            k=i;
	
	            if (f_t[i].data=="^") {
	                if (i>= n-1) continue;
	                if (f_t[i+1].type==0 && f_t[i+1].data!="t" && isinteger(f_t[i+1].data)) continue;
	            }	
	            splitfunction(f_t,k);
	            tok.del();
	            tok.reserve(k-i);
	            for (int j=i;j<k;++j) tok.push_back(f_t[j]);
	            l=toklist.length()+1;
	            toklist.push_back(std::move(tok));
	            //Replace operators from i to k-1 by xl:
	            f_t[i] = token("x"+val::ToString(l),1);
	            for (int j=k;j<n;++j) f_t[j-k+i+1]=std::move(f_t[j]);
	            n-=k-i-1;
	            f_t.resize(n);
	        }
	    }
	    //std::cout<<"\n Nach äußere Subst. f_t : ";
	    //for (auto& value : f_t) std::cout<<value.data<<" ";
	    simplifypolynomial(f_t);
	
	    back_subst(f_t,toklist,nx);
	    n=f_t.length();
	    //std::cout<<"\n after back-substitution f_t : ";
	    //for (auto& value : f_t) std::cout<<value.data<<" ";
		int isequal;
	    do {
	        found=0;
	        for (i=0;i<n;++i) {
	            if (f_t[i].data=="*" || f_t[i].data=="/") {
	                k=i+1;
	                h_t=splitfunction(f_t,k);
	                tok=splitfunction(f_t,k);
	                if (h_t[0].data=="sqrt" && tok[0].data=="sqrt") {
	                    //std::cout<<"\nHere!";
	                    m1=h_t.length(),m2=tok.length(),m=m1+m2-1;
	                    isequal = 1;
	                    if (m1 != m2) isequal = 0;
	                    else {
							for (int j = 1; j < m1; j++) {
								if (h_t[j].data != tok[j].data || h_t[j].type != tok[j].type) {
									isequal = 0;
									break;
								}
							}
						}	 
	                    h1.del();
	                    if (isequal) {
							if (f_t[i].data == "*") {
								h1.reserve(m1-1);
								for (int j = 1; j < m1; ++j) h1.push_back(tok[j]);
							}
							else {
								h1.reserve(1); h1.push_back(token("1",0));
							}
							squeeze(f_t,h1,i,k);
							n = f_t.length();
							break;
						}
	                    h1.reserve(m);
	                    if (f_t[i].data=="*") h1[0] = token("*",2);
	                    else h1[0] = token("/",2);
	                    for (l=1;l<m1;++l) h1[l] = h_t[l];
	                    for (l=1;l<m2;++l) h1[l+m1-1] = tok[l];
	                    f_t[i].data="sqrt";
	                    squeeze(f_t,h1,i+1,k);
	                    n=f_t.length();
	                    found=1;
	                    //std::cout<<"\nexp simplified f_t = ";
	                    //for (const auto& value : f_t) std::cout<<value.data + " ";
	                    break;
	                }
	                if (f_t[i].data=="*" &&  h_t[0].data=="sqrt" && tok[0].data=="*" && tok[1].data=="sqrt") {
						m1 = h_t.length();
						l = 1;
						tok11 = splitfunction(tok,l);
						m2 = tok11.length();
						isequal = 1;
						if (m1 != m2) isequal = 0;
						else {
							for (int j = 1; j < m1; ++j) {
								if (h_t[j].data != tok11[j].data || h_t[j].type != tok11[j].type) {
									isequal = 0;
									break;
								}
							}
						}
						if (isequal) {
							tok12 = splitfunction(tok,l);
							m2 = h_t.length() + tok12.length() -1;
							h1.del();
							h1.reserve(m2);
							for (int j = 1; j < m1; ++j) h1.push_back(h_t[j]);
							for (int j = 0; j < tok12.length(); ++j) h1.push_back(tok12[j]);
							squeeze(f_t,h1,i+1,k);
							n = f_t.length();
							break;
						}
						
						l=1;
						m2=tok.length();
						h1.del();
						h1.reserve(m=m1+m2-1);
						tok1 = splitfunction(tok,l);
						tok2= splitfunction(tok,l);
						h1.push_back(h_t[0]); h1.push_back(token("*",2));
						for (l=1;l<m1;++l) h1.push_back(h_t[l]);
						for (l=1;l<tok1.length();++l) h1.push_back(tok1[l]);
						for (l=0;l<tok2.length();++l) h1.push_back(tok2[l]);
						squeeze(f_t,h1,i+1,k);
						n=f_t.length();
						found=1;
						break;
					}
					if (f_t[i].data == "/" && h_t[0].data=="sqrt" && tok[0].data=="*" && tok[1].data=="sqrt") {
						m1 = h_t.length();
						l = 1;
						tok11 = splitfunction(tok,l);
						m2 = tok11.length();
						isequal = 1;
						if (m1 != m2) isequal = 0;
						else {
							for (int j = 1; j < m1; ++j) {
								if (h_t[j].data != tok11[j].data || h_t[j].type != tok11[j].type) {
									isequal = 0;
									break;
								}
							}
						}
						if (isequal) {
							tok12 = splitfunction(tok,l);
							m2 = tok12.length();
							h1.del();
							h1.reserve(m2);
							for (int j = 0; j < tok12.length(); ++j) h1.push_back(tok12[j]);
							squeeze(f_t,h1,i,k);
							n = f_t.length();
							break;
						}

						m = h_t.length() + tok.length() -2;
						f_t[i].data = "*";
						h1.del();
						h1.reserve(m);
						h1.push_back(token("sqrt",2)); h1.push_back(token("/",2));
						for (l = 1; l < h_t.length(); ++l) h1.push_back(h_t[l]);
						for (l = 2; l < tok.length(); ++l) h1.push_back(tok[l]);
						squeeze(f_t,h1,i+1,k);
						n=f_t.length();
						found=1;
						break;
					}						
	            }
	        }
	    }
	    while (found);	
	    //std::cout<<"\n after product f_t : ";
	    //for (auto& value : f_t) std::cout<<value.data<<" ";
	}

    
    // sqrt (h^n)
    for (i = 0; i < n - 2; ++i) {
        if (f_t[i].data == "sqrt" && f_t[i+1].data == "^" && f_t[i+2].type == 0) {
            if (!val::isinteger(f_t[i+2].data)) continue;
            p = val::FromString<int>(f_t[i+2].data);
            if (p < 3) continue;
            if (p % 2 == 0) {
                q = (p - 2)/2;
                p = 2;
            }
            else {
                q = (p-1)/2;
                p = 1;
            }
            k = i+2;
            h_t=splitfunction(f_t,k);
            tok=splitfunction(f_t,k);
            m = tok.length();
            if (p == 1) {
                if (q == 1) {
                    tok1.reserve(2*m+2);
                    tok1.push_back(token("*",2)); tok1.push_back(token("sqrt",2));
                    tok1.append(tok); tok1.append(tok);
                }
                else {
                    tok1.reserve(2*m+4);
                    tok1.push_back(token("*",2)); tok1.push_back(token("sqrt",2));
                    tok1.append(tok); tok1.push_back(token("^",2)); tok1.push_back(token(val::ToString(q),0));
                    tok1.append(tok);
                }
            }
            else {
                if (q == 1) {
                    tok1.reserve(2*m+4);
                    tok1.push_back(token("*",2)); tok1.push_back(token("sqrt",2));
                    tok1.push_back(token("^",2)); tok1.push_back(token("2",0));
                    tok1.append(tok); tok1.append(tok);
                }
                else {
                    tok1.reserve(2*m+6);
                    tok1.push_back(token("*",2)); tok1.push_back(token("sqrt",2));
                    tok1.push_back(token("^",2)); tok1.push_back(token("2",0)); tok1.append(tok);
                    tok1.push_back(token("^",2)); tok1.push_back(token(val::ToString(q),0));
                    tok1.append(tok);
                }
            }
            squeeze(f_t,tok1,i,k);
            i = k - 1;
            n = f_t.length();
            tok1.del();
        }
    }
    
    // sqrt(h)^n
    Glist<token> G;
    std::string sf; 
    for (i = 0; i < n - 3; ++i) {
		if (f_t[i].data != "^") continue;
		k = i + 1;
		tok = splitfunction(f_t,k);
		h_t = splitfunction(f_t,k);
		if (h_t[0].data != "sqrt") continue;
		//std::cout<<"\n h_t:\n";
		//for (int l = 0; l < h_t.length(); ++l) std::cout<<h_t[l].data<<"  ";
		G.dellist();
		for (const auto& v : tok) G.push(v);
		sf = get_infix(G);
		if (!val::isinteger(sf)) continue;
		//std::cout<<"\n sf = "<<sf;
		p = FromString<int>(sf);
		q = 0;
		if (p%2) continue;//q = 1;
		p /= 2;
		m = h_t.length();
		if (q) ++m;
		if (p < 0) ++m;
		tok1.reserve(m);
		if (p <  0) tok1.push_back(token("m",2));
		sf = ToString(abs(p));
		tok1.push_back(token(sf,0)); 
		if (q) tok1.push_back(token("sqrt",2));
		for (int j = 1; j < h_t.length(); ++j) tok1.push_back(h_t[j]);
		squeeze(f_t,tok1,i+1,k);
		n = f_t.length();
		tok1.del();
	}
}





void valfunction::simplify_qsin(d_array<token> &f_t)
{
	int i,k,n=f_t.length(),l;
	d_array<token> h1,h2,h;

	for (i=0;i<n-3;++i) {
		k=i+1;
		if (!(f_t[i].data=="^" && f_t[i+1].data=="2" && f_t[i+2].data=="sin")) continue;
		splitfunction(f_t,k);
		h1=splitfunction(f_t,k);
		l=1;
		h2=splitfunction(h1,l);
		//std::cout<<"\nh2 = ";
        //for (const auto& value : h2) std::cout<<value.data + " ";
		h.del();
		h.reserve(h2.length()+5);
		h.push_back(token("-",2)); h.push_back(token("^",2));h.push_back(token("2",0)); h.push_back(token("cos",2));
		for (l=0;l<h2.length();++l) h.push_back(h2[l]);
		h.push_back(token("1",0));
		squeeze(f_t,h,i,k);
		n=f_t.length();
	}
	//std::cout<<"\nNach qsin f = ";
    //for (const auto& value : f_t) std::cout<<value.data + " ";
	//std::cout<<std::endl;
}

valfunction::valfunction(const valfunction& f)
{
    val::GlistIterator<token> it;
    for (it=f.Gdat;it;it++) Gdat.inserttoend(it());
    t=f.t;
    s_infix=f.s_infix;
    nvar=f.nvar;
}




pol<rational> valfunction::getpolynomial() const
{
    pol<rational> v2,value;
    if (Gdat.isempty()) return value;

    GlistIterator<valfunction::token> iT;
    Glist<pol<rational>> G;


    for (iT=Gdat;iT;iT++) {
        G.resetactual();
        if (iT().type==0) {
            if (iT().data=="t") G.inserttohead(pol<rational>(rational(t),0));
            else G.inserttohead(pol<rational>(val::FromString<rational>(iT().data),0));
        }
        else if (iT().type==1) {G.inserttohead(pol<rational>(1,1));}
        else {
            value=pol<rational>();
            if (iT().data=="+") {   //case "+":
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value+=G.actualvalue();G.skiphead();}
                G.inserttohead(value);
            }
            else if (iT().data=="-") {  // case "-":
                if (!G.isempty()) {v2=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();value-=v2;}
                G.inserttohead(value);
            }
            else if (iT().data=="m") {
                if (!G.isempty()) G.actualvalue()=-G.actualvalue();
            }
            else if (iT().data=="*") {  //case "*":
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value*=G.actualvalue();G.skiphead();}
                G.inserttohead(value);
            }
            else if (iT().data=="/") {  //case "/":
                if (!G.isempty()) {v2=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();value/=v2;}
                G.inserttohead(value);
            }
            else if (iT().data=="^") { //case "^":
                if (!G.isempty()) {v2=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {
                    value=G.actualvalue();
                    G.skiphead();
                    value=val::power(value,int(val::integer(v2.LC())));
                }
                G.inserttohead(value);
            }
        }
    }
    G.resetactual();
    value=pol<rational>();
    if (!G.isempty()) value=G.actualvalue();
    return value;
}



n_polynom<rational> valfunction::getn_polynom() const
{
    n_polynom<rational> v2,value;
    if (Gdat.isempty()) return value;

    GlistIterator<valfunction::token> iT;
    Glist<n_polynom<rational>> G;
    n_expo X(nvar);
    int i,j,n;


    for (iT=Gdat;iT;iT++) {
        G.resetactual();
        for (i=0;i<nvar;++i) X[i]=0;
        if (iT().type==0) {
            if (iT().data=="t") G.inserttohead(n_polynom<rational>(rational(t),X));
            else G.inserttohead(n_polynom<rational>(val::FromString<rational>(iT().data),X));
        }
        else if (iT().type==1) {
                j=1;
                if (iT().data=="x" || iT().data=="x1") X[0]=1;
                else {
                    n=val::FromString<int>(fparser::findnumber(iT().data,j));
                    X[n-1]=1;
                }
                G.inserttohead(n_polynom<rational>(1,X));
        }
        else {
            value=n_polynom<rational>();
            if (iT().data=="+") {   //case "+":
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value+=G.actualvalue();G.skiphead();}
                G.inserttohead(value);
            }
            else if (iT().data=="-") {  // case "-":
                if (!G.isempty()) {v2=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();value-=v2;}
                G.inserttohead(value);
            }
            else if (iT().data=="m") {
                if (!G.isempty()) G.actualvalue()=-G.actualvalue();
            }
            else if (iT().data=="*") {  //case "*":
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value*=G.actualvalue();G.skiphead();}
                G.inserttohead(value);
            }
            else if (iT().data=="/") {  //case "/":
                if (!G.isempty()) {v2=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();value/=v2;}
                G.inserttohead(value);
            }
            else if (iT().data=="^") { //case "^":
                if (!G.isempty()) {v2=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {
                    value=G.actualvalue();
                    G.skiphead();
                    value=val::power(value,int(val::integer(v2.LC())));
                }
                G.inserttohead(value);
            }
        }
    }
    G.resetactual();
    value=n_polynom<rational>();
    if (!G.isempty()) value=G.actualvalue();
    return value;
}

s_polynom<rational> valfunction::gets_polynom() const
{
    return val::To_s_polynom(getn_polynom());
}


rationalfunction valfunction::getrationalfunction() const
{
    rationalfunction v2,value=val::zero_element<rationalfunction>();
    if (Gdat.isempty()) return value;

    GlistIterator<valfunction::token> iT;
    Glist<rationalfunction> G;

    for (iT=Gdat;iT;iT++) {
        G.resetactual();
        if (iT().type==0) {
            if (iT().data=="t") G.inserttohead(rationalfunction(pol<rational>(rational(t),0)));
            else G.inserttohead(rationalfunction(pol<rational>(val::FromString<rational>(iT().data),0)));
        }
        else if (iT().type==1) {G.inserttohead(rationalfunction(pol<rational>(1,1)));} //std::cout<<"  variable ";}
        else {
            value=val::zero_element<rationalfunction>();
            if (iT().data=="+") {   //case "+":
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value+=G.actualvalue();G.skiphead();}
                G.inserttohead(value);
            }
            else if (iT().data=="-") {  // case "-":
                if (!G.isempty()) {v2=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();value-=v2;}
                G.inserttohead(value);
            }
            else if (iT().data=="m") {
                if (!G.isempty()) G.actualvalue()=-G.actualvalue();
            }
            else if (iT().data=="*") {  //case "*":
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value*=G.actualvalue();G.skiphead();}
                G.inserttohead(value);
            }
            else if (iT().data=="/") {  //case "/":
                if (!G.isempty()) {v2=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();value/=v2;}
                G.inserttohead(value);
            }
            else if (iT().data=="^") { //case "^":
                if (!G.isempty()) {v2=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {
                    value=G.actualvalue();
                    G.skiphead();
                    value=val::power(value,int(val::integer(v2.nominator().LC())));
                }
                G.inserttohead(value);
            }
        }
    }
    G.resetactual();
    value=val::zero_element<rationalfunction>();
    if (!G.isempty()) value=G.actualvalue();
    return value;
}


double valfunction::operator() (const double& x) const
{
    using namespace val;
    double value=0.0,v2;
    if (Gdat.isempty()) return 0.0;

    GlistIterator<valfunction::token> iT;
    Glist<double> G;
    int i;

    for (iT=Gdat;iT;iT++) {
        G.resetactual();
        if (iT().type==0) {
            if (iT().data=="PI") G.inserttohead(val::PI);
            else if (iT().data=="t") G.inserttohead(t);
            else G.inserttohead(val::FromString<double>(iT().data));
        }
        else if (iT().type==1) {
			if (iT().data=="x" || iT().data=="x1") G.inserttohead(x);
			else G.inserttohead(0);
		}
        else {
            value=0.0;
            if (iT().data=="+") {   //case "+":
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value+=G.actualvalue();G.skiphead();}
                G.inserttohead(value);
            }
            else if (iT().data=="-") {  // case "-":
                if (!G.isempty()) {v2=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();value-=v2;}
                G.inserttohead(value);
            }
            else if (iT().data=="m") {
                if (!G.isempty()) G.actualvalue()*=-1.0;
            }
            else if (iT().data=="*") {  //case "*":
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value*=G.actualvalue();G.skiphead();}
                G.inserttohead(value);
            }
            else if (iT().data=="/") {  //case "/":
                if (!G.isempty()) {v2=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();value/=v2;}
                G.inserttohead(value);
            }
            else if (iT().data=="^") { //case "^":
                if (!G.isempty()) {v2=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {
                    value=G.actualvalue();
                    G.skiphead();
                    if (value<0) {
                        if (fparser::isinteger(v2)) value=val::power(value,int(val::integer(val::rational(v2))));
                        else value=val::exp(value,v2);
                    }
                    else value=val::exp(value,v2);
                }
                G.inserttohead(value);
            }
            else if ((i = fparser::getindexoffunction(iT().data)) != -1) {
				if (!G.isempty()) G.actualvalue()=fparser::functionpairs[i].y(G.actualvalue());
			}
            
            /*
            else if (iT().data=="sqrt") {
                if (!G.isempty()) G.actualvalue()=val::sqrt(G.actualvalue());
            }
            else if (iT().data=="exp") {
                if (!G.isempty()) G.actualvalue()=val::exp(G.actualvalue());
            }
            else if (iT().data=="log") {
                if (!G.isempty()) G.actualvalue()=val::log(G.actualvalue());
            }
            else if (iT().data=="sin") {
                if (!G.isempty()) G.actualvalue()=val::sin(G.actualvalue());
            }            //default: break;
            else if (iT().data=="cos") {
                if (!G.isempty()) G.actualvalue()=val::cos(G.actualvalue());
            }
            else if (iT().data=="tan") {
                if (!G.isempty()) G.actualvalue()=val::tan(G.actualvalue());
            }
            else if (iT().data=="abs") {
                if (!G.isempty()) G.actualvalue()=val::abs(G.actualvalue());
            }
            else if (iT().data=="arcsin") {
                if (!G.isempty()) G.actualvalue()=val::arcsin(G.actualvalue());
            }
            else if (iT().data=="arccos") {
                if (!G.isempty()) G.actualvalue()=val::arccos(G.actualvalue());
            }
            else if (iT().data=="arctan") {
                if (!G.isempty()) G.actualvalue()=val::arctan(G.actualvalue());
            }

            else if (iT().data=="sinh") {
                if (!G.isempty()) G.actualvalue()=val::sinh(G.actualvalue());
            }            //default: break;
            else if (iT().data=="cosh") {
                if (!G.isempty()) G.actualvalue()=val::cosh(G.actualvalue());
            }
            else if (iT().data=="tanh") {
                if (!G.isempty()) G.actualvalue()=val::tanh(G.actualvalue());
            }
            else if (iT().data=="arsinh") {
                if (!G.isempty()) G.actualvalue()=val::arsinh(G.actualvalue());
            }
            else if (iT().data=="arcosh") {
                if (!G.isempty()) G.actualvalue()=val::arcosh(G.actualvalue());
            }
            else if (iT().data=="artanh") {
                if (!G.isempty()) G.actualvalue()=val::artanh(G.actualvalue());
            }
            */
        }
    }
    G.resetactual();
    value=0.0;
    if (!G.isempty()) value=G.actualvalue();
    return value;
}

double valfunction::operator() (const vector<double>& x) const
{
    using namespace val;
    double value=0.0,v2;
    if (Gdat.isempty()) return 0.0;

    GlistIterator<valfunction::token> iT;
    Glist<double> G;
    int i,k;

    for (iT=Gdat;iT;iT++) {
        G.resetactual();
        if (iT().type==0) {
            if (iT().data=="PI") G.inserttohead(val::PI);
            else if (iT().data=="t") G.inserttohead(t);
            else G.inserttohead(val::FromString<double>(iT().data));
        }
        else if (iT().type==1) {
			if (iT().data=="x" || iT().data=="x1") G.inserttohead(x(0));
			else {
				i=1;
				k=val::FromString<int>(fparser::findnumber(iT().data,i));
				G.inserttohead(x(k-1));
			}
		}
        else {
            value=0.0;
            if (iT().data=="+") {   //case "+":
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value+=G.actualvalue();G.skiphead();}
                G.inserttohead(value);
            }
            else if (iT().data=="-") {  // case "-":
                if (!G.isempty()) {v2=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();value-=v2;}
                G.inserttohead(value);
            }
            else if (iT().data=="m") {
                if (!G.isempty()) G.actualvalue()*=-1.0;
            }
            else if (iT().data=="*") {  //case "*":
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value*=G.actualvalue();G.skiphead();}
                G.inserttohead(value);
            }
            else if (iT().data=="/") {  //case "/":
                if (!G.isempty()) {v2=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();value/=v2;}
                G.inserttohead(value);
            }
            else if (iT().data=="^") { //case "^":
                if (!G.isempty()) {v2=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {
                    value=G.actualvalue();
                    G.skiphead();
                    if (value<0) {
                        if (fparser::isinteger(v2)) value=val::power(value,int(val::integer(val::rational(v2))));
                        else value=val::exp(value,v2);
                    }
                    else value=val::exp(value,v2);
                }
                G.inserttohead(value);
            }
            else if ((i = fparser::getindexoffunction(iT().data)) != -1) {
				if (!G.isempty()) G.actualvalue() = fparser::functionpairs[i].y(G.actualvalue());
			}
            /*
            else if (iT().data=="sqrt") {
                if (!G.isempty()) G.actualvalue()=val::sqrt(G.actualvalue());
            }
            else if (iT().data=="exp") {
                if (!G.isempty()) G.actualvalue()=val::exp(G.actualvalue());
            }
            else if (iT().data=="log") {
                if (!G.isempty()) G.actualvalue()=val::log(G.actualvalue());
            }
            else if (iT().data=="sin") {
                if (!G.isempty()) G.actualvalue()=val::sin(G.actualvalue());
            }            //default: break;
            else if (iT().data=="cos") {
                if (!G.isempty()) G.actualvalue()=val::cos(G.actualvalue());
            }
            else if (iT().data=="tan") {
                if (!G.isempty()) G.actualvalue()=val::tan(G.actualvalue());
            }
            else if (iT().data=="abs") {
                if (!G.isempty()) G.actualvalue()=val::abs(G.actualvalue());
            }
            else if (iT().data=="arcsin") {
                if (!G.isempty()) G.actualvalue()=val::arcsin(G.actualvalue());
            }
            else if (iT().data=="arccos") {
                if (!G.isempty()) G.actualvalue()=val::arccos(G.actualvalue());
            }
            else if (iT().data=="arctan") {
                if (!G.isempty()) G.actualvalue()=val::arctan(G.actualvalue());
            }
            */
        }
    }
    G.resetactual();
    value=0.0;
    if (!G.isempty()) value=G.actualvalue();
    return value;
}


valfunction valfunction::operator() (const valfunction& f) const
{
    if (is_zero()) return valfunction();
    std::string s,sf;
    int n=s_infix.length();

    if (f.is_zero()) sf = "0";
    else sf= "(" + f.s_infix + ")";

    for (int i=0;i<n;++i) {
        if (i>0 && i<n-1 && s_infix[i-1]=='e' && s_infix[i] == 'x' && s_infix[i+1] == 'p') s+=s_infix[i];
        else if (s_infix[i]=='x') s+=sf;
        else s+= s_infix[i];
    }

    return valfunction(s);
}

valfunction valfunction::operator() (const vector<valfunction> &F) const
{
	if (is_zero()) return valfunction();
	int i,n= s_infix.length(),k,dim = val::Min(nvar,F.dimension());
	std::string s;
	vector<std::string> sf("0",nvar);

	for (i=0;i<dim;++i) {
        sf(i) = "(" + F(i).s_infix + ")";
	}

	for (i=0;i<n;++i) {
        if (s_infix[i]=='z') s+=sf(2);
        else if (s_infix[i]=='y') s+=sf(1);
        else if (s_infix[i]=='x') {
            if (i<n-1 && s_infix[i+1]>='1' && s_infix[i+1]<='9') {
                ++i;
                k=val::FromString<int>(fparser::findnumber(s_infix,i));
                s+=sf(k-1);
                --i;
            }
            else s+=sf(0);
        }
        else s+= s_infix[i];
	}

	return valfunction (s);
}


// Implementation of Shunting - Yard - algorithm:
const valfunction& valfunction::infix_to_postfix(const std::string &s)
{
    using namespace val;
    int i=0,n=s.size(),nG=0,failed=1,vnumber;//, j, syntax=1;
    std::string out="",s_number, sf;
    Glist<s_stack> G,OP;
    s_stack t,*tlast;

    if (!Gdat.isempty()) Gdat.dellist();
    s_infix="";

    if (!fparser::checkbrackets(s)) return *this;//out;  // Brackets not set correctly => *this = 0-function
    s_infix=s;
    // Set G in infix-notation;

    while (i<n) {
        failed=0;
        //syntax=1;
        if (s[i]=='+') {
            t = s_stack("+",s_stack::OPERATOR,2);
            i++;
        }
        else if (s[i]=='-') {
            if (G.isempty()) t = s_stack("m",s_stack::OPERATOR,2);
            else if (G[nG-1].data==")" || G[nG-1].type==0 || G[nG-1].type==1) t = s_stack("-",s_stack::OPERATOR,2);
            else t = s_stack("m",s_stack::OPERATOR,4);
            i++;
        }
        else if (s[i]=='*' || s[i]=='/') {
            out="";
            out+=s[i];
            t = s_stack(out,s_stack::OPERATOR,3);
            i++;
        }
        else if (s[i]=='^') {
            t = s_stack("^",s_stack::OPERATOR,5,0);
            i++;
        }
        else if (s[i]=='(') {
            t = s_stack("(",s_stack::LBRACKET);
            i++;
        }
        else if (s[i]==')') {
            t = s_stack(")",s_stack::RBRACKET);
            i++;
        }
        else if (s[i]>='0' && s[i] <='9') {
            t = s_stack(fparser::findnumber(s,i),s_stack::NUMBER);
        }
        else if (s[i]=='x') {
            out="x";
            if (i!=n-1 && s[i+1]>='0' && s[i+1]<='9') {
                ++i;
                s_number = fparser::findnumber(s,i);
                out+=s_number;
                vnumber = val::FromString<int>(s_number);
                nvar = val::Max(nvar,vnumber);
            }
            else {++i;}
            t = s_stack(out,s_stack::VARIABLE);
        }
        else if (s[i]=='y') {
			t = s_stack("x2",s_stack::VARIABLE);
			if (nvar<2) nvar=2;
			++i;
		}
        else if (s[i]=='z') {
			t = s_stack("x3",s_stack::VARIABLE);
			if (nvar<3) nvar=3;
			++i;
		}
        else if (s[i] == 't' && i < n-1 && s[i+1] != 'a') {
			t = s_stack("t",s_stack::NUMBER);
			++i;
		}
		else if ((sf = fparser::getstringfunction(s,i)) != "") {
			t = s_stack(sf,s_stack::OPERATOR,5);
			i+= sf.length();
		}
		/*
        else if (s[i]=='e') {
            if (fparser::foundpattern(s,"exp",i)) {
                t = s_stack("exp",s_stack::OPERATOR,5);
            }
            else syntax=0;
        }
        else if (s[i]=='l') {
            if (fparser::foundpattern(s,"log",i)) {
                t = s_stack("log",s_stack::OPERATOR,5);
            }
            else syntax=0;
        }
        else if (s[i]=='s') {
            j=i;
            if (fparser::foundpattern(s,"sin",i)) {
                t = s_stack("sin",s_stack::OPERATOR,5);
            }
            else {
                i=j;
                if (fparser::foundpattern(s,"sqrt",i)) {
                    t = s_stack("sqrt",s_stack::OPERATOR,5);
                }
                else syntax=0;
            }
        }
        else if (s[i]=='c') {
            if (fparser::foundpattern(s,"cos",i)) {
                t = s_stack("cos",s_stack::OPERATOR,5);
            }
            else syntax=0;
        }
        else if (s[i]=='t') {
            j=i;
            if (fparser::foundpattern(s,"tan",i)) {
                t = s_stack("tan",s_stack::OPERATOR,5);
            }
            else {
                t = s_stack("t",s_stack::NUMBER);
                i=j+1;
            }
        }
        else if (s[i]=='P') {
            if (fparser::foundpattern(s,"PI",i)) {
                t = s_stack("PI",s_stack::NUMBER);
            }
            else syntax=0;
        }
        else if (s[i]=='a') {
            if (i!=n-1 && s[i+1]=='b') {
                if (fparser::foundpattern(s,"abs",i)) {
                    t = s_stack("abs",s_stack::OPERATOR,5);
                }
                else syntax=0;
            }
            else if (fparser::foundpattern(s,"arc",i)) {
                syntax=0;
                if (s[i]=='s') {
                    if (fparser::foundpattern(s,"sin",i)) {
                        t = s_stack("arcsin",s_stack::OPERATOR,5);
                        syntax=1;
                    }
                }
                else if (s[i]=='c') {
                    if (fparser::foundpattern(s,"cos",i)) {
                        t = s_stack("arccos",s_stack::OPERATOR,5);
                        syntax=1;
                    }
                }
                else if (s[i]=='t') {
                    if (fparser::foundpattern(s,"tan",i)) {
                        t = s_stack("arctan",s_stack::OPERATOR,5);
                        syntax=1;
                    }
                }
            }
            else syntax=0;
        }
        */
        else {failed=1;++i;}
        
        /*
        if (!syntax) {
            s_infix="";
            return *this;
        }
        */
        
        if (!failed) {
            if (!G.isempty()) {
                tlast = &(G[nG-1]);
                if (t.data=="(") {
                    if (tlast->type==s_stack::NUMBER || tlast->data==")" || tlast->type==s_stack::VARIABLE) {
                        G.inserttoend(s_stack("*",s_stack::OPERATOR,3));
                        nG++;
                    }
                }
                else if (t.type==s_stack::OPERATOR && t.precedence==5 && t.data!="^") {
                    if (tlast->type==s_stack::NUMBER || tlast->data==")" || tlast->type==s_stack::VARIABLE) {
                        G.inserttoend(s_stack("*",s_stack::OPERATOR,3));
                        nG++;
                    }
                }
                else if (t.type==s_stack::VARIABLE) {
                    if (tlast->type==s_stack::NUMBER || tlast->data==")" || tlast->type==s_stack::VARIABLE) {
                        G.inserttoend(s_stack("*",s_stack::OPERATOR,3));
                        nG++;
                    }
                }
                else if (t.type==s_stack::NUMBER) {
                    if (tlast->type==s_stack::NUMBER || tlast->data==")" || tlast->type==s_stack::VARIABLE) {
                        G.inserttoend(s_stack("*",s_stack::OPERATOR,3));
                        nG++;
                    }
                }
            }
            G.inserttoend(t);
            nG++;
        }
    }

    // Apply now Hunting-Yard-algorithm:

    G.resetactual();
    while (!G.isempty()) {
         t = G.actualvalue(); G.skiphead();G.resetactual();
         if (t.type==s_stack::NUMBER || t.type==s_stack::VARIABLE) Gdat.inserttoend(valfunction::token(t.data,int(t.type)));
         else if (t.type==s_stack::OPERATOR) {
            while (!OP.isempty() && OP.actualvalue().type==s_stack::OPERATOR) {
                if (OP.actualvalue().precedence>t.precedence  || (OP.actualvalue().precedence==t.precedence && OP.actualvalue().leftassociativ)) {
                    Gdat.inserttoend(valfunction::token(OP.actualvalue().data,2));
                    OP.skiphead(); OP.resetactual();
                }
                else break;
            }
            OP.inserttohead(t); OP.resetactual();
         }
         else if (t.type==s_stack::LBRACKET) {OP.inserttohead(t);OP.resetactual();}
         else if (t.type==s_stack::RBRACKET) {
            while (!OP.isempty() && OP.actualvalue().type!=s_stack::LBRACKET) {
                Gdat.inserttoend(valfunction::token(OP.actualvalue().data,int(OP.actualvalue().type)));
                OP.skiphead();OP.resetactual();
            }
            if (!OP.isempty()) {OP.skiphead();OP.resetactual();}
         }
    }

    while (!OP.isempty()) {
         Gdat.inserttoend(valfunction::token(OP.actualvalue().data,int(OP.actualvalue().type)));
        OP.skiphead();OP.resetactual();
    }

    return *this;
}



void valfunction::print() const
{
    val::GlistIterator<valfunction::token> It;
    for (It=Gdat;It;It++) std::cout<<It().data<<" ";
}


int valfunction::isconst() const
{
	for (const auto &value : Gdat) {
		if (value.data[0]=='x' || value.data[0]=='y' || value.data[0]=='z') return 0;
	}
	return 1;
}

int valfunction::isrationalfunction() const
{
    int i=0;

    for (const auto& value : Gdat) {
        if (value.type==2) {
            if (value.data!="m" && value.data!="+" && value.data!="-" && value.data!="*" && value.data!="/" && value.data!="^") return 0;
            if (value.data=="^") {
                if (i==0) return 0;
                if (Gdat[i-1].type!=0) return 0;
                const std::string &s=Gdat[i-1].data;
                if (s=="t") return 0;
                if (!isinteger(s)) return 0;
            }
        }
        else if (value.type==0 && (value.data=="t" || value.data=="PI")) return 0;
        ++i;
    }

    return 1;
}


int valfunction::isdifferentiable() const
{
	for (const auto& value : Gdat) {
		if (value.data=="abs") return 0;
	}
	return 1;
}

int valfunction::islinearfunction() const
{
    for (const auto& value : Gdat) {
        if (value.type==2) {
            if (value.data!="m" && value.data!="+" && value.data!="-" && value.data!="*" && value.data!="/") return 0;
        }
        //else if (value.type==0 && (value.data=="t" || value.data=="PI")) return 0;
    }

    // f is recognized as rational and max degree of variables is 1. Check if f is polynomial:

    int n = s_infix.length(),nklammerauf=0,nklammerzu=0;

    for (int i=0;i<n;++i) {
        if (s_infix[i]=='/') {
            if (i==n-1) return 1;
            ++i;
            for (;i<n;++i) {
                if (s_infix[i]=='(') nklammerauf++;
                if (s_infix[i]=='x' || s_infix[i]=='y' || s_infix[i]=='z') return 0;
                if (s_infix[i]==')') nklammerzu++;
                if (nklammerauf==nklammerzu) break;
            }
        }
    }

    return 1;
}


int valfunction::ispolynomialfunction() const
{
    int n = s_infix.length(),nklammerauf=0,nklammerzu=0;

    if (n==0) return 1;
    if (!isrationalfunction()) return 0;

    for (int i=0;i<n;++i) {
        if (s_infix[i]=='/') {
            if (i==n-1) return 1;
            ++i;
            for (;i<n;++i) {
                if (s_infix[i]=='(') nklammerauf++;
                if (s_infix[i]=='x' || s_infix[i]=='y' || s_infix[i]=='z') return 0;
                if (s_infix[i]==')') nklammerzu++;
                if (nklammerauf==nklammerzu) break;
            }
        }
    }
    return 1;
}

void valfunction::print_prefix() const
{
    int n = Gdat.length();
    if (!n) return;
    val::d_array<std::string> prefix(n);
    auto it=Gdat.begin();

    for (int i=n-1;i>=0;--i,it++) prefix[i] = it().data;
    for (const auto& t : prefix) std::cout<<t<<" ";
}




d_array<std::string> valfunction::get_prefix() const
{
    int n = Gdat.length();
    if (!n) return d_array<std::string>();
    d_array<std::string> prefix(n);
    auto it=Gdat.begin();

    for (int i=n-1;i>=0;--i,it++) prefix[i] = it().data;
    return prefix;

}


std::string valfunction::get_infix(const Glist<valfunction::token>& Gdat,int nvar)
{
    Glist<std::string> G,Gtoken;
    std::string op1,op2,tok1,tok2;
    //int first=1;

    for (const auto& value : Gdat)
    {
        //std::cout<<"\n G = ";
        //for (const auto& x : G) std::cout<<x<<"  ";
        // Push operands
        if (value.type<=1) {
           if (nvar<=3) {
			   if (value.data=="x1") {G.push("x");Gtoken.push("x");}
			   else if (value.data=="x2") {G.push("y");Gtoken.push("y");}
			   else if (value.data=="x3") {G.push("z");Gtoken.push("z");}
			   else {G.push(value.data);Gtoken.push(value.data);}
		   }
           else {G.push(value.data); Gtoken.push(value.data);}
           G.resetactual();
           Gtoken.resetactual();
        }

        // We assume that input is
        // a valid postfix and expect
        // an operator.
        else if (value.data=="+" || value.data=="-" || value.data=="*" || value.data=="/" || value.data=="^") {
            if (G.length()<2) return "";
            op1 = G.getelement();
            G.skiphead();
            op2 = G.getelement();
            G.skiphead();
            tok1 = Gtoken.getelement();
            Gtoken.skiphead();
            tok2= Gtoken.getelement();
            Gtoken.skiphead();
            Gtoken.push(value.data + " " + tok1+ " " +tok2);
            Gtoken.resetactual();
            if (value.data=="-") {
                if (fparser::is_sum_operator(tok1)) {
                    op1 = "(" + op1 + ")";
                }
            }
            else if (value.data=="*") {
                if (fparser::is_sum_operator(tok2) || tok2[0]=='/') op2 = "(" + op2+ ")";
                if (fparser::is_sum_operator(tok1) || tok1[0]=='m') op1 = "(" + op1 + ")";
            }
            else if (value.data=="/") {
                if (fparser::is_sum_operator(tok2)) op2 = "(" + op2+ ")";
                if (fparser::is_arithmetic_operator(tok1)) op1 = "(" + op1 + ")";
            }
            else if (value.data=="^") {
                if (fparser::is_arithmetic_operator(tok2)) op2 = "(" + op2+ ")";
                if (fparser::is_arithmetic_operator(tok1) || tok1[0]=='m') op1 = "(" + op1 + ")";
            }

            G.push(op2 + value.data + op1);
            G.resetactual();
        }
        else {
            if (G.isempty()) return "";
            op1 = G.getelement();G.skiphead();
            op2 = value.data;
            tok1=Gtoken.getelement();Gtoken.skiphead();
            Gtoken.push(value.data + tok1);Gtoken.resetactual();
            if (op2=="m") {
                op2 = "-";
                if (fparser::is_sum_operator(tok1)) op1="(" + op1 + ")";
            }
            else op1 = "(" + op1 + ")";
            G.push(op2 + op1);
            G.resetactual();
        }
    }
    //std::cout<<"\n G = ";
    //for (const auto& x : G) std::cout<<x<<"  ";

    // There must be a single element
    // in stack now which is the required
    // infix.
    if (!G.isempty()) return G.getelement();
    else return std::string("0");
}


std::string valfunction::get_infix() const
{
    return get_infix(Gdat,nvar);
}


d_array<valfunction::token> valfunction::splitfunction(const val::d_array<valfunction::token>& f,int& i)
{
    int n=f.length();
    if (!n) return d_array<valfunction::token>();
    Glist<int> is,soll;
    d_array<valfunction::token> g;
    g.reserve(n);

    for (;i<n;++i) {
        if (f[i].type<=1) {
            if (!is.isempty()) is[0]++;
        }
        else if (f[i].data=="+" || f[i].data=="-" || f[i].data=="*" || f[i].data=="/" || f[i].data=="^") {
            is.push(0);
            soll.push(2);
        }
        else {
            is.push(0);
            soll.push(1);
        }
        g.push_back(f[i]);
        while (!soll.isempty() && soll[0]==is[0]) {
            soll.skiphead();is.skiphead();
            if (!is.isempty()) is[0]++;
        }
        if (is.isempty()) {
            ++i;
            break;
        }
    }
    return g;
}


void valfunction::simplifypolynomial(d_array<token> &f_t)
{
    int n=f_t.length(),j,k,dim;
    fraction<n_polynom<rational>> f,g;
    n_polynom<rational> f_nom,f_denom;
    Glist<fraction<n_polynom<rational>>> G;
    rational r;
    n_expo X;

    for (int i=n-1;i>=0;--i) {
        if (f_t[i].type==0) { // constants
            r=val::FromString<rational>(f_t[i].data);
            G.push(fraction<n_polynom<rational>>(n_polynom<rational>(r)));
        }
        else if (f_t[i].data[0]=='x') { // variables
            j=1;k=val::FromString<int>(fparser::findnumber(f_t[i].data,j));
            if (k<=0) return;
            X=n_expo(0,k);X[k-1]=1;
            r=rational(1);
            G.push(fraction<n_polynom<rational>>(n_polynom<rational>(r,X)));
        }
        else if (f_t[i].data=="+" && G.length()>=2) {
            G.resetactual();
            g=G.actualvalue();G.skiphead();f=G.actualvalue();G.skiphead();
            G.push(f+g);
        }
        else if (f_t[i].data=="-" && G.length()>=2) {
            G.resetactual();
            g=G.actualvalue();G.skiphead();f=G.actualvalue();G.skiphead();
            G.push(f-g);
        }
        else if (f_t[i].data=="*" && G.length()>=2) {
            G.resetactual();
            g=G.actualvalue();G.skiphead();f=G.actualvalue();G.skiphead();
            G.push(f*g);
        }
        else if (f_t[i].data=="/" && G.length()>=2) {
            G.resetactual();
            g=G.actualvalue();G.skiphead();f=G.actualvalue();G.skiphead();
            G.push(f/g);
        }
        else if (f_t[i].data=="^" && G.length()>=2) {
            G.resetactual();
            g=G.actualvalue();G.skiphead();f=G.actualvalue();G.skiphead();
            G.push(val::power(f,int(nominator(g.nominator().LC()))));
        }
        else if (f_t[i].data=="m" && G.length()>=1) {
            G.resetactual();
            G.actualvalue() = -G.actualvalue();
        }
    }
    f = fraction<n_polynom<rational>>();
    G.resetactual();
    if (!G.isempty()) f= std::move(G.actualvalue());

    // Decomposition of numerator and denominator (split x_i-powers):
    dim=0;
    f_nom=f.nominator();
    f_denom = f.denominator();
    n_polynom<rational> hgcd = ::gcd(f_nom,f_denom);
    f_nom/=hgcd;
    f_denom/=hgcd;
    for (const auto& monom : f_nom) {
        dim = val::Max(dim,monom.actualterm().dimension());
    }
    for (const auto& monom : f_denom) {
        dim = val::Max(dim,monom.actualterm().dimension());
    }

    // Case: univariate polynomials:
    if (dim<=1) {
        //std::cout<<"\n OK!";
        pol<rational> uf,ug;
        std::string sF;
        valfunction fz,fn;

        for (const auto& value : f_nom) uf.insert(value.actualcoef(),value.actualterm()[0]);
        for (const auto& value : f_denom) ug.insert(value.actualcoef(),value.actualterm()[0]);
        rationalfunction F(uf,ug);
        ug = F.denominator(); uf = F.nominator();
        if (ug.degree()==0) {
            r = ug.LC(); r = rational(1)/r;
            uf*= r;
            sF = val::PolToString(uf);
            fz.infix_to_postfix(sF);
        }
        else {
            pol<integer> iuf,iug;
            rational cuf,cug;

            primitivpart(uf,iuf,cuf);
            primitivpart(ug,iug,cug);
            r=cuf/cug;
            iuf*=nominator(r);
            iug*=denominator(r);
            sF = val::PolToString(iuf);
            fz.infix_to_postfix(sF);
            sF = val::PolToString(iug);
            fn.infix_to_postfix(sF);
            fz.Gdat.append(std::move(fn.Gdat));
            fz.Gdat.push_back(token("/",2));
        }
        f_t.del();
        f_t.reserve(n=fz.Gdat.length());
        auto it=fz.Gdat.begin();
        for (j=n-1;j>=0;--j,it++) {
            f_t[j] =it();
            if (f_t[j].data=="x") f_t[j].data="x1";
        }
        //std::cout<<"\nf_t = ";
        //for (j=0;j<n;++j) std::cout<<f_t[j].data<<" ";
        return;
    }
    //
    if (totaldegree(f_denom)==0) {
        if (f_denom.LC() != rational(1)) {
            r=f_denom.LC();
            r = rational(1)/r;
            f_nom*=r;
            f_denom*=r;
        }
    }
    else {
        rational rn,rz;
        rz = content(f_nom);
        f_nom*=rational(1)/rz;
        rn= content(f_denom);
        f_denom*=rational(1)/rn;
        r = rz/rn;
        f_nom*=rational(nominator(r));
        f_denom*=rational(denominator(r));
    }

    d_array<int> num_pot(10000000,dim),denum_pot(10000000,dim);
    for (int i=0;i<dim;++i) {
        if (f_nom.iszero()) num_pot[i]=0;
        for (const auto& monom : f_nom) {
            num_pot[i] = val::Min(num_pot[i],monom.actualterm()[i]);
        }
        if (f_denom.iszero()) denum_pot[i] =0;
        for (const auto& monom : f_denom) {
            denum_pot[i] = val::Min(denum_pot[i],monom.actualterm()[i]);
        }
    }
    //std::cout<<"\nnum_pot   = ";for (const auto& value : num_pot) std::cout<<value<<" , ";
    //std::cout<<"\ndenum_pot = ";for (const auto& value : denum_pot) std::cout<<value<<" , ";
    n_expo Y(dim);
    for (int i=0;i<dim;++i) Y[i]=num_pot[i];
    f_nom/=Y;
    for (int i=0;i<dim;++i) Y[i]=denum_pot[i];
    f_denom/=Y;
    //std::cout<<"\nf_nom = "<<MPolToString(f_nom);
    //std::cout<<"\nf_denom = "<<MPolToString(f_denom);
    for (int i=0;i<dim;++i) num_pot[i]-=denum_pot[i];
    //std::cout<<"\nnum_pot   = ";for (const auto& value : num_pot) std::cout<<value<<" , ";

    // Convert to infix:
    std::string s_fnom = MPolToString(f_nom), s_fdenom=MPolToString(f_denom),mult_nom="",mult_denom="";
    int n_nom=0,n_denom=0;
    for (int i=0;i<dim;++i) {
        if (num_pot[i]>0) {
            n_nom++;
            mult_nom+="x" + ToString(i+1);
            if (num_pot[i]>1) mult_nom+="^" + ToString(num_pot[i]);
            mult_nom+="*";
        }
        else if (num_pot[i]<0) {
            n_denom++;
            mult_denom+="x" + ToString(i+1);
            if (num_pot[i]<-1) mult_denom+="^" + ToString(val::abs(num_pot[i]));
            mult_denom+="*";
        }
    }
    if (mult_nom!="") mult_nom.resize(mult_nom.length()-1);
    if (mult_denom!="") mult_denom.resize(mult_denom.length()-1);
    //std::cout<<"\nf_nom = "<<s_fnom<<" ;  mult_nom = "<<mult_nom;
    //std::cout<<"\nf_denom = "<<s_fdenom<<"  ; mult_denom = "<<mult_denom;

    // Convert to postfix
    valfunction fz,mz,fn,mn;
    if (!f_nom.iszero()) {
        fz.infix_to_postfix(s_fnom);
        if (mult_nom!="") {
            mz.infix_to_postfix(mult_nom);
            if (s_fnom=="1") fz.Gdat = std::move(mz.Gdat);
            else if (s_fnom=="-1") {
                fz.Gdat= std::move(mz.Gdat);
                fz.Gdat.push_back(token("m",2));
            }
            else {
                fz.Gdat.appendtoend(mz.Gdat);
                fz.Gdat.push_back(token("*",2));
            }
        }
    }

    if (s_fdenom=="1") {
        if (mult_denom!="") {
            mn.infix_to_postfix(mult_denom);
            fz.Gdat.appendtoend(mn.Gdat);
            fz.Gdat.push_back(token("/",2));
        }
    }
    else {
        fn.infix_to_postfix(s_fdenom);
        if (mult_denom!="") {
            mn.infix_to_postfix(mult_denom);
            fn.Gdat.appendtoend(mn.Gdat);
            fn.Gdat.push_back(token("*",2));
        }
        fz.Gdat.appendtoend(fn.Gdat);
        fz.Gdat.push_back(token("/",2));
    }

    // Convert to prefix
    f_t.del();
    f_t.reserve(dim=fz.Gdat.length());
    auto it=fz.Gdat.begin();
    for (int i=dim-1;i>=0;--i,it++) f_t[i] = it();
}


void valfunction::simplify(int extended)
{
    if (is_zero() || Gdat.isempty() ) return;
    int n=Gdat.length(),nx=0,i,k,found,l,d,ispot=0,isdouble=0,qsin=0,qcos=0;//simplifyagain=0;
    d_array<token> f_t(n),tok,h_t;
    d_array<d_array<token>> toklist;
    auto it=Gdat.begin();
    valfunction h;

    h.nvar=nvar;

    h.s_infix="1";

    for (i=n-1;i>=0;--i,it++) f_t[i] = it();

    //std::cout<<"\n Start: f_t = ";
    //for (const auto& value : f_t) std::cout<<value.data + " ";
    //std::cout<<"\n nvar = "<<nvar;
    //std::cout<<std::endl;

    qsin = qcos = 0;
    for (i=0;i<n-3;++i) {
		if (f_t[i].data=="^" && f_t[i+1].data=="2") {
			if (f_t[i+2].data=="sin") qsin=1;
			if (f_t[i+2].data=="cos") qcos=1;
		}
    }

    if (qsin && qcos) {
        simplify_qsin(f_t);         // trigonometric pythagoras
		n=f_t.length();
	}



    // exp - rules
    if (has_operator(f_t,"exp")) {
        simplify_exp(f_t,nvar,1);      // simplify also exp - products
        n=f_t.length();
    }

    if (has_operator(f_t,"sqrt")) {
        simplify_sqrt(f_t,nvar,1);
        n = f_t.length();
    }
    
    //std::cout<<"\n After first simplifications: f_t = ";
    //for (const auto& value : f_t) std::cout<<value.data + " ";


    // Inner functions:
    for (i=0;i<n;) {
        if (f_t[i].type==2 && f_t[i].data!="+" && f_t[i].data!="-" && f_t[i].data!="*" && f_t[i].data!="/" && f_t[i].data!="m") {
            if (f_t[i].data=="^") ispot=1;
            else ispot=0;
            //std::cout<<"\n f_t = ";
            //for (const auto& value : f_t) std::cout<<value.data + " ";
            //std::cout<<"\n i = "<<i;
            h.Gdat.dellist();
            k=i;
            ++i;
            h_t = splitfunction(f_t,i);
            //std::cout<<"\n h_t = ";
            //for (const auto& value : h_t) std::cout<<value.data + " ";
            l=i;
            i=k+1;
            for (auto& value : h_t) h.Gdat.push(std::move(value));
            h.simplify(0);
            d=h.Gdat.length();
            h_t.del();
            h_t.reserve(d);
            for (k=d-1,it=h.Gdat.begin();k>=0;it++,--k) h_t[k] = it();
            if (d == 0) {
                h_t.reserve(1); h_t.push_back(token("0",0));
            }
            squeeze(f_t,h_t,i,l);
            i+=d;
            n=f_t.length();

            if (ispot) {
                h.Gdat.dellist();
                k=i;
                h_t = splitfunction(f_t,i);
                l=i;
                i=k;
                for (auto& value : h_t) h.Gdat.push(std::move(value));
                h.simplify(0);
                d=h.Gdat.length();
                h_t.del();
                h_t.reserve(d);
                for (k=d-1,it=h.Gdat.begin();k>=0;it++,--k) h_t[k] = it();
                squeeze(f_t,h_t,i,l);
                i+=d;
                n=f_t.length();
            }

            //std::cout<<"\n f_t = ";
            //for (const auto& value : f_t) std::cout<<value.data + " ";
            //std::cout<<"\n i = "<<i;
        }
        else ++i;
    }

    //std::cout<<"\n After inner substitutions, f_t = ";
    //for (auto& value : f_t) std::cout<<value.data<<" ";

    // outer - functions: ------------------------------------------------

    // exp - rules
    if (has_operator(f_t,"exp")) {
        simplify_exp(f_t,nvar,0);       // don't simplify products
        n=f_t.length();
    }


    // log - rules
    if (has_operator(f_t,"log")) {
        simplify_log(f_t,extended);
        n = f_t.length();
    }
    // sqrt - rules
    if (has_operator(f_t,"sqrt")) {
        simplify_sqrt(f_t);
        n = f_t.length();
    }

    // Create toklist
    subst_var_t_pi(f_t,toklist,nx,nvar);
    //std::cout<<"\n Nach innere Subst. f_t nun: ";
    //for (auto& value : f_t) std::cout<<value.data<<" ";
    //std::cout<<"\n toklist.length() = "<<toklist.length();
    //std::cout<<"\n nvar = "<<nvar;

    // Replace operators:
    for (i=n-1;i>=0;--i) {
        if (f_t[i].type==2 && f_t[i].data!="+" && f_t[i].data!="-" && f_t[i].data!="*" && f_t[i].data!="/" && f_t[i].data!="m") {
            k=i;

            if (f_t[i].data=="^") {
                if (i>= n-1) continue;
                if (f_t[i+1].type==0 && f_t[i+1].data!="t" && isinteger(f_t[i+1].data)) continue;
            }

            //std::cout<<"\n Operator: "<<f_t[k].data;
            splitfunction(f_t,k);
            tok.del();
            tok.reserve(k-i);
            for (int j=i;j<k;++j) tok.push_back(f_t[j]);
            //std::cout<<"\n tok = ";
            //for (const auto&  value : tok) std::cout<<value.data<<" ";

            // check if tok is already in toklist.
            l=0;
            found=0;
            //std::cout<<std::endl;
            //for (const auto &t : tok) std::cout<<t.data<<" ";
            for (const auto& value : toklist) {
                //std::cout<<std::endl;
                found=0;
                //for (const auto &t : value) std::cout<<t.data<<" ";
                if (value.length() != tok.length()) {++l;continue;}
                found=1;
                for (int j=0;j<value.length();++j) {
                    //std::cout<<tok[j].data<<"  "<<value[j].data;
                    if (tok[j].data!=value[j].data) {
                        found=0;break;
                    }
                }
                if (found) {break;}
                ++l;
            }
            if (!found) {
                l=toklist.length()+1;
                toklist.push_back(std::move(tok));
            }
            else ++l;
            // Replace operators from i to k-1 by xl:
            f_t[i] = token("x"+val::ToString(l),1);
            for (int j=k;j<n;++j) f_t[j-k+i+1]=std::move(f_t[j]);
            n-=k-i-1;
            f_t.resize(n);
            //std::cout<<"\n Substitution i , k , n: "<<i<<" , "<<k<<" , "<<n<<"  ,f_t = ";
            //for (auto& value : f_t) std::cout<<value.data<<" ";
        }
    }
    //std::cout<<"\n After exterior substitution, f_t =  ";
    //for (auto& value : f_t) std::cout<<value.data<<" ";
    //std::cout<<"\n";
    for (int i=0;i<n;++i) {
        if (f_t[i].type==0 && fparser::has_point(f_t[i].data)) {
            isdouble=1;
            break;
        }
    }
    simplifypolynomial(f_t);
    //std::cout<<"\n Nach simplifypol f_t: ";
    //for (auto& value : f_t) std::cout<<value.data<<" ";

    // Back-substitution:
    back_subst(f_t,toklist,nx);
    //
    if (isdouble) to_double(f_t);

    Gdat.dellist();
    //std::cout<<"\n f_t nach Rücksubstitution: ";
    //for (auto& value : f_t) std::cout<<value.data<<" ";
    std::string s_number;


	nvar=1;
    for (auto& value : f_t) {
		if (value.data != "x" && value.data[0]=='x') {
			i=1;
			s_number = fparser::findnumber(value.data,i);
			nvar=val::Max(nvar,val::FromString<int>(s_number));
		}
		Gdat.push(value);
	}
    s_infix=get_infix(Gdat,nvar);


    //std::cout<<"\n s_infix at end of simplify = "<<s_infix;

    //if (simplifyagain) simplify(extended);
    //std::cout<<"\n f_t = "<<s_infix;
}


valfunction valfunction::operator+ (const valfunction& g) const
{
    if (g.is_zero()) return *this;
    if (is_zero()) return g;
    return valfunction(s_infix + "+" + g.s_infix);
}


valfunction valfunction::operator- (const valfunction& g) const
{
    if (g.is_zero()) return *this;
    std::string s_g = g.get_prefix()[0];
    return valfunction(s_infix + "-(" + g.s_infix + ")");
}


valfunction valfunction::operator -() const
{
    if (is_zero()) return *this;
    return valfunction("-(" + s_infix +")");
}


valfunction valfunction::operator* (const valfunction& g) const
{
    if (g.is_zero()) return g;
    if (is_zero()) return *this;
    if (s_infix=="1") return g;
    if (g.s_infix=="1") return *this;
    if (s_infix=="-1") return -g;
    if (g.s_infix=="-1") return -(*this);
    std::string ps = get_prefix()[0], ps_g = g.get_prefix()[0],s=s_infix,s_g=g.s_infix;

    s = "(" + s + ")";
    s_g = "(" + s_g + ")";
    return valfunction(s + "*" + s_g);
}


valfunction valfunction::operator/ (const valfunction& g) const
{
    if (g.is_zero()) Error::error("\nval::valfunction operator /: Division by zero!!!");
    if (g.s_infix=="1") return *this;
    if (is_zero()) return *this;
    std::string ps = get_prefix()[0], ps_g = g.get_prefix()[0],s=s_infix,s_g=g.s_infix;

    s = "(" + s + ")";
    s_g = "(" + s_g + ")";

    return valfunction(s + "/" + s_g);
}


valfunction valfunction::operator^ (int n) const
{
    if (n==0) return valfunction("1");
    if (is_zero()) return *this;
    if (n==1) return *this;
    std::string ps = get_prefix()[0],s=s_infix,s_n=val::ToString(n);

    s = "(" + s + ")";
    s_n= "(" + s_n + ")";

    return valfunction(s+ "^" + s_n);
}


valfunction valfunction::derive(int k) const
{
    if (Gdat.isempty() || is_zero() || k<=0 || k>nvar) return valfunction();
    int n=Gdat.length(),j=1;
    auto it = Gdat.begin();
    d_array<token> f_t(n),g_t,h_t;
    valfunction g,h,g1,g2,g3;
    std::string svar="x" + val::ToString(k);

    g.nvar = h.nvar = g1.nvar = g2.nvar = g3.nvar = nvar;

    // Case: rational function:
    if (nvar==1 && isrationalfunction()) {
        rationalfunction F=getrationalfunction();
        pol<rational> f_denom = F.denominator();

        if (f_denom.degree()==0) {
            //std::cout<<"\n Polynom!\n";
            rational r = f_denom.LC(); r = rational(1)/r;
            pol<rational> g = r*F.nominator();
            g = g.derive();
            return valfunction(PolToString(g));
        }
        //std::cout<<"\nRationale Funktion!\n";
        rational cont,c1;
        pol<integer> fnom,fdnom;
        valfunction f,g;

        F=fparser::derive(F);
        primitivpart(F.nominator(),fnom,cont);
        primitivpart(F.denominator(),fdnom,c1);
        cont/=c1;
        fnom*=val::nominator(cont);
        fdnom*=val::denominator(cont);
        f.infix_to_postfix(PolToString(toRationalPolynom(fnom))); g.infix_to_postfix(PolToString(toRationalPolynom(fdnom)));
        //std::cout<<"\n fdnom = "<<PolToString(fdnom) + "\n";
        //f.print();std::cout<<std::endl;
        //g.print();std::cout<<std::endl;
        f.Gdat.append(std::move(g.Gdat));
        f.Gdat.push_back(token("/",2));
        f.s_infix=get_infix(f.Gdat,nvar);
        //f.print();std::cout<<std::endl;
        return f;
    }

    for (int i=n-1;i>=0;--i,it++) f_t[i] = it();
    g_t=splitfunction(f_t,j);
    for (const auto& value: g_t) g.Gdat.push(value);
    g.s_infix = get_infix(g.Gdat,nvar);

    if (f_t[0].type==0) return valfunction(); // Constant derivation
    else if (k==1 && (f_t[0].data=="x" || f_t[0].data=="x1")) return valfunction("1");
    else if (k!=1 && f_t[0].data==svar) return valfunction("1");
    else if (f_t[0].data=="+" || f_t[0].data=="-" || f_t[0].data=="*" || f_t[0].data=="/" || f_t[0].data=="^") {
        h_t=splitfunction(f_t,j);
        for (const auto& value: h_t) h.Gdat.push(value);
        h.s_infix = get_infix(h.Gdat,nvar);

        if (f_t[0].data=="+") return (h.derive(k) + g.derive(k));
        if (f_t[0].data=="-") return (h.derive(k) - g.derive(k));
        if (f_t[0].data=="*") {
            return (h.derive(k)*g + h*g.derive(k));
        }
        if (f_t[0].data=="/") {
            return (h.derive(k)*g - h*g.derive(k))/(g^2);
        }
        if (f_t[0].data=="^") {
            if (!has_variable(g_t)) {
                g3=g; g3.Gdat.push_back(token("1",0)); g3.Gdat.push_back(token("-",2));
                g3.simplify(1);
                g2 = h;
                g2.Gdat.append(std::move(g3.Gdat));
                g2.Gdat.push_back(token("^",2));
                g2.simplify(1);
                g2.s_infix=get_infix(g2.Gdat,nvar);
                return g*h.derive(k)*g2;
            }
            else {
                g1=h; g1.Gdat.push_back(token("log",2)); g1.s_infix=get_infix(g1.Gdat,nvar); // g1 = log(h)
                g2=g*g1; g2.s_infix=get_infix(g2.Gdat,nvar);
                g2=g2.derive(k);
                return g2*(*this);
            }
        }
    }
    else if (f_t[0].data=="m") return -g.derive(k);
    else if (f_t[0].data=="exp") return g.derive(k)*(*this);
    else if (f_t[0].data=="log" && f_t[1].data == "abs") {
		h.Gdat.push(f_t[0]);
		for (int i = 2; i < n; ++i) h.Gdat.push(f_t[i]);
		h.s_infix = get_infix(h.Gdat,nvar);
		return(h.derive(k));
	}
    else if (f_t[0].data=="log") return g.derive(k)/g;
    else if (f_t[0].data=="sqrt") {
        return g.derive(k)/(valfunction("2")*(*this));
    }
    else if (f_t[0].data=="sin") {
        h.Gdat = g.Gdat;
        h.Gdat.push_back(token("cos",2));
        //std::cout<<"\nOK: ";
        //for (const auto& value : h.Gdat) std::cout<<value.data<<"  ";
        h.s_infix = get_infix(h.Gdat,nvar);
        return g.derive(k)*h;
    }
    else if (f_t[0].data=="cos") {
        h.Gdat = g.Gdat;
        h.Gdat.push_back(token("sin",2));
        h.s_infix = get_infix(h.Gdat,nvar);
        return -g.derive(k)*h;
    }
    else if (f_t[0].data=="tan") {
        h.Gdat=Gdat;
        h.Gdat.push_back(token("2",0));h.Gdat.push_back(token("^",2));
        h.Gdat.push_back(token("1",0));h.Gdat.push_back(token("+",2));
        h.s_infix=get_infix(h.Gdat,nvar);
        // h = 1 +(*this)^2;
        return g.derive(k)*h;
    }
    else if (f_t[0].data=="arcsin") {
        h.Gdat=g.Gdat;
        h.Gdat.push_back(token("2",0));h.Gdat.push_back(token("^",2));
        h.Gdat.push(token("1",0));h.Gdat.push_back(token("-",2));h.Gdat.push_back(token("sqrt",2));
        // now: h = sqrt(1-g^2);
        h.s_infix=get_infix(h.Gdat,nvar);
        return g.derive(k)/h;
    }
    else if (f_t[0].data=="arccos") {
        h.Gdat=g.Gdat;
        h.Gdat.push_back(token("2",0));h.Gdat.push_back(token("^",2));
        h.Gdat.push(token("1",0));h.Gdat.push_back(token("-",2));h.Gdat.push_back(token("sqrt",2));
        // now: h = sqrt(1-g^2);
        h.s_infix=get_infix(h.Gdat,nvar);
        return -g.derive(k)/h;
    }
    else if (f_t[0].data=="arctan") {
        h.Gdat=g.Gdat;
        h.Gdat.push_back(token("2",0));h.Gdat.push_back(token("^",2));
        h.Gdat.push_back(token("1",0));h.Gdat.push_back(token("+",2));
        // now: h = x^2 +1;
        h.s_infix=get_infix(h.Gdat,nvar);
        return g.derive(k)/h;
    }
    else if (f_t[0].data=="sinh") {
        h.Gdat = g.Gdat;
        h.Gdat.push_back(token("cosh",2));
        //std::cout<<"\nOK: ";
        //for (const auto& value : h.Gdat) std::cout<<value.data<<"  ";
        h.s_infix = get_infix(h.Gdat,nvar);
        return g.derive(k)*h;
    }
    else if (f_t[0].data=="cosh") {
        h.Gdat = g.Gdat;
        h.Gdat.push_back(token("sinh",2));
        h.s_infix = get_infix(h.Gdat,nvar);
        return g.derive(k)*h;
    }
    else if (f_t[0].data=="tanh") {
        h.Gdat=Gdat;
        h.Gdat.push_back(token("2",0));h.Gdat.push_back(token("^",2));
        h.Gdat.push_back(token("1",0));h.Gdat.push_back(token("-",2));
        h.s_infix=get_infix(h.Gdat,nvar);
        // h = 1 +(*this)^2;
        return -g.derive(k)*h;
    }
    else if (f_t[0].data == "arsinh") {
		std::string sg = g.getinfixnotation(), sF = "1/sqrt((" + sg + ")^2 + 1)";
		return g.derive() * valfunction(sF); 
	}
    else if (f_t[0].data == "arcosh") {
		std::string sg = g.getinfixnotation(), sF = "1/sqrt((" + sg + ")^2 - 1)";
		return g.derive() * valfunction(sF); 
	}
    else if (f_t[0].data == "artanh") {
		std::string sg = g.getinfixnotation(), sF = "1/(1 -(" + sg + ")^2)";
		return g.derive() * valfunction(sF); 
	}
		
    else return valfunction();

    return g;
}


valfunction valfunction::getfirstargument() const
{
    valfunction f;
    if (is_zero()) return f;

    int n=Gdat.length(),j=1;
    d_array<token> f_t(n),g_t;
    auto it = Gdat.begin();

    for (int i=n-1;i>=0;--i,it++) f_t[i] = it();
    g_t=splitfunction(f_t,j);
    if (j<n) g_t=splitfunction(f_t,j);
    for (const auto& value: g_t) f.Gdat.push(value);
    f.s_infix = get_infix(f.Gdat);
    return f;
}

valfunction valfunction::getsecondargument() const
{
    valfunction f;
    if (is_zero()) return f;

    int n=Gdat.length(),j=1;
    d_array<token> f_t(n),g_t;
    auto it = Gdat.begin();

    for (int i=n-1;i>=0;--i,it++) f_t[i] = it();
    g_t=splitfunction(f_t,j);
    for (const auto& value: g_t) f.Gdat.push(value);
    f.s_infix = get_infix(f.Gdat);
    return f;
}



Glist<double> valfunction::double_roots(const double &x1,const double &x2,int iterations,const double &epsilon,int methoditerations) const
{
	Glist<double> d_roots;

	if (is_zero() || nvar>1) return d_roots;

	if (isrationalfunction()) {
        rationalfunction F = getrationalfunction();
        pol<rational> f=F.nominator();
        vector<double> realroots;
        double z;

        f /= val::gcd(f,f.derive());
        pol<double> df = ToDoublePolynom(f);
        realRoots(df,realroots,epsilon);
        for (int i=0;i<realroots.dimension();++i) {
            z=realroots[i];
            if (abs(z)<epsilon) z=0.0;
            d_roots.push_back(z);
        }
        d_roots.sort();
        return d_roots;
	}

    int n=Gdat.length(),j=1;
    auto it = Gdat.begin();
    d_array<token> f_t(n),g_t,h_t;
    Glist<double> g_roots,h_roots;
    valfunction g,h;
    std::string oper;

    g.nvar=h.nvar=1;

    for (int i=n-1;i>=0;--i,it++) f_t[i] = it();
    oper = f_t[0].data;

    if (oper=="+" || oper =="-" || oper == "^" || oper =="sin" || oper =="cos" || oper =="tan" || oper =="arcsin" || oper =="arccos" || oper =="arctan") {
        int i,n;
        if (x2<=x1) return d_roots;
        double x,delta = (x2-x1)/double(iterations),y1,xd,xs1,xs2;

        for (i=0,x=x1;i<iterations-1;++i,x+=delta) {
            if (abs(y1=operator()(x))<=epsilon) {
                if (!d_roots.isempty()) {
                    n = d_roots.length();
                    if (fparser::issimilar(x,0.0,epsilon)) xd=0.0;
                    else xd=x;
                    if (!fparser::issimilar(d_roots[n-1],xd,epsilon)) d_roots.push_back(xd);
                }
                else d_roots.push_back(x);
                continue;
            }
            if (isNaN(y1)) continue;
            xd=x+delta;
            xs1=x; xs2=xd;
            if (SecantMethod(*this,xs1,xs2,epsilon,methoditerations)>=0)  {
                if (fparser::issimilar(xs2,0.0,epsilon)) xs2=0.0;
                if (xs2>=x && xs2<=xd) d_roots.push_back(xs2);
            }
        }
        if (abs(operator()(x2))<=epsilon) {
            if (!d_roots.isempty()) {
                n = d_roots.length();
                if (!fparser::issimilar(d_roots[n-1],x2,epsilon)) d_roots.push_back(x2);
            }
            else d_roots.push_back(x2);
        }

        return d_roots;
    }
    else if (oper =="exp") return d_roots;

    g_t=splitfunction(f_t,j);
    for (const auto& value: g_t) g.Gdat.push(value);
    g.s_infix = get_infix(g.Gdat,nvar);

    if (oper == "sqrt" || oper == "m" || oper=="log" || oper == "abs") {
        if (oper=="log") g-= valfunction("1");
        return g.double_roots(x1,x2,iterations,epsilon);
    }

    h_t=splitfunction(f_t,j);
    for (const auto& value: h_t) h.Gdat.push(value);
    h.s_infix = get_infix(h.Gdat,nvar);

    if (oper=="*") {
        g_roots=g.double_roots(x1,x2,iterations,epsilon);
        h_roots=h.double_roots(x1,x2,iterations,epsilon);

        d_roots=fparser::unify(g_roots,h_roots,epsilon);
        int i=0;
        double value;
        n=d_roots.length();
        while (i<n) {
			value = operator()(d_roots[0]);
            if (val::isNaN(value)) {d_roots.skiphead();n--;}
            ++i;
        }

        return d_roots;
    }
    if (oper=="/") {
		d_roots=h.double_roots(x1,x2,iterations,epsilon);
        int i=0;
        double value;
        n=d_roots.length();
        while (i<n) {
			value = g(d_roots[0]);
            if (val::isNaN(value) || value == val::Inf || value == -val::Inf || value == 0.0)  {d_roots.skiphead();n--;}
            ++i;
        }
	}


	return d_roots;
}

Glist<GPair<double>> valfunction::get_undefined_intervals(const double &x1,const double &x2,int iterations,const double &epsilon,int methoditerations) const
{
	Glist<GPair<double>> intervals;
	if (is_zero() || nvar>1) return intervals;

	if (isrationalfunction()) {
        rationalfunction F = getrationalfunction();
        pol<rational> f= F.denominator();
        vector<double> realroots;
        double z;

        f /= val::gcd(f,f.derive());
        //std::cout<<"\n f = \n"<<f;
        pol<double> df = ToDoublePolynom(f);
        realRoots(df,realroots,epsilon);
        realroots.sort();
        for (int i=0;i<realroots.dimension();++i)  {
            z = realroots(i);
            if (abs(z)<epsilon) z = 0.0;
            intervals.push_back(GPair<double>(z,z));
        }
        return intervals;
	}

    int n=Gdat.length(),i,j=1;
    auto it = Gdat.begin();
    d_array<token> f_t(n),g_t,h_t;
    Glist<GPair<double>> g_interval,h_interval;
    valfunction g,h;
    std::string oper;

    g.nvar=h.nvar=1;

    for (i=n-1;i>=0;--i,it++) f_t[i] = it();
    oper = f_t[0].data;

    g_t = splitfunction(f_t,j);
    for (const auto& value: g_t) g.Gdat.push(value);
    g.s_infix = get_infix(g.Gdat,nvar);

    if (oper=="abs" || oper== "exp" || oper=="sin" || oper=="cos" || oper=="m") {
        return g.get_undefined_intervals(x1,x2,iterations,epsilon,methoditerations);
    }
    if (oper=="log" || oper=="sqrt") {
        Glist<double> roots = g.double_roots(x1,x2,iterations,epsilon,methoditerations);
        int israt=g.isrationalfunction(),n=roots.length();
        double x,y;

        for (auto & z : roots)  {
            if (abs(z)<epsilon) z = 0.0;
            g_interval.push_back(GPair<double>(z,z));
        }

        h_interval = g.get_undefined_intervals(x1,x2,iterations,epsilon,methoditerations);
        g_interval = fparser::unify(g_interval,h_interval,epsilon);
        n=g_interval.length();
        if (n==0) {
            if (israt) x = 0.0;
            else x = (x1+x2)/2.0;
            if ((y=g(x))<0.0 || (abs(y)<epsilon && oper=="log")) {
                if (israt) g_interval.push_back(GPair<double>(-val::Inf,val::Inf));
                else g_interval.push_back(GPair<double>(x1,x2));
            }
            return g_interval;
        }
        h_interval.dellist();

        if (israt) x=g_interval[0].x-1.0;
        else x = (x1+g_interval[0].x)/2.0;
        if ((y=g(x))<0.0 || (abs(y)<epsilon && oper=="log")) {
            if (israt) g_interval[0].x=-val::Inf;
            else g_interval[0].x = x1;
        }
        h_interval.push_back(g_interval[0]);
        j=0;
        for (i=1;i<n;++i) {
            x = (h_interval[j].y + g_interval[i].x)/2.0;
            if ((y=g(x))<0.0 || (abs(y)<epsilon && oper=="log")) {
                h_interval[j].y = g_interval[i].y;
            }
            else {
                h_interval.push_back(g_interval[i]);
                ++j;
            }
        }
        if (israt) x=h_interval[j].y+1.0;
        else x = (h_interval[j].y + x2)/2.0;
        if ((y=g(x))<0.0 || (abs(y)<epsilon && oper=="log")) {
            if (israt) h_interval[j].y=val::Inf;
            else h_interval[j].y=x2;
        }
        return h_interval;
    }

    if (oper=="tan") {
        h=valfunction("cos(" + g.getinfixnotation()+")");
        g_interval = g.get_undefined_intervals(x1,x2,iterations,epsilon,methoditerations);
        Glist<double> roots = h.double_roots(x1,x2,iterations,epsilon,methoditerations);

        for (auto &z : roots) {
            if (abs(z)<epsilon) z=0;
            h_interval.push_back(GPair<double>(z,z));
        }
        return fparser::unify(g_interval,h_interval,epsilon);
    }

    h_t = splitfunction(f_t,j);
    for (const auto& value: h_t) h.Gdat.push(value);
    h.s_infix = get_infix(h.Gdat,nvar);

    if (oper=="*" || oper=="-" || oper=="+") {
        g_interval = g.get_undefined_intervals(x1,x2,iterations,epsilon,methoditerations);
        h_interval = h.get_undefined_intervals(x1,x2,iterations,epsilon,methoditerations);
        return fparser::unify(g_interval,h_interval,epsilon);
    }

    if (oper=="/") {
        Glist<double> roots=g.double_roots(x1,x2,iterations,epsilon,methoditerations);
        g_interval = g.get_undefined_intervals(x1,x2,iterations,epsilon,methoditerations);
        h_interval = h.get_undefined_intervals(x1,x2,iterations,epsilon,methoditerations);
        h_interval=fparser::unify(g_interval,h_interval,epsilon);
        g_interval.dellist();
        for (auto &z : roots) {
            if (abs(z)<epsilon) z =0.0;
            g_interval.push_back(GPair<double>(z,z));
        }
        return fparser::unify(g_interval,h_interval,epsilon);
    }

    if (oper=="^") {
        valfunction f("log("+h.getinfixnotation()+")*(" + g.getinfixnotation()+")");
        return f.get_undefined_intervals(x1,x2,iterations,epsilon,methoditerations);
    }

	return intervals;
}

} //end namespace val
