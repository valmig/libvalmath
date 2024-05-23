#ifndef FUNCTION_PARSER_H_INCLUDED
#define FUNCTION_PARSER_H_INCLUDED

#include <string>
#include <Glist.h>
#include <d_array.h>
#include <pol.h>

namespace fparser
{
std::string findnumber(const std::string &s, int &i);
}


namespace val
{

template <class T> class fraction;
//template <class T> class pol;
template <class T> class vector;
template <class T> class n_polynom;
template <class T> class s_polynom;
class rational;

typedef fraction<pol<rational>> rationalfunction;

class DLL_PUBLIC valfunction
{
private:
    struct token {
        std::string data="";
        int type=0;
        token() = default;
        token(const std::string &s,const int t) : data(s),type(t) {}
        //token(token&& t) :data(std::move(t.data)),type(t.type) {}
        token(const token& t) : data(t.data), type(t.type) {}
        //token& operator= (const token&)=default;
    };
    struct s_stack {
        std::string data="";
        enum data_type {NUMBER,VARIABLE,OPERATOR,LBRACKET,RBRACKET};
        data_type type=NUMBER;
        int precedence=0;
        int leftassociativ=0;
        s_stack() = default;
        s_stack(const std::string& s,data_type t,int p=0,int l=1) : data(s), type(t), precedence(p),leftassociativ(l) {}
    };
    val::Glist<token> Gdat;
    double t=0;
    std::string s_infix="";
    int nvar=1;
    //
    //static int infix_to_postfix(const std::string &s,Glist<token> &Gdat);
    static val::d_array<token> splitfunction(const val::d_array<token>& f,int& i);
    static void squeeze(d_array<token> &f,const d_array<token> &g,int a,int b);
    static void to_double(d_array<token> &f);
    static int has_variable(const d_array<token> &f);
    static int has_operator(const d_array<token> &f,const std::string &op);
    static void simplify_exp(d_array<token> &f,int nvar=1,int prod = 0);
    static void simplify_log(d_array<token> &f,int extended);
    static void simplify_sqrt(d_array<token> &f,int nvar=1,int prod = 0);
    static void simplify_qsin(d_array<token> &f);
    static void subst_var_t_pi(d_array<token> &f_t,d_array<d_array<token>> &toklist,int &nx,int nvar=1);
    static void back_subst(d_array<token> & f_t,const d_array<d_array<token>> &toklist,int nx);
    static std::string get_infix(const Glist<token>& Gdat,int nvar=1);
    static const std::string zero_string;
public:
    valfunction() = default;
    explicit valfunction(const std::string &s,int simp=1) {infix_to_postfix(s);if (simp) {simplify();} else {s_infix = get_infix(Gdat,nvar);}}
    valfunction(const valfunction&);
    valfunction(valfunction &&f) : Gdat(std::move(f.Gdat)), t(std::move(f.t)), s_infix(std::move(f.s_infix)) {nvar=f.nvar;}//{Gdat=std::move(f.Gdat);t=f.t;}
    valfunction& operator= (valfunction&& f) {Gdat=std::move(f.Gdat);t=std::move(f.t);s_infix=std::move(f.s_infix);nvar=f.nvar;return *this;}
    valfunction& operator= (const valfunction&) = default;
    //valfunction& operator= (valfunction f) {Gdat=std::move(f.Gdat);t=std::move(f.t);s_infix=std::move(f.s_infix);nvar=f.nvar;return *this;}
    const valfunction& infix_to_postfix(const std::string &s);
    double operator() (const double&) const;
    double operator() (const vector<double>&) const;
    valfunction operator() (const valfunction&) const;
    valfunction operator() (const vector<valfunction>&) const;
    template <class T> T rationaleval(const T& r) const;
    template <class T> T rationaleval(const vector<T> &v) const;
    pol<rational> getpolynomial() const;
    n_polynom<rational> getn_polynom() const;
    s_polynom<rational> gets_polynom() const;
    template <class T> val::pol<T> getunivarpol(const T& a,const std::string var="x") const; // polynomial p(var) = f(var,a,a,...,a);
    rationalfunction getrationalfunction(int reduced = 1) const;
    void setparameter(const double &a) {t=a;}
    const double& getparameter() const {return t;}
    const std::string& getinfixnotation() const;
    std::string getfirstoperator() const;
    int numberofvariables() const {return nvar;}
    int isconst() const;
    int isconst(int k) const;
    int isrationalfunction() const;
    int ispolynomialfunction() const;
    int isdifferentiable() const;
    int islinearfunction() const;
    void print() const;
    void print_prefix() const;
    d_array<std::string> get_prefix() const;
    std::string get_infix() const;
    //
    int is_zero() const {return (s_infix=="" || s_infix=="0" || Gdat.isempty());}
    //
    valfunction operator +(const valfunction &g) const;
    const valfunction& operator +=(const valfunction& g) {*this = *this + g;return *this;}
    const valfunction& operator -=(const valfunction& g) {*this = *this - g;return *this;}
    valfunction operator -(const valfunction &g) const;
    valfunction operator -() const;
    valfunction operator *(const valfunction &g) const;
    valfunction operator /(const valfunction &g) const;
    valfunction operator ^(int n) const;
    //
    valfunction derive(int k=1) const;
    valfunction getfirstargument() const;
    valfunction getsecondargument() const;
    // zeros in [x1,x2]
    Glist<double> double_roots(const double &x1, const double &x2,int iterations,const double &epsilon=1e-9,int methoditerations=20) const;
    // intervals where function is undefined in [x1,x21
    Glist<GPair<double>> get_undefined_intervals(const double &x1,const double &x2,int iterations,const double &epsilon=1e-9,int methoditerations=20) const;
    static void simplifypolynomial(d_array<token> &f_t);
    void simplify(int extended = 0);
};


template <class T>
T valfunction::rationaleval(const T& x) const
{
    T zero=val::zero_element<T>(),value,v2;
    if (Gdat.isempty()) return zero;

    GlistIterator<valfunction::token> iT;
    Glist<T> G;
    int i=0,exp;


    for (iT=Gdat;iT;++i,iT++) {
        G.resetactual();
        if (iT().type==0) {
            if (iT().data!="t" && iT().data!="PI")  {
				G.inserttohead(val::FromString<T>(iT().data));
			}
        }
        else if (iT().type==1) {
			if (iT().data=="x" || iT().data=="x1") G.inserttohead(x);
			else G.inserttohead(zero);
		} //std::cout<<"  variable ";}
        else {
            value=zero;
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
            else if (iT().data=="^" && i>0) { //case "^":
                if (!G.isempty()) {
					G.skiphead();
					if (Gdat[i-1].data=="m" && i-1>0) {
						exp = -FromString<int>(Gdat[i-2].data);
					}
					else exp = FromString<int>(Gdat[i-1].data);
				}
                if (!G.isempty()) {
                    value=G.actualvalue();
                    G.skiphead();
                    //std::cout<<"\n value = "<<value;
                    //std::cout<<"\n exp = "<<exp<<std::endl;
                    value=val::power(value,exp);
                }
                G.inserttohead(value);
            }
		}
	}
    G.resetactual();
    value=zero;
    if (!G.isempty()) value=G.actualvalue();
    return value;
}



template <class T>
T valfunction::rationaleval(const vector<T>& x) const
{
    T zero=val::zero_element<T>(),value,v2;
    if (Gdat.isempty()) return zero;

    GlistIterator<valfunction::token> iT;
    Glist<T> G;
    int i=0,exp,j,k;


    for (iT=Gdat;iT;++i,iT++) {
        G.resetactual();
        if (iT().type==0) {
            if (iT().data!="t" && iT().data!="PI")  {
				G.inserttohead(val::FromString<T>(iT().data));
			}
        }
        else if (iT().type==1) {
			if (iT().data=="x" || iT().data=="x1") G.inserttohead(x(0));
			else {
				j=1;
				k=val::FromString<int>(fparser::findnumber(iT().data,j));
				G.inserttohead(x(k-1));
			}
		} //std::cout<<"  variable ";}
        else {
            value=zero;
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
            else if (iT().data=="^" && i>0) { //case "^":
                if (!G.isempty()) {
					G.skiphead();
					if (Gdat[i-1].data=="m" && i-1>0) {
						exp = -FromString<int>(Gdat[i-2].data);
					}
					else exp = FromString<int>(Gdat[i-1].data);
				}
                if (!G.isempty()) {
                    value=G.actualvalue();
                    G.skiphead();
                    //std::cout<<"\n value = "<<value;
                    //std::cout<<"\n exp = "<<exp<<std::endl;
                    value=val::power(value,exp);
                }
                G.inserttohead(value);
            }
		}
	}
    G.resetactual();
    value=zero;
    if (!G.isempty()) value=G.actualvalue();
    return value;
}

template <class T>
val::pol<T> valfunction::getunivarpol(const T& a,const std::string var) const
{
    using namespace val;
    pol<T> Y(T(1),1),value,v2;  //
    int exponent;

    if (Gdat.isempty()) return value;

    GlistIterator<valfunction::token> iT;
    Glist<pol<T>> G;


    for (iT=Gdat;iT;iT++) {
        G.resetactual();
        if (iT().type==0) {
            if (iT().data=="t") G.inserttohead(pol<T>(T(t),0));
            else G.inserttohead(pol<T>(val::FromString<T>(iT().data),0));
        }
        else if (iT().type==1) {
            if (iT().data != var) G.inserttohead(pol<T>(a,0));
            else G.inserttohead(Y);
        } //std::cout<<"  variable ";
        else {
            value.del();  // value=0
            if (iT().data=="+") {   //case "+":
                if (!G.isempty()) {value=std::move(G.actualvalue());G.skiphead();}
                if (!G.isempty()) {value+=G.actualvalue();G.skiphead();}
                G.inserttohead(value);
            }
            else if (iT().data=="-") {  // case "-":
                if (!G.isempty()) {v2=G.actualvalue();G.skiphead();}
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();value-=v2;}
                G.inserttohead(std::move(value));
            }
            else if (iT().data=="m") {
                if (!G.isempty()) G.actualvalue()*=T(-1);
            }
            else if (iT().data=="*") {  //case "*":
                if (!G.isempty()) {value=std::move(G.actualvalue());G.skiphead();}
                if (!G.isempty()) {value*=G.actualvalue();G.skiphead();}
                G.inserttohead(std::move(value));
            }
            else if (iT().data=="/") {  //case "/":
                if (!G.isempty()) {v2=std::move(G.actualvalue());G.skiphead();}
                if (!G.isempty()) {value=G.actualvalue();G.skiphead();value/=v2;}
                G.inserttohead(std::move(value));
            }
            else if (iT().data=="^") { //case "^":
                if (!G.isempty()) {v2=std::move(G.actualvalue());G.skiphead();}
                exponent=int(v2.leader());
                if (!G.isempty()) {
                    value=G.actualvalue();
                    G.skiphead();
                    value=val::power(value,exponent);
                }
                G.inserttohead(std::move(value));
            }
        }
    }
    G.resetactual();
    value.del();
    if (!G.isempty()) value=std::move(G.actualvalue());
    return value;
}



} //end namespace val


#endif // FUNCTION_PARSER_H_INCLUDED
