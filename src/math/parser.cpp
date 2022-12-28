#include <parser.h>
#include <rational.h>
#include <Glist.h>
#include <n_polynom.h>
#include <pol.h>
#include <val_utils.h>



namespace parser
{


template <class T>
val::n_polynom<T> power(const val::n_polynom<T> &a,int n)
{
    using namespace val;

    n_polynom<T> y,x(T(1));

    if (n<=0) return x;

    if (a.iszero()) return y;

    y =a;

    while (n!=0) {
       if (n%2!=0) {
           x*=y;
       }
       y*=y;
       n=n/2;
    }
    return x;
}


struct token {
    char op;
    val::n_polynom<val::rational> p;
    int isoperator;           //token = operator => oprator =1 , sonst 0;
    token (char o,const val::n_polynom<val::rational> &f,int is) :op(o),p(f),isoperator(is) {}
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



val::n_polynom<val::rational> pop(val::Glist<val::n_polynom<val::rational> > &G)
{
    val::n_polynom<val::rational> f;
    if (G.isempty()) return f;

    if (!G.actualvalid()) return f;

    f = std::move(G.actualvalue());
    G.skiphead();
    G.resetactual();
    return f;
}


// push f to G;
void push(val::n_polynom<val::rational> &f,val::Glist<val::n_polynom<val::rational> > &G)
{
    G.inserttohead(std::move(f));
    G.resetactual();
}

val::n_polynom<val::rational> parse_from_polish_notation(val::Glist<token> &T)
{
    using namespace val;

    val::n_polynom<val::rational> f,h;
    GlistManipulator<token> iT;
    Glist<n_polynom<rational> > G;


    for (iT=T;iT;iT++) {
        if (!iT().isoperator) {
            push(iT().p,G);
            //std::cout<<" \n G.length() = "<<G.length();
        }
        else {
            switch (iT().op) {
            case '+':
                f = pop(G) + pop(G);
                push(f,G);
                break;
            case '-':
                h=pop(G);f=pop(G);
                f-=h;
                //std::cout<<"\n f = "<<f;
                push(f,G);
                break;
            case 'm':
                G.resetactual();
                if (!G.isempty()) G.actualvalue()*=rational(-1);
                break;
            case '*':
                f = pop(G) * pop(G);
                push(f,G);
                break;
            case '/':
                h=pop(G);f=pop(G);
                f.edivby(h.LC());
                push(f,G);
                break;
            case '^':
                h=pop(G);f=pop(G);
                f = power(f,int(nominator(h.LC())));
                push(f,G);
                break;
            default: break;
            }
            //std::cout<<" \n G.length() = "<<G.length();
        }
    }

    f.del();
    if (!G.isempty()) f= std::move(G.actualvalue());
    return f;
}




int is_number(const val::n_polynom<val::rational> &f)
{
    if (f.iszero()) return 1;
    if (f.length() > 1) return 0;
    if (f.LT().iszero()) return 1;
    else return 0;
}


void set_token_list(const std::string & s,val::Glist<token> &T)
{
    using namespace val;
    int i,n=s.size(),k;
    n_polynom<rational> f;
    n_expo X;

    if (!T.isempty()) T.dellist();

    for (i=0;i<n;i++) {
        if (s[i]==' ' || s[i]=='\n') continue;
        else if (s[i]=='+') {
            T.inserttoend(token('+',f,1));
        }
        else if (((s[i]=='-' && (i==n-1)) || (s[i]=='-' && s[i+1]==' '))) {
            T.inserttoend(token('-',f,1));
        }
        else if (s[i]=='m') {
            T.inserttoend(token('m',f,1));
        }
        else if (s[i]=='*') {
            T.inserttoend(token('*',f,1));
        }
        else if (s[i]=='/') {
            T.inserttoend(token('/',f,1));
        }
        else if (s[i]=='^') {
            T.inserttoend(token('^',f,1));
        }
        else {
            if (s[i]=='x') {
                if (i==n-1) X=n_expo(1,1);
                else {
                    k=int(s[i+1]-48);
                    X=n_expo(0,k);
                    X(k-1) =1;
                    i++;
                }
                f=n_polynom<rational>(rational(1),X);
                T.inserttoend(token(' ',f,0));
                f.del();
            }
            else {
                std::string s1;
                while ((i<n) && (s[i] !=' ' || s[i]=='\n')) {
                    s1+=s[i];
                    i++;
                }
                f=n_polynom<rational>(FromString<rational>(s1),n_expo());
                T.inserttoend(token(' ',f,0));
                f.del();
            }
        }

    }
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
        else if (s[i]=='e') {
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

std::string findvariable(const std::string &s, int &i)
{
    int n=s.size();
    std::string out;

    if (s[i]=='y') {i++;return "x2";}
    if (s[i]=='z') {i++;return "x3";}
    if (s[i]=='x') {
        if (i==n-1) {i++;return "x1";}
        else if (s[i+1] >= '1' && s[i+1] <='9') {
            out="x";
            out+=s[i+1];
            i+=2;
        }
        else {out="x1";i++;}
    }
    else i++;
    return out;
}


// Anwendung des shunting - yard - Algortihmus:
std::string infix_to_postfix(const std::string &s)
{
    using namespace val;
    int i=0,n=s.size(),nG=0;
    std::string out;
    Glist<s_stack> G,OP;
    GlistIterator<s_stack> IG;
    s_stack t;

    if (!checkbrackets(s)) return out;  //Fehler bei Klammerung;

    // Setze G in infix -Notation;

    while (i<n) {
        if (s[i]=='+') {//|| s[i]=='-') {
            //out="";
            //out+=s[i];
            //std::cout<<" \n out =   "<<out;
            G.inserttoend(s_stack("+",s_stack::OPERATOR,2));
            i++;
            nG++;
        }
        else if (s[i]=='-') {
            if (G.isempty()) G.inserttoend(s_stack("m",s_stack::OPERATOR,2));
            else if (G[nG-1].data==")" || G[nG-1].type==0 || G[nG-1].type==1) G.inserttoend(s_stack("-",s_stack::OPERATOR,2));
            else G.inserttoend(s_stack("m",s_stack::OPERATOR,3));
            nG++;
            i++;
        }
        else if (s[i]=='*' || s[i]=='/') {
            out="";
            out+=s[i];
            G.inserttoend(s_stack(out,s_stack::OPERATOR,3));
            nG++;
            i++;
        }
        else if (s[i]=='^') {
            G.inserttoend(s_stack("^",s_stack::OPERATOR,4,0));
            nG++;
            i++;
        }
        else if (s[i]=='(') {
            G.inserttoend(s_stack("(",s_stack::LBRACKET));
            nG++;
            i++;
        }
        else if (s[i]==')') {
            G.inserttoend(s_stack(")",s_stack::RBRACKET));
            nG++;
            i++;
        }
        else if (s[i]>='0' && s[i] <='9') {
            G.inserttoend(s_stack(findnumber(s,i),s_stack::NUMBER));
            nG++;
        }
        else if (s[i]=='x' || s[i]=='y' || s[i]=='z') {
            G.inserttoend(s_stack(findvariable(s,i),s_stack::VARIABLE));
            nG++;
        }
        else i++;
    }

    out="";
    //std::cout<<std::endl;
    //for (IG=G;IG;IG++) std::cout<<IG().data<<" ";

    // Nun hunting-yard-algorithmus:

    G.resetactual();
    while (!G.isempty()) {
         t = G.actualvalue(); G.skiphead();G.resetactual();
         if (t.type==s_stack::NUMBER || t.type==s_stack::VARIABLE) out+= t.data + " ";
         else if (t.type==s_stack::OPERATOR) {
            while (!OP.isempty() && OP.actualvalue().type==s_stack::OPERATOR) {
                if (OP.actualvalue().precedence>t.precedence  || (OP.actualvalue().precedence==t.precedence && OP.actualvalue().leftassociativ)) {
                    out+=OP.actualvalue().data + " ";
                    OP.skiphead(); OP.resetactual();
                }
                else break;
            }
            OP.inserttohead(t); OP.resetactual();
         }
         else if (t.type==s_stack::LBRACKET) {OP.inserttohead(t);OP.resetactual();}
         else if (t.type==s_stack::RBRACKET) {
            while (!OP.isempty() && OP.actualvalue().type!=s_stack::LBRACKET) {
                out+= OP.actualvalue().data + " "; OP.skiphead();OP.resetactual();
            }
            if (!OP.isempty()) {OP.skiphead();OP.resetactual();}
         }
    }

    while (!OP.isempty()) {
        out+=OP.actualvalue().data + " ";
        OP.skiphead();OP.resetactual();
    }

    return out;
}

} // end namespace parser



namespace val
{
val::n_polynom<val::rational> parse_n_polynomial(const std::string &s)
{
    std::string t = parser::infix_to_postfix(s);

    Glist<parser::token> T;

    parser::set_token_list(t,T);

    return parser::parse_from_polish_notation(T);
}


val::pol<val::rational> parse_u_polynomial(const std::string &s)
{
    n_polynom<rational> g;
    n_polynomIterator<rational> pg;
    pol<rational> f;

    g=val::parse_n_polynomial(s);
    if (g.iszero()) return f;

    for (pg=g;pg;pg++) {
        f.insert(pg.actualcoef(),pg.actualterm()(0));
    }
    return f;
}


} //end namespace val




