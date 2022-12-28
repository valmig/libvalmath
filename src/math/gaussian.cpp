
#include <gaussian.h>
#include <val_utils.h>

namespace val
{

gaussian gaussian::operator +(const gaussian &x) const
{
   gaussian y(re+x.re,im+x.im);
   return y;
}


gaussian gaussian::operator -() const
{
   gaussian y(re,im);
   y.re.changesign(); y.im.changesign();
   return y;
}


gaussian gaussian::operator -(const gaussian &x) const
{
   gaussian y(re-x.re,im-x.im);
   return y;
}


gaussian gaussian::operator *(const gaussian &x) const
{
  gaussian y;
  y.re= re*x.re - im*x.im;
  y.im= re*x.im + im*x.re;
  return y;
}


gaussian gaussian::operator /(const gaussian &x) const
{
  gaussian y;
  y.re= (re*x.re + im*x.im)/(x.re*x.re + x.im*x.im);
  y.im= (im*x.re - re*x.im)/(x.re*x.re + x.im*x.im);
  return y;
}


const gaussian& gaussian::operator =(const gaussian &x)
{
  re = x.re;
  im = x.im;
  return *this;
}


const gaussian& gaussian::operator =(gaussian &&x)
{
  re = std::move(x.re);
  im = std::move(x.im);
  return *this;
}


const gaussian& gaussian::operator +=(const gaussian &x)
{
   *this = *this +x;
   return *this;
}

const gaussian& gaussian::operator -=(const gaussian &x)
{
   *this = *this - x;
   return *this;
}

const gaussian& gaussian::operator *=(const gaussian &x)
{
   *this = *this * x;
   return *this;
}

const gaussian& gaussian::operator /=(const gaussian &x)
{
   *this = *this/x;
   return *this;
}


gaussian gaussian::conjugate() const
{
    gaussian y(re,im);
    y.im.changesign();
    return y;
}


const gaussian& gaussian::operator =(const rational &x)
{
   re = x;
   im = rational(0);
   return *this;
}


const gaussian& gaussian::operator=(rational &&x)
{
    re = std::move(x);
    im = rational(0);
    return *this;
}


int gaussian::operator ==(const gaussian &x) const
{
  return ((re==x.re) && (im==x.im));
}


int gaussian::operator !=(const gaussian &x) const
{
  return ((re!=x.re) || (im!=x.im));
}


int gaussian::operator ==(const rational &x) const
{
  return ((im.iszero()) && (re==x));
}

int gaussian::operator !=(const rational &x) const
{
  return ((!im.iszero()) || (re!=x));
}


gaussian operator *(const rational &x,const gaussian &z)
{
   gaussian y(x*z.re,x*z.im);
   //y.re= x*z.re;
   //y.im= x*z.im;
   return y;
}




std::ostream& operator <<(std::ostream& os,const gaussian &x)
{

  if (x.re.iszero() && x.im.iszero()) {
        os<<'0';
        return os;
  }
  if (!x.re.iszero()) os<<x.re;
  if (!x.im.iszero()) {
    if (x.im.signum()>0 && !x.re.iszero()) os<<'+';
    os<<x.im<<'i';
  }
  return os;
}


// Vorlaeufig:
std::istream& operator >>(std::istream& is,gaussian& x)
{

  x=gaussian();
  char *buf,*hbuf;
  int i,l,k=1,wert=0,sign=1;
  std::string is_re,is_im;
  //std::istream is_re,is_im;
  //integer ten(10);

  buf = new char[1000];



  l=0;
  //is.seekg(1l);
  while ((wert!=10) && (wert!=32) && is) {
	wert=is.get();
	//cout<<wert;
	if ((wert==10) || (wert==32)){
		if (l!=0) buf[l]='\0';
		else wert=0;
	}
	else buf[l++]=wert;
	if (l>=(k*1000)) {
        k++;
        hbuf= new char[k*1000];
        for (i=0;i<((k-1)*1000);i++) hbuf[i]=buf[i];
        delete[] buf;
        buf=hbuf;
        hbuf=NULL;
	}

	if (!is) l--;
  }

  //suche i: nach i werden alle Zeichen ignoriert:
  for (i=0;i<l;i++) {
       if (buf[i]=='i') break;
  }
  if (i==l) {// Imaginär-Teil = 0
      is_re = std::string(buf,i);
      x.re=FromString<rational>(is_re);
      //x.re=gaussian::char_to_double(buf,i);
      delete[] buf;
      //std::cout<<"\n is (re) ready: "<<x<<std::endl;;
      return is;
  }
  // Nun ist Imginär-Teil !=0:
  l=i;
  // suche nach Trennendem + oder -
  i=0;
  while (i<l && (buf[i]=='+' || buf[i]=='-')) {
        if (buf[i]=='-') sign*=-1;
        i++; //+ und Minus werden zunächst ignoriert.
  }

  if (i==l) {
      x.im=rational(sign*1);
      delete[] buf;
      //std::cout<<"\n is (im) ready";
      return is;
  }
  //Nun ist i<l und buf[i]!=+ -
  for (;i<l;i++) {
       if (buf[i]=='+' || buf[i]=='-') {
           if (buf[i-1]!='e') break;
       }
  }
  if (i==l) {  // Zahl ist rein -imaginär
      is_im = std::string(buf,i);
      x.im = FromString<rational>(is_im);
      //x.im=gaussian::char_to_double(buf,i);
      delete[] buf;
      return is;
  }
  //x.zaehler=integer::char_to_integer(buf,i);
  //Zahl hat im und re -Teil != 0
  is_re = std::string(buf,i);
  x.re = FromString<rational>(is_re);
  //x.re=gaussian::char_to_double(buf,i);
  if (i<l-1) {
    hbuf=&buf[i+1];
    //x.nenner=integer::char_to_integer(hbuf,l-1-i);
    is_im = std::string(hbuf,l-i-1);
    x.im = FromString<rational>(is_im);
    //x.im=gaussian::char_to_double(hbuf,l-1-i);
    if (buf[i]=='-') x.im.changesign();
    if (x.im.iszero()) {
        sign=1;
        x.im=rational(1);
        for (;i<l;i++) {
            if (buf[i]=='-') sign*=-1;
            else {
                x.im=rational(0);
                break;
            }
        }
        if (sign<0) x.im.changesign();
    }
    hbuf=NULL;
  }
  else {
    if (buf[i]=='-') x.im =rational(-1);
    else x.im=rational(1);
  }

  //if (x.nenner==0) hilfratiofktn::errornenner();

  //if (rational::gekuerzt) x.kuerze();

  delete[] buf;
  return is;
}








} // end namespace val
