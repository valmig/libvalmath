#ifndef PARSER_H_INCLUDED
#define PARSER_H_INCLUDED

#include <string>
#include <val_basics.h>

namespace val
{

template <class T> class n_polynom;
template <class T> class pol;
class rational;


DLL_PUBLIC val::n_polynom<val::rational> parse_n_polynomial(const std::string &s);
DLL_PUBLIC val::pol<val::rational> parse_u_polynomial(const std::string &s);



} // end namespace val

#endif // PARSER_H_INCLUDED
