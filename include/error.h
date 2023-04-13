#ifndef ERROR_H_INCLUDED
#define ERROR_H_INCLUDED

#include <val_basics.h>
#include <val_utils.h>

namespace val
{
typedef void (char_function) (const char*);
typedef int (char1_function) (const char*);

DLL_PUBLIC void SetErrorMessage(char_function &);
DLL_PUBLIC void SetWarningMessage(char1_function &);

namespace Error
{
DLL_PUBLIC void ErrorMessage(const char*);
DLL_PUBLIC void error(const char*);
extern DLL_PUBLIC void (*errormessage)(const char*);
DLL_PUBLIC int WarningMessage(const char*);
DLL_PUBLIC void warning(const char*);
extern DLL_PUBLIC int (*warningmessage)(const char*);

template <class T>
void error(const std::string &message,const T& t)
{
    std::string msg = "\n" + val::gettypename(t) + message;
    error(msg.c_str());
}

} // end namespace Error


} // end namespace val

#endif // ERROR_H_INCLUDED
