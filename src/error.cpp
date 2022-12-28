#include <error.h>
#include <iostream>
#include <cstdlib>

namespace val
{

void SetErrorMessage(char_function &F)
{
    Error::errormessage=&F;
}

void SetWarningMessage(char1_function &F)
{
    Error::warningmessage=&F;
}


namespace Error
{
void ErrorMessage(const char* c) {std::cout<<c;}

void error(const char* c) {errormessage(c);exit(-1);}

void (*errormessage)(const char*) = &ErrorMessage;

int WarningMessage(const char* c)
{
    int a;

    std::cout<<c;
    std::cout<<"\nInput 1 to continue or 0 to abort!";
    std::cin>>a;
    if (a==1) return 1;
    else return 0;
}

void warning(const char* c)
{
    if (warningmessage(c)) return;
    else exit(-1);
}

int (*warningmessage)(const char*) = &WarningMessage;

}

}  //End namespace val
