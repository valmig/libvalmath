#ifndef VAL_UTILS_CPP_INCLUDED
#define VAL_UTILS_CPP_INCLUDED

#include <val_utils.h>
#include <cstdlib>

#ifdef __WIN32
#include <Windows.h>
#include <cstdio>
#include <cstring>
#endif

//#include <rational.h>
//#include <complex.h>
//#include <sstream>

namespace val
{

/*
std::string ToString(const int& a)
{
   std::stringstream ss;//create a stringstream
   ss << a;//add number to the stream
   return ss.str();
}


std::string ToString(const double& a)
{
   std::stringstream ss;//create a stringstream
   ss <<a;//add number to the stream
   return ss.str();
}
*/

char* StringToChar(const std::string& s)
{
   int n=0,i;
   char *c=NULL;


    n=s.length();
    c = new char[n+1];
    for (i=0;i<n;i++) c[i]=s[i];
    if (n>0) c[n]='\0';
    return c;
}

#ifdef __WIN32
int system(const std::string &command)
{
    char            *tmp_command, *cmd_exe_path;
    int         ret_val;
    size_t          len;
    PROCESS_INFORMATION process_info = {0};
    STARTUPINFOA        startup_info = {0};


    len = strlen(command.c_str());
    tmp_command = malloc(len + 4);
    tmp_command[0] = 0x2F; // '/'
    tmp_command[1] = 0x63; // 'c'
    tmp_command[2] = 0x20; // <space>;
    memcpy(tmp_command + 3, command.c_str(), len + 1);

    startup_info.cb = sizeof(STARTUPINFOA);
    cmd_exe_path = getenv("COMSPEC");
    _flushall();  // required for Windows system() calls, probably a good idea here too

    if (CreateProcessA(cmd_exe_path, tmp_command, NULL, NULL, 0, CREATE_NO_WINDOW, NULL, NULL, &startup_info, &process_info)) {
        WaitForSingleObject(process_info.hProcess, INFINITE);
        GetExitCodeProcess(process_info.hProcess, &ret_val);
        CloseHandle(process_info.hProcess);
        CloseHandle(process_info.hThread);
    }

    free((void *) tmp_command);

    return(ret_val);
}

#else
int system(const std::string &command, int silent)
{
    std::string c = command, quiet = " > /dev/null";
    if (silent) c += quiet;
    return std::system(c.c_str());
}
#endif

/*
int isinteger(const std::string &s)
{
    int i;
    if (s=="") return 0;
    for (i=0;s[i]!='\0';i++) {
        if (s[i]=='-' && i==0 ) continue;
        else if (s[i]<48 || s[i] > 57) return 0;
    }
    return 1;
}
*/

/*
#ifdef RATION_H

rational char_to_rational(const char* s)
{
   std::stringstream ss;
   rational x;

   ss.str(s);
   ss>>x;
   return x;
}

rational string_to_rational(const std::string& s)
{
   std::stringstream ss;
   rational x;

   ss.str(s);
   ss>>x;
   return x;
}

#endif

#ifdef COMPLEX_H_INCLUDED

val::complex string_to_complex(const std::string& s)
{
   std::stringstream ss;
   val::complex x;

   ss.str(s);
   ss>>x;
   return x;
}

#endif
*/


} // End namespace val
#endif // VAL_UTILS_CPP_INCLUDED
