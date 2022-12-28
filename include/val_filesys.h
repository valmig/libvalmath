#ifndef VAL_FILESYS_H_INCLUDED
#define VAL_FILESYS_H_INCLUDED

#include <string>
#include <Glist.h>



namespace val
{


namespace filesys
{
    extern DLL_PUBLIC const int _ISFILE,_ISDIR,_NHIDDEN,_HIDDEN,_ALL,_ALLNHIDDEN,_ALLHIDDEN;
}


// Dir -path + filename of current exe
DLL_PUBLIC std::string GetExeFileName();
// dir path of current exe
DLL_PUBLIC std::string GetExeDir();

DLL_PUBLIC int GetExeFileName(std::string &path,std::string &name);

DLL_PUBLIC std::string GetCurrentDir();
DLL_PUBLIC std::string CurrentUser();
DLL_PUBLIC std::string CurrentHomeDir();


DLL_PUBLIC val::Glist<std::string> MakeFileList(const std::string &dirname,const int FLAG=filesys::_ALL,const std::string &wildcard="*");


DLL_PUBLIC int DirExists(const std::string&);

DLL_PUBLIC int CreateDir(const std::string&);

// Deletes empty directory
DLL_PUBLIC int DeleteDir(const std::string&);

DLL_PUBLIC int ChangeDir(const std::string&);

DLL_PUBLIC int FileExists(const std::string&);

DLL_PUBLIC int Copy_File(const std::string& oldfile,const std::string& newfile,int overwrite=0);

DLL_PUBLIC int Delete_File(const std::string&);



} // end namespace val



#endif // VAL_FILESYS_H_INCLUDED
