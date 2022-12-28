
#include <val_filesys.h>

#ifdef _WIN32

#include <windows.h>
#include <Lmcons.h>
#include <Shlobj.h>


#else    // Linux:

#include <unistd.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <dirent.h>
#include <fstream>
#include <pwd.h>

#endif // _WIN32

namespace val
{


namespace filesys
{
    const int _ISFILE=1;
    const int _ISDIR=2;
    const int _NHIDDEN=4;
    const int _HIDDEN=8;
    const int _ALL=15;
    const int _ALLNHIDDEN=7;
    const int _ALLHIDDEN=11;
}




// check if wildcard coincides with the beginning of s.
int is_at_top(const std::string &s,const std::string &wildcard);

// check if wildcard coincides with the end of s,
int is_at_bottom(const std::string &s,const std::string &wildcard);

// check if wildcard coincides anywhere inside s.
int is_anywhere(const std::string &s,const std::string &wildcard);

int wildcardmatch(const std::string &s,const std::string &wildcard);



//

int is_at_top(const std::string &s,const std::string &wildcard)
{
    if (wildcard=="") return 1;
    int lw=wildcard.length(),ls=s.length();
    if (lw>ls) return 0;
    for (int i=0;i<lw;i++) if (s[i]!=wildcard[i]) return 0;
    return 1;
}

int is_at_bottom(const std::string &s,const std::string &wildcard)
{
    if (wildcard=="") return 1;
    int lw=wildcard.length(),ls=s.length();
    if (lw>ls) return 0;
    for (int i=0;i<lw;i++) if (s[ls-i-1]!=wildcard[lw-i-1]) return 0;
    return 1;
}

int is_anywhere(const std::string &s,const std::string &wildcard)
{
    if (wildcard=="") return 1;
    int i=s.find(wildcard);

    if (i==-1) return 0;
    else return 1;
}


int wildcardmatch(const std::string &s,const std::string &wildcard)
{
    if (wildcard=="") return s==wildcard;
    if (wildcard=="*") return 1;

    int l=wildcard.length();
    std::string wcard;

    if (l==1) return wildcard==s;

    if (wildcard[0]=='*' && wildcard[l-1]=='*' && l>2) {
        wcard=std::string(wildcard,1,l-2);
        return is_anywhere(s,wcard);
    }
    else if (wildcard[0]=='*') {
        wcard=std::string(wildcard,1,l-1);
        return is_at_bottom(s,wcard);
    }
    else if (wildcard[l-1]=='*') {
        wcard=std::string(wildcard,0,l-2);
        return is_at_top(s,wcard);
    }
    else return wildcard==s;
}



#ifdef _WIN32

std::string CurrentUser()
{
    char username[UNLEN+1];
    DWORD username_len = UNLEN+1;
    GetUserName(username, &username_len);
    return std::string(username);
}

std::string CurrentHomeDir()
{
    char path [MAX_PATH +1];
    if (SHGetSpecialFolderPathA(HWND_DESKTOP, path,CSIDL_PROFILE, FALSE)) return std::string(path);
    else return "";
}

std::string GetExeFileName()
{
  char buffer[MAX_PATH];
  GetModuleFileName( NULL, buffer, MAX_PATH );
  return std::string(buffer);
}

const char dirchar='\\';


std::string GetCurrentDir()
{
    char NPath[MAX_PATH];
    GetCurrentDirectory(MAX_PATH, NPath);
    return std::string(NPath);
}

val::Glist<std::string> MakeFileList(const std::string &dirname,const int FLAG,const std::string &wildcard)
{
    val::Glist<std::string> G;
    std::string wcard="",filename,dir,fullname;
    HANDLE hFind;
    WIN32_FIND_DATA FindData;
    DWORD FileAt;

    dir=dirname+"\\*.*";
    hFind = FindFirstFile(dir.c_str(), &FindData);
    if (hFind == INVALID_HANDLE_VALUE) {
        //std::cout<<"\nCannot find directory "<<dirname;
        return G;
    }
    filename=std::string(FindData.cFileName);
    if (wildcardmatch(filename,wildcard)) {
        fullname=dirname+dirchar+filename;
        FileAt=GetFileAttributes(fullname.c_str());
        if ((FileAt & FILE_ATTRIBUTE_HIDDEN) && (FLAG & filesys::_NHIDDEN)) ;
        if (FLAG==filesys::_ALL) G.sinsert(filename);
        else if ((FLAG & filesys::_HIDDEN) && (!(FileAt & FILE_ATTRIBUTE_HIDDEN))) ;
        else if ((FileAt & FILE_ATTRIBUTE_DIRECTORY) && (FLAG & filesys::_ISDIR)) G.sinsert(filename);
        else if ((FileAt & FILE_ATTRIBUTE_ARCHIVE) && (FLAG & filesys::_ISFILE)) G.sinsert(filename);
    }
    while (FindNextFile(hFind, &FindData)) {
        filename=std::string(FindData.cFileName);
        if (!wildcardmatch(filename,wildcard)) continue;
        fullname=dirname+dirchar+filename;
        FileAt=GetFileAttributes(fullname.c_str());
        //std::cout<<"\n"<<filename<<"  "<<(FileAt & FILE_ATTRIBUTE_HIDDEN)<<"  "<<(FLAG & filesys::_ISHIDE);
        if ((FileAt & FILE_ATTRIBUTE_HIDDEN) && (FLAG & filesys::_NHIDDEN)) continue;
        if (FLAG==filesys::_ALL) G.sinsert(filename);
        else if ((FLAG & filesys::_HIDDEN) && (!(FileAt & FILE_ATTRIBUTE_HIDDEN))) continue;
        else if ((FileAt & FILE_ATTRIBUTE_DIRECTORY) && (FLAG & filesys::_ISDIR)) G.sinsert(filename);
        else if ((FileAt & FILE_ATTRIBUTE_ARCHIVE) && (FLAG & filesys::_ISFILE)) G.sinsert(filename);
    }

    return G;
}


int DirExists(const std::string& dirname)
{
    DWORD ftyp = GetFileAttributesA(dirname.c_str());
    if (ftyp == INVALID_FILE_ATTRIBUTES)
        return 0;  //something is wrong with your path!

    if (ftyp & FILE_ATTRIBUTE_DIRECTORY) return 1;   // this is a directory!

    return 0;    // this is not a directory!
}

int CreateDir(const std::string& dirname)
{
    return CreateDirectory(dirname.c_str(),NULL);
}

int ChangeDir(const std::string& dirname)
{
    return SetCurrentDirectory(dirname.c_str());
}

int DeleteDir(const std::string& dirname)
{
    return RemoveDirectory(dirname.c_str());
}

int FileExists(const std::string& dirname)
{
    DWORD ftyp = GetFileAttributesA(dirname.c_str());
    if (ftyp == INVALID_FILE_ATTRIBUTES)
        return 0;  //something is wrong with your path!

    if (ftyp & FILE_ATTRIBUTE_DIRECTORY) return 0;   // this is a directory!

    return 1;    // this is not a directory!
}

int Copy_File(const std::string& oldfile,const std::string& newfile,int overwrite)
{
    return ::CopyFile(oldfile.c_str(),newfile.c_str(),!overwrite);
}

int Delete_File(const std::string& filename)
{
    return ::DeleteFile(filename.c_str());
}



#else    // Linux:


std::string GetExeFileName()
{
	char result[ PATH_MAX ];
	ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
	//found = string(result).find_last_of("/");
	//return(string(result).substr(0,found) + "/");
	return std::string( result, (count > 0) ? count : 0 );
}

const char dirchar='/';



std::string GetCurrentDir()
{
    char *buffer=getcwd( 0, 0 ) ;
    std::string dir(buffer);
    free(buffer);
    return std::string(dir);
}


std::string CurrentUser()
{
    struct passwd *pwd = getpwuid(getuid());
    if (pwd)
        return pwd->pw_name;
    else
        return "";
}

std::string CurrentHomeDir()
{
    struct passwd *pwd = getpwuid(getuid());
    if (pwd)
        return pwd->pw_dir;
    else
        return "";
}


val::Glist<std::string> MakeFileList(const std::string &dirname,const int FLAG,const std::string &wildcard)
{
    DIR *workdir=NULL;
    struct dirent *workdirent=NULL;
    val::Glist<std::string> G;
    std::string wcard="",filename;


    workdir = opendir(dirname.c_str());
    if (workdir==NULL) {
        //std::cout<<"\nCannot find directory "<<dirname;
        return G;
    }
    while((workdirent=readdir(workdir)) !=NULL )  {
            filename=std::string(workdirent->d_name);
            if (!wildcardmatch(filename,wildcard)) continue;
            if (FLAG==filesys::_ALL) G.sinsert(filename);
            else if (filename!="" && filename[0]== '.' && (FLAG & filesys::_NHIDDEN)) continue;
            else if ((FLAG & filesys::_HIDDEN) && filename[0]!='.') continue;
            else if (workdirent->d_type==DT_REG && (FLAG & filesys::_ISFILE)) G.sinsert(filename);
            else if (workdirent->d_type==DT_DIR && (FLAG & filesys::_ISDIR)) G.sinsert(filename);
    }
    closedir(workdir);

    //for (G.resetactual();G.actualvalid();G.moveactual()) std::cout<<G.actualvalue()<<std::endl;
    return G;
}


int DirExists(const std::string& dirname)
{
    struct stat statbuf;

    if (stat(dirname.c_str(), &statbuf) != -1) {
        if (S_ISDIR(statbuf.st_mode)) return 1;
        else return 0;
    }
    else return 0;
}

int CreateDir(const std::string& dirname)
{
    int i= mkdir(dirname.c_str(),0755);
    if (i==0) return 1;
    else return 0;
}

int DeleteDir(const std::string &dirname)
{
    int i= rmdir(dirname.c_str());
    if (i==0) return 1;
    else return 0;
}

int ChangeDir(const std::string& dirname)
{
    int i= chdir(dirname.c_str());
    if (i==0) return 1;
    else return 0;
}

int FileExists(const std::string& filename)
{
    struct stat statbuf;

    if (stat(filename.c_str(), &statbuf) != -1) {
        if (S_ISREG(statbuf.st_mode)) return 1;
        else return 0;
    }
    else return 0;
}

int Copy_File(const std::string& oldfile,const std::string& newfile,int overwrite)
{
    if (!overwrite && val::FileExists(newfile)) return 0;
    std::ifstream source(oldfile,std::ios::binary);
    if (!source) return 0;
    std::ofstream target(newfile,std::ios::binary);
    if (!target) return 0;
    target<< source.rdbuf();
    return 1;
}

int Delete_File(const std::string& filename)
{
    if (!unlink(filename.c_str())) return 1;
    else return 0;
}



#endif // _Win32



std::string GetExeDir()
{
    std::string filename=GetExeFileName();
    return filename.substr(0,filename.find_last_of(dirchar));
}


int GetExeFileName(std::string &path,std::string &name)
{
    std::string filename=GetExeFileName();
    int pos;
    if (filename=="") return 0;
    path=filename.substr(0,pos=filename.find_last_of(dirchar));
    name=filename.substr(pos+1);
    return 1;
}






}  // end namespace val;
