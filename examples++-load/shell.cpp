// Example C++ function "myfunction", dynamically loaded into "load.edp"
// ---------------------------------------------------------------------
// $Id$
#include <ff++.hpp>
#include <unistd.h>
#include <sys/stat.h>
#include <dirent.h> 
#include <strings.h>
#ifdef _WIN32
#include <windows.h>
#endif
#ifdef _WIN32
const char sepdir='\\';
#else
const char sepdir='/';
#endif
extern const char *  prognamearg;
extern const char *  edpfilenamearg;

//#include <tr1/regex>
using namespace Fem2D;
long ff_chdir(string * c) { return chdir(c->c_str());}
long ff_rmdir(string * c) { return rmdir(c->c_str());}
long ff_unlink(string * c) { return unlink(c->c_str());}
#ifndef _WIN32
long ff_mkdir(string * c,long mm) {mode_t m=mm;cout << " mkdir " << *c << "mode =" << m << endl;  return mkdir(c->c_str(),m);}
#endif
long ff_chmod(string * c,long mm) {mode_t m=mm;cout << " mkdir " << *c << "mode =" << m << endl;  return chmod(c->c_str(),m);}
long ff_mkdir(string * c) {mode_t m=07777;
#ifdef _WIN32
  return mkdir(c->c_str());
#else
  return mkdir(c->c_str(),m);
#endif
  }

long ff_stat(string * c) {
struct stat buff;
return stat(c->c_str(),&buff);}
long   ff_isdir(string * c) {
struct stat buff;
 if( 0==  stat(c->c_str(),&buff))
   return  ( buff.st_mode  & S_IFDIR) ? 1 : 0 ;
 else return  -1; // err 
}

DIR ** OpenDir(DIR **pp,string *n)
{
  *pp=opendir(n->c_str());
  return pp; 
}
string * ReadDir(Stack s,DIR ** const &dirpp){
  if(*dirpp) 
    {
      struct dirent * dp =   readdir(*dirpp);
      if(dp)
	return Add2StackOfPtr2Free(s,new string(dp->d_name));
    }
  return Add2StackOfPtr2Free(s,new string(""));
} 
 inline AnyType CloseDir(Stack,const AnyType &x){
  DIR ** dirpp=PGetAny<DIR*>(x);
  if(*dirpp) (void)closedir(*dirpp);
  *dirpp=0;
  return  Nothing;
}
#ifdef _WIN32 
string * ffgetenv(Stack s,string * const & k)
{
    const int LEN = 4096;
    char envv[LEN];

     GetEnvironmentVariable(k->c_str(),envv,LEN);
    
    return Add2StackOfPtr2Free(s,new string(envv));
}

long  ffsetenv(string * const & k,string * const & v)
{
 
    char * vv = strcpy((char*)malloc(v->size()+2),v->c_str());
    char * kk = strcpy((char*)malloc(k->size()+2),k->c_str());
    SetEnvironmentVariable(vv,kk);
    return SetEnvironmentVariable(vv,kk) ;
}
long  ffunsetenv(string * const & k)
{
    SetEnvironmentVariable(k->c_str(),0);
    return 0 ;
}
#else
string * ffgetenv(Stack s,string * const & k)
{
    const char *env = getenv(k->c_str());
    if(!env) env ="";
    return Add2StackOfPtr2Free(s,new string(env));
}

long  ffsetenv(string * const & k,string * const & v)
{
    char * vv = strcpy((char*)malloc(v->size()+2),v->c_str());
    char * kk = strcpy((char*)malloc(k->size()+2),k->c_str());
    long r= setenv(vv,kk,1);
    return r ;
}
long  ffunsetenv(string * const & k)
{
    long r= unsetenv(k->c_str());
    return r ;
}
#endif
string  dirname(const string * ppath)
{
    const string & path= *ppath;
    int i,l=path.length();
    for(i=l-1;i>=0;i--)
        if(path[i]==sepdir) break;
    if(i==0) return ".";
    else if(i==1) return "/";
    else return path.substr(0,i-1);
}
string  basename(const string * ppath)
{
    const string & path= *ppath;
    int i,l=path.length();
    for(i=l-1;i>=0;i--)
        if(path[i]==sepdir) {i++;break;}
    if(i<0) i=0;
    return path.substr(i);
}


string * ff_dirname(Stack s,string * const &path)
{
    return Add2StackOfPtr2Free(s,new string(dirname(path)));
}
string * ff_basename(Stack s,string * const & path)
{
    return Add2StackOfPtr2Free(s,new string(basename(path)));
}
long copyfile(string * const & filecp, string * const & target)
{
    int tagetisdir = ff_isdir(target);
    string filein = *filecp;
    string filetarget = *target;
    if( verbosity>9)
        cout << "  cpfile :"<< filein << "-> " <<filetarget << " "<< tagetisdir << endl;
    
    if(tagetisdir ==1 )
    {
        int i,l=filein.length();
        for(i=l-1;i>=0;i--)
            if(filein[i]==sepdir) break;
        if(i<0) i=0;
        // cout << filein << " " << i << " " << l << endl;
        filetarget +=sepdir;
         filetarget +=filein.substr(i);
    }
    
    FILE* source = fopen(filein.c_str(), "rb");
    FILE* dest = fopen(filetarget.c_str(), "wb");
    if( verbosity>1)
        cout << "  cpfile :"<< filein << "-> " <<filetarget << endl;
    if(!source && !dest)
    {
        cout << " erreur copy file form " << endl;
        cout << " file in    : " << filein    << " " <<source << endl;
        cout << " file taget : " << filetarget << " " <<dest << endl;
        ffassert(0);
        return -1; // erreur
        
    }
    char buf[BUFSIZ];
    size_t size;

    while ((size = fread(buf, 1, BUFSIZ, source))) {
        fwrite(buf, 1, size, dest);
    }
    
    fclose(source);
    fclose(dest);
    return 0;
    
    
}



extern  mylex *zzzfff;

static void init(){
  Dcl_Type< DIR **  > (0,CloseDir,0);
  zzzfff->Add("Directory",atype<DIR ** >());
  TheOperators->Add("<-",   new OneOperator2<DIR **,DIR **,string *  >(OpenDir));
  Global.Add("readdir","(",new OneOperator1s_<string*,DIR**>(ReadDir));
  Global.New("modeRWXu long",CConstant<bool>(0700));
  Global.New("modeRWXg long",CConstant<bool>(070));
  Global.New("modeRWXo long",CConstant<bool>(07));

  Global.New("modeRu long",CConstant<bool>(0400));
  Global.New("modeRg long",CConstant<bool>(040));
  Global.New("modeRo long",CConstant<bool>(04));

  Global.New("modeWu long",CConstant<bool>(0200));
  Global.New("modeWg long",CConstant<bool>(020));
  Global.New("modeWo long",CConstant<bool>(02));

  Global.New("modeXu long",CConstant<bool>(0100));
  Global.New("modeXg long",CConstant<bool>(010));
  Global.New("modeXo long",CConstant<bool>(01));

  Global.Add("unlink","(",new OneOperator1<long,string*>(ff_unlink));
  Global.Add("rmdir","(",new OneOperator1<long,string*>(ff_rmdir));
  Global.Add("cddir","(",new OneOperator1<long,string*>(ff_chdir));
    Global.Add("chdir","(",new OneOperator1<long,string*>(ff_chdir));
    Global.Add("basename","(",new OneOperator1s_<string*,string*>(ff_basename));
  Global.Add("dirname","(",new OneOperator1s_<string*,string*>(ff_dirname));
  #ifndef _WIN32
  Global.Add("mkdir","(",new OneOperator2<long,string*,long>(ff_mkdir));
  #endif
  Global.Add("chmod","(",new OneOperator2<long,string*,long>(ff_chmod));
  Global.Add("mkdir","(",new OneOperator1<long,string*>(ff_mkdir));
  Global.Add("cpfile","(",new OneOperator2_<long,string*,string*>(copyfile));
  Global.Add("stat","(",new OneOperator1<long,string*>(ff_stat));
  Global.Add("isdir","(",new OneOperator1<long,string*>(ff_isdir));
  Global.Add("getenv","(",new OneOperator1s_<string*,string*>(ffgetenv));
  Global.Add("setenv","(",new OneOperator2_<long,string*,string*>(ffsetenv));
  Global.Add("unsetenv","(",new OneOperator1_<long,string*>(ffunsetenv));
    static string edpfilenameargstr=edpfilenamearg;
    Global.New("edpfilenamearg",CConstant<string *>(&edpfilenameargstr));
}

LOADFUNC(init);
