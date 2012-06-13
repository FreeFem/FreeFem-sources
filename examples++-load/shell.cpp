// Example C++ function "myfunction", dynamically loaded into "load.edp"
// ---------------------------------------------------------------------
// $Id$
#include <ff++.hpp>
#include <unistd.h>
#include <sys/stat.h>
#include <dirent.h> 
using namespace Fem2D;
long ff_chdir(string * c) { return chdir(c->c_str());}
long ff_rmdir(string * c) { return rmdir(c->c_str());}
long ff_unlink(string * c) { return unlink(c->c_str());}
long ff_mkdir(string * c,long mm) {mode_t m=mm;cout << " mkdir " << *c << "mode =" << m << endl;  return mkdir(c->c_str(),m);}
long ff_chmod(string * c,long mm) {mode_t m=mm;cout << " mkdir " << *c << "mode =" << m << endl;  return chmod(c->c_str(),m);}
long ff_mkdir(string * c) {mode_t m=07777;cout << " mkdir " << *c <<" mode =" << m << endl;  return mkdir(c->c_str(),m);}
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


extern  mylex *zzzfff;

void init(){
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
  Global.Add("mkdir","(",new OneOperator2<long,string*,long>(ff_mkdir));
  Global.Add("chmod","(",new OneOperator2<long,string*,long>(ff_chmod));
  Global.Add("mkdir","(",new OneOperator1<long,string*>(ff_mkdir));
  Global.Add("stat","(",new OneOperator1<long,string*>(ff_stat));
  Global.Add("isdir","(",new OneOperator1<long,string*>(ff_isdir));
}

LOADFUNC(init);
