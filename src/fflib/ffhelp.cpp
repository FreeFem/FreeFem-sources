/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : online help in a FreeFEM script (in the kernel: no load)
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE

// Usage in a script:
//   help(plot);         // a function: its help (if written) and its signatures
//   help(medit);        // after load "medit"
//   help("fvscheme");   // a keyword (a type, solve, ...) or a function given by its name
//   help();             // the list of the topics
// the written help of a topic is the text file help/<name>.txt of a directory of the include path (idp/help in the
// sources, installed with the idp files): not in the C++ sources; a file whose first line is @other gives the help
// of other; the signatures, the named parameters and the members of the types are found in the kernel

#include "ff++.hpp"
#include "lex.hpp"
#include <sstream>
#include <set>
#include <fstream>
#include <dirent.h>
#include <cxxabi.h>
#include <csignal>
#include <csetjmp>

using namespace std;

// ---------------------------------------------------------------------------
//  the written help: the files help/<name>.txt of the include path
// ---------------------------------------------------------------------------
static vector< string > HelpDirs( ) {
  vector< string > d;
  if (zzzfff)
    for (mylex::ICffincludedir i = zzzfff->ffincludedir.begin( ); i != zzzfff->ffincludedir.end( ); ++i) {
      string s = *i;
      if (s.size( ) && s[s.size( ) - 1] != '/' && s[s.size( ) - 1] != '\\') s += '/';
      d.push_back(s + "help/");
    }
  return d;
}

// the written help of name in text, false if there is none (a name, not a path)
static bool FindTopic(const string &name, string &text, int depth = 0) {
  if (name.empty( ) || name.find_first_of("/\\.") != string::npos || depth > 4) return false;
  vector< string > dirs = HelpDirs( );
  for (size_t k = 0; k < dirs.size( ); k++) {
    ifstream f((dirs[k] + name + ".txt").c_str( ));
    if (!f) continue;
    ostringstream o;
    o << f.rdbuf( );
    text = o.str( );
    while (text.size( ) && (text[text.size( ) - 1] == '\n' || text[text.size( ) - 1] == '\r')) text.erase(text.size( ) - 1);
    if (text.size( ) > 1 && text[0] == '@') {// an alias: the help of another topic
      string other = text.substr(1, text.find_first_of("\r\n") - 1);
      return FindTopic(other, text, depth + 1);
    }
    return true;
  }
  return false;
}

// the names of the topics with a written help (the files *.txt of the help directories)
static set< string > HelpTopicNames( ) {
  set< string > names;
  vector< string > dirs = HelpDirs( );
  for (size_t k = 0; k < dirs.size( ); k++) {
    DIR *d = opendir(dirs[k].c_str( ));
    if (!d) continue;
    for (struct dirent *e = readdir(d); e; e = readdir(d)) {
      string n = e->d_name;
      if (n.size( ) > 4 && n.compare(n.size( ) - 4, 4, ".txt") == 0) names.insert(n.substr(0, n.size( ) - 4));
    }
    closedir(d);
  }
  return names;
}

// ---------------------------------------------------------------------------
//  the functions of FreeFEM: name <-> Polymorphic in the tables of identifiers
// ---------------------------------------------------------------------------
static const char *FindName(const Polymorphic *p) {
  for (ListOfTOfId::const_iterator t = tables_of_identifier.begin( ); t != tables_of_identifier.end( );
       ++t)
    for (TableOfIdentifier::const_iterator i = (*t)->m.begin( ); i != (*t)->m.end( ); ++i)
      if (i->second.second == p) return i->first;
  return 0;
}

static const Polymorphic *FindPolymorphic(const string &name) {
  for (ListOfTOfId::const_iterator t = tables_of_identifier.begin( ); t != tables_of_identifier.end( );
       ++t) {
    TableOfIdentifier::const_iterator i = (*t)->m.find(name.c_str( ));
    if (i != (*t)->m.end( )) return dynamic_cast< const Polymorphic * >(i->second.second);
  }
  return 0;
}

// ---------------------------------------------------------------------------
//  the signatures: FreeFEM prints the C++ types ( <N5Fem2D4MeshE> ), we print the FreeFEM types
// ---------------------------------------------------------------------------
static void ReplaceAll(string &s, const string &a, const string &b) {
  for (size_t i = 0; (i = s.find(a, i)) != string::npos; i += b.size( )) s.replace(i, a.size( ), b);
}

static string TypeName(const string &mangled) {
  if (mangled == typeid(void *).name( )) return "mpiComm";// pcommworld (comm=): the only void* parameter
  int status = 0;
  char *d = abi::__cxa_demangle(mangled.c_str( ), 0, 0, &status);
  string t = (status == 0 && d) ? d : mangled;
  free(d);
  ReplaceAll(t, "std::__1::", "");
  ReplaceAll(t, " >", ">");// the demangler of MinGW: KN_<complex > (Windows)
  ReplaceAll(t, "__cxx11::", "");// libstdc++ (Linux, Windows): std::__cxx11::string
  ReplaceAll(t, "std::", "");
  ReplaceAll(t, "basic_string<char, char_traits<char>, allocator<char> >", "string");
  ReplaceAll(t, "basic_string<char, char_traits<char>, allocator<char>>", "string");
  ReplaceAll(t, "Fem2D::", "");
  ReplaceAll(t, "complex<double>", "complex");
  ReplaceAll(t, " const", "");
  while (t.size( ) && (t[t.size( ) - 1] == '*' || t[t.size( ) - 1] == '&' || t[t.size( ) - 1] == ' '))
    t.erase(t.size( ) - 1);
  static const char *const tr[][2] = {
    {"long", "int"},     {"double", "real"},           {"Mesh", "mesh"},
    {"Mesh3", "mesh3"},  {"MeshS", "meshS"},           {"MeshL", "meshL"},
    {"E_Array", "[...]"}, {"KN<double>", "real[int]"}, {"KN_<double>", "real[int]"},
    {"KN<long>", "int[int]"}, {"KN_<long>", "int[int]"}, {"KN<complex>", "complex[int]"},
    {"KN_<complex>", "complex[int]"}, {"KNM<double>", "real[int,int]"},
    {"Polymorphic", "func"}, {"CDomainOfIntegration", "integration domain"},
    {"NothingType", "void"}, {"ostream", "ofstream"}};
  for (size_t i = 0; i < sizeof(tr) / sizeof(tr[0]); i++)
    if (t == tr[i][0]) return tr[i][1];
  if (t.compare(0, 19, "GQuadratureFormular") == 0) return "quadrature formula";// qft=, qfe=, ...
  // the intermediate types of the members of the meshes and matrices (help("mesh"))
  ReplaceAll(t, "(anonymous namespace)::", "");
  if (t.find("::BE") != string::npos) return "boundary elements";
  if (t.find("::Adj") != string::npos) return "adjacent elements";
  if (t.compare(0, 18, "GlgBoundaryElement") == 0 || t == "lgBoundaryEdge") return "boundary element";
  if (t.compare(0, 10, "GlgElement") == 0 || t == "lgElement") return "element";
  if (t.compare(0, 9, "GlgVertex") == 0 || t == "lgVertex") return "vertex";
  if (t == "MeshPoint") return "point";
  if (t.compare(0, 18, "newpMatrice_Creuse") == 0 || t.compare(0, 13, "SparseMat_Add") == 0) return "matrix";
  if (t.compare(0, 10, "TheCoefMat") == 0) return "real[int] (the coefficients)";
  if (t.compare(0, 10, "TheDiagMat") == 0) return "real[int] (the diagonal)";
  if (t == "Inv_KN_long") return "int[int]^-1";
  if (t == "v_fes" || t == "v_fesL" || t == "v_fesS" || t == "v_fes3") return "fespace";// varf a; a(Vh, Vh)
  if (t == "vect_generic_v_fes") return "[fespace,...]";
  if (t.compare(0, 26, "Call_CompositeFormBilinear") == 0) return "matrix (blocks)";
  if (t.compare(0, 17, "Call_FormBilinear") == 0) return "matrix";
  if (t.compare(0, 15, "Call_FormLinear") == 0) return "real[int]";
  if (t == "SubArray") return "range a:b";// s(a:b)
  if (t == "SubString") return "string";
  if (t.compare(0, 14, "Matrice_Creuse") == 0) return "matrix";
  if (t.compare(0, 15, "map<string, str") == 0) return "string[string]";
  // the expressions on arrays (ex. sin(a))
  if (t.compare(0, 12, "F_KN_<double") == 0) return "real[int]";
  if (t.compare(0, 13, "F_KN_<complex") == 0) return "complex[int]";
  return t;
}

// the text of "o << *p":  "\t (\t  <R> :   <A1>, <A2>... )"  ->  "    R name(A1, A2, ...)"
static string Signatures(const Polymorphic *p, const string &name) {
  ostringstream o;
  o << *p;
  istringstream in(o.str( ));
  string line, res;
  while (getline(in, line)) {
    size_t c = line.find(" : ");
    size_t a = line.find('(');
    if (a == string::npos || c == string::npos || c < a) continue;
    string r, args;
    size_t b = line.find('<', a), e = line.find('>', b);
    if (b == string::npos || e == string::npos || e > c) continue;
    r = TypeName(line.substr(b + 1, e - b - 1));
    for (size_t i = c; (b = line.find('<', i)) != string::npos; i = e + 1) {
      e = line.find('>', b);
      if (e == string::npos) break;
      args += (args.size( ) ? ", " : "") + TypeName(line.substr(b + 1, e - b - 1));
    }
    if (line.find("...", c) != string::npos) args += args.size( ) ? ", ..." : "...";
    res += "    " + r + " " + name + "(" + args + ")\n";
  }
  return res;
}

// ---------------------------------------------------------------------------
//  the named parameters: the code of each operator of f called with dummy arguments, SetNameParam throws
//  its list (ffhelp_probe, AFunction2.cpp); nothing is executed; an operator which looks at its arguments
//  before SetNameParam may fail: its exception, or a crash (signal), only stops its probe
// ---------------------------------------------------------------------------
bool ffhelp_probe = false;

class E_HelpDummy : public E_F0 {// a dummy argument: an evaluation stops the probe of this operator
 public:
  AnyType operator( )(Stack) const { throw ffHelpNamedParams{0, nullptr}; }
};

#ifdef _WIN32
static jmp_buf help_jmp;
#define HELP_SETJMP(b) setjmp(b)
#define HELP_LONGJMP(b) longjmp(b, 1)
#else
static sigjmp_buf help_jmp;
#define HELP_SETJMP(b) sigsetjmp(b, 1)
#define HELP_LONGJMP(b) siglongjmp(b, 1)
#endif
static void HelpSignal(int) { HELP_LONGJMP(help_jmp); }

class NullBuffer : public streambuf {// the messages of the probe are not printed
 protected:
  int overflow(int c) { return c; }
};

// the named parameters of one operator: (name, type), empty if none or if the probe failed
static vector< pair< string, string > > ProbeOperator(const OneOperator *oo) {
  vector< pair< string, string > > res;
  int n = oo->NbTypes( );
  C_F0 *a = new C_F0[n > 0 ? n : 1];
  for (int i = 0; i < n; i++) a[i] = C_F0(Type_Expr(oo->TypeAt(i), new E_HelpDummy));
  basicAC_F0 args;
  args = pair< int, C_F0 * >(n, a);
  ffHelpNamedParams r{0, nullptr};
  static const int sigs[] = {SIGSEGV, SIGFPE, SIGILL, SIGABRT
#ifdef SIGBUS
                             , SIGBUS
#endif
  };
  const int nsigs = sizeof(sigs) / sizeof(sigs[0]);
  void (*old[8])(int);
  NullBuffer nb;
  streambuf *bout = cout.rdbuf(&nb);// (in FreeFEM cerr is cout: throwassert.hpp)
  for (int k = 0; k < nsigs; k++) old[k] = signal(sigs[k], HelpSignal);
  ffhelp_probe = true;
  if (HELP_SETJMP(help_jmp) == 0) {
    try {
      oo->code(args);
    } catch (ffHelpNamedParams &e) {
      r = e;
    } catch (...) {
    }
  }
  ffhelp_probe = false;
  for (int k = 0; k < nsigs; k++) signal(sigs[k], old[k]);
  cout.rdbuf(bout);
  delete[] a;
  for (int i = 0; i < r.n; i++) res.push_back(make_pair(string(r.l[i].name), TypeName(r.l[i].type->name( ))));
  return res;
}

// the named parameters of all the operators of the function, each name once (its types joined by " or ")
static string NamedParameters(const Polymorphic *p) {
  vector< pair< string, string > > all;
  for (const OneOperator *oo = p->FindOp("("); oo; oo = oo->Next( )) {
    vector< pair< string, string > > l = ProbeOperator(oo);
    for (size_t i = 0; i < l.size( ); i++) {
      size_t j = 0;
      while (j < all.size( ) && all[j].first != l[i].first) j++;
      if (j == all.size( ))
        all.push_back(l[i]);
      else if (("|" + all[j].second + "|").find("|" + l[i].second + "|") == string::npos)
        all[j].second += "|" + l[i].second;
    }
  }
  string res;
  for (size_t j = 0; j < all.size( ); j++) {
    string t = all[j].second;
    ReplaceAll(t, "|", " or ");
    res += "    " + all[j].first + "= " + t + "\n";
  }
  return res;
}

// the type of the keyword name (mesh, meshL, matrix, ...), 0 if name is not a type
static aType FindType(const string &name) {
  aType t = 0;
  int r = 0;
  if (zzzfff && zzzfff->InMotClef(name.c_str( ), t, r)) return t;
  return 0;
}

// the members of a type, read in its table of identifiers: x.name (op "."), x(...) and x[...] (keys "(" and
// "["), and the members of their results (depth 2: Th[k][j].x), one line "access : type" each
static void MemberLines(aType t, const string &pre, int depth, set< aType > &path, vector< string > &out,
                        set< string > &uniq) {
  aType ts[2] = {t, t->right( )};
  for (int k = 0; k < 2; k++) {
    if (!ts[k] || (k == 1 && ts[1] == ts[0])) continue;
    const TableOfIdentifier &ti = ts[k]->Members( );
    for (TableOfIdentifier::const_iterator i = ti.m.begin( ); i != ti.m.end( ); ++i) {
      const Polymorphic *p = dynamic_cast< const Polymorphic * >(i->second.second);
      if (!p) continue;
      string key = i->first;
      bool call = key == "(" || key == "[";
      for (const OneOperator *oo = p->FindOp(call ? "" : "."); oo; oo = oo->Next( )) {
        ostringstream o;
        o << *oo;// "<R> :  <A0>, <A1>, ..." (A0: the object)
        string sg = o.str( );
        vector< string > tk;
        for (size_t a = 0, b; (a = sg.find('<', a)) != string::npos; a = b + 1) {
          b = sg.find('>', a);
          if (b == string::npos) break;
          tk.push_back(sg.substr(a + 1, b - a - 1));
        }
        if (tk.empty( )) continue;
        string acc;
        if (call) {
          string args;
          for (size_t a = 2; a < tk.size( ); a++) args += (args.size( ) ? "," : "") + TypeName(tk[a]);
          acc = pre + (key == "(" ? "(" + args + ")" : "[" + args + "]");
        } else
          acc = pre + "." + key;
        string rt = TypeName(tk[0]), line = acc + " : " + rt;
        size_t nout = out.size( );
        if (uniq.insert(line).second) out.push_back(line);
        // the members of the result, not those of a scalar (a real is also a function of (x,y,z))
        static const char *const scal[] = {"real", "int", "bool", "complex", "string"};
        bool scalar = false;
        for (int q = 0; q < 5; q++) scalar = scalar || rt == scal[q];
        map< const string, basicForEachType * >::const_iterator r = map_type.find(tk[0]);
        if (!scalar && depth < 2 && r != map_type.end( ) && r->second && !path.count(r->second) && out.size( ) < 300) {
          path.insert(r->second);
          MemberLines(r->second, acc, depth + 1, path, out, uniq);
          path.erase(r->second);
          // an intermediate result used only with arguments (Th.be(i), A.resize(n,m)): no line Th.be
          bool onlycall = out.size( ) > nout + 1;
          for (size_t q = nout + 1; q < out.size( ) && onlycall; q++)
            onlycall = out[q].compare(0, acc.size( ) + 1, acc + "(") == 0;
          if (onlycall && out.size( ) > nout && out[nout] == line) out.erase(out.begin( ) + nout);
        }
      }
    }
  }
}

static string Members(aType t, const string &name) {
  vector< string > out;
  set< string > uniq;
  set< aType > path;
  path.insert(t);
  string x = name.compare(0, 4, "mesh") == 0 ? "Th" : "x";
  MemberLines(t, x, 0, path, out, uniq);
  string res;
  for (size_t i = 0; i < out.size( ); i++) res += "    " + out[i] + "\n";
  return res;
}

static void ffShowHelp(const string &name, const Polymorphic *p) {
  if (mpirank) return;
  string ht;
  bool h = FindTopic(name, ht);
  cout << "\n ---------------- help(" << name << ") ----------------\n";
  if (h) cout << ht << "\n";
  if (p) {
    string np = NamedParameters(p);
    cout << (h ? "\n" : "") << "  Signatures:\n" << Signatures(p, name);
    if (np.size( )) cout << "\n  Named parameters:\n" << np;
  }
  aType ty = p ? 0 : FindType(name);
  if (ty) {
    string mb = Members(ty, name);
    cout << (h ? "\n" : "") << "  " << name << " is a type";
    if (mb.size( )) cout << "; its members (" << (name.compare(0, 4, "mesh") == 0 ? "Th" : "x") << " a " << name
                        << "):\n" << mb;
    else cout << ".\n";
  }
  if (!h && !p && !ty) {
    cout << "  No help for \"" << name << "\" (neither a function nor a known topic).\n"
         << "  help() gives the list of the topics.";
    if (name == "medit" || name == "savesol") cout << " For medit: load \"medit\".";
    cout << "\n";
  }
  cout << " -----------------------------------------------\n" << endl;
}

static void ffShowTopics( ) {
  if (mpirank) return;
  cout << "\n ---------------- help ----------------\n"
       << "  help(function) or help(\"keyword\"). Topics with a written help:\n  ";
  set< string > names = HelpTopicNames( );
  size_t l = 2, k = 0, nt = names.size( );
  if (!nt) cout << "(none: no directory help/ in the include path)\n";
  for (set< string >::const_iterator i = names.begin( ); i != names.end( ); ++i) {
    const string &n = *i;
    if (l + n.size( ) + 2 > 78) cout << "\n  ", l = 2;
    cout << n << (++k < nt ? ", " : "\n");
    l += n.size( ) + 2;
  }
  cout << "  For any other function, help(f) prints its signatures and named parameters.\n"
       << " --------------------------------------\n"
       << endl;
}

// ---------------------------------------------------------------------------
//  the operators help(f), help("name"), help()
// ---------------------------------------------------------------------------
class E_Help : public E_F0 {
 public:
  string name;              // the name of the function, help(f)
  const Polymorphic *poly;  // the function, help(f)
  Expression ename;         // the name, help("name")
  bool topics;              // help()

  E_Help( ) : poly(0), ename(0), topics(true) {}
  E_Help(const string &n, const Polymorphic *p) : name(n), poly(p), ename(0), topics(false) {}
  E_Help(Expression e) : poly(0), ename(e), topics(false) {}

  AnyType operator( )(Stack stack) const {
    if (topics)
      ffShowTopics( );
    else if (ename) {
      string n = *GetAny< string * >((*ename)(stack));
      ffShowHelp(n, FindPolymorphic(n));
    } else
      ffShowHelp(name, poly);
    return SetAny< long >(0L);
  }
  operator aType( ) const { return atype< long >( ); }
};

class HelpFunction : public OneOperator {
 public:
  HelpFunction( ) : OneOperator(atype< long >( ), atype< Polymorphic * >( )) {}
  E_F0 *code(const basicAC_F0 &args) const {
    const Polymorphic *p = dynamic_cast< const Polymorphic * >(args[0].LeftValue( ));
    if (!p) CompileError("help(f): f must be a function");
    const char *n = FindName(p);
    return new E_Help(n ? n : "?", p);
  }
};

class HelpName : public OneOperator {
 public:
  HelpName( ) : OneOperator(atype< long >( ), atype< string * >( )) {}
  E_F0 *code(const basicAC_F0 &args) const { return new E_Help(CastTo< string * >(args[0])); }
};

class HelpTopics : public OneOperator {
 public:
  HelpTopics( ) : OneOperator(atype< long >( )) {}
  E_F0 *code(const basicAC_F0 &args) const { return new E_Help( ); }
};

// called by Init_map_type (AFunction.cpp): help in the kernel, no load
void Load_Init_help( ) {
  Global.Add("help", "(", new HelpFunction, new HelpName, new HelpTopics);
}
