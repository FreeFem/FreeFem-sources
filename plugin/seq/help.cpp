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
// SUMMARY : online help in a FreeFEM script
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

// Usage:
//   load "help"
//   help(plot);       // a function: its help (if written below) and its signatures
//   help(medit);      // after load "medit"
//   help("solve");    // a keyword or a function given by its name
//   help();           // the list of the topics

#include "ff++.hpp"
#include <sstream>
#include <cxxabi.h>

using namespace std;

// ---------------------------------------------------------------------------
//  the written help: name, text
// ---------------------------------------------------------------------------
struct HelpTopic {
  const char *name;
  const char *text;
};

static const HelpTopic helptopics[] = {
  {"plot",
   "plot(objet1, objet2, ..., parametres) : trace dans la fenetre ffglut.\n"
   "\n"
   "  Objets : un maillage (mesh, mesh3, meshS, meshL), une fonction element fini u,\n"
   "  un champ de vecteurs [ux,uy] ou [ux,uy,uz], une courbe [x,y] (tableaux real[int]).\n"
   "\n"
   "  Parametres principaux :\n"
   "    wait=true          attend une touche avant de continuer le script\n"
   "    cmm=\"texte\"        titre (voir les options du graphe 1D plus bas)\n"
   "    fill=true          isovaleurs remplies\n"
   "    value=true         legende des valeurs des couleurs\n"
   "    nbiso=n            nombre d'isovaleurs (20 par defaut)\n"
   "    viso=tab           isovaleurs imposees (real[int])\n"
   "    nbarrow=n, varrow=tab   idem pour les fleches\n"
   "    coef=c             longueur des fleches\n"
   "    dim=3              vue 3D d'une solution 2D\n"
   "    bb=[[xmin,ymin],[xmax,ymax]]   cadre (axes fixes d'un graphe 1D)\n"
   "    aspectratio=true   meme echelle en x et en y\n"
   "    boundary=false     sans le bord du maillage\n"
   "    bw=true, grey=true noir et blanc, niveaux de gris\n"
   "    hsv=tab            palette de couleurs (triplets h,s,v)\n"
   "    ps=\"f.eps\", pdf=\"f.pdf\", svg=\"f.svg\"   sauvegarde dans un fichier\n"
   "    prev=true          garde le point de vue du plot precedent\n"
   "\n"
   "  Graphe 1D (courbes [x,y] ou fonction sur un meshL) : options dans cmm,\n"
   "  separees par | :\n"
   "    cmm=\"titre | xlabel=... | ylabel=... | legend=a,b | logx | logy | loglog\"\n"
   "\n"
   "  Touches dans la fenetre : ? aide, entree plot suivant, p plot precedent,\n"
   "  f remplissage, m maillage, v valeurs, 3 vue 2D/3D, k lignes/points ou partie\n"
   "  complexe, G grille, x y echelle log (graphe 1D), s surface/isosurfaces (3D),\n"
   "  F1 plan de coupe (3D), F2 le deplacer a la souris, F3 le figer sur l'objet,\n"
   "  w capture PPM, W capture PNG, Echap quitter.\n"
   "  Rendu soigne 2D et 3D (par defaut) : R active/desactive, C change de palette\n"
   "  (turbo, viridis, cool to warm, arc-en-ciel, gris), v barre de couleurs,\n"
   "  n N moins/plus d'isocontours (2D) ou d'isosurfaces (3D), o isocontours 2D\n"
   "  affiches ou masques, P valeur de la fonction sous la souris (2D et 3D).\n"
   "  Pour un rendu ParaView : include \"ffparaview.idp\" puis paraview(\"u\", Th, u);\n"
   "\n"
   "  Exemples :\n"
   "    plot(Th, wait=1);\n"
   "    plot(u, fill=1, value=1, cmm=\"solution\");\n"
   "    plot([x,y1], [x,y2], cmm=\"Titre | xlabel=t | legend=a,b\");"},

  {"medit",
   "medit(\"noms\", Th, sol1, sol2, ..., parametres) : trace dans la fenetre ffmedit\n"
   "  (load \"medit\").\n"
   "\n"
   "  \"noms\" : un nom par solution, separes par des espaces. Les options du\n"
   "  graphe 1D se mettent apres une barre verticale :\n"
   "    \"u w | titre | xlabel=... | ylabel=... | legend=a,b | logx | logy | loglog\n"
   "     | xrange=a:b | yrange=a:b\"\n"
   "  Une solution vectorielle s'ecrit entre crochets : [vx,vy] ou [vx,vy,vz].\n"
   "\n"
   "  Parametres :\n"
   "    wait=true          le script attend entree dans la fenetre ;\n"
   "                       wait=false : le script continue (animation)\n"
   "    order=0            solution donnee par element (P0) ; 1 par defaut (sommets)\n"
   "    bb=[[xmin,ymin],[xmax,ymax]]   axes fixes du graphe 1D\n"
   "    meditff=\"prog\"     programme ffmedit a lancer\n"
   "    save=\"f.sol\"       sauvegarde la solution dans un fichier\n"
   "    bin=1              transmission binaire (defaut)\n"
   "\n"
   "  Touches (les memes que ffglut) : entree plot suivant, p precedent, Echap quitter,\n"
   "  ? aide, f remplissage, o isolignes, n N moins/plus d'isolignes, m maillage,\n"
   "  B bord, l eclairage, v palette, P valeur sous la souris, j axes, 3 relief 2D,\n"
   "  h H hauteur du relief, k partie complexe, s isosurfaces, V fleches,\n"
   "  a A longueur des fleches, e E vue eclatee, T normales, U superposition,\n"
   "  F1 plan de coupe, F2 le deplacer, F3 le figer sur l'objet, = vue initiale,\n"
   "  z Z zoom, w capture PPM, b fond noir/blanc. Le reste : bouton droit.\n"
   "  Graphe 1D : k lignes/points, G grille, x y log, A figer les axes, i axes auto.\n"
   "\n"
   "  Exemples :\n"
   "    medit(\"u\", Th, u, wait=1);\n"
   "    medit(\"u grad\", Th3, u, [gx,gy,gz]);\n"
   "    medit(\"p\", Th, p, order=0);\n"
   "    medit(\"u | Titre | xlabel=x\", ThL, u, bb=[[0,-1],[1,1]]);"},

  {"savesol",
   "savesol(\"fichier.sol\", Th, sol1, sol2, ..., order=1) : ecrit des solutions au\n"
   "  format medit (load \"medit\"). order=0 pour une solution par element.\n"
   "  Exemple : savesol(\"u.sol\", Th, u);"},

  {"square",
   "square(nx, ny, [fx,fy], flags=, label=, region=) : maillage 2D du carre unite\n"
   "  ]0,1[^2 avec nx x ny intervalles, transforme par (x,y) -> (fx,fy).\n"
   "  Labels du bord : 1 bas, 2 droite, 3 haut, 4 gauche.\n"
   "  Exemples :\n"
   "    mesh Th = square(10, 10);\n"
   "    mesh Th = square(20, 10, [2*x, y]);"},

  {"buildmesh",
   "buildmesh(b1(n1) + b2(n2) + ...) : maillage 2D du domaine limite par des bords\n"
   "  (border), parcourus dans le sens direct ; n = nombre de points sur le bord\n"
   "  (negatif pour inverser le sens).\n"
   "  Exemple :\n"
   "    border C(t=0, 2*pi) { x=cos(t); y=sin(t); label=1; }\n"
   "    mesh Th = buildmesh(C(50));"},

  {"border",
   "border nom(t=a, b) { x=...; y=...; label=n; } : bord parametre utilise par\n"
   "  buildmesh. Voir help(buildmesh)."},

  {"movemesh",
   "movemesh(Th, [fx, fy]) : deplace les sommets du maillage par (x,y) -> (fx,fy).\n"
   "  En 3D : movemesh(Th3, [fx, fy, fz]).\n"
   "  Exemple : Th = movemesh(Th, [x + 0.1*y, y]);"},

  {"adaptmesh",
   "adaptmesh(Th, u, ...) : adapte le maillage a la fonction u (ou a plusieurs).\n"
   "  Parametres courants : err= (erreur visee), hmin=, hmax= (tailles d'arete),\n"
   "  nbvx= (nombre maximal de sommets), iso=true (maillage isotrope),\n"
   "  ratio= (rapport de taille entre aretes voisines).\n"
   "  Exemple : Th = adaptmesh(Th, u, err=0.01, hmin=1e-3);"},

  {"cube",
   "cube(nx, ny, nz, [fx,fy,fz], flags=, label=, region=) : maillage 3D du cube\n"
   "  unite (load \"msh3\" si necessaire).\n"
   "  Exemple : mesh3 Th = cube(10, 10, 10);"},

  {"segment",
   "segment(n, [fx,fy,fz]) : maillage 1D (meshL) du segment [0,1] avec n\n"
   "  intervalles, transforme par x -> (fx,fy,fz) (load \"msh3\" si necessaire).\n"
   "  Exemple : meshL ThL = segment(20, [4*x, 0, 0]);"},

  {"savemesh",
   "savemesh(Th, \"fichier.msh\") : sauvegarde un maillage ; l'extension choisit le\n"
   "  format (.msh FreeFEM, .mesh medit, ...).\n"
   "  Exemple : savemesh(Th, \"Th.msh\");"},

  {"readmesh",
   "readmesh(\"fichier.msh\") : lit un maillage 2D (readmesh3 pour un maillage 3D).\n"
   "  Exemple : mesh Th = readmesh(\"Th.msh\");"},

  {"int1d",
   "int1d(Th, labels)(f) : integrale de f sur les bords de labels donnes.\n"
   "  Exemple : real l = int1d(Th, 1, 3)(1.);"},

  {"int2d",
   "int2d(Th)(f) : integrale de f sur le maillage 2D (ou sur une region :\n"
   "  int2d(Th, region)(f)).\n"
   "  Exemple : real m = int2d(Th)(u);"},

  {"int3d",
   "int3d(Th)(f) : integrale de f sur le maillage 3D.\n"
   "  Exemple : real m = int3d(Th3)(u);"},

  {"on",
   "on(labels, u=g) : condition de Dirichlet u = g sur les bords de labels donnes,\n"
   "  dans solve, problem ou varf.\n"
   "  Exemple : + on(1, 2, 3, 4, u=0)"},

  {"solve",
   "solve nom(u, v) = forme bilineaire + forme lineaire + conditions ;\n"
   "  definit et resout immediatement le probleme. Les termes de la forme lineaire\n"
   "  s'ecrivent avec un signe moins.\n"
   "  Exemple :\n"
   "    fespace Vh(Th, P1); Vh u, v;\n"
   "    solve Poisson(u, v) = int2d(Th)(dx(u)*dx(v) + dy(u)*dy(v))\n"
   "                        - int2d(Th)(1*v) + on(1, 2, 3, 4, u=0);"},

  {"problem",
   "problem nom(u, v) = ... ; meme syntaxe que solve, mais le probleme est resolu\n"
   "  seulement quand on ecrit nom; (utile dans une boucle). Voir help(\"solve\")."},

  {"varf",
   "varf a(u, v) = forme ; definit une forme variationnelle pour construire\n"
   "  matrices et seconds membres :\n"
   "    matrix A = a(Vh, Vh);  real[int] b = l(0, Vh);"},

  {"fespace",
   "fespace Vh(Th, P1); : espace d'elements finis sur le maillage Th.\n"
   "  Elements courants : P0, P1, P2, P1dc, P1b, RT0, P1nc ; espace produit :\n"
   "  fespace Wh(Th, [P2, P2, P1]);"},

  {"interpolate",
   "interpolate(Wh, Vh) : matrice d'interpolation de Vh vers Wh.\n"
   "  Pour interpoler une fonction il suffit d'ecrire : Wh w = u;"},

  {"help",
   "help(fonction) ou help(\"mot-cle\") : aide d'une fonction ou d'un mot-cle\n"
   "  (load \"help\"). help() donne la liste des sujets. Pour une fonction sans aide\n"
   "  ecrite, help affiche ses signatures (types des arguments)."},
};

static const int nhelptopics = sizeof(helptopics) / sizeof(HelpTopic);

static const HelpTopic *FindTopic(const string &name) {
  for (int i = 0; i < nhelptopics; i++)
    if (name == helptopics[i].name) return &helptopics[i];
  return 0;
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
  int status = 0;
  char *d = abi::__cxa_demangle(mangled.c_str( ), 0, 0, &status);
  string t = (status == 0 && d) ? d : mangled;
  free(d);
  ReplaceAll(t, "std::__1::", "");
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
    {"Polymorphic", "func"}, {"CDomainOfIntegration", "domaine d'integration"},
    {"NothingType", "void"}, {"ostream", "ofstream"}};
  for (size_t i = 0; i < sizeof(tr) / sizeof(tr[0]); i++)
    if (t == tr[i][0]) return tr[i][1];
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

static void ffShowHelp(const string &name, const Polymorphic *p) {
  const HelpTopic *h = FindTopic(name);
  cout << "\n ---------------- help(" << name << ") ----------------\n";
  if (h) cout << h->text << "\n";
  if (p) cout << (h ? "\n" : "") << "  Signatures :\n" << Signatures(p, name);
  if (!h && !p) {
    cout << "  Pas d'aide pour \"" << name << "\" (ni fonction, ni sujet connu).\n"
         << "  help() donne la liste des sujets.";
    if (name == "medit" || name == "savesol") cout << " Pour medit : load \"medit\".";
    cout << "\n";
  }
  cout << " -----------------------------------------------\n" << endl;
}

static void ffShowTopics( ) {
  cout << "\n ---------------- help ----------------\n"
       << "  help(fonction) ou help(\"mot-cle\"). Sujets avec une aide ecrite :\n  ";
  int l = 2;
  for (int i = 0; i < nhelptopics; i++) {
    string n = helptopics[i].name;
    if (l + n.size( ) + 2 > 78) cout << "\n  ", l = 2;
    cout << n << (i + 1 < nhelptopics ? ", " : "\n");
    l += n.size( ) + 2;
  }
  cout << "  Pour toute autre fonction, help(f) affiche ses signatures.\n"
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

static void Load_Init( ) {
  Global.Add("help", "(", new HelpFunction, new HelpName, new HelpTopics);
}

LOADFUNC(Load_Init)
