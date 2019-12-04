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
// SUMMARY : FastMarchingBundle
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Jean-Marie Mirebeau
// E-MAIL  : jean-marie.mirebeau@math.u-psud.fr

#ifndef SORTED_LIST_H
#define SORTED_LIST_H

#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include "float.h"
#include "BasicMath.h"
using namespace std;

/*
 *
 * The template class Tab implements arrays of arbitrary length, which does not need to be fixed in
 * advance.
 *
 * The template class SortedList implements sorted lists of arbitrary length, which does not need to
 * be fixed in advance either. Common usages, such as insersion and deletion, access to minimal and
 * maximal element, have the complexity log(n).
 *
 * The template class RBTree implements Balanced Trees using the Red/Black labels approach. It
 * should not be used directly, but only through SortedList.
 */

template< class TabElement >
class Tab;
template< class TreeLabel >
class SortedList;
template< class TreeLabel >
class RBTree;

/****************************** Arrays of arbitrary length : Tab
 * *******************************************/

template< class TabElement >
class Tab {
 public:
  Tab( ) : cardMax(startCard), growIndex(0), max_accessed_pos(-1) {
    elements[growIndex++].resize(startCard);
  };
  Tab(const Tab< TabElement > &tab) : max_accessed_pos(tab.max_accessed_pos) {
    cout << "Tab constructor Warning : copying Tab of cardinality " << tab.cardMax
         << "; max accessed pos : " << tab.max_accessed_pos << endl;
    elements[growIndex++].resize(startCard);

    for (int i = 0; i < tab.cardMax; i++) {
      Element(i) = tab[i];
    }
  }

  TabElement &operator[](int pos) {    // TabElement & Element(int pos)
    assert_msg(pos >= 0, "Tab::Element Error : Negative index " << pos);
    if (pos >= cardMax) {
      const bool hasGrown = grow( );    // contents of an assertion are not executed if NDEBUG
      assert_msg(hasGrown, "Tab::Element Error : Maximum array size excessed. " << pos);
      return operator[](pos);
    }

    max_accessed_pos = max(pos, max_accessed_pos);
    if (pos < startCard) {
      return elements[0][pos];
    }

    int i, pow;

    for (i = growIndex - 1, pow = cardMax / 2; pos < pow; i--, pow /= 2) {
    }

    ;
    // log(cardMax) complexity in worst case. Unit in average thanks to reversed loop
    return elements[i][pos - pow];
  };
  const TabElement &operator[](int pos) const {
    assert_msg(pos >= 0, "Tab::Element Error : Negative index " << pos);
    assert_msg(pos <= max_accessed_pos,
               "Tab::Element const Error : max_accessed_pos exceeded. " << pos);
    if (pos < startCard) {
      return elements[0][pos];
    }

    int i, pow;

    for (i = growIndex - 1, pow = cardMax / 2; pos < pow; i--, pow /= 2) {
    }

    ;
    return elements[i][pos - pow];
  };
  void export_content(const char *filename, Format_Math format = Mathematica,
                      bool one_per_line = false) const {
    ofstream data_out;

    data_out.open(filename);
    print_array(data_out << format, *this, one_per_line);
    data_out.close( );
  }

  void sort( );
  int max_accessed_pos;
  int card( ) const { return max_accessed_pos + 1; }

  TabElement *next( ) { return &operator[](max_accessed_pos + 1); }

  int index(TabElement const *ptr) const {
    const int j0 = int(ptr - &elements[0][0]);

    if (0 <= j0 && j0 < startCard) {
      return j0;
    }

    int i, pow;

    for (i = growIndex - 1, pow = cardMax / 2; i >= 1; i--, pow /= 2) {
      const int j = int(ptr - &elements[i][0]);
      if (0 <= j && j < pow) {
        return pow + j;
      }
    }

    cout << "Tab::index error : element does not belong to tab" << endl;
    return -1;
  }

 private:
  int cardMax;
  int growIndex;
  const static int startCard = 4;
  const static int growMax = 30;

  vector< TabElement > elements[growMax];

  bool grow( ) {
    if (growIndex == growMax) {
      return false;
    }

    elements[growIndex++].resize(cardMax);
    cardMax *= 2;
    return true;
  }

  TabElement &Element(int pos) { return operator[](pos); }

  const TabElement &Element(int pos) const { return operator[](pos); }
};

template< class TabElement >
void Tab< TabElement >::sort( ) {
  if (max_accessed_pos > 50) {    // tri n ln(n).
    SortedList< TabElement > list;

    for (int i = 0; i <= max_accessed_pos; ++i) {
      list.insert(Element(i));
    }

    for (int i = 0; i <= max_accessed_pos; ++i) {
      Element(i) = list.pop( );
    }

    return;
  }

  TabElement swap;

  for (int i = 0; i < max_accessed_pos; ++i) {
    if (operator[](i + 1) < operator[](i)) {
      swap = Element(i + 1);
      Element(i + 1) = Element(i);
      int j;

      for (j = i - 1; j >= 0 && Element(j) > swap; --j) {
        Element(j + 1) = Element(j);
      }

      Element(j + 1) = swap;
    }
  }
}

template< class E >
void print_array(ostream &f, const Tab< E > &tab, bool one_per_line = false) {
  const int N = tab.max_accessed_pos + 1;

  if (one_per_line) {
    for (int i = 0; i < N; i++) {
      f << tab[i] << endl;
    }
  } else {
    for (int i = 0; i < N; i++) {
      f << tab[i] << " ";
    }
  }
}

template< class E >
void print_array(ostream_math f, const Tab< E > &tab, bool one_per_line = false) {
  if (f.format == Mathematica) {
    const int N = tab.max_accessed_pos + 1;
    if (N <= 0) {
      f << "{}";
      return;
    }

    f << "{";

    for (int i = 0; i < N; i++) {
      f << tab[i];
      if (i < N - 1) {
        f << ",";
      }
    }

    f << "}";
  } else {
    print_array(f.os, tab, one_per_line);
    return;
  }
}

template< class TabElement >
ostream &operator<<(ostream &f, const Tab< TabElement > &tab) {
  print_array(f, tab);
  return f;
}

template< class TabElement >
ostream_math operator<<(ostream_math f, const Tab< TabElement > &tab) {
  print_array(f, tab);
  return f;
}

// ********************** Reservoir *********************
// Reservoir of copies of a given object.  References are consistent over time.
template< class E >
class Reservoir {
  Tab< E > reserve;
  mutable vector< E * > unused;

 public:
  int card( ) const { return int(reserve.card( ) - unused.size( )); }

  E *next( ) {
    if (unused.empty( )) {
      return reserve.next( );
    }

    E *e = unused.back( );
    unused.pop_back( );
    return e;
  }    // bizarre that pop is void

  bool free(E *e) {
    if (reserve.index(e) < 0) {
      return false;
    }

    unused.push_back(e);
    return true;
  }

  void enumerate(vector< E * > &elems);
  void enumerate(vector< const E * > &elems) const;
};

template< class E >
void Reservoir< E >::enumerate(vector< E * > &elems) {
  sort(unused.begin( ), unused.end( ));
  int u = 0;    // vector<E*>::const_iterator u does not work for some reason

  for (int i = 0; i < reserve.card( ); ++i) {
    E *const e = &reserve[i];
    if (e == unused[u]) {
      ++u;
    } else {
      elems.push_back(e);
    }
  }
}

template< class E >
void Reservoir< E >::enumerate(vector< const E * > &elems) const {
  sort(unused.begin( ), unused.end( ));
  int u = 0;    // vector<E*>::const_iterator u does not work for some reason

  for (int i = 0; i < reserve.card( ); ++i) {
    const E *const e = &reserve[i];
    if (u < unused.size( ) && e == unused[u]) {
      ++u;
    } else {
      elems.push_back(e);
    }
  }
}

/************************ Safe Vector (at debug time) *********************/
// same as standard library's vector, but with additional checks at DEBUG time. Should be zero
// overhead at non debug time

template< class E >
class safe_vector : public vector< E > {
 public:
  E &operator[](size_t n) {
    assert(0 <= n && n < this->size( ));
    return this->vector< E >::operator[](n);
  }

  const E &operator[](size_t n) const {
    assert(0 <= n && n < this->size( ));
    return this->vector< E >::operator[](n);
  }

  E &front( ) {
    assert(0 < this->size( ));
    return this->vector< E >::front( );
  }

  const E &front( ) const {
    assert(0 < this->size( ));
    return this->vector< E >::front( );
  }

  E &back( ) {
    assert(0 < this->size( ));
    return this->vector< E >::back( );
  }

  const E &back( ) const {
    assert(0 < this->size( ));
    return this->vector< E >::back( );
  }
};

/************************ Sorted List *********************/

// new implementation : a simple interface with std::set
template< class E >
class SortedList {
  set< E > s;

 public:
  int Card( ) { return int(s.size( )); }

  bool contains(E e) { return s.count(e); }

  E min( ) { return *s.begin( ); }    // empty cases ?

  E max( ) { return *--s.end( ); }

  bool insert(E e) { return s.insert(e).second; }

  bool remove(E e) { return s.erase(e); }

  E pop( ) {
    E e = min( );
    remove(e);
    return e;
  }

  void print(ostream &f) {
    Tab< E > tab;
    enumerate(tab);
    f << tab;
  }

  // int enumerate(Tab<E> &tab){set<E>::const_iterator it; int i; for(it=s.begin(), i=0;
  // it!=s.end(); ++it, ++i){tab[i]=*it;} return i;} int enumerate(Tab<E> &tab){for(int i=0;
  // i<Card(); ++i){tab[i]=s[i];} return Card();}
  int enumerate(Tab< E > &tab);
  void clear( ) { s.clear( ); }
};

template<>
inline int SortedList< RZ >::enumerate(Tab< RZ > &tab) {
  set< RZ >::const_iterator it;
  int i;

  for (it = s.begin( ), i = 0; it != s.end( ); ++it, ++i) {
    tab[i] = *it;
  }

  return i;
}

// old implementation based on personal construction of red black trees. Works fine, but the
// standard library might be safer and/or faster (?)
/*
 * template <class TreeLabel>
 * class SortedList {
 * public:
 *  const int & Card;
 *
 *  SortedList() : newNode(0), deleteNode(0), card(0), Card(card), cardMax(startCard), growIndex(1){
 *      nodes[0]    = new RBTree<TreeLabel>    [2*cardMax];
 *      pNodes      = new RBTree<TreeLabel> *  [2*cardMax];
 *      for(int i=0; i<2*cardMax; i++) {pNodes[i] = nodes[0]+i;}
 *  }
 *  SortedList(const SortedList & list){
 *      cout << "Sorted list constructor warning : copying list of cardinal " << list.Card << endl;
 *      Tab<TreeLabel> tab; list.enumerate(tab); for(int i=0; i<list.Card; i++) insert(tab[i]);
 *  }
 *  ~SortedList(){
 *      for(int i=0; i<growIndex; i++) delete nodes[i];
 *      delete pNodes;
 *  }
 *  bool contains(TreeLabel m){return root.contains(m);}
 *  TreeLabel min(){return root.min();}
 *  TreeLabel max(){return root.max();}
 *  bool insert(TreeLabel m){
 *      if(root.insert(m, pNodes[newNode], pNodes[newNode+1])){
 *          card++; newNode = (newNode+2)%(2*cardMax);
 *          if(card==cardMax) grow();
 *          return true;}
 *      return false;
 *  }
 *  bool remove(TreeLabel m){
 *      if(root.remove(m, pNodes[deleteNode], pNodes[deleteNode+1])){card--; deleteNode =
 * (deleteNode+2)%(2*cardMax); return true;} return false;
 *  }
 *  TreeLabel pop();
 *  void print(ostream& f) const {f << "{ "; root.printInOrder(f); f << "}";}
 *  void printTree(ostream& f) const {f << root;}
 *  int enumerate(Tab<TreeLabel> & tab){int counter=0; root.enumerate(tab, counter); return
 * counter;} void clear(){while(Card>0) pop();} private: const static int startCard = 64; const
 * static int growMax   = 25; int growIndex;
 *
 *  RBTree<TreeLabel> root;
 *  RBTree<TreeLabel> * nodes[growMax];
 *  RBTree<TreeLabel> * * pNodes;
 *  int newNode;
 *  int deleteNode;
 *
 *  int card;   //attention : le nombre de noeuds est le double du cardinal
 *  int cardMax;
 *
 *  bool grow(){
 *      //lorsque grow est appelé, on doit avoir card==cardMax, et donc newNode == deleteNode
 *      if(growIndex == growMax) return false;
 *      nodes[growIndex] = new RBTree<TreeLabel> [2*cardMax]; //création des nouveaux noeuds. même
 * nombre que tous ceux créés jusque alors cardMax*=2; delete pNodes; pNodes = new RBTree<TreeLabel>
 * * [2*cardMax]; //création des nouveaux pointeurs, autant que de noeuds au total for(int i=0;
 * i<cardMax; i++) pNodes[i] = nodes[growIndex]+i; newNode = 0;    deleteNode = cardMax;
 * growIndex++; return true;
 *  }
 * };
 *
 *
 * template<class TreeLabel>
 * TreeLabel SortedList<TreeLabel>::pop(){
 *  TreeLabel ans;
 *  if(card<=0) return root.ELEMENT_MAX;
 *  ans = root.pop(pNodes[deleteNode], pNodes[deleteNode+1]);
 *  card--; deleteNode = (deleteNode+2)%(2*cardMax);
 *  return ans;
 * }
 *
 * template<class TreeLabel> ostream& operator <<(ostream& f, const SortedList<TreeLabel> & list){
 *  list.print(f); return f;}
 * template<class TreeLabel> ostream_math operator <<(ostream_math f, const SortedList<TreeLabel> &
 * list){ if(f.format==Mathematica) {Tab<TreeLabel> tab; list.enumerate(tab); f<<tab;} else f.os <<
 * list; return f;}
 *
 * // ********************** RBTree *********************
 *
 * // Note : si le besoin s'en fait sentir, il est envisageable de diviser par 2 l'occupation
 * mémoire,
 * // en faisant porter un noeud non trivial aux feuilles.
 *
 * enum RBL {Red, Black, Leaf};
 * enum AP {Absorbed, Propagated};
 * enum APb {Absorbed_true, Absorbed_false, Propagated_true};
 *
 * inline bool APb2bool(APb a){return a!=Absorbed_false;}
 * inline AP APb2AP(APb a){return a==Propagated_true ? Propagated : Absorbed;}
 * inline APb AP_b2APb(AP a, bool b){return a==Propagated ? Propagated_true : (b==true ?
 * Absorbed_true : Absorbed_false);}
 *
 * template <class TreeLabel> class RBTree {
 * public:
 *  RBTree () :color(Leaf),n(),left(NULL),right(NULL){}
 *  bool contains (TreeLabel m);
 *  bool checkColor();
 *  TreeLabel min();
 *  TreeLabel max();
 *  int black_height();
 *  bool check();
 *  bool insert(TreeLabel m,    RBTree<TreeLabel> * leftTree,       RBTree<TreeLabel> * rightTree);
 *  bool remove(TreeLabel m,    RBTree<TreeLabel> * & leftTree,     RBTree<TreeLabel> * &
 * rightTree); TreeLabel pop(              RBTree<TreeLabel> * & leftTree,     RBTree<TreeLabel> * &
 * rightTree); //renvoie le plus petit élément, et le supprime void print(ostream& f) const; void
 * printInOrder(ostream& f) const; bool isEmpty(){return color==Leaf;} void enumerate(Tab<TreeLabel>
 * & tab, int & counter); void reset(); private: friend class SortedList<TreeLabel>; RBL color;
 *  TreeLabel n;
 *  RBTree * left;
 *  RBTree * right;
 *
 *  void conflict();
 *  bool rec_insert(TreeLabel m, RBTree<TreeLabel> * leftTree, RBTree<TreeLabel> * rightTree);
 *  APb rec_remove(TreeLabel m, RBTree<TreeLabel> * & leftTree, RBTree<TreeLabel> * & rightTree);
 *  AP unbalanced_right();
 *  AP unbalanced_left();
 *
 *  static const TreeLabel ELEMENT_MIN;
 *  static const TreeLabel ELEMENT_MAX;
 * };
 *
 * template<class TreeLabel>
 * inline ostream& operator <<(ostream& f, const RBTree<TreeLabel> & tree ){
 *  tree.print(f);
 *  return f;
 * }
 *
 * template <class TreeLabel>
 * void RBTree<TreeLabel>::print(std::ostream& f) const {
 *  if(color==Leaf) {f << "Leaf"; return;}
 *  f << n << " ";
 *  if(color==Red) {f << "Red";} else f << "Black";
 *  f << "( " << *left << ", " << *right << ")";
 * }
 *
 * template <class TreeLabel>
 * void RBTree<TreeLabel>::printInOrder(std::ostream& f) const {
 *  if(color==Leaf) return;
 *  left->printInOrder(f);
 *  if(!left->isEmpty()) f << ", ";
 *  f << n;
 *  if(!right->isEmpty()) f << ", ";
 *  right->printInOrder(f);
 * }
 *
 * template <class TreeLabel>
 * bool RBTree<TreeLabel>::contains (TreeLabel m){return n==m || (color != Leaf && (m<n ?
 * left->contains(m) : right->contains(m)) ); }
 *
 * template <class TreeLabel>
 * bool RBTree<TreeLabel>::checkColor(){
 *  return color==Leaf || (left->checkColor() && right->checkColor() && (color != Red ||
 * (left->color != Red && right->color != Red)));
 * }
 *
 * template <class TreeLabel>
 * TreeLabel RBTree<TreeLabel>::min()  {return color==Leaf ? ELEMENT_MAX : (left->color==Leaf ? n :
 * left->min() );}
 *
 * template <class TreeLabel>
 * TreeLabel RBTree<TreeLabel>::max()  {return color==Leaf ? ELEMENT_MIN : (right->color==Leaf ? n :
 * right->max() );}
 *
 * template <class TreeLabel>
 * int RBTree<TreeLabel>::black_height(){
 *  if(color==Leaf) return 0;
 *  int lh = left->black_height();
 *  int rh = right->black_height();
 *  return lh==rh ? (color == Black) + lh : INT_MIN;
 * }
 * // par construction, la hauteur noire à gauche et à droite doivent être égales. En cas de
 * différence, la valeur reçue est négative.
 *
 * template <class TreeLabel>
 * bool RBTree<TreeLabel>::check(){return checkColor() && black_height() >= 0;}
 *
 *
 * template <class TreeLabel>
 * void RBTree<TreeLabel>::conflict(){
 *  if(color!=Black) return;
 *
 *  if(left->color==Red && left->left->color==Red){
 *      RBTree * OldLeft = left;
 *      RBTree * OldLeftLeft = left->left;
 *      RBTree * t1 = OldLeftLeft->left;
 *      RBTree * t2 = OldLeftLeft->right;
 *      RBTree * t3 = OldLeft->right;
 *      RBTree * t4 = right;
 *
 *      const TreeLabel n1 = OldLeftLeft->n;
 *      const TreeLabel n2 = OldLeft->n;
 *      const TreeLabel n3 = n;
 *
 *      color = Red;
 *      left  = OldLeftLeft;    left->color = Black;
 *      right = OldLeft;        right->color = Black;
 *
 *      left->left = t1;      left->right = t2;
 *      right->left = t3;     right->right = t4;
 *
 *      left->n   = n1;
 *      n         = n2;
 *      right->n  = n3;
 *      return;
 *  }
 *  if(left->color == Red && left->right->color == Red){
 *      //RBTree * OldLeft = left;
 *      RBTree * OldLeftRight = left->right;
 *      //RBTree * t1 = OldLeft->left;
 *      RBTree * t2 = OldLeftRight->left;
 *      RBTree * t3 = OldLeftRight->right;
 *      RBTree * t4 = right;
 *
 *      //const TreeLabel n1 = OldLeft->n;
 *      const TreeLabel n2 = OldLeftRight->n;
 *      const TreeLabel n3 = n;
 *
 *      color = Red;
 *      //left = OldLeft;       //inchangé
 *      left->color = Black;
 *      right = OldLeftRight;
 *      right->color = Black;
 *
 *      //left->left = t1;    //inchangé
 *      left->right  = t2;
 *      right->left  = t3;
 *      right->right = t4;
 *
 *      //left->n   = n1;
 *      n           = n2;
 *      right->n  = n3;
 *      return;
 *  }
 *  if(right->color == Red && right->right->color == Red){
 *      RBTree * OldRight = right;
 *      RBTree * OldRightRight = right->right;
 *      RBTree * t1 = left;
 *      RBTree * t2 = OldRight->left;
 *      RBTree * t3 = OldRightRight->left;
 *      RBTree * t4 = OldRightRight->right;
 *
 *      const TreeLabel n1 = n;
 *      const TreeLabel n2 = OldRight->n;
 *      const TreeLabel n3 = OldRightRight->n;
 *
 *      color = Red;
 *      left    = OldRight;         left->color = Black;
 *      right   = OldRightRight;    right->color = Black;
 *
 *      left->left = t1;      left->right = t2;
 *      right->left = t3;     right->right = t4;
 *
 *      left->n   = n1;
 *      n           = n2;
 *      right->n  = n3;
 *      return;
 *  }
 *  if(right->color == Red && right->left->color == Red){
 *      //RBTree * OldRight = right;
 *      RBTree * OldRightLeft = right->left;
 *      RBTree * t1 = left;
 *      RBTree * t2 = OldRightLeft->left;
 *      RBTree * t3 = OldRightLeft->right;
 *      //RBTree * t4 = OldRight->right;
 *
 *      const TreeLabel n1 = n;
 *      const TreeLabel n2 = OldRightLeft->n;
 *      //            const TreeLabel n3 = OldRight->n;
 *
 *      color = Red;
 *      left = OldRightLeft;
 *      left->color = Black;
 *      //right = OldRight;      //inchangé
 *      right->color = Black;
 *
 *      left->left = t1;      left->right = t2;
 *      right->left = t3;     //right->right = t4; //inchangé
 *
 *      left->n   = n1;
 *      n         = n2;
 *      //right->n  = n3;
 *      return;
 *  }
 * }
 *
 * template <class TreeLabel>
 * bool RBTree<TreeLabel>::rec_insert(TreeLabel m, RBTree<TreeLabel> * leftTree, RBTree<TreeLabel> *
 * rightTree){ if(color==Leaf){color=Red; n=m; left = leftTree; right = rightTree; return true;}
 *  if(n==m) return false;
 *  const bool ans = m < n ? left->rec_insert(m, leftTree, rightTree) : right->rec_insert(m,
 * leftTree, rightTree); conflict(); return ans;
 * }
 *
 * template <class TreeLabel>
 * bool RBTree<TreeLabel>::insert(TreeLabel m, RBTree<TreeLabel> * leftTree, RBTree<TreeLabel> *
 * rightTree){ const bool ans = rec_insert(m, leftTree, rightTree); color = Black; return ans;
 * }
 *
 * template <class TreeLabel>
 * AP RBTree<TreeLabel>::unbalanced_right(){ //branche de droite plus légère que celle de gauche
 * (suite à une délétion) if(color == Red && left->color == Black){ color = Black; left->color =
 * Red; conflict(); return Absorbed;
 *  }
 *  if(color == Black && left->color == Red){
 *      RBTree * OldLeft = left;
 *      RBTree * t1 = left->left;
 *      RBTree * t2 = left->right;
 *      RBTree * t3 = right;
 *      const TreeLabel n1 = left->n;
 *      const TreeLabel n2 = n;
 *
 *      left = t1;
 *      right = OldLeft;
 *      right->left = t2;
 *      right->right = t3;
 *
 *      right->color = Black;
 *      t2->color = Red;  //par construction t2->color == Black initialement
 *
 *      n = n1; right->n = n2;
 *      right->conflict();
 *      return Absorbed;
 *
 *  }
 *  if(color == Black && left->color == Black){
 *      left->color = Red;
 *      conflict();
 *      return Propagated; //l'arbre a été allégé, le déséquilibre est propagé
 *  }
 *  return Absorbed; // on n'est pas censé en arriver là
 * }
 *
 * template <class TreeLabel>
 * AP RBTree<TreeLabel>::unbalanced_left(){ //branche de gauche plus légère que celle de droite
 * (suite à une délétion) if(color == Red && right->color == Black){ color = Black;
 *      right->color=Red;
 *      conflict();
 *      return Absorbed;
 *  }
 *  if(color == Black && right->color == Red){
 *      RBTree * OldRight = right;
 *      RBTree * t1 = left;
 *      RBTree * t2 = right->left;
 *      RBTree * t3 = right->right;
 *      const TreeLabel n1 = n;
 *      const TreeLabel n2 = right->n;
 *
 *      left = OldRight;
 *      left->left = t1;
 *      left->right = t2;
 *      right = t3;
 *
 *      left->color = Black;
 *      t2->color = Red;  //par construction t2->color == Black initialement
 *
 *      left->n = n1; n = n2;
 *      left->conflict();
 *      return Absorbed;
 *
 *  }
 *  if(color == Black && right->color == Black){
 *      right->color = Red;
 *      conflict();
 *      return Propagated; //l'arbre a été allégé, le déséquilibre est propagé
 *  }
 *  return Absorbed; // on n'est pas censé en arriver là
 * }
 *
 * template <class TreeLabel>
 * APb RBTree<TreeLabel>::rec_remove(TreeLabel m, RBTree<TreeLabel> * & leftTree, RBTree<TreeLabel>
 * * & rightTree){ if(color == Leaf) return Absorbed_false; //rien ne se passe if(m<n){ const APb
 * ans = left->rec_remove(m, leftTree, rightTree); conflict(); if(ans != Propagated_true) return
 * ans; return AP_b2APb(unbalanced_left(), APb2bool(ans));
 *  }
 *  if(m>n){
 *      const APb ans = right->rec_remove(m, leftTree, rightTree);
 *      conflict();
 *      if(ans != Propagated_true) return ans;
 *      return AP_b2APb(unbalanced_right(), APb2bool(ans));
 *  }
 *  //cas m==n
 *  if(left->color != Leaf){
 *      const TreeLabel maxLeft = left->max();
 *      n = maxLeft;
 *      const APb ans = left->rec_remove(maxLeft, leftTree, rightTree);
 *      conflict();
 *      if(ans != Propagated_true) return ans;
 *      return AP_b2APb(unbalanced_left(), APb2bool(ans));
 *  }
 *  if(right->color!=Leaf){
 *      const TreeLabel minRight = right->min();
 *      n = minRight;
 *      const APb ans = right->rec_remove(minRight, leftTree, rightTree);
 *      conflict();
 *      if(ans != Propagated_true) return ans;
 *      return AP_b2APb(unbalanced_right(), APb2bool(ans));
 *  }
 *  RBL OldColor = color;
 *  leftTree = left;
 *  rightTree = right;
 *  reset();
 *  if(OldColor == Black) return Propagated_true;
 *  return Absorbed_true;
 * }
 *
 * template<class TreeLabel> void RBTree<TreeLabel>::reset() {color = Leaf; n=ELEMENT_MIN; left =
 * NULL; right = NULL;}
 *
 * template <class TreeLabel>
 * bool RBTree<TreeLabel>::remove(TreeLabel m, RBTree<TreeLabel> * & leftTree, RBTree<TreeLabel> * &
 * rightTree){ const APb ans = rec_remove(m, leftTree, rightTree); if(color == Red) color = Black;
 *  return APb2bool(ans);
 * }
 *
 * template <class TreeLabel>
 * TreeLabel RBTree<TreeLabel>::pop(RBTree<TreeLabel> * & leftTree, RBTree<TreeLabel> * &
 * rightTree){ const TreeLabel m=min(); remove(m, leftTree, rightTree); return m;
 * }
 *
 * template <class TreeLabel>
 * void RBTree<TreeLabel>::enumerate(Tab<TreeLabel> & tab, int & counter){
 *  if(color==Leaf) return;
 *  left->enumerate(tab, counter);
 *  tab[counter++] = n; //RZ(n.distance, n.number);
 *  right->enumerate(tab,counter);
 * }
 *
 * template<>                const double    RBTree<double>   ::ELEMENT_MIN = DBL_MIN;
 * template<>                const int       RBTree<int>      ::ELEMENT_MIN = INT_MIN;
 * template<class TreeLabel> const TreeLabel RBTree<TreeLabel>::ELEMENT_MIN = TreeLabel().MIN();
 * template<>                const double    RBTree<double>   ::ELEMENT_MAX = DBL_MAX;
 * template<>                const int       RBTree<int>      ::ELEMENT_MAX = INT_MAX;
 * template<class TreeLabel> const TreeLabel RBTree<TreeLabel>::ELEMENT_MAX = TreeLabel().MAX();
 */

#endif
