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
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

#ifndef RENUMB_HPP
#define RENUMB_HPP

// Licensing:
//
// This code is distributed under the GNU LGPL license.
//
// Modified:
//
// 08 March 2013
//
// Author:
//
// John Burkardt
//

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstring>

namespace renumb {
  void i4vec_print(int n, int a[], string title);
  void adj_print(int node_num, int adj_num, int adj_row[], int adj[], string title);
  void adj_print_some(int node_num, int node_lo, int node_hi, int adj_num, int adj_row[], int adj[],
                      string title);
  int adj_bandwidth(int node_num, int adj_num, int adj_row[], int adj[]);
  int adj_perm_bandwidth(int node_num, int adj_num, int adj_row[], int adj[], int perm[],
                         int perm_inv[]);
  int *genrcm(int node_num, int adj_num, int adj_row[], int adj[]);
  void rcm(int root, int adj_num, int adj_row[], int adj[], int mask[], int perm[], int *iccsze,
           int node_num);
  void root_find(int *root, int adj_num, int adj_row[], int adj[], int mask[], int *level_num,
                 int level_row[], int level[], int node_num);
  void i4vec_reverse(int n, int a[]);
  void level_set(int root, int adj_num, int adj_row[], int adj[], int mask[], int *level_num,
                 int level_row[], int level[], int node_num);
  void degree(int root, int adj_num, int adj_row[], int adj[], int mask[], int deg[], int *iccsze,
              int ls[], int node_num);
  int *perm_inverse3(int n, int perm[]);

  void i4vec_print(int n, int a[], string title) {
    // ****************************************************************************80
    //
    // Purpose:
    //
    // I4VEC_PRINT prints an I4VEC.
    //
    // Licensing:
    //
    // This code is distributed under the GNU LGPL license.
    //
    // Modified:
    //
    // 25 February 2003
    //
    // Author:
    //
    // John Burkardt
    //
    // Parameters:
    //
    // Input, int N, the number of components of the vector.
    //
    // Input, int A[N], the vector to be printed.
    //
    // Input, string TITLE, a title.
    //
    int i;

    cout << "\n";
    cout << title << "\n";
    cout << "\n";

    for (i = 0; i < n; i++) {
      cout << "  " << setw(8) << i << "  " << setw(8) << a[i] << "\n";
    }

    return;
  }

  void adj_print(int node_num, int adj_num, int adj_row[], int adj[], string title) {
    // ****************************************************************************80
    //
    // Purpose:
    //
    // ADJ_PRINT prints adjacency information.
    //
    // Discussion:
    //
    // The list has the form:
    //
    // Row   Nonzeros
    //
    // 1       2   5   9
    // 2       7   8   9   15   78   79   81  86  91  99
    // 100 103
    // 3      48  49  53
    //
    // Licensing:
    //
    // This code is distributed under the GNU LGPL license.
    //
    // Modified:
    //
    // 18 December 2002
    //
    // Author:
    //
    // John Burkardt
    //
    // Parameters:
    //
    // Input, int NODE_NUM, the number of nodes.
    //
    // Input, int ADJ_NUM, the number of adjacency entries.
    //
    // Input, int ADJ_ROW[NODE_NUM+1], organizes the adjacency entries
    // into rows.  The entries for row I are in entries ADJ_ROW(I)
    // through ADJ_ROW(I+1)-1.
    //
    // Input, int ADJ[ADJ_NUM], the adjacency structure, which contains,
    // for each row, the column indices of the nonzero entries.
    //
    // Input, string TITLE, a title.
    //
    adj_print_some(node_num, 0, node_num - 1, adj_num, adj_row, adj, title);

    return;
  }

  int adj_bandwidth(int node_num, int adj_num, int adj_row[], int adj[]) {
    // ****************************************************************************80
    //
    // Purpose:
    //
    // ADJ_BANDWIDTH computes the bandwidth of an adjacency matrix.
    //
    // Licensing:
    //
    // This code is distributed under the GNU LGPL license.
    //
    // Modified:
    //
    // 03 January 2007
    //
    // Author:
    //
    // John Burkardt
    //
    // Reference:
    //
    // Alan George, Joseph Liu,
    // Computer Solution of Large Sparse Positive Definite Systems,
    // Prentice Hall, 1981.
    //
    // Parameters:
    //
    // Input, int NODE_NUM, the number of nodes.
    //
    // Input, int ADJ_NUM, the number of adjacency entries.
    //
    // Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
    // in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
    //
    // Input, int ADJ[ADJ_NUM], the adjacency structure.
    // For each row, it contains the column indices of the nonzero entries.
    //
    // Output, int ADJ_BANDWIDTH, the bandwidth of the adjacency
    // matrix.
    //
    int band_hi;
    int band_lo;
    int col;
    int i;
    int j;
    int value;

    band_lo = 0;
    band_hi = 0;

    for (i = 0; i < node_num; i++) {
      for (j = adj_row[i]; j <= adj_row[i + 1] - 1; j++) {
        col = adj[j];
        band_lo = std::max(band_lo, i - col);
        band_hi = std::max(band_hi, col - i);
      }
    }

    value = band_lo + 1 + band_hi;

    return value;
  }

  int adj_perm_bandwidth(int node_num, int adj_num, int adj_row[], int adj[], int perm[],
                         int perm_inv[]) {
    // ****************************************************************************80
    //
    // Purpose:
    //
    // ADJ_PERM_BANDWIDTH computes the bandwidth of a permuted adjacency matrix.
    //
    // Discussion:
    //
    // The matrix is defined by the adjacency information and a permutation.
    //
    // The routine also computes the bandwidth and the size of the envelope.
    //
    // Licensing:
    //
    // This code is distributed under the GNU LGPL license.
    //
    // Modified:
    //
    // 05 January 2007
    //
    // Author:
    //
    // John Burkardt
    //
    // Reference:
    //
    // Alan George, Joseph Liu,
    // Computer Solution of Large Sparse Positive Definite Systems,
    // Prentice Hall, 1981.
    //
    // Parameters:
    //
    // Input, int NODE_NUM, the number of nodes.
    //
    // Input, int ADJ_NUM, the number of adjacency entries.
    //
    // Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
    // in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
    //
    // Input, int ADJ[ADJ_NUM], the adjacency structure.
    // For each row, it contains the column indices of the nonzero entries.
    //
    // Input, int PERM[NODE_NUM], PERM_INV(NODE_NUM), the permutation
    // and inverse permutation.
    //
    // Output, int ADJ_PERM_BANDWIDTH, the bandwidth of the permuted
    // adjacency matrix.
    //
    int band_hi;
    int band_lo;
    int bandwidth;
    int col;
    int i;
    int j;

    band_lo = 0;
    band_hi = 0;

    for (i = 0; i < node_num; i++) {
      for (j = adj_row[perm[i]]; j <= adj_row[perm[i] + 1] - 1; j++) {
        col = perm_inv[adj[j]];
        band_lo = std::max(band_lo, i - col);
        band_hi = std::max(band_hi, col - i);
      }
    }

    bandwidth = band_lo + 1 + band_hi;

    return bandwidth;
  }

  int *genrcm(int node_num, int adj_num, int adj_row[], int adj[]) {
    // ****************************************************************************80
    //
    // Purpose:
    //
    // GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
    //
    // Discussion:
    //
    // For each connected component in the graph, the routine obtains
    // an ordering by calling RCM.
    //
    // Licensing:
    //
    // This code is distributed under the GNU LGPL license.
    //
    // Modified:
    //
    // 14 May 2011
    //
    // Author:
    //
    // Original FORTRAN77 version by Alan George, Joseph Liu.
    // C++ version by John Burkardt.
    //
    // Reference:
    //
    // Alan George, Joseph Liu,
    // Computer Solution of Large Sparse Positive Definite Systems,
    // Prentice Hall, 1981.
    //
    // Parameters:
    //
    // Input, int NODE_NUM, the number of nodes.
    //
    // Input, int ADJ_NUM, the number of adjacency entries.
    //
    // Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
    // in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
    //
    // Input, int  ADJ[ADJ_NUM], the adjacency structure.
    // For each row, it contains the column indices of the nonzero entries.
    //
    // Output, int GENRCM[NODE_NUM], the RCM ordering.
    //
    // Local Parameters:
    //
    // Local, int  LEVEL_ROW[NODE_NUM+1], the index vector for a level
    // structure.  The level structure is stored in the currently unused
    // spaces in the permutation vector PERM.
    //
    // Local, int MASK[NODE_NUM], marks variables that have been numbered.
    //
    int i;
    int iccsze;
    int level_num;
    int *level_row;
    int *mask;
    int num;
    int *perm;
    int root;

    //
    // Assuming the input dat is 0 based, add 1 to ADJ_ROW and ADJ,
    // because GENRCM uses 1-based indexing!
    //
    for (i = 0; i < node_num + 1; i++) {
      adj_row[i] = adj_row[i] + 1;
    }

    for (i = 0; i < adj_num; i++) {
      adj[i] = adj[i] + 1;
    }

    perm = new int[node_num];
    level_row = new int[node_num + 1];
    mask = new int[node_num];

    for (i = 0; i < node_num; i++) {
      mask[i] = 1;
    }

    num = 1;

    for (i = 0; i < node_num; i++) {
      //
      // For each masked connected component...
      //
      if (mask[i] != 0) {
        root = i + 1;
        //
        // Find a pseudo-peripheral node ROOT.  The level structure found by
        // ROOT_FIND is stored starting at PERM(NUM).
        //
        root_find(&root, adj_num, adj_row, adj, mask, &level_num, level_row, perm + num - 1,
                  node_num);
        //
        // RCM orders the component using ROOT as the starting node.
        //
        rcm(root, adj_num, adj_row, adj, mask, perm + num - 1, &iccsze, node_num);

        num = num + iccsze;
      }

      //
      // We can stop once every node is in one of the connected components.
      //
      if (node_num < num) {
        break;
      }
    }

    delete[] level_row;
    delete[] mask;

    //
    // PERM is computed as a 1-based vector.
    // Rewrite it as a 0-based vector.
    //
    for (i = 0; i < node_num; i++) {
      perm[i] = perm[i] - 1;
    }

    //
    // Subtract 1 from ADJ_ROW and ADJ because GENRCM used 1-based indexing!
    //
    for (i = 0; i < node_num + 1; i++) {
      adj_row[i] = adj_row[i] - 1;
    }

    for (i = 0; i < adj_num; i++) {
      adj[i] = adj[i] - 1;
    }

    return perm;
  }

  void adj_print_some(int node_num, int node_lo, int node_hi, int adj_num, int adj_row[], int adj[],
                      string title) {
    // ****************************************************************************80
    //
    // Purpose:
    //
    // ADJ_PRINT_SOME prints some adjacency information.
    //
    // Discussion:
    //
    // The list has the form:
    //
    // Row   Nonzeros
    //
    // 1       2   5   9
    // 2       7   8   9   15   78   79   81  86  91  99
    // 100 103
    // 3      48  49  53
    //
    // Licensing:
    //
    // This code is distributed under the GNU LGPL license.
    //
    // Modified:
    //
    // 04 January 2007
    //
    // Author:
    //
    // John Burkardt
    //
    // Parameters:
    //
    // Input, int NODE_NUM, the number of nodes.
    //
    // Input, int NODE_LO, NODE_HI, the first and last nodes for
    // which the adjacency information is to be printed.
    //
    // Input, int ADJ_NUM, the number of adjacency entries.
    //
    // Input, int ADJ_ROW[NODE_NUM+1], organizes the adjacency entries
    // into rows.  The entries for row I are in entries ADJ_ROW(I)
    // through ADJ_ROW(I+1)-1.
    //
    // Input, int ADJ[ADJ_NUM], the adjacency structure, which contains,
    // for each row, the column indices of the nonzero entries.
    //
    // Input, string TITLE, a title to be printed.
    //
    int i;
    int j;
    int jhi;
    int jlo;
    int jmax;
    int jmin;

    cout << "\n";
    cout << title << "\n";
    cout << "  Sparse adjacency structure:\n";
    cout << "\n";
    cout << "  Number of nodes       = " << node_num << "\n";
    ;
    cout << "  Number of adjacencies = " << adj_num << "\n";
    cout << "\n";
    cout << "  Node   Min   Max          Nonzeros \n";
    cout << "\n";

    for (i = node_lo; i <= node_hi; i++) {
      jmin = adj_row[i];
      jmax = adj_row[i + 1] - 1;

      if (jmax < jmin) {
        cout << "  " << setw(4) << i << "  " << setw(4) << jmin << "  " << setw(4) << jmax << "\n";
      } else {
        for (jlo = jmin; jlo <= jmax; jlo = jlo + 5) {
          jhi = std::min(jlo + 4, jmax);

          if (jlo == jmin) {
            cout << "  " << setw(4) << i << "  " << setw(4) << jmin << "  " << setw(4) << jmax
                 << "   ";

            for (j = jlo; j <= jhi; j++) {
              cout << setw(8) << adj[j];
            }

            cout << "\n";
          } else {
            cout << "                     ";

            for (j = jlo; j <= jhi; j++) {
              cout << setw(8) << adj[j];
            }

            cout << "\n";
          }
        }
      }
    }

    return;
  }

  void rcm(int root, int adj_num, int adj_row[], int adj[], int mask[], int perm[], int *iccsze,
           int node_num) {
    // ****************************************************************************80
    //
    // Purpose:
    //
    // RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
    //
    // Discussion:
    //
    // The connected component is specified by a node ROOT and a mask.
    // The numbering starts at the root node.
    //
    // An outline of the algorithm is as follows:
    //
    // X(1) = ROOT.
    //
    // for ( I = 1 to N-1)
    // Find all unlabeled neighbors of X(I),
    // assign them the next available labels, in order of increasing degree.
    //
    // When done, reverse the ordering.
    //
    // Licensing:
    //
    // This code is distributed under the GNU LGPL license.
    //
    // Modified:
    //
    // 05 January 2007
    //
    // Author:
    //
    // Original FORTRAN77 version by Alan George, Joseph Liu.
    // C++ version by John Burkardt.
    //
    // Reference:
    //
    // Alan George, Joseph Liu,
    // Computer Solution of Large Sparse Positive Definite Systems,
    // Prentice Hall, 1981.
    //
    // Parameters:
    //
    // Input, int ROOT, the node that defines the connected component.
    // It is used as the starting point for the RCM ordering.
    //
    // Input, int ADJ_NUM, the number of adjacency entries.
    //
    // Input, int ADJ_ROW(NODE_NUM+1).  Information about row I is stored
    // in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
    //
    // Input, int ADJ(ADJ_NUM), the adjacency structure.
    // For each row, it contains the column indices of the nonzero entries.
    //
    // Input/output, int MASK(NODE_NUM), a mask for the nodes.  Only
    // those nodes with nonzero input mask values are considered by the
    // routine.  The nodes numbered by RCM will have their mask values
    // set to zero.
    //
    // Output, int PERM(NODE_NUM), the RCM ordering.
    //
    // Output, int ICCSZE, the size of the connected component
    // that has been numbered.
    //
    // Input, int NODE_NUM, the number of nodes.
    //
    // Local Parameters:
    //
    // Workspace, int DEG[NODE_NUM], a temporary vector used to hold
    // the degree of the nodes in the section graph specified by mask and root.
    //
    int *deg;
    int fnbr;
    int i;
    int j;
    int jstop;
    int jstrt;
    int k;
    int l;
    int lbegin;
    int lnbr;
    int lperm;
    int lvlend;
    int nbr;
    int node;

    //
    // Find the degrees of the nodes in the component specified by MASK and ROOT.
    //
    deg = new int[node_num];

    degree(root, adj_num, adj_row, adj, mask, deg, iccsze, perm, node_num);

    mask[root - 1] = 0;

    if (*iccsze <= 1) {
      delete[] deg;
      return;
    }

    lvlend = 0;
    lnbr = 1;

    //
    // LBEGIN and LVLEND point to the beginning and
    // the end of the current level respectively.
    //
    while (lvlend < lnbr) {
      lbegin = lvlend + 1;
      lvlend = lnbr;

      for (i = lbegin; i <= lvlend; i++) {
        //
        // For each node in the current level...
        //
        node = perm[i - 1];
        jstrt = adj_row[node - 1];
        jstop = adj_row[node] - 1;
        //
        // Find the unnumbered neighbors of NODE.
        //
        // FNBR and LNBR point to the first and last neighbors
        // of the current node in PERM.
        //
        fnbr = lnbr + 1;

        for (j = jstrt; j <= jstop; j++) {
          nbr = adj[j - 1];

          if (mask[nbr - 1] != 0) {
            lnbr = lnbr + 1;
            mask[nbr - 1] = 0;
            perm[lnbr - 1] = nbr;
          }
        }

        //
        // If no neighbors, skip to next node in this level.
        //
        if (lnbr <= fnbr) {
          continue;
        }

        //
        // Sort the neighbors of NODE in increasing order by degree.
        // Linear insertion is used.
        //
        k = fnbr;

        while (k < lnbr) {
          l = k;
          k = k + 1;
          nbr = perm[k - 1];

          while (fnbr < l) {
            lperm = perm[l - 1];

            if (deg[lperm - 1] <= deg[nbr - 1]) {
              break;
            }

            perm[l] = lperm;
            l = l - 1;
          }

          perm[l] = nbr;
        }
      }
    }

    //
    // We now have the Cuthill-McKee ordering.  Reverse it.
    //
    i4vec_reverse(*iccsze, perm);

    delete[] deg;

    return;
  }

  void root_find(int *root, int adj_num, int adj_row[], int adj[], int mask[], int *level_num,
                 int level_row[], int level[], int node_num) {
    // ****************************************************************************80
    //
    // Purpose:
    //
    // ROOT_FIND finds a pseudo-peripheral node.
    //
    // Discussion:
    //
    // The diameter of a graph is the maximum distance (number of edges)
    // between any two nodes of the graph.
    //
    // The eccentricity of a node is the maximum distance between that
    // node and any other node of the graph.
    //
    // A peripheral node is a node whose eccentricity equals the
    // diameter of the graph.
    //
    // A pseudo-peripheral node is an approximation to a peripheral node;
    // it may be a peripheral node, but all we know is that we tried our
    // best.
    //
    // The routine is given a graph, and seeks pseudo-peripheral nodes,
    // using a modified version of the scheme of Gibbs, Poole and
    // Stockmeyer.  It determines such a node for the section subgraph
    // specified by MASK and ROOT.
    //
    // The routine also determines the level structure associated with
    // the given pseudo-peripheral node; that is, how far each node
    // is from the pseudo-peripheral node.  The level structure is
    // returned as a list of nodes LS, and pointers to the beginning
    // of the list of nodes that are at a distance of 0, 1, 2, ...,
    // NODE_NUM-1 from the pseudo-peripheral node.
    //
    // Licensing:
    //
    // This code is distributed under the GNU LGPL license.
    //
    // Modified:
    //
    // 05 January 2007
    //
    // Author:
    //
    // Original FORTRAN77 version by Alan George, Joseph Liu.
    // C++ version by John Burkardt.
    //
    // Reference:
    //
    // Alan George, Joseph Liu,
    // Computer Solution of Large Sparse Positive Definite Systems,
    // Prentice Hall, 1981.
    //
    // Norman Gibbs, William Poole, Paul Stockmeyer,
    // An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix,
    // SIAM Journal on Numerical Analysis,
    // Volume 13, pages 236-250, 1976.
    //
    // Norman Gibbs,
    // Algorithm 509: A Hybrid Profile Reduction Algorithm,
    // ACM Transactions on Mathematical Software,
    // Volume 2, pages 378-387, 1976.
    //
    // Parameters:
    //
    // Input/output, int *ROOT.  On input, ROOT is a node in the
    // the component of the graph for which a pseudo-peripheral node is
    // sought.  On output, ROOT is the pseudo-peripheral node obtained.
    //
    // Input, int ADJ_NUM, the number of adjacency entries.
    //
    // Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
    // in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
    //
    // Input, int ADJ[ADJ_NUM], the adjacency structure.
    // For each row, it contains the column indices of the nonzero entries.
    //
    // Input, int MASK[NODE_NUM], specifies a section subgraph.  Nodes
    // for which MASK is zero are ignored by FNROOT.
    //
    // Output, int *LEVEL_NUM, is the number of levels in the level structure
    // rooted at the node ROOT.
    //
    // Output, int LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the
    // level structure array pair containing the level structure found.
    //
    // Input, int NODE_NUM, the number of nodes.
    //
    int iccsze;
    int j;
    int jstrt;
    int k;
    int kstop;
    int kstrt;
    int level_num2;
    int mindeg;
    int nabor;
    int ndeg;
    int node;

    //
    // Determine the level structure rooted at ROOT.
    //
    level_set(*root, adj_num, adj_row, adj, mask, level_num, level_row, level, node_num);
    //
    // Count the number of nodes in this level structure.
    //
    iccsze = level_row[*level_num] - 1;
    //
    // Extreme case:
    // A complete graph has a level set of only a single level.
    // Every node is equally good (or bad).
    //
    if (*level_num == 1) {
      return;
    }

    //
    // Extreme case:
    // A "line graph" 0--0--0--0--0 has every node in its only level.
    // By chance, we've stumbled on the ideal root.
    //
    if (*level_num == iccsze) {
      return;
    }

    //
    // Pick any node from the last level that has minimum degree
    // as the starting point to generate a new level set.
    //
    for (;;) {
      mindeg = iccsze;

      jstrt = level_row[*level_num - 1];
      *root = level[jstrt - 1];

      if (jstrt < iccsze) {
        for (j = jstrt; j <= iccsze; j++) {
          node = level[j - 1];
          ndeg = 0;
          kstrt = adj_row[node - 1];
          kstop = adj_row[node] - 1;

          for (k = kstrt; k <= kstop; k++) {
            nabor = adj[k - 1];
            if (0 < mask[nabor - 1]) {
              ndeg = ndeg + 1;
            }
          }

          if (ndeg < mindeg) {
            *root = node;
            mindeg = ndeg;
          }
        }
      }

      //
      // Generate the rooted level structure associated with this node.
      //
      level_set(*root, adj_num, adj_row, adj, mask, &level_num2, level_row, level, node_num);
      //
      // If the number of levels did not increase, accept the new ROOT.
      //
      if (level_num2 <= *level_num) {
        break;
      }

      *level_num = level_num2;
      //
      // In the unlikely case that ROOT is one endpoint of a line graph,
      // we can exit now.
      //
      if (iccsze <= *level_num) {
        break;
      }
    }

    return;
  }

  void i4vec_reverse(int n, int a[]) {
    // ****************************************************************************80
    //
    // Purpose:
    //
    // I4VEC_REVERSE reverses the elements of an I4VEC.
    //
    // Example:
    //
    // Input:
    //
    // N = 5,
    // A = ( 11, 12, 13, 14, 15 ).
    //
    // Output:
    //
    // A = ( 15, 14, 13, 12, 11 ).
    //
    // Licensing:
    //
    // This code is distributed under the GNU LGPL license.
    //
    // Modified:
    //
    // 22 September 2005
    //
    // Author:
    //
    // John Burkardt
    //
    // Parameters:
    //
    // Input, int N, the number of entries in the array.
    //
    // Input/output, int A(N), the array to be reversed.
    //
    int i;
    int j;

    for (i = 0; i < n / 2; i++) {
      j = a[i];
      a[i] = a[n - 1 - i];
      a[n - 1 - i] = j;
    }

    return;
  }

  void level_set(int root, int adj_num, int adj_row[], int adj[], int mask[], int *level_num,
                 int level_row[], int level[], int node_num) {
    // ****************************************************************************80
    //
    // Purpose:
    //
    // LEVEL_SET generates the connected level structure rooted at a given node.
    //
    // Discussion:
    //
    // Only nodes for which MASK is nonzero will be considered.
    //
    // The root node chosen by the user is assigned level 1, and masked.
    // All (unmasked) nodes reachable from a node in level 1 are
    // assigned level 2 and masked.  The process continues until there
    // are no unmasked nodes adjacent to any node in the current level.
    // The number of levels may vary between 2 and NODE_NUM.
    //
    // Licensing:
    //
    // This code is distributed under the GNU LGPL license.
    //
    // Modified:
    //
    // 05 January 2007
    //
    // Author:
    //
    // Original FORTRAN77 version by Alan George, Joseph Liu.
    // C++ version by John Burkardt.
    //
    // Reference:
    //
    // Alan George, Joseph Liu,
    // Computer Solution of Large Sparse Positive Definite Systems,
    // Prentice Hall, 1981.
    //
    // Parameters:
    //
    // Input, int ROOT, the node at which the level structure
    // is to be rooted.
    //
    // Input, int ADJ_NUM, the number of adjacency entries.
    //
    // Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
    // in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
    //
    // Input, int ADJ[ADJ_NUM], the adjacency structure.
    // For each row, it contains the column indices of the nonzero entries.
    //
    // Input/output, int MASK[NODE_NUM].  On input, only nodes with nonzero
    // MASK are to be processed.  On output, those nodes which were included
    // in the level set have MASK set to 1.
    //
    // Output, int *LEVEL_NUM, the number of levels in the level
    // structure.  ROOT is in level 1.  The neighbors of ROOT
    // are in level 2, and so on.
    //
    // Output, int LEVEL_ROW[NODE_NUM+1], LEVEL[NODE_NUM], the rooted
    // level structure.
    //
    // Input, int NODE_NUM, the number of nodes.
    //
    int i;
    int iccsze;
    int j;
    int jstop;
    int jstrt;
    int lbegin;
    int lvlend;
    int lvsize;
    int nbr;
    int node;

    mask[root - 1] = 0;
    level[0] = root;
    *level_num = 0;
    lvlend = 0;
    iccsze = 1;

    //
    // LBEGIN is the pointer to the beginning of the current level, and
    // LVLEND points to the end of this level.
    //
    for (;;) {
      lbegin = lvlend + 1;
      lvlend = iccsze;
      *level_num = *level_num + 1;
      level_row[*level_num - 1] = lbegin;

      //
      // Generate the next level by finding all the masked neighbors of nodes
      // in the current level.
      //
      for (i = lbegin; i <= lvlend; i++) {
        node = level[i - 1];
        jstrt = adj_row[node - 1];
        jstop = adj_row[node] - 1;

        for (j = jstrt; j <= jstop; j++) {
          nbr = adj[j - 1];

          if (mask[nbr - 1] != 0) {
            iccsze = iccsze + 1;
            level[iccsze - 1] = nbr;
            mask[nbr - 1] = 0;
          }
        }
      }

      //
      // Compute the current level width (the number of nodes encountered.)
      // If it is positive, generate the next level.
      //
      lvsize = iccsze - lvlend;

      if (lvsize <= 0) {
        break;
      }
    }

    level_row[*level_num] = lvlend + 1;

    //
    // Reset MASK to 1 for the nodes in the level structure.
    //
    for (i = 0; i < iccsze; i++) {
      mask[level[i] - 1] = 1;
    }

    return;
  }

  void degree(int root, int adj_num, int adj_row[], int adj[], int mask[], int deg[], int *iccsze,
              int ls[], int node_num) {
    // ****************************************************************************80
    //
    // Purpose:
    //
    // DEGREE computes the degrees of the nodes in the connected component.
    //
    // Discussion:
    //
    // The connected component is specified by MASK and ROOT.
    // Nodes for which MASK is zero are ignored.
    //
    // Licensing:
    //
    // This code is distributed under the GNU LGPL license.
    //
    // Modified:
    //
    // 05 January 2007
    //
    // Author:
    //
    // Original FORTRAN77 version by Alan George, Joseph Liu.
    // C++ version by John Burkardt.
    //
    // Reference:
    //
    // Alan George, Joseph Liu,
    // Computer Solution of Large Sparse Positive Definite Systems,
    // Prentice Hall, 1981.
    //
    // Parameters:
    //
    // Input, int ROOT, the node that defines the connected component.
    //
    // Input, int ADJ_NUM, the number of adjacency entries.
    //
    // Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
    // in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
    //
    // Input, int ADJ[ADJ_NUM], the adjacency structure.
    // For each row, it contains the column indices of the nonzero entries.
    //
    // Input, int MASK[NODE_NUM], is nonzero for those nodes which are
    // to be considered.
    //
    // Output, int DEG[NODE_NUM], contains, for each  node in the connected
    // component, its degree.
    //
    // Output, int *ICCSIZE, the number of nodes in the connected component.
    //
    // Output, int LS[NODE_NUM], stores in entries 1 through ICCSIZE the nodes
    // in the connected component, starting with ROOT, and proceeding
    // by levels.
    //
    // Input, int NODE_NUM, the number of nodes.
    //
    int i;
    int ideg;
    int j;
    int jstop;
    int jstrt;
    int lbegin;
    int lvlend;
    int lvsize;
    int nbr;
    int node;

    //
    // The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
    //
    ls[0] = root;
    adj_row[root - 1] = -adj_row[root - 1];
    lvlend = 0;
    *iccsze = 1;

    //
    // LBEGIN is the pointer to the beginning of the current level, and
    // LVLEND points to the end of this level.
    //
    for (;;) {
      lbegin = lvlend + 1;
      lvlend = *iccsze;

      //
      // Find the degrees of nodes in the current level,
      // and at the same time, generate the next level.
      //
      for (i = lbegin; i <= lvlend; i++) {
        node = ls[i - 1];
        jstrt = -adj_row[node - 1];
        jstop = abs(adj_row[node]) - 1;
        ideg = 0;

        for (j = jstrt; j <= jstop; j++) {
          nbr = adj[j - 1];

          if (mask[nbr - 1] != 0) {
            ideg = ideg + 1;

            if (0 <= adj_row[nbr - 1]) {
              adj_row[nbr - 1] = -adj_row[nbr - 1];
              *iccsze = *iccsze + 1;
              ls[*iccsze - 1] = nbr;
            }
          }
        }

        deg[node - 1] = ideg;
      }

      //
      // Compute the current level width.
      //
      lvsize = *iccsze - lvlend;
      //
      // If the current level width is nonzero, generate another level.
      //
      if (lvsize == 0) {
        break;
      }
    }

    //
    // Reset ADJ_ROW to its correct sign and return.
    //
    for (i = 0; i < *iccsze; i++) {
      node = ls[i] - 1;
      adj_row[node] = -adj_row[node];
    }

    return;
  }

  int *perm_inverse3(int n, int perm[]) {
    // ****************************************************************************80
    //
    // Purpose:
    //
    // PERM_INVERSE3 produces the inverse of a given permutation.
    //
    // Licensing:
    //
    // This code is distributed under the GNU LGPL license.
    //
    // Modified:
    //
    // 14 May 2011
    //
    // Author:
    //
    // John Burkardt
    //
    // Parameters:
    //
    // Input, int N, the number of items permuted.
    //
    // Input, int PERM[N], a permutation.
    //
    // Output, int PERM_INVERSE3[N], the inverse permutation.
    //
    int i;
    int *perm_inv;

    perm_inv = new int[n];

    for (i = 0; i < n; i++) {
      perm_inv[perm[i]] = i;
    }

    return perm_inv;
  }
}    // namespace renumb

#endif
