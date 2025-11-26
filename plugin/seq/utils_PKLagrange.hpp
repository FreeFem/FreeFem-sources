#include "ff++.hpp"
#include <cmath>
#include <iostream>
#include <map>
#include <vector>
using namespace std;

void NumInternPoint(int k, vector<R2> &Pt) {
    /*Numbering of interior points
   index i is associated with the coordinate lambda3 and j with lambda2
   therefore the coordinate for lambda1 is k − (i + j). We iterate only
   between 1 and k − 2 to avoid the edge nodes
   NB: We start at the vertex node, then move down to the next row,
   numbering from right to left, and continue this pattern for the following rows
   */
    for (int i = k - 2; i >= 1; i--) {
        for (int j = 1; j <= k - 2 - i + 1; j++) {
            Pt.push_back(R2(double(k - 1 - i - (j - 1)) / k, double(i) / k));
            // cout<<"("<<j<<","<<k-1-i-(j-1)<<","<<i<<")";
        }
    }
}

//Construction of point coordinates
vector<R2> PtConstruction(int k) {
    int NptPerV = k - 1;
    vector<R2> Pt; //(Ndof);
    // vertices
    Pt.push_back(R2(0, 0));
    Pt.push_back(R2(1., 0));
    Pt.push_back(R2(0, 1.));
    // bottom edge
    for (int i = 0; i < NptPerV; i++)
        Pt.push_back(R2(double(k - i - 1) / k, double(i + 1) / k));
    // right edge
    for (int i = 0; i < NptPerV; i++)
        Pt.push_back(R2(0 / k, double(k - i - 1) / k));
    // left edge
    for (int i = 0; i < NptPerV; i++)
        Pt.push_back(R2(double(i + 1) / k, 0. / k));
    // Interior of the element
    NumInternPoint(k, Pt);
    return Pt;
}

void FillDataLagrange(int k, vector<int> &Data) {
    /* Fill Data array in PKlagrange element (see P3 or P4
    for more details)
    */

    int ndof = (k + 1) * (k + 2) * 0.5;
    Data.resize(5 * ndof + 3, 0);

    // First row(the support number  of the node of the df)
    // Data[0]=0;
    Data[1] = 1;
    Data[2] = 2;
    // Second row(the number of the df on  the node)
    /*Data[ndof]=0;
    Data[ndof+1]=0;
    Data[ndof+2]=0;*/
    // 3rd row
    Data[2 * ndof] = 0;
    Data[2 * ndof + 1] = 1;
    Data[2 * ndof + 2] = 2;
    // 4th row
    /*Data[3*ndof]=0;
    Data[3*ndof+1]=0;
    Data[3*ndof+2]=0;*/
    // 5th row
    // Data[4*ndof]=0;
    Data[4 * ndof + 1] = 1;
    Data[4 * ndof + 2] = 2;
    for (int i = 3; i < ndof; i++) {
        if (i <= k + 1) {
            Data[i] = 3;
            Data[ndof + i] = i - 3;

        } else {
            if (i <= 2 * k) {
                Data[i] = 4;
                Data[ndof + i] = i - 3 - (k) + 1;
            } else if (i <= 3 * k - 1) {
                Data[i] = 5;
                Data[ndof + i] = i - 3 - 2 * (k - 1);
            } else {
                Data[i] = 6;
                Data[ndof + i] = i - 3 - 3 * (k - 1);
            }
        }
        // 3row
        Data[2 * ndof + i] = Data[i];
        // 4 th row
        // Data[3*ndof+i]=0;
        // 5th row
        Data[4 * ndof + i] = i;
    }
    Data[5 * ndof + 2] = ndof;
    /*Data[5*ndof+1]=0;
    Data[5*ndof]=0;*/
    // Pi_h_coef.resize(ndof,1.);
}

vector<pair<int, int>> ExchangeidxVector(int k) {
    /* Determine the necessary permutations
    for the internal points located on the edges of an element*/

    int nbrPerm = (k % 2) ? (k - 1) : (k - 2); // Number of permutation per axis
    vector<pair<int, int>> permut;
    // Axis 0
    for (int i = 0; i < nbrPerm / 2; i++) {
        permut.push_back(make_pair(3 + i, 3 + (k - 2) - i));
    }
    // Axis 1
    for (int i = 0; i < nbrPerm / 2; i++) {
        permut.push_back(make_pair(3 + k - 1 + i, 3 + (2 * k - 3) - i));
    }
    // Axis 2
    for (int i = 0; i < nbrPerm / 2; i++) {
        permut.push_back(make_pair(3 + 2 * k - 2 + i, 3 + (3 * k - 4) - i));
    }
    return permut;
}

void FillOther(vector<int> &other, int PK, int Ndof) {
    /*Used to fill Other table in PK
    The function creates a permutation array other
     that symmetrically renumbers the internal points on the edges for a
     mesh element of order PK(see Element_P3 or P4 files for more details)
    */
    for (int i = 0; i < Ndof; i++)
        other[i] = i;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < ((PK - 1) / 2); j++) {
            other[3 + i * (PK - 1) + j] = 3 + (i + 1) * (PK - 1) - j - 1;
            other[3 + (i + 1) * (PK - 1) - j - 1] = 3 + i * (PK - 1) + j;
        }
    }
}