#include <iostream>
#include "ff++.hpp"
#include "AddNewFE.h"
#include "utils_PKLagrange.hpp"


long factorial(long n) {
    if (n <= 1)
        return 1;
    long result = 1;
    for (long i = 2; i <= n; ++i) {
        result *= i;
    }
    return result;
}


vector<vector<int>> generate_barycoordinates(int Pk){
    int ndof= (Pk+1)*(Pk+2)*(Pk+3)/6;
    vector<vector<int>> pts(ndof,vector<int>(4,0));
    // Vertices
    pts[0]={Pk,0,0,0};
    pts[1]={0,Pk,0,0};
    pts[2]={0,0,Pk,0};
    pts[3]={0,0,0,Pk};
    
    // Edges 
    int edgePairs[6][2] = {
        {0, 1},
        {0, 2},
        {0, 3},
        {1, 2},
        {1, 3},
        {2, 3}
    };
    int c = 4; // 
    for (int e = 0; e < 6; ++e) {
        int i = edgePairs[e][0];
        int j = edgePairs[e][1];
        for (int k = 1; k < Pk; ++k) {
            // barycentriques: j = Pk-k, i = k, for the others = 0
            pts[c][j] = Pk - k;
            pts[c][i] = k;
            ++c;
        }
    }
    assert((c-4)==6*(Pk-1));
    // faces
    int faces[4][3]=  {{3,2,1}, {0,2,3},{ 3,1,0},{ 0,1,2}};

     c = 4 + 6*(Pk-1); // shift edges*(Pk-1)
    for(int f=0; f<4; f++){

        int a = faces[f][0];
        int b = faces[f][1];
        int d = faces[f][2];

        for(int i=1; i<=Pk-2; i++){
            for(int j=1; j<=Pk-1-i; j++){

                int k = Pk - i - j;

                if(k >= 1){

                    pts[c][a] = k;
                    pts[c][b] = j;
                    pts[c][d] = i;

                    c++;
                }
            }
        }
        }
    assert((c-(4 + 6*(Pk-1)))==4*(Pk-1)*(Pk-2)/2);
    // Volume 
    for(int i=1; i<=Pk-3; i++){
        for(int j=1; j<=Pk-2-i; j++){
            for(int k=1; k<=Pk-1-i-j; k++){

                int l = Pk - i - j - k;

                if(l >= 1){
                    pts[c][0] = i;
                    pts[c][1] = j;
                    pts[c][2] = k;
                    pts[c][3] = l;
                    c++;
                }
            }
        }
    }
    assert(c==ndof);
    /*// affihage
    for (int i=0;i<pts.size();i++) cout<<"{"<<pts[i][0]<<" , "<<pts[i][1]<<" , "<<pts[i][2]<<" , "<<pts[i][3]<<"} ";
    cout<<endl;*/
    return pts;
}



void BasisFctPK3D(int K, vector<vector<int>> &coef, vector<vector<int>> &shift, vector<int> &ff,vector<vector<int>> & partitionK) {
    // Silverster formula ~1968 for 3D lagrange elements
    // Compute Phi_{i,j,k,l}(lambda1,lambda2,lambda3,lambda4) basis functions; while i+j+k+l=K and lambda are the barycentric
    // coordinates
    //vector<vector<int>> partitionK =generate_barycoordinates(K);
    partitionK =generate_barycoordinates(K);
    int idx = 0;
    for (auto &partition : partitionK) {
        int i = partition[0];
        int j = partition[1];
        int k = partition[2];
        int l = partition[3];

        if (i + j + k +l  == K) {
            int ID = 0;
            int denom = factorial(i) * factorial(j) * factorial(k)* factorial(l);
            ff[idx] = denom;
            // cout<<denom<<" Point -->"<< idx<<"( "<<i<<" "<<j<<" "<<k<<" )," << "denom= "<< tgamma(i+1)<<"
            // "<<tgamma(j+1)<<" "<<tgamma(k+1)<<" ff= " <<ff[idx]<< endl;
            if (i > 0) {
                for (int ii = 0; ii <= i - 1; ii++) {

                    coef[idx][ID] = 0;
                    shift[idx][ID] = ii;
                    ID++;
                }
            }
            if (j > 0) {
                for (int jj = 0; jj <= j - 1; jj++) {
                    coef[idx][ID] = 1;
                    shift[idx][ID] = jj;
                    ID++;
                }
            }
            if (k > 0) {
                for (int kk = 0; kk <= k - 1; kk++) {
                    coef[idx][ID] = 2;
                    shift[idx][ID] = kk;
                    ID++;
                }
            }
            if (l > 0) {
                for (int lc = 0; lc <= l - 1; lc++) {
                    coef[idx][ID] = 3;
                    shift[idx][ID] = lc;
                    ID++;
                }
            }
            idx++;
        }
    }
    // affihage
    /*cout<<"Polynome:\n";
    for (int i=0;i<coef.size();i++) cout<<"{"<<coef[i][0]<<" , "<<coef[i][1]<<" , "<<coef[i][2]<<" , "<<coef[i][3]<<"} ";
    cout<<endl;*/
}

vector<pair<int, int>> ExchangeidxVector3D(int k) {
    int nbrPerm = (k % 2) ? (k - 1) : (k - 2);
    vector<pair<int, int>> permut;

    int half = nbrPerm / 2;

    for (int axis = 0; axis < 6; axis++) {
        int base = 4 + axis * (k - 1);
        int mirror = base + (k - 2);

        for (int i = 0; i < half; i++) {
            permut.emplace_back(base + i, mirror - i);
        }
    }

    return permut;
}

/*std::vector<std::vector<std::vector<int>>> permutationslist (int kp){
    std::vector<std::vector<std::vector<int>>>  permutlist;
    for (int idx=0;idx<kp-4;idx++){
        int nbrperm=int((idx+2)/2.);
        for(int p=0;p<nbrperm; p++){
            //0102
            int FirstPts0102=(idx+1)+p*(kp-3)-p*(p-1)/2;
            int LastPts0102=1+(idx+1)*(kp-3) -(idx*(idx-1)/2) -(p*(kp-3-idx)+p*(p-1)/2);
            //0112
            int FirstPts0112=(kp-3)-(idx+1)+p*(kp-2)-p*(p-1)/2;
            int LastPts0112=(kp-3)+(idx+1)*(kp-3) -(idx*(idx+1)/2) -(p*(kp-2-idx)+p*(p-1)/2);
            //0212
            int FirstPts0212= ((kp-2)*(kp-1)/2-1)-((idx+2)*(idx+1)/2)-(idx+1)+(p);
            int LastPts0212=((kp-2)*(kp-1)/2-1)-((idx+2)*(idx+1)/2)-p;
            //if (p!=nbrperm-1) 
                permutlist.push_back({{-(1),FirstPts0102,LastPts0102},{FirstPts0112,-(1),LastPts0112},{FirstPts0212,LastPts0212,-(1)}});
              

        }
    }

    return permutlist; 
}*/

inline int numerotationPermutation(int kp, std::vector<int> barcoord,std::vector<int> permu){
    int i = barcoord[permu[0]];
    int j = barcoord[permu[1]];
    int k = barcoord[permu[2]];

    //std::cout<<"\t old nbr: "<<double((barcoord[0])*((kp-3)+1))-(barcoord[0]*(barcoord[0]-1))/2.+barcoord[1]<<" new one: "<< double((i)*((kp-3)+1))-(i*(i-1))/2.+j<<std::endl;
    //return double((j)*((kp-3)+1))-(j*(j-1))/2.+k;
    return double((k)*((kp-3)+1))-(k*(k-1))/2.+j;
    
    //return double(j*((kp-3)+1))-(j*(j-1))/2.+k;

}

inline int numerotationLocIntern(int kp, int i , int j, int k){
    assert(i+j+k==(kp-3));
    //return double(j*((kp-3)+1))-(j*(j-1))/2.+k;
    return double((k)*((kp-3)+1))-(k*(k-1))/2.+j;

}