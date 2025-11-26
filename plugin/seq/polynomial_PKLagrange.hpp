#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

using namespace std;
typedef double REAL;

long factorial(long n) {
    if (n <= 1)
        return 1;
    long result = 1;
    for (long i = 2; i <= n; ++i) {
        result *= i;
    }
    return result;
}

vector<vector<long>> generate_partitions(int sum, int nparts) {
    /*Generates all compositions (or ordered partitions) of an integer sum into nparts parts.*/
    vector<vector<long>> partitions;

    if (nparts == 1) {
        partitions.push_back({int(sum)});
        return partitions;
    }

    for (int first = 0; first <= sum; first++) {
        vector<vector<long>> rest_partitions = generate_partitions(sum - first, nparts - 1);
        for (auto &rest : rest_partitions) {
            vector<long> comp;
            comp.push_back(first);
            comp.insert(comp.end(), rest.begin(), rest.end());
            partitions.push_back(comp);
        }
    }
    return partitions;
}

vector<vector<long>> sort_partitions_geometric(vector<vector<long>> partitions, int K) {
    /*Function to sort the partitions according to the geometric order P_K*/
    vector<vector<long>> sorted;
    // 1. Vertices (i+j+k=K )
    sorted.push_back({K, 0, 0});
    sorted.push_back({0, K, 0});
    sorted.push_back({0, 0, K});
    // 2. Edges (i+j+k=K one coordiante = 0, the others > 0)
    // Order: Edge k=0 , Edge j=0, then i=0

    vector<vector<long>> edge1, edge2, edge3;
    for (auto &p : partitions) {
        if (p[2] == 0 && p[0] > 0 && p[1] > 0)
            edge1.push_back(p);
        if (p[1] == 0 && p[0] > 0 && p[2] > 0)
            edge2.push_back(p);
        if (p[0] == 0 && p[1] > 0 && p[2] > 0)
            edge3.push_back(p);
    }

    // Sort Edge 0
    sort(edge1.begin(), edge1.end(), [](const vector<long> &a, const vector<long> &b) { return a[0] > b[0]; });

    // Sort Edge 1
    sort(edge2.begin(), edge2.end(), [](const vector<long> &a, const vector<long> &b) { return a[0] < b[0]; });

    // Sort Edge 2
    sort(edge3.begin(), edge3.end(), [](const vector<long> &a, const vector<long> &b) { return a[1] > b[1]; });

    sorted.insert(sorted.end(), edge3.begin(), edge3.end());
    sorted.insert(sorted.end(), edge2.begin(), edge2.end());
    sorted.insert(sorted.end(), edge1.begin(), edge1.end());

    // 3. Interior point (i,j,k  > 0)
    vector<vector<long>> interior;
    for (auto &p : partitions) {
        if (p[0] > 0 && p[1] > 0 && p[2] > 0) {
            interior.push_back(p);
        }
    }
    // Sort internal points to match PK order points
    sort(interior.begin(), interior.end(), [](const vector<long> &a, const vector<long> &b) {
        if (a[2] != b[2])
            return a[2] > b[2];
        if (a[0] != b[0])
            return a[0] < b[0];
        return a[1] > b[1];
    });

    sorted.insert(sorted.end(), interior.begin(), interior.end());

    return sorted;
}

void BasisFctPK(int K, vector<vector<long>> &coef, vector<vector<long>> &shift, vector<long> &ff, vector<long> &il,
                vector<long> &jl, vector<long> &kl) {
    // Silverster formula ~1968
    // Compute Phi_{i,j,k}(lambda1,lambda2,lambda3) basis functions; while i+j+k=K and lambda are the barycentric
    // coordinates

    vector<vector<long>> unsortedpartitionK = generate_partitions(K, 3);
    vector<vector<long>> partitionK = sort_partitions_geometric(unsortedpartitionK, K);
    int idx = 0;
    for (auto &partition : partitionK) {
        long i = partition[0];
        long j = partition[1];
        long k = partition[2];
        if (i + j + k == K) {
            int ID = 0;
            long denom = factorial(i) * factorial(j) * factorial(k);
            ff[idx] = denom;
            il[idx] = i;
            jl[idx] = j;
            kl[idx] = k;
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
            idx++;
        }
    }
}
