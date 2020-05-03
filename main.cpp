#include "nwta.h"
#include "dnq.h"
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <chrono>

using namespace std;

vector<vector<double>> getDistMatrix(string& name){
    fstream problem(R"(C:\Users\sesip\source\repos\tsp_hopfield_network\tsp\)" + name +".txt");
    int n;
    problem >> n;
    vector<vector<double>> dist(n, vector<double>(n, 0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            problem >> dist[i][j];
        }
    }

    problem.close();
    return dist;
}

vector<vector<double>> getDistMatrixWithNalog(vector<vector<double>>& dist, int iStart, double p) {
    vector<vector<double>> newDist(dist.size(), vector<double>(dist.size(), 0));

    for (int i = 0; i < dist.size(); i++) {
        for (int j = 0; j < dist.size(); j++) {
            newDist[i][j] = dist[i][j];

            if (i == j) {
                newDist[i][j] = pow(10, 5);
            }

            if (j == iStart) {
                newDist[i][j] *= p;
            }
        }
    }

    return newDist;
}

int main(int argc, char** argv)
{
    double beta = 0.1;
    double eta = 10;
    double lambda = 1;
    double tau = pow(10, 3);
    double eps = pow(10, -5);
    double p = pow(10,6);
    int iStart = 0;
    double deltaT = 0.025;

    //string name = "a280";
    //vector<vector<double>> dist = getDistMatrix(name);
    //vector<vector<double>> distWithNalog = getDistMatrixWithNalog(dist, iStart, p);

    //int n = dist.size();
    int n = 10;

    time_t t;

    srand((unsigned) time(&t));
    vector<vector<double>> tmp = {
        {10, 0.43, 0.82, 0.85, 0.74, 0.1, 0.22, 0.25, 0.37, 0.35},
        {0.43, 10, 0.8, 0.75, 0.4, 0.34, 0.2, 0.1, 0.8, 0.6},
        {0.82, 0.8, 10, 0.2, 0.6, 0.7, 0.75, 0.7, 0.1, 0.70},
        {0.85, 0.75, 0.2, 10, 0.5, 0.8, 0.85, 0.8, 0.1, 0.8},
        {0.74,0.4, 0.6, 0.5, 10, 0.75, 0.7, 0.64, 0.68, 9.99},
        {0.1, 0.34, 0.7, 0.8, 0.75, 10, 0.1, 0.25, 0.8, 0.4},
        {0.22, 0.2, 0.75, 0.85, 0.7, 0.1, 10, 0.14, 9.99, 0.5},
        {0.25, 0.1, 0.7, 0.8, 0.64, 0.25, 0.14, 10, 0.9, 0.6},
        {0.37, 0.8, 0.1, 0.1, 0.68, 0.8, 9.99, 0.9, 10, 0.8},
        {0.35, 0.6, 0.7, 0.8, 9.99, 0.4, 0.5, 0.6, 0.8, 10}
    };

    vector<vector<double>> dist(closestNeigboursMatrix(tmp, 3));

    for (int i = 0; i < 10; i++) {
        vector<vector<double>> u(generateRandMatr(n));
        vector<vector<double>> res(solveNwta(u, dist, n, beta, eta, lambda, tau, eps, deltaT));
        vector<vector<int>> chains(wtaWithCycles(res, n));

        if (chains.size() > 1) {
            res = solveSecondPhase(tmp, chains);
            return 0;
        }
        else {
            cout << "Chains: " << endl;
            for (int i = 0; i < chains[0].size(); i++) {
                cout << chains[0][i] << endl;
            }

            cout << endl;
            cout << endl;
            cout << endl;
            cout << endl;
            cout << endl;

        }
       /* cout << endl;
        cout << "Res: " << endl;
        printMatrix(res);
        cout << endl;*/
    }

    //vector<vector<int>> testChains = {
    //    {9},
    //    {2, 8, 3},
    //    {0, 5, 6, 7, 1},
    //    {4}
    //};

    //vector<vector<double>> res = solveSecondPhase(tmp, testChains);
    //double K = 0;
    //for (int i = 0; i < testChains.size(); i++) {
    //    if (testChains[i].size() > 1) {
    //        K++;
    //    }
    //}
    //if (testChains.size() == 1) {
    //    cout << "SOLVED" << endl;
    //    for (int i = 0; i < testChains[0].size(); i++) {
    //        cout << testChains[0][i] << ' ';
    //    }
    //    cout << endl;
    //} else {
    //    vector<vector<double>> newD = getDistanceMatrixForSecondPhase(tmp, testChains);
    //    printMatrix(newD);
    //    double minD = INT_MAX;
    //    double maxD = 0;
    //    
    //    for(int i = 0; i < newD.size(); i++) {
    //        for (int j = 0; j < newD[0].size(); j++) {
    //            if (newD[i][j] > maxD) {
    //                maxD = newD[i][j];
    //            } else if (newD[i][j] != 0 && newD[i][j] < minD) {
    //                minD = newD[i][j];
    //            }
    //        }
    //    }
    //    
    //    double C = 0.001;
    //    double N = newD.size() + K + 3/C;
    //    double A = 3 + C;
    //    double B = A + minD/maxD;
    //    double D = 1 / maxD;
    //    vector<vector<double>> v = chnSimulation(newD, A, B, C, D, N, K, testChains);
    //    printMatrix(v);
    //}

    return 0;
    //for (unsigned int i = 0; i < 10; i++) {
    //    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    //    double pathLen = solve(dist, distWithNalog, n, beta, eta, lambda, tau, eps, deltaT, iStart);

    //    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " s" << std::endl;
    //    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms" << std::endl;

    //    cout << pathLen << endl;
    //}

    //return 0;
}
