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
        cout << "Chains: " << endl;
        vector<vector<int>> wtaWithCycles(res, n);
        cout << endl;
        cout << "Res: " << endl;
        printMatrix(res);
        cout << endl;
    }


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
