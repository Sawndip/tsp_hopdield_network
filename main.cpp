#include "nwta.h"
#include "dnq.h"
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include <armadillo>
using namespace arma;

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

    string name = "gr17";
    vector<vector<double>> dist = getDistMatrix(name);
    printMatrix(dist);
    vector<vector<double>> distWithNalog = getDistMatrixWithNalog(dist, iStart, p);
    printMatrix(distWithNalog);
    vector<vector<double>> closeNeighborsDist = closestNeigboursMatrix(dist, 3);
    printMatrix(closeNeighborsDist);

    int n = dist.size();

    time_t t;
    srand((unsigned)time(&t));

    for (int i = 0; i < 10; i++) {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        vector<vector<double>> u(generateRandMatr(n));
        vector<vector<double>> res(solveNwta(u, closeNeighborsDist, n, beta, eta, lambda, tau, eps, deltaT));
        vector<vector<int>> chains(wtaWithCycles(res, n));
        cout << "CHAINS: " << endl;
        printMatrix<int>(chains);
        int sum = 0;
        for (int i = 0; i < chains.size(); i++) {
            if (chains[i].size() == 1) {
                sum++;
            }
        };
        if (sum == n) {
            cout << "BROKEN" << endl;
        } else if (chains.size() == 1) {
            cout << "SOLVED" << endl;
            for (int i = 0; i < chains[0].size(); i++) {
                cout << chains[0][i] << ' ';
            }
            cout << endl;
        } else {
            mat v = fromMatlab(chains, dist, 1e-05);
            cout << v << endl;
        }

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " s" << std::endl;
        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms" << std::endl;
    }
    return 0;
}
