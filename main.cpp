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
    time_t t;
    srand((unsigned)time(&t));

    vector<string> names = { "gr17", "gr24", "fri26", "eil51", "ch130", "d198", "a280", "pcb442", "pr1002" };
    //vector<string> names = { "gr17", "fri26", "eil51" };
    for (int idx = 0; idx < names.size(); idx++) {
        int iters = 10;
        double meanLen = 0;
        double minLen = pow(10, 10);
        double maxLen = 0;

        double minTime = pow(10, 10);
        double maxTIme = 0;
        double meanTime = 0;

        string name = names[idx];
        vector<vector<double>> dist = getDistMatrix(name);
        vector<vector<double>> distWithNalog = getDistMatrixWithNalog(dist, iStart, p);
        vector<vector<double>> closeNeighborsDist = closestNeigboursMatrix(dist, 3);
        int n = dist.size();

        for (int iter = 0; iter < iters; iter++) {
            std::clock_t c_start = std::clock();

            vector<vector<double>> u(generateRandMatr(n));
            vector<vector<double>> res(solveNwta(u, closeNeighborsDist, n, beta, eta, lambda, tau, eps, deltaT));
            vector<vector<int>> chains(wtaWithCycles(res, n));
            double len = 0;
            if (chains.size() == 1) {
                len = getLenByPath(dist, chains[0]);
            } else {
                std::pair<vec, double> result = secondPhaseSim(chains, dist, 1e-05);
                len = result.second;
            }
            std::clock_t c_end = std::clock();
            double t = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;

            meanLen += len;
            if (minLen > len) {
                minLen = len;
            }
            meanTime += t;
            if (minTime > t) {
                minTime = t;
            }
        }

        std::cout << "Name: " << name << endl << "Min len: " << minLen << "; Mean len: " << meanLen/iters << "; Min time: " << minTime << "; Mean time: " << meanTime / iters << endl;
    }

    return 0;
}
