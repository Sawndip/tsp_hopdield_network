// awtf.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "nwta.h"
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <chrono>


using namespace std;
vector<vector<double>> getDistMatrix(string& name){
    fstream problem(R"(C:\Users\sesip\source\repos\helloworl\tsp\)" + name +".txt");
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

    string name = "a280";
    vector<vector<double>> dist = getDistMatrix(name);
    vector<vector<double>> distWithNalog = getDistMatrixWithNalog(dist, iStart, p);

    int n = dist.size();
    time_t t;

    srand((unsigned) time(&t));

    for (unsigned int i = 0; i < 10; i++) {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        double pathLen = solve(dist, distWithNalog, n, beta, eta, lambda, tau, eps, deltaT, iStart);

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " s" << std::endl;
        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms" << std::endl;

        cout << pathLen << endl;
    }

    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
