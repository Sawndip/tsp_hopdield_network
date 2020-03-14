// awtf.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "nwta.h"
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <random>
using namespace std;
vector<vector<double>> getDistMatrix(string& name, double p, int iStart){
    fstream problem(R"(C:\Users\sesip\source\repos\helloworl\tsp\)" + name +".txt");
    int n;
    problem >> n;
    vector<vector<double>> dist(n, vector<double>(n, 0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double d;
            problem >> d;
            if (i == iStart) {
                dist[i][j] = p;
            } else {
                dist[i][j] = d;
            }
        }
    }

    problem.close();
    return dist;
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
    double deltaT = pow(10, -5);

    string name = "gr17";
    vector<vector<double>> dist = getDistMatrix(name, p, iStart);
    int n = dist.size();
    minstd_rand simple_rand;
    simple_rand.seed(42);

    vector<vector<double>> startVal(n, vector<double>(n, 0));
    ofstream f("C:\\Users\\sesip\\source\\repos\\helloworl\\a.txt");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            startVal[i][j] = double(simple_rand() % 100) / 100 - (1. / 2.);
            f << startVal[i][j] << ' ';
        }
    }


    vector<vector<double>> nwtaRes = solveNwta(startVal, dist, n, beta,  eta, lambda, tau, eps, deltaT);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << nwtaRes[i][j] << ' ';
        }

        cout << endl;
    }
    cout << endl;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << endl;

    wta(nwtaRes, n, iStart);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << nwtaRes[i][j] << ' ';
        }

        cout << endl;
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
