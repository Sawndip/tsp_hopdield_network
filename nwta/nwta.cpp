#include <iostream>
#include <vector>
#include <cmath>
#include <2opt/2opt.h>
#include <fstream>
#include <chrono>

using namespace std;

void printMatrix(vector<vector<double>> &m) {
    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m.size(); j++) {
            cout << m[i][j] << ' ';
        }

        cout << endl;
    }
}
template<class T>
void printVec(vector<T>& vec) {
    for (int i=0; i<n; i++) {
        cout << vec[i] << ' ';
    }
    cout << endl;
}

double fScal(double u, double beta) {
    return 1 / (1 + exp(-beta * u));
}

void f(vector<vector<double>> &v, vector<vector<double>> &u, int n, double beta) {
    for (int i = 0; i < n; i++) {
        for (int j  = 0; j < n; j++) {
            v[i][j] = fScal(u[i][j], beta);
        }
    }
}

double sumByRows(vector<vector<double>> &v, int i) {
    double sum = 0;

    for (int k = 0; k < v.size(); k++) {
        sum += v[i][k];
    }

    return sum;
}

double sumByCols(vector<vector<double>> &v, int i) {
    double sum = 0;

    for (int k = 0; k < v.size(); k++) {
        sum += v[k][i];
    }

    return sum;
}

int getMaxIdx(vector<double > &v) {
    double max = 0;
    int maxIdx = 0;
    for (int i = 0; i < v.size(); i++) {
        if (v[i] > max) {
            max = v[i];
            maxIdx = i;
        }
    }

    return maxIdx;
}

vector<vector<double>> solveNwta(
        vector<vector<double>> &u,
        vector<vector<double>> &d,
        int n,
        double beta,
        double eta,
        double lmbd,
        double tau,
        double eps,
        double deltaT
) {
    vector<vector<double>> v(n, vector<double>(n, 0.));
    f(v, u, n, beta);
    int t = 0;
    bool cond = false;
    int tmpCond = 0;
    
    while (!cond && t < 100) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double v1 = sumByRows(v, i);
                double v2 = sumByCols(v, j);

                double  partSum = v1 + v2 - 2;

                if (partSum < eps) {
                    return v;
                }
                
                u[i][j] = u[i][j] + deltaT * ((-eta * partSum) - lmbd * d[i][j] * exp(-t / tau));
                v[i][j] = fScal(u[i][j], beta);
            }
        }


        t += 1;

        if (t % 1000 == 0) {
            cout << t << endl;
        }

        cond = tmpCond == 0;
        tmpCond = 0;
    }

    return v;
}

void fillColWithZeros(vector<vector<double>> &v, int colIdx) {
    for (int i = 0; i < v.size(); i++) {
        v[i][colIdx] = 0;
    }
}

int wta(vector<vector<double>> &v, int n, int iStart) {
    int i = iStart;
    int count = 0;

    while (count < n) {
        int jMax = getMaxIdx(v[i]);
        v[i] = vector<double>(n, 0.);
        fillColWithZeros(v, jMax);
        v[i][jMax] = 1;
        i = jMax;

        if (i == iStart) {
            break;
        }

        count++;
    }

    return count == n - 1;
}

vector<int> getPathByMatrix(vector<vector<double>>& m, int iStart) {
    vector<int> path = vector<int>(m.size(), 0);
    int count = 0;
    int i = iStart;
    while (count < m.size()){
        int jMax = getMaxIdx(m[i]);
        path[count] = jMax;
        i = jMax;
        count++;
    }

    return path;
}

double getLenByPath(vector<vector<double>>& d, vector<int> path) {
    int count = 0;
    double len = 0;
    int i = 0;

    while (count < d.size()) {
        len += d[i][path[count]];
        i = path[count];
        count += 1;
    }

    return len;
}

double solve(
        vector<vector<double>>& dist,
        vector<vector<double>>& distWithNalog,
        int n,
        double beta,
        double eta,
        double lambda,
        double tau,
        double eps,
        double deltaT,
        int iStart
) {
    bool noCycles = false;
    vector<vector<double>> nwtaRes;
    fstream startData(R"(C:\Users\sesip\source\repos\helloworl\nwta\2opt\startData.txt)");

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    while (!noCycles) {
        vector<vector<double>> startVal(n, vector<double>(n, 0));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                startVal[i][j] = (float) rand() / RAND_MAX - 0.5;
            }
        }

        nwtaRes = solveNwta(startVal, distWithNalog, n, beta,  eta, lambda, tau, eps, deltaT);
        noCycles = wta(nwtaRes, n, iStart);
        cout << noCycles << '\n';
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time nwta difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms" << std::endl;

    vector<int> path = getPathByMatrix(nwtaRes, iStart);

    cout << getLenByPath(dist, path);
    cout << endl;

    path = tsp2opt(distWithNalog, path);
    return getLenByPath(dist, path);
}