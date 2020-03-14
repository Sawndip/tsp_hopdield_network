#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

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

    return  maxIdx;
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
    vector<vector<double>> v(n, vector<double>(n, 0));
    f(v, u, n, beta);
    int t = 0;
    bool cond = false;
    int tmpCond = 0;
    
    while (!cond && t < 100000) {
        double v1 = sumByRows(v, i);
        double v2 = sumByCols(v, j);

        double  partSum = v1 + v2 - 2;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {


                if (partSum > eps) {
                    tmpCond += 1;
                }
                
                u[i][j] = u[i][j] + deltaT * ((-eta * partSum) - lmbd * d[i][j] * exp(-t / tau));
                v[i][j] = fScal(u[i][j], beta);
            }
        }
        if (t % 1000 == 0) {
            cout << t << endl;
        }
        t += 1;
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

void wta(vector<vector<double>> &v, int n, int iStart) {
    int i = INT_MAX;

    while (i != iStart) {
        int count = 0;
        i = iStart;
        while (count < n) {
            int jMax = getMaxIdx(v[i]);
            v[i] = vector<double>(n, 0);
            fillColWithZeros(v, jMax);
            v[i][jMax] = 1;
            i = jMax;

            if (i == iStart) {
                break;
            }

            count++;
        }
    }
}


/*
 *
def wta(v, n):
    count = 0
    i = 10000000
    while i != iStart:
        count = 0
        i = iStart
        while count < n:
            jMax = np.argmax(v[i])
            v[i] = np.zeros(n)
            v[:, jMax] = np.zeros(n)
            v[i, jMax] = 1
            i = jMax
            if i == iStart:
                break
            count += 1

    return count == n - 1*/