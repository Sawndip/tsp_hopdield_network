#include <vector>
using namespace std;
vector<vector<double>> solveNwta(
        vector<vector<double>> &u ,
        vector<vector<double>> &d,
        int n,
        double beta,
        double eta,
        double lmbd,
        double tau,
        double eps,
        double deltaT
);

void wta(vector<vector<double>> &v, int n, int iStart);
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
);

template<class T>
void printMatrix(vector<vector<T>>& m) {
    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m[i].size(); j++) {
            cout << m[i][j] << ' ';
        }

        cout << endl;
    }
}
vector<vector<double>> generateRandMatr(int n);
int getMaxIdx(vector<double >& v);
void fillColWithZeros(vector<vector<double>>& v, int colIdx);
double getLenByPath(vector<vector<double>>& d, vector<int> path);