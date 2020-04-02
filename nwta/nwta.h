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

void printMatrix(vector<vector<double>> &m);

template<typename  T>
void printVec(vector<T>& vec);