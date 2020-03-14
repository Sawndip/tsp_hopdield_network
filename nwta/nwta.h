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
