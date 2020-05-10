#include <vector>
#include <armadillo>
using namespace arma; 
using namespace std;

vector<vector<double>> closestNeigboursMatrix(vector<vector<double>>& d, int tau);
vector<vector<int>> wtaWithCycles(vector<vector<double>>& v, int n);
vector<vector<double>> getDistanceMatrixForSecondPhase(
	vector<vector<double>>& dist,
	vector<vector<int>>& chains
);

vector<vector<double>> chnSimulation(
	vector<vector<double>>& d,
	double A,
	double B,
	double C,
	double D,
	double N,
	int n,
	int K,
	vector<vector<int>> &chains
);

vector<vector<double>> solveSecondPhase(vector<vector<double>>& originalDist, vector<vector<int>>& chains);
std::pair<vec, double> secondPhaseSim(vector<vector<int>>& chains, vector<vector<double>>& dist, double C);
