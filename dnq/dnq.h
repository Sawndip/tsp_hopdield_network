#include <vector>

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
	vector<vector<int>> &chains
);