#include "dnq.h"
#include "nwta.h"
#include <vector>
#include <queue>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <algorithm>

using namespace std;

vector<int> findBestIndices(vector<double>& vec, const int& N)
{
	vector<int> indices(vec.size());
	std::iota(indices.begin(), indices.end(), 0); // fill with 0,1,2,...

	std::partial_sort(indices.begin(), indices.begin() + N, indices.end(),
		[&vec](int i, int j) {return vec[i] < vec[j]; });

	return vector<int>(indices.begin(), indices.begin() + N);
}


vector<vector<double>> closestNeigboursMatrix(vector<vector<double>>& d, int tau)
{
	int n = d.size();
	vector<int> temp(tau);
	iota(begin(temp), end(temp), 0);
	vector<vector<int>> closeNeighborIdx(n, vector<int>(tau, 0));
	vector<vector<double>> closeNeighbor(n, vector<double>(n, 100.));
	
	for (int i = 0; i < n; i++) {
		closeNeighborIdx[i] = findBestIndices(d[i], tau);
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < tau; j++) {

			int idx = closeNeighborIdx[i][j];
			closeNeighbor[i][idx] = d[i][idx];
			closeNeighbor[idx][i] = d[i][idx];
		}
		
	}
	
	return closeNeighbor;
}

vector<vector<int>> wtaWithCycles(vector<vector<double>>& v, int n) {
	vector<vector<int>> chains;
	unordered_map<int, bool> umap;

	for (int i = 0; i < n; i++) {
		int iStart = i;
		int iTemp = i;

		if (umap.find(i) == umap.end()) {
			umap[i] = true;
			chains.push_back(vector<int>(1, iStart));

			do {
				int jMax = getMaxIdx(v[iTemp]);
				v[iTemp] = vector<double>(n, 0.);
				fillColWithZeros(v, jMax);
				v[iTemp][jMax] = 1;
				iTemp = jMax;
				umap[iTemp] = true;

				if (iTemp != iStart) {
					chains[chains.size() - 1].push_back(iTemp);
				}
			} while (iTemp != iStart);
		}
		else {
			continue;
		}
	}

	for (int i = 0; i < chains.size(); i++) {
		for (int j = 0; j < chains[i].size(); j++) {
			cout << chains[i][j] << ' ';
		}

		cout << endl;
	}

	return chains;
}

vector<int> getNewProblemSize(vector<vector<int>>& chains) {
	vector<int> citiesVector;
	for (int i = 0; i < chains.size(); i++) {
		citiesVector.push_back(chains[i][0]);
		if (chains[i].size() > 1) {
			citiesVector.push_back(chains[i][chains[i].size() - 1]);
		}
	}

	return citiesVector;
}

unordered_map<int, int> getNeghiborsMap(vector<vector<int>>& chains) {
	unordered_map<int, int> umap;

	for (int i = 0; i < chains.size(); i++) {
		if (chains[i].size() > 1) {
			umap[chains[i][0]] = chains[i][chains[i].size() - 1];
			umap[chains[i][chains[i].size() - 1]] = chains[i][0];
		}
	}

	return umap;
}
vector<vector<double>> getDistanceMatrixForSecondPhase(
	vector<vector<double>>& dist,
	vector<vector<int>>& chains
) {
	vector<int> citiesVec(getNewProblemSize(chains));
	unordered_map<int, int> neighborsMap = getNeghiborsMap(chains);

	int n = citiesVec.size();
	vector<vector<double>> newDist(n, vector<double>(n, 0.));
	
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			if (i == j || neighborsMap[citiesVec[i]] == citiesVec[j]) {
				newDist[i][j] = 0;
			}
			else {
				newDist[i][j] = dist[citiesVec[i]][citiesVec[j]];
				newDist[j][i] = dist[citiesVec[i]][citiesVec[j]];
			}
		}
	}

	return newDist;
}

int kd(int i, int j) {
	if (i == j)
		return 1;

	return 0;
}
int mod(int a, int b) { return (a % b + b) % b; };

vector<vector<double>> generateRandMatr(int n, int m) {
	vector<vector<double>> startVal(n, vector<double>(m, 0.));

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {

			startVal[i][j] = 0.5 + pow(10, -7) * ((double)rand() / RAND_MAX - 0.5);
		}
	}

	return startVal;
}
double Txiyj(
	vector<vector<double>>& d,
	int x, int i, int y, int j,
	double A, double B, double C, double D,
	unordered_map<int, int>& neighbours
) {
	int m = 4;
    int xc;
    int yc;

    if (neighbours.find(x) == neighbours.end()) {
        xc = x;
	} else {
        xc = neighbours[x];
    }
    if (neighbours.find(y) == neighbours.end()) {
        yc = x;
    } else {
        yc = neighbours[y];
    }

	return -(
		A * kd(x, y) * (1 + kd(x, xc)) * (1 - kd(i, j))
		+ B * (1 - kd(x, y)) * kd(i, j)
		+ C
		+ D * (kd(i, mod(j - 1, m)) * d[yc][x] + kd(i, mod(j + 1, m) * d[xc][y]))
	);
}

vector<vector<double>> chnSimulation(vector<vector<double>>& d, double A, double B, double C, double D, double N, vector<vector<int>> &chains) {
	double u_zero = 1;
	int n = d.size();
	int m = 4;
	double eps = 0.000001;
	unordered_map<int, int> neigboursMap = getNeghiborsMap(chains);
	vector<vector<double>> v = generateRandMatr(n, m);
	vector<vector<double>> dv = generateRandMatr(n, m);

	int t = 0;
	int iter = 0;
	bool cont = true;
	int maxIter = 1000;
	double dt = 0.001;

	while (iter < maxIter && cont) {
		double firstSum = 0;
		double secondSum = 0;

		double  minK = -INT_MAX;
		double  maxK = INT_MAX;
		for (int x = 0; x < n; x++) {
			for (int i = 0; i < m; i++) {
				double du = 0;
				
				for (int y = 0; y < n; y++) {
					for (int j = 0; j < m; j++) {
						du += (Txiyj(d, x, i, y, j, A, B, C, D, neigboursMap)*v[y][j]);
					}
				}

				du += C * N;

				dv[x][i] = 2 / u_zero * v[x][i] * (1 - v[x][i]) * du;
				if (dv[x][i] != 0) {
					double left = -v[x][i] / dv[x][i];
					double right = (1 - v[x][i]) / dv[x][i];

					double a = v[x][i] + left * dv[x][i];
					double b = v[x][i] + right * dv[x][i];

					if (left > right) {
						if (left < maxK) {
							maxK = left;
						}
						if (right < 0) {
							right = 0;
						}
						if (right > minK) {
							minK = right;
						}
					}
					else {
						if (right < maxK) {
							maxK = right;
						}
						if (left < 0) {
							left = 0;
						}
						if (left > minK) {
							minK = left;
						}
					}
				}

				firstSum += du * dv[x][i];
			}
		}

        for (int x = 0; x < n; x++) {
            for (int i = 0; i < m; i++) {
                for (int y = 0; y < n; y++) {
                    for (int j = 0; j < m; j++) {
                        secondSum += dv[x][i] * Txiyj(d, x, i, y, j, A, B, C, D, neigboursMap) * dv[y][j];
                    }
                }
            }
        }
		double tpm = firstSum / secondSum;
        if (secondSum < 0) {
            dt = min(maxK, -firstSum / secondSum) ;
        } else {
            dt = maxK;
        }

        double maxDif = 0.;

        for (int x = 0; x < n; x++) {
            for (int i = 0; i < m; i++) {
                double newV = v[x][i] + dv[x][i] * dt;
                double dif = abs(v[x][i] - newV);
				v[x][i] = newV;
	
                if (dif > maxDif) {
                    maxDif = dif;
                }
            }
        }

        cont = maxDif > eps;
		iter++;
	}

	return v;
}