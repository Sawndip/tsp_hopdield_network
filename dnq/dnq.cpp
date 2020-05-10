#include "dnq.h"
#include "nwta.h"
#include <vector>
#include <queue>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <algorithm>
#include <armadillo>
using namespace arma;


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
	vector<vector<double>> closeNeighbor(n, vector<double>(n, pow(10, 5)));
	
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

	//for (int i = 0; i < chains.size(); i++) {
	//	for (int j = 0; j < chains[i].size(); j++) {
	//		cout << chains[i][j] << ' ';
	//	}

	//	cout << endl;
	//}

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
		umap[chains[i][0]] = chains[i][chains[i].size() - 1];
		umap[chains[i][chains[i].size() - 1]] = chains[i][0];
	}

	return umap;
}
vector<vector<double>> getDistanceMatrixForSecondPhase(
	vector<vector<double>>& dist,
	vector<vector<int>>& chains
) {
	vector<int> citiesVec(getNewProblemSize(chains));
	unordered_map<int, int> neighborsMap = getNeghiborsMap(chains);

	int n = dist.size();
	vector<vector<double>> newDist(n, vector<double>(n, 0.));
	
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			if (i == j || neighborsMap.find(i) == neighborsMap.end()|| neighborsMap.find(j) == neighborsMap.end() || neighborsMap[i] == j) {
				newDist[i][j] = 0;
			}
			else {
				newDist[i][j] = dist[i][j];
				newDist[j][i] = dist[i][j];
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
	double A, double B, double C, double D, int m,
	unordered_map<int, int>& neighbours
) {
    int xc;
    int yc;

    if (neighbours.find(x) == neighbours.end()) {
        xc = x;
	} else {
        xc = neighbours[x];
    }
    if (neighbours.find(y) == neighbours.end()) {
        yc = y;
    } else {
        yc = neighbours[y];
    }

	double a = d[yc][x];
	double b = d[xc][y];
	return -(
		A * kd(x, y) * (1 + kd(x, xc)) * (1 - kd(i, j))
		+ B * (1 - kd(x, y)) * kd(i, j)
		+ C
		+ D * (kd(i, mod(j - 1, m)) * d[yc][x] + kd(i, mod(j + 1, m) * d[xc][y]))
	);
}

vector<vector<double>> chnSimulation(vector<vector<double>>& d, double A, double B, double C, double D, double N, int n, int K, vector<vector<int>> &chains) {
	double u_zero = 0.1;
	int m = n - K;
	double eps = 0.000001;
	unordered_map<int, int> neigboursMap = getNeghiborsMap(chains);
	vector<vector<double>> v = generateRandMatr(n, m);
	vector<vector<double>> dv = generateRandMatr(n, m);

	int t = 0;
	int iter = 0;
	bool cont = true;
	int maxIter = 1000;
	double dt = 0.00001;
    vector<int> cities = getNewProblemSize(chains);

	while (iter < maxIter && cont) {
		double firstSum = 0;
		double secondSum = 0;

		double  minK = -INT_MAX;
		double  maxK = INT_MAX;
		for (int x = 0; x < n; x++) {
		    int x_c = cities[x];
			for (int i = 0; i < m; i++) {
				double du = 0;
				
				for (int y = 0; y < n; y++) {
				    int y_c = cities[y];
					for (int j = 0; j < m; j++) {
						du += (Txiyj(d, x_c, i, y_c, j, A, B, C, D, m, neigboursMap)*v[y][j]);
					}
				}

				du += C * N;

				dv[x][i] = 2 / u_zero * v[x][i] * (1 - v[x][i]) * du;
				if (dv[x][i] != 0) {
					double left = -v[x][i] / dv[x][i];
					double right = (1 - v[x][i]) / dv[x][i];

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
            int x_c = cities[x];
            for (int i = 0; i < m; i++) {
                for (int y = 0; y < n; y++) {
                    int y_c = cities[y];
                    for (int j = 0; j < m; j++) {
                        secondSum += dv[x][i] * Txiyj(d, x_c, i, y_c, j, A, B, C, D, m, neigboursMap) * dv[y][j];
                    }
                }
            }
        }

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
	std::cout << "ITER:: " << iter << endl;
	return v;
}

mat invTransormFunc(double u_zero, mat& v) {
	return u_zero / 2 * log(v / (1 - v));
}

double invTransormFunc(double u_zero, double v) {
	return u_zero / 2 * log(v / (1 - v));
}

mat transferFcn(double u_zero, mat v) {
	return 0.5 * (1 + tanh(v / u_zero));	
}

mat weightMatrixTimesVvectorized(colvec& sumVrow, rowvec& sumVcol, double sumV, mat& v, int n, int K, mat& dist, double A, double B, double C, double D) {
	uvec deltaPrev(n - K);
	uvec deltaNext(n - K);

	deltaPrev(0) = n - K - 1;
	deltaNext(n - K - 1) = 0;

	for (int i = 1; i < n - K; i++) {
		deltaPrev(i) = i - 1;
		deltaNext(i-1) = i;
	}
	uvec tmp = linspace<uvec>(0, 2 * K - 1, 2 * K);
	uvec reOrder(reshape(flipud(reshape(tmp, 2, K)), 2 * K, 1));

	mat termA(repmat(sumVrow, 1, size(v, 1)) - v);

	termA.rows(0, 2 * K - 1) = termA.rows(0, 2 * K - 1) + termA.rows(reOrder);
	if (2 * K < termA.n_rows) {
		termA.rows(2 * K, termA.n_rows - 1) = 2 * termA.rows(2 * K, termA.n_rows - 1);
	}

	mat termB = repmat(sumVcol, size(v, 0), 1) - v;

	mat dycx(dist);
	mat dxcy(dist);

	dycx.cols(0, 2 * K - 1) = dycx.cols(reOrder);

	dxcy.rows(0, 2 * K - 1) = dxcy.rows(reOrder);

	mat termD = dxcy * v.cols(deltaNext) + dycx * v.cols(deltaPrev);

	return -A * termA - B * termB - C * sumV - D * termD;
}

std::pair<mat, double> computeTour(double dUaux, mat& distanceMat, mat& v, double E, int n, int K) {
	mat d = dUaux * distanceMat;
	v(find(v > 1 - pow(10, -1 * E))).fill(1.);
	v(find(v < pow(10, -1 * E))).fill(0.);

	if (all(any(v == 1, 0)) && arma::sum(arma::sum(v)) == n - K) {
		uvec provisionalOrder(n-K);
		for (int i = 0; i < n - K; i++) {
			provisionalOrder(i) = v.col(i).index_max();
		}
		mat visitOrder = zeros(1, n);
		int iNew = 0;
		int iOld = 0;

		while (iNew < n) {
			visitOrder(iNew) = provisionalOrder(iOld);
			if (visitOrder(iNew) < 2 * K) {
				if (!mod(visitOrder(iNew), 2)) {
					visitOrder(iNew + 1) = visitOrder(iNew) + 1;
					iNew++;
				}
				else {
					visitOrder(iNew + 1) = visitOrder(iNew) - 1;
					iNew++;
				}
			}
			iNew++;
			iOld++;
		}
		double sum = 0;
		for (int i = 1; i < size(visitOrder); i++) {
			sum += d(visitOrder(i - 1), visitOrder(i));
		}
		sum += d(visitOrder(0), visitOrder(size(visitOrder) - 1));
		return std::make_pair(visitOrder, sum);
	}
	else {
		return std::make_pair(mat(1, 1), -1);
	}
}

std::pair<vec, double> secondPhaseSim(vector<vector<int>>& chains, vector<vector<double>>& dist, double C) {
	C = 1e-5;
	double energy = 0.;
	int R_iter = 20;
	double Q = 0.8;
	double E = 13;
	double u_zero = 0.3;
	double inChainSum = 0;
	vector<double> chainsEndpoint;
	vector<double> singleCities;
	vector<vector<int>> trueChains;
	int n = 0;
	int K = 0;

	for (int i = 0; i < chains.size(); i++) {
		n++;
		if (chains[i].size() > 1) {
			chainsEndpoint.push_back(chains[i][0]);
			chainsEndpoint.push_back(chains[i][chains[i].size() - 1]);
			trueChains.push_back(chains[i]);
			for (int j = 1; j < chains[i].size(); j++) {
				inChainSum += dist[chains[i][j-1]][chains[i][j]];
			}
			inChainSum += dist[chains[i][0]][chains[i][chains[i].size() - 1]];

			K++;
			n++;
		}
		else {
			singleCities.push_back(chains[i][0]);
		}
	}

	chainsEndpoint.insert(chainsEndpoint.end(), singleCities.begin(), singleCities.end());

	vec newCities(chainsEndpoint);
	mat newDistMat(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				newDistMat(i, i) = 0;
			}
			else {
				newDistMat(i, j) = dist[newCities(i)][newCities(j)];
			}
		}
	}


	for (int x = 0; x < 2 * K; x += 2) {
		newDistMat(x + 1, x) = 0;
		newDistMat(x, x + 1) = 0;
	}

	mat distanceToIgnore = ones(n, n);

	for (int i = 1; i < K; i++) {
		distanceToIgnore.submat(2 * (i - 1), 2 * (i - 1), 2 * i - 1, 2 * i - 1) = zeros(2, 2);
	}

	for (int i = K; i < n; i++) {
		distanceToIgnore(i, i) = 0.;
	}


	double dUpper = newDistMat.max();
	newDistMat /= dUpper;
	double dUaux = dUpper;
	dUpper = 1;

	double dL = newDistMat(find(newDistMat > 0)).min();

	double D = 1 / dUpper;
	double rho = dL / dUpper;
	double Np = n - K + 3 / C;
	double A = 3 + C;
	double B = A + rho;
	
	mat One = ones(n, n);
	mat I = eye(n, n);

	mat JK(I);
	for (int i = 1; i <= K; i++) {
		JK.submat(2 * (i - 1), 2 * (i - 1), 2 * i - 1, 2 * i - 1) = ones(2,2) - eye(2, 2);
	}

	mat DK(JK * newDistMat);

	DK.diag().zeros();

	mat IK(JK + I);

	vec vxi_col = solve(B + (n - K) * C * One + (n - K - 1) * A * IK - B * I + D * (DK + DK.t()), ones(n)) * C * Np;

	mat v(repmat(vxi_col, 1, n - K));

	v += (randu(n, n - K) - 0.5) * 10e-10;

	mat u(invTransormFunc(u_zero, v));

	double stopC1 = pow(10, -1 * E);
	double stopC2 = pow(10, -1.5 * E);
	double maxDiffV = 1;
	bool unstable = false;
	double u_e = invTransormFunc(u_zero, stopC1);
	double ib = C * Np;

	rowvec sumVcol = arma::sum(v, 0);
	colvec  sumVrow = arma::sum(v, 1);

	double sumV = arma::sum(sumVcol);
	mat dU(zeros(n, n - K));

	int iter = 0;
	int maxIter = 1000;

	while (iter < maxIter && (maxDiffV > stopC1 || (maxDiffV > stopC1 && unstable))) {
		double dt = 10e100;
		mat TV = weightMatrixTimesVvectorized(sumVrow, sumVcol, sumV, v, n, K, newDistMat, A, B, C, D);
		dU = TV + ib;
		mat dV = 2 / u_zero * v % (1 - v) % dU;
		umat interiorV = u > u_e && u < -u_e;
		umat criteria1 = interiorV && dV < 0;
		if (any(any(criteria1))) {
			mat VdV = -v / dV;
			dt = std::min(dt, VdV(find(criteria1)).min());
		}

		umat criteria2 = interiorV && dV > 0;

		if (any(any(criteria2))) {
			mat antVdV = (1 - v) / dV;
			dt = std::min(dt, antVdV(find(criteria2)).min());
		}

		umat criteria3 = (u <= u_e && dU > 0) || (u >= -u_e && dU < 0);

		unstable = any(any(criteria3));

		if (unstable) {
			dt = min(dt, min(min((abs(u(find(criteria3))) - u_e) / abs(dU(find(criteria3))))));
		}

		double S1 = arma::sum(arma::sum(dV % dU));
		
		sumVcol = arma::sum(dV, 0);
		sumVrow = arma::sum(dV, 1);
		sumV = arma::sum(sumVcol);

		mat TdV = weightMatrixTimesVvectorized(sumVrow, sumVcol, sumV, dV, n, K, newDistMat, A, B, C, D);

		double S2 = -arma::sum(arma::sum(dV % TdV));

		bool sw_optimal = false;
		if (S2 > 0) {
			dt = min(dt, S1 / S2);
			sw_optimal = true;
		}

		if (iter < R_iter && !sw_optimal) {
			dt = Q * dt;
		}

		mat vPrev(v);

		umat borderV = (u <= u_e) || (u >= -u_e);

		uvec uBorder(find(borderV));

		u(uBorder) = u(uBorder) + dU(uBorder) * dt;

		v(uBorder) = transferFcn(u_zero, u(uBorder));
		u(find(borderV && u <= u_e)).fill(u_e);
		v(find(borderV && u <= u_e)).fill(0); // Stays in the border
		u(find(borderV && u >= -u_e)).fill(-u_e);
		v(find(borderV && u >= -u_e)).fill(1); // Stays in the border


		uvec VupdateInInterior = find(interiorV && (vPrev + dV * dt > stopC1) && (vPrev + dV * dt < 1 - stopC1));
		v(VupdateInInterior) = vPrev(VupdateInInterior) + dV(VupdateInInterior) * dt; // Stays in the interior

		uvec VupdateInBorder0 = find(interiorV && (vPrev + dV*dt <= stopC1));
		uvec VupdateInBorder1 = find(interiorV && (vPrev + dV*dt >= 1 - stopC1));

		u(VupdateInBorder0).fill(u_e);
		u(VupdateInBorder1).fill(-u_e);
		v(VupdateInBorder0).fill(0);
		v(VupdateInBorder1).fill(1);

		maxDiffV = max(max(abs(vPrev - v)));

		sumVcol = arma::sum(v, 0);
		sumVrow = arma::sum(v, 1);
		sumV = arma::sum(sumVcol);

		iter = iter + 1;
	}
	std::pair<mat, double> res = computeTour(dUaux, newDistMat, v, E, n, K);
	if (res.second == -1) {
		return res;
	}
	mat chainsOrder = res.first;
	vec finalTour = zeros(dist.size());
	int iNew = 0;

	for (int i = 0; i < size(res.first); i++) {
		finalTour(iNew) = chainsOrder(i);
		iNew++;

		if (chainsOrder(i) < 2 * K) {
			for (int j = 0; j < trueChains.size(); j++) {
				if (trueChains[j][0] == chainsOrder(i)) {
					for (int k = 1; k < trueChains[j].size() - 1; k++) {
						finalTour(iNew) = trueChains[j][k];
						iNew++;

					};
				}
				else if (trueChains[j][trueChains[j].size() - 1] == chainsOrder(i)) {
					for (int k = trueChains[j].size() - 2; k <= 0; k++) {
						finalTour(iNew) = trueChains[j][k];
						iNew++;
					};
				}
			}
		}
	}
	return std::make_pair(finalTour, res.second + inChainSum);
}


vector<vector<double>> solveSecondPhase(vector<vector<double>>& originalDist, vector<vector<int>>& chains) {
    vector<vector<double>> newD = getDistanceMatrixForSecondPhase(originalDist, chains);

	printMatrix<double>(newD);


	int K = 0;
	int n = 0;
	for(int i = 0; i < chains.size(); i++) {
		n++;
		if (chains[i].size() > 1) {
			K++;
			n++;
		}
	}
	double minD = INT_MAX;
    double maxD = 0;
    
    for(int i = 0; i < newD.size(); i++) {
        for (int j = 0; j < newD[0].size(); j++) {
            if (newD[i][j] > maxD) {
                maxD = newD[i][j];
            } else if (newD[i][j] != 0 && newD[i][j] < minD) {
                minD = newD[i][j];
            }
        }
    }
     
    double C = 0.00001;
    double N = double(n - K) + 3/C;
    double A = 3 + C;
    double B = A + minD/maxD;
    double D = 1 / maxD;

    vector<vector<double>> v = chnSimulation(newD, A, B, C, D, N, n, K, chains);
	unordered_map<int, vector<int>> chainsMap;

	for (int i = 0; i < chains.size(); i++) {
		chainsMap[chains[i][0]] = chains[i];
		chainsMap[chains[i][chains[i].size() - 1]] = chains[i];
		reverse(begin(chainsMap[chains[i][chains[i].size() - 1]]), end(chainsMap[chains[i][chains[i].size() - 1]]));
		int a = 0;
	}

	vector<int> asd = getNewProblemSize(chains);
	for (int i = 0; i < v[0].size(); i++) {
		for (int x = 0; x < v.size(); x++) {
			if (v[x][i] == 1) {
				vector<int> chain = chainsMap[asd[x]];
				for (int k = 0; k < chain.size(); k++) {
					std::cout << chain[k] << endl;
				}
			}
		}
	}
	return v;
}