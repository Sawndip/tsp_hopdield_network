#include "dnq.h"
#include "nwta.h"
#include <vector>
#include <queue>
#include <iostream>
#include <numeric>
#include <unordered_map>
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
			closeNeigbor[i][idx] = d[i][idx];
			closeNeigbor[idx][i] = d[i][idx];
		}
		
	}
	
	return closeNeigbor;
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
