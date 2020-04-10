#include <vector>
#include <iostream>
using namespace std;

double routeCost(vector<vector<double>>& d, vector<int> path) {
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

template<typename T>
void printVec(vector<T>& vec) {
    for (int i = 0; i < vec.size(); i++) {
        cout << vec[i] << ' ';
    }

    cout << endl;
}

vector<int> swap2opt(vector<int>& route, int i, int k) {
    vector<int> newRoute(route);
    int dec = 0;

    for (int c = i; c < k + 1; c++) {
        newRoute[c] = route[k - dec];
        dec++;
    }

    return newRoute;
}


vector<int> tsp2opt(vector<vector<double>>& graph, vector<int>& route) {
    int size = route.size();

    int improve = 0;
    vector<int> bestFoundRoute = vector<int>(route);
    bool improved = true;

    while (improved) {
        improved = false;
        double best_distance = routeCost(graph, bestFoundRoute);

        for (int i = 0; i < size - 1; i++) {
            for (int k = i + 1; k < size; k++) {
                vector<int> newRoute = swap2opt(route, i, k);

                double new_distance = routeCost(graph, newRoute);

                if (new_distance < best_distance) {
                    // Improvement found so reset
                    bestFoundRoute = newRoute;
                    best_distance = new_distance;
                    break;
                }
            }
        }
    }
    return bestFoundRoute;

}
