#include <vector>
#include <iostream>
using namespace std;

//double routeCost(vector<vector<double>>& graph, vector<int>& path) {
//    double cost = 0;
//
//    for (int i = 0; i < path.size() - 1; i++) {
//        cost += graph[path[i]][path[i+1]];
//    }
//
//    cost += graph[path[path.size() - 1]][path[0]];
//
//    return cost;
//}

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
    // Get tour size
    int size = route.size();

    // repeat until no improvement is made
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
/*vector<int> tsp2opt(vector<vector<double>>& graph, vector<int>& route) {
    bool improved = true;
    vector<int> bestFoundRoute = vector<int>(route);
    double bestFoundRouteCost = routeCost(graph, bestFoundRoute);
    cout << "Fitst route " <<  bestFoundRouteCost << endl;
    int iter = 0;

    while (improved) {
        improved = false;

        for (int i = 1; i < bestFoundRoute.size() - 1; i++) {
            for (int k = i + 1; k < bestFoundRoute.size() - 1; k++) {
                vector<int> newRoute(swap2opt(bestFoundRoute, i, k));
                double  newRouteCost = routeCost(graph, newRoute);

                if (newRouteCost < bestFoundRouteCost) {
                    bestFoundRoute = newRoute;
                    bestFoundRouteCost = newRouteCost;
                    improved = true;
                }

            }

            if (improved) {
                break;
            }
        }

        iter++;
    }

    cout << iter << endl;
    return bestFoundRoute;
}*/

//        def tsp_2_opt(graph, route):
//improved = True
//best_found_route = route
//best_found_route_cost = route_cost(graph, best_found_route)
//iter = 0
//while improved and iter < 1000:
//improved = False
//for i in range(1, len(best_found_route) - 1):
//for k in range(i + 1, len(best_found_route) - 1):
//new_route = _swap_2opt(best_found_route, i, k)
//new_route_cost = route_cost(graph, new_route)
//if new_route_cost < best_found_route_cost:
//best_found_route_cost = new_route_cost
//best_found_route = new_route
//improved = True
//break
//if improved:
//break
//iter += 1
//return best_found_route