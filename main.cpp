#include <iostream>
#include <vector>
#include <random>
#include <limits>
#include <queue>
#include <stack>
#include <ctime>

using namespace std;

// Function to find the index of the vertex with the minimum combined value (distance + traffic lights)
int minCombinedValue(const vector<double>& distances, const vector<int>& trafficSignals, const vector<bool>& visited) {
    double minCombined = numeric_limits<double>::infinity();
    int minIndex = -1;

    for (int i = 0; i < distances.size(); ++i) {
        if (!visited[i] && (distances[i] + trafficSignals[i]) < minCombined) {
            minCombined = distances[i] + trafficSignals[i];
            minIndex = i;
        }
    }

    return minIndex;
}

// Function to perform Dijkstra's algorithm to find the shortest combined value path
void dijkstra(const vector<vector<double>>& distances, const vector<vector<int>>& trafficSignals, const vector<vector<bool>>& trafficDensity, int source, vector<double>& shortestDistances, vector<int>& previousVertices) {
    int numVertices = distances.size();
    vector<bool> visited(numVertices, false);

    shortestDistances.assign(numVertices, numeric_limits<double>::infinity());
    shortestDistances[source] = 0.0;
    previousVertices.assign(numVertices, -1);

    for (int count = 0; count < numVertices - 1; ++count) {
        int u = minCombinedValue(shortestDistances, trafficSignals[source], visited);
        visited[u] = true;

    for (int v = 0; v < numVertices; ++v) {
        if (!visited[v] && distances[u][v] != numeric_limits<double>::infinity() && (shortestDistances[u] + distances[u][v]) < shortestDistances[v]) {
            double extraDistance = 0.0;

            // Add 30 meters if traffic density is true
            if (trafficDensity[u][v]) {
            extraDistance += 30.0;
            }

        // Add 10 meters for each traffic signal
            extraDistance += trafficSignals[u][v] * 10.0;

            shortestDistances[v] = shortestDistances[u] + distances[u][v] + extraDistance;
            previousVertices[v] = u; // Store the previous vertex for reconstructing the path
    }
}

    }
}

// Function to calculate the total traffic signals along the path
int calculateTotalTrafficSignals(const vector<int>& previousVertices, int source, int destination, const vector<vector<int>>& trafficSignals) {
    int totalSignals = 0;
    int currentVertex = destination;

    while (currentVertex != source) {
        int previousVertex = previousVertices[currentVertex];
        totalSignals += trafficSignals[currentVertex][previousVertex];
        currentVertex = previousVertex;
    }

    return totalSignals;
}

// Function to print the optimal path
void printPath(const vector<int>& previousVertices, int source, int destination) {
    stack<int> path;
    int currentVertex = destination;

    while (currentVertex != source) {
        path.push(currentVertex);
        currentVertex = previousVertices[currentVertex];
    }

    cout << source << " -> ";
    while (!path.empty()) {
        cout << path.top();
        path.pop();
        if (!path.empty()) {
            cout << " -> ";
        }
    }
    cout << endl;
}

int main() {
    int numCoordinates;
    cout << "Enter number of locations: ";
    cin >> numCoordinates;

    // Create a random number generator
    random_device rd;
   mt19937 gen(static_cast<unsigned int>(time(0)));
    uniform_real_distribution<double> dist(1.0, 100.0);  // Adjust the range of distances as needed

    // Create the adjacency matrix filled with infinity (unconnected by default)
    vector<vector<double>> adjacencyMatrix(numCoordinates, vector<double>(numCoordinates, numeric_limits<double>::infinity()));
    vector<vector<int>> trafficSignalsMatrix(numCoordinates, vector<int>(numCoordinates, 0));
    vector<vector<bool>> trafficDensityMatrix(numCoordinates, vector<bool>(numCoordinates, false));

    // Generate random distances for connected coordinates, random traffic signals, and random traffic density
    for (int i = 0; i < numCoordinates; ++i) {
        for (int j = 0; j < numCoordinates; ++j) {
            if (i == j) {
                adjacencyMatrix[i][j] = 0.0;  // Set distance to zero when i == j
            } else if (dist(gen) < 60.0) {  // Adjust the probability as needed for the level of connectivity
                double distance = dist(gen);
                int numTrafficSignals = rand() % 5 + 1; // Random number of traffic signals (1 to 5)
                bool highTrafficDensity = rand() % 2 == 1; // Randomly set high traffic density to true or false
                adjacencyMatrix[i][j] = distance;
                adjacencyMatrix[j][i] = distance;  // Since it's an undirected graph
                trafficSignalsMatrix[i][j] = numTrafficSignals;
                trafficSignalsMatrix[j][i] = numTrafficSignals; // Traffic signals are symmetric
                trafficDensityMatrix[i][j] = highTrafficDensity;
                trafficDensityMatrix[j][i] = highTrafficDensity; // Traffic density is symmetric
            }
        }
    }

    int source, destination;
    cout << "Enter source coordinate: ";
    cin >> source;
    cout << "Enter destination coordinate: ";
    cin >> destination;

    // Print the adjacency matrix, traffic signals matrix, and traffic density matrix
    cout << "Adjacency Matrix (Distances):" << endl;
    for (int i = 0; i < numCoordinates; ++i) {
        for (int j = 0; j < numCoordinates; ++j) {
            if (adjacencyMatrix[i][j] == numeric_limits<double>::infinity()) {
                cout << "inf\t";
            } else {
                cout << adjacencyMatrix[i][j] << "\t";
            }
        }
        cout << endl;
    }

    cout << "Traffic Signals Matrix:" << endl;
    for (int i = 0; i < numCoordinates; ++i) {
        for (int j = 0; j < numCoordinates; ++j) {
            cout << trafficSignalsMatrix[i][j] << "\t";
        }
        cout << endl;
    }

    cout << "Traffic Density Matrix:" << endl;
    for (int i = 0; i < numCoordinates; ++i) {
        for (int j = 0; j < numCoordinates; ++j) {
            cout << (trafficDensityMatrix[i][j] ? "true" : "false") << "\t";
        }
        cout << endl;
    }
    

    vector<double> shortestDistances;
    vector<int> previousVertices;
    dijkstra(adjacencyMatrix, trafficSignalsMatrix, trafficDensityMatrix, source, shortestDistances, previousVertices);
     int highTrafficDensitySignals = 0;
        int currentVertex = destination;

        while (currentVertex != source) {
            int previousVertex = previousVertices[currentVertex];
            if (trafficDensityMatrix[currentVertex][previousVertex]) {
                highTrafficDensitySignals++;
            }
            currentVertex = previousVertex;
        }
    
    if (shortestDistances[destination] != numeric_limits<double>::infinity()) {
        cout << "Shortest combined value from " << source << " to " << destination << " is: " << shortestDistances[destination] << " meters" << endl;
                cout << "Shortest combined value from " << source << " to " << destination << " is: " << shortestDistances[destination] << " meters" << endl;
        cout << "Optimal path: ";
        printPath(previousVertices, source, destination);
        int totalTrafficSignals = calculateTotalTrafficSignals(previousVertices, source, destination, trafficSignalsMatrix);
        cout << "Total traffic signals on the path: " << totalTrafficSignals << endl;
        cout << "Number of traffic signals with high traffic density on the path: " << highTrafficDensitySignals << endl;
    } else {
        cout << "There is no path from " << source << " to " << destination << endl;
    }
    

    return 0;
}
