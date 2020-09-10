// ProjectCode.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <queue>
#include "fiboqueue.h"
#include <list>
#include <limits>
#include <vector>
#include <utility> // for pair
#include <chrono>
#include <tuple>
#include <random>

using namespace std;

// define chrono clock
typedef std::chrono::high_resolution_clock Clock;



// GRAPH OBJECTS - distance list and distance matrix

// structures to contain edge information
template <class T>
struct Edge {
	int source;
	pair<int, T> dest; // first is destination vertex, second is edge weight
};

// standard graph object - a weighted directed graph by default - distance list representation
template <class T>
class StandardGraph
{
public:
	vector<vector<pair<int, T>>> list; 	// a vector of vectors to represent an adjacency list
	int size; // number of vertices
	bool directed; // whether the graph is directed or undirected

	// construct the graph - passed a list of edges and number of vertices - directed unless otherwise specified
	StandardGraph(vector<Edge<T>> const &edges, int N, bool dir = true)
	{
		size = N;
		directed = dir;

		// resize the vector to N elements of type vector<int>
		list.resize(N);

		// add edges to the directed graph
		for (const Edge<T> &edge : edges)
		{
			// insert edge at the end of source vertex list
			list[edge.source].push_back(edge.dest);

			// if undirected, add reverse edge
			if (!directed)
				list[edge.dest.first].push_back(make_pair(edge.source, edge.dest.second));
		}
	}

	// print the contents of this graph to console
	void printGraph()
	{
		// iterate over every vertex
		for (int i = 0; i < size; i++)
		{
			// print current vertex number
			std::cout << " " << i << " --> ";

			// print all neighbouring vertices of vertex i
			for (const pair<int, T> &neighbour : list[i])
				std::cout << neighbour.first << "/" << neighbour.second << " ";
			std::cout << "\n";
		}
		std::cout << "\n";
	}

	// returns a graph with the same edges as this one
	StandardGraph<T> deepcopy()
	{
		vector<Edge<T>> edges;

		for (int i = 0; i < size; i++) {
			for (const pair<int, T> &edge : list[i]) {
				edges.push_back({ i, make_pair(edge.first, edge.second) });
			}
		}

		StandardGraph<T> copy(edges, size);

		return copy;
	}

	// adds an empty vertex of label N + 1
	void addVertex()
	{
		// add a new list to the parent list
		vector<pair<int, T>> vertex;
		list.push_back(vertex);

		size = size + 1;
	}

	// removes all edges associated with a vertex
	// does not actually remove the vertex label unless the vertex is the last vertex - this is because the vertex label is implied from the index
	void removeVertex(int label)
	{
		for (int i = 0; i < size; i++) { // iterate over every vertex
			vector<pair<int, T>> vertex = list[i]; // list[i] must be copied to avoid changing container size as edges are deleted
			for (pair<int, T> edge : vertex) { // and every edge leaving that vertex

				// remove all edges with this vertex as a source
				if (i == label) {
					removeEdge(label, edge.first);
				}

				// remove all edges with this vertex as a destination
				if (edge.first == label) {
					removeEdge(i, label);
				}
			}
		}

		// resize list if vertex is last vertex
		if (label + 1 == size) {
			size = size - 1;
			list.erase(list.begin() + label);
		}

	}

	// add an edge from source to dest with a weight
	void addEdge(int source, int dest, T weight)
	{
		list[source].push_back(make_pair(dest, weight));

		// if graph is undirected add the opposite edge too
		if (!directed)
			list[dest].push_back(make_pair(source, weight));
	}

	// remove edge from source to dest
	void removeEdge(int source, int dest)
	{
		list[source].erase(remove_if(
			list[source].begin(),
			list[source].end(),
			[dest](pair<int, T> edge) {return edge.first == dest; }
		));

		// if graph is undirected remove the opposite edge too
		if (!directed) {
			list[dest].erase(remove_if(
				list[dest].begin(),
				list[dest].end(),
				[source](pair<int, T> edge) {return edge.first == source; }
			));
		}
		
	}

	// reweight edge from source to dest 
	void reweightEdge(int source, int dest, T weight)
	{
		for (pair<int, T> &edge : list[source]) {
			if (edge.first == dest) {
				edge.second = weight;
			}
		}

		// if undirected reweight edge from dest to source too
		if (!directed) {
			for (pair<int, T> &edge : list[dest]) {
				if (edge.first == source) {
					edge.second = weight;
				}
			}
		}
	}
};

// object to represent a distance matrix
template <class T>
class DistanceMatrix
{
public: 
	vector<T> matrix; // 1D vector to store the 2D matrix
	int size; // number of vertices
	
	// construct a distance matrix of size N with all elements intilised to positive infinity
	DistanceMatrix(int N)
	{
		size = N;

		matrix.resize(N*N);

		// iterate over every source destination pair and set it to positive infinity
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				setEdge(i, j, std::numeric_limits<T>::max());
			}
		}
	}

	// construct a distance matrix of size N using the given list of edges - a distance matrix representing a graph
	DistanceMatrix(vector<Edge<T>> const &edges, int N)
	{
		size = N;

		matrix.resize(N*N);

		// iterate over every source destination pair and set it to positive infinity
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				setEdge(i, j, std::numeric_limits<T>::max());
			}
		}

		// add edges from the list to the matrix
		for (const Edge<T> &edge : edges) {
			setEdge(edge.source, edge.dest.first, edge.dest.second);
		}
	}

	// set element in index equivalent to [source][dest] to weight
	void setEdge(int source, int dest, T weight) {
		matrix[(source * size) + dest] = weight;
	}

	// return element in index equivalent to [source][dest]
	T getEdge(int source, int dest) const {
		return matrix[(source * size) + dest];
	}

	// copy a graph in distance list format into distance matrix format
	void fill_with_graph(const StandardGraph<T> &graph)
	{
		for (int i = 0; i < size; i++) {
			for (const pair<int, T> &dest : graph.list[i]) {
				setEdge(i, dest.first, dest.second);
			}
		}
	}

	// print the contents of this distance matrix to console - formating breaks when using elements above 2 characters
	void printMatrix()
	{
		// print header
		std::cout << "   | ";
		for (int i = 0; i < size; i++)
			std::cout << i << "  ";
		std::cout << "\n";
		std::cout << " --|-";
		for (int i = 0; i < size; i++)
			std::cout << "---";
		std::cout << "\n";

		for (int i = 0; i < size; ++i) {
			std::cout << " " << i << " | ";
			for (int j = 0; j < size; ++j) {
				// if element is infinity print a null symbol
				if (getEdge(i, j) == std::numeric_limits<T>::max()) {
					std::cout << "-  ";
				}
				// else print the element value
				else {
					std::cout << getEdge(i, j) << " ";
					// if element is only one character print another space
					if (getEdge(i, j) < 10) {
						std::cout << " ";
					}
				}
			}
			std::cout << "\n";
		}
		std::cout << "\n";
	}
};




// UTILITY FUNCTIONS - loading graphs and solutions, loading graphs from files and randomly generating graphs

// returns a pre-defined int graph with 8 vertices and 11 edges
StandardGraph<int> loadKnownGraph()
{
	// Starting graph:
	// - 7 - - - - - -
	// - - 3 3 - - - -
	// 1 - - - 5 - - -
	// - - - - - - - -
	// - 4 - 8 - - - 5
	// - - 1 - - - - -
	// - - - - - - - 3
	// - - - - - 2 - -

	vector<Edge<int>> edges =
	{
		{ 0, make_pair(1, 7) },
		{ 1, make_pair(2, 3) }, { 1, make_pair(3, 3) },
		{ 2, make_pair(0, 1) }, { 2, make_pair(4, 5) },

		{ 4, make_pair(1, 4) }, { 4, make_pair(3, 8) },	{ 4, make_pair(7, 5) },
		{ 5, make_pair(2, 1) },
		{ 6, make_pair(7, 3) },
		{ 7, make_pair(5, 2) }
	};

	// number of vertices in the graph
	int N = 8;

	// construct graph
	StandardGraph<int> graph(edges, N);

	return graph;
}

// returns a pre-defined int graph with 4 vertices and 7 edges including negative weights
StandardGraph<int> loadKnownGraphNegatives()
{
	// Starting graph:
	// -  -  -  2
	// 6  -  3  -
	// 4  -  -  5 
	// -  -7 -3 - 

	vector<Edge<int>> edges =
	{
		{ 0, make_pair(3, 2) },
		{ 1, make_pair(0, 6) }, { 1, make_pair(2, 3) },
		{ 2, make_pair(0, 4) }, { 2, make_pair(3, 5) },
		{ 3, make_pair(1, -7) }, { 3, make_pair(2, -3) }
	};

	// number of vertices in the graph
	int N = 4;

	// construct graph
	StandardGraph<int> graph(edges, N);

	return graph;
}

// returns a pre-defined int graph with 8 vertices and 11 edges
StandardGraph<int> loadKnownGraphUnweighted()
{
	// Starting graph:
	// - 1 - - - - - -
	// - - 1 1 - - - -
	// 1 - - - 1 - - -
	// - - - - - - - -
	// - 1 - 1 - - - 1
	// - - 1 - - - - -
	// - - - - - - - 1
	// - - - - - 1 - -

	vector<Edge<int>> edges =
	{
		{ 0, make_pair(1, 1) },
		{ 1, make_pair(2, 1) }, { 1, make_pair(3, 1) },
		{ 2, make_pair(0, 1) }, { 2, make_pair(4, 1) },

		{ 4, make_pair(1, 1) }, { 4, make_pair(3, 1) },	{ 4, make_pair(7, 1) },
		{ 5, make_pair(2, 1) },
		{ 6, make_pair(7, 1) },
		{ 7, make_pair(5, 1) }
	};

	// number of vertices in the graph
	int N = 8;

	// construct graph
	StandardGraph<int> graph(edges, N);

	return graph;
}

// returns the known APSP solution of graph from loadKnownGraph()
StandardGraph<int> loadKnownSolution()
{

	// Solution:
	// 0  7  10 10 15 22 - 20
	// 4  0  3  3  8  15 - 13
	// 1  8  0  11 5  12 - 10
	// -  -  -  0  -  -  -  -
	// 8  4  7  7  0  7  -  5
	// 2  9  1  12 6  0  -  11
	// 7  14 6  17 11 5  0  3
	// 4  11 3  14 8  2  -  0

	vector<Edge<int>> edges =
	{
		{ 0, make_pair(0, 0) }, { 0, make_pair(1, 7) }, { 0, make_pair(2, 10) }, { 0, make_pair(3, 10) },
		{ 0, make_pair(4, 15) }, { 0, make_pair(5, 22) }, { 0, make_pair(7, 20) },

		{ 1, make_pair(0, 4) }, { 1, make_pair(1, 0) }, { 1, make_pair(2, 3) }, { 1, make_pair(3, 3) },
		{ 1, make_pair(4, 8) }, { 1, make_pair(5, 15) }, { 1, make_pair(7, 13) },

		{ 2, make_pair(0, 1) }, { 2, make_pair(1, 8) }, { 2, make_pair(2, 0) }, { 2, make_pair(3, 11) },
		{ 2, make_pair(4, 5) }, { 2, make_pair(5, 12) }, { 2, make_pair(7, 13) },

		{ 3, make_pair(3, 0) },

		{ 4, make_pair(0, 8) }, { 4, make_pair(1, 4) }, { 4, make_pair(2, 7) }, { 4, make_pair(3, 7) },
		{ 4, make_pair(4, 0) }, { 4, make_pair(5, 7) }, { 4, make_pair(7, 5) },

		{ 5, make_pair(0, 2) }, { 5, make_pair(1, 9) }, { 5, make_pair(2, 1) }, { 5, make_pair(3, 12) },
		{ 5, make_pair(4, 6) }, { 5, make_pair(5, 0) }, { 5, make_pair(7, 11) },

		{ 6, make_pair(0, 7) }, { 6, make_pair(1, 14) }, { 6, make_pair(2, 6) }, { 6, make_pair(3, 17) },
		{ 6, make_pair(4, 11) }, { 6, make_pair(5, 5) }, { 6, make_pair(6, 0) }, { 6, make_pair(7, 3) },

		{ 7, make_pair(0, 4) }, { 7, make_pair(1, 11) }, { 7, make_pair(2, 3) }, { 7, make_pair(3, 14) },
		{ 7, make_pair(4, 8) }, { 7, make_pair(5, 2) }, { 7, make_pair(7, 0) }
	};

	// number of vertices in the graph
	int N = 8;

	// construct graph
	StandardGraph<int> graph(edges, N);

	return graph;
}

// returns the known APSP solution of graph from loadKnownGraph()
StandardGraph<int> loadKnownSolutionUnweighted()
{

	// Solution:
	// 0  1  2  2  3  5  -  4
	// 2  0  1  1  2  4  -  3
	// 1  2  0  2  1  3  -  2
	// -  -  -  0  -  -  -  -
	// 3  1  2  1  0  2  -  1
	// 2  3  1  3  2  0  -  3
	// 4  5  3  5  4  2  0  1
	// 3  4  2  4  3  1  -  0

	vector<Edge<int>> edges =
	{
		{ 0, make_pair(0, 0) }, { 0, make_pair(1, 1) }, { 0, make_pair(2, 2) }, { 0, make_pair(3, 2) },
		{ 0, make_pair(4, 3) }, { 0, make_pair(5, 5) }, { 0, make_pair(7, 4) },

		{ 1, make_pair(0, 2) }, { 1, make_pair(1, 0) }, { 1, make_pair(2, 1) }, { 1, make_pair(3, 1) },
		{ 1, make_pair(4, 2) }, { 1, make_pair(5, 4) }, { 1, make_pair(7, 3) },

		{ 2, make_pair(0, 1) }, { 2, make_pair(1, 2) }, { 2, make_pair(2, 0) }, { 2, make_pair(3, 2) },
		{ 2, make_pair(4, 1) }, { 2, make_pair(5, 3) }, { 2, make_pair(7, 2) },

		{ 3, make_pair(3, 0) },

		{ 4, make_pair(0, 3) }, { 4, make_pair(1, 1) }, { 4, make_pair(2, 2) }, { 4, make_pair(3, 1) },
		{ 4, make_pair(4, 0) }, { 4, make_pair(5, 2) }, { 4, make_pair(7, 1) },

		{ 5, make_pair(0, 2) }, { 5, make_pair(1, 3) }, { 5, make_pair(2, 1) }, { 5, make_pair(3, 3) },
		{ 5, make_pair(4, 2) }, { 5, make_pair(5, 0) }, { 5, make_pair(7, 3) },

		{ 6, make_pair(0, 4) }, { 6, make_pair(1, 5) }, { 6, make_pair(2, 3) }, { 6, make_pair(3, 5) },
		{ 6, make_pair(4, 4) }, { 6, make_pair(5, 2) }, { 6, make_pair(6, 0) }, { 6, make_pair(7, 1) },

		{ 7, make_pair(0, 3) }, { 7, make_pair(1, 4) }, { 7, make_pair(2, 2) }, { 7, make_pair(3, 4) },
		{ 7, make_pair(4, 3) }, { 7, make_pair(5, 1) }, { 7, make_pair(7, 0) }
	};

	// number of vertices in the graph
	int N = 8;

	// construct graph
	StandardGraph<int> graph(edges, N);

	return graph;
}

// returns the known APSP solution of graph from loadKnownGraphNegatives()
StandardGraph<int> loadKnownSolutionNegatives()
{

	// Solution:
	// 0  -5 -2  2
	// 6  0  3   8
	// 4  -2 0   5 
	// -1 -7 -4  0  

	vector<Edge<int>> edges =
	{
		{ 0, make_pair(0, 0) }, { 1, make_pair(1, 0) }, { 2, make_pair(2, 0) }, { 3, make_pair(3, 0) },

		{ 0, make_pair(1, -5) }, { 0, make_pair(2, -2) }, { 0, make_pair(3, 2) },
		{ 1, make_pair(0, 6) }, { 1, make_pair(2, 3) }, { 1, make_pair(3, 8) },
		{ 2, make_pair(0, 4) }, { 2, make_pair(1, -2) }, { 2, make_pair(3, 5) },
		{ 3, make_pair(0, -1) }, { 3, make_pair(1, -7) }, { 3, make_pair(2, -4) },
	};

	// number of vertices in the graph
	int N = 4;

	// construct graph
	StandardGraph<int> graph(edges, N);

	return graph;
}

// builds and returns a graph contained in a file - takes a file path - for file format see writeup
StandardGraph<double> loadRealGraph(string filename)
{
	vector<Edge<double>> edges;

	int N = 0;
	double maxWeight = std::numeric_limits<double>::min();
	double minWeight = std::numeric_limits<double>::max();
	
	string line;
	ifstream file(filename);
	if (file.is_open())
	{
		while (getline(file, line))
		{
			istringstream iss(line);
			vector<string> tokens{ istream_iterator<string>{iss}, istream_iterator<string>{} };
			int source = stoi(tokens[0]) - 1; // - 1 because vertices in the source file are labelled 1 to N but we need 0 to N - 1
			int dest = stoi(tokens[1]) - 1;
			double weight = stoi(tokens[2]);
			edges.push_back({ source, make_pair(dest, weight) });

			// find the max vertex label
			if (source + 1 > N)
				N = source + 1;
			if (dest + 1 > N)
				N = dest + 1;

			// find the minimum edge weight
			if (weight < minWeight)
				minWeight = weight;
			// find the maximum edge weight
			if (weight > maxWeight)
				maxWeight = weight;
		}
		file.close();
	}
	else
		std::cout << " Unable to open file";

	// construct graph
	StandardGraph<double> graph(edges, N);

	std::cout << "\n Graph loaded:";
	std::cout << "\n  - Number of vertices: " << N;
	std::cout << "\n  - Number of edges: " << edges.size();
	std::cout << "\n  - Minimum edge weight: " << minWeight;
	std::cout << "\n  - Maximum edge weight: " << maxWeight << "\n";

	return graph;
}

// builds and returns a half the graph contained in a file - used to test structure of graphs too large, does not find the full APSP solution
// requires knowning the max vertex label, if lower is true we add every vertex < max_vertex/2
StandardGraph<double> loadRealGraphSplit(string filename, int max_vertex, bool lower)
{
	vector<Edge<double>> edges;

	int middle = max_vertex / 2;
	int N = 0; // - 1 because vertices in the source file are labelled 1 to N but we need 0 to N - 1
	double maxWeight = std::numeric_limits<double>::min();
	double minWeight = std::numeric_limits<double>::max();

	string line;
	ifstream file(filename);
	if (file.is_open())
	{
		while (getline(file, line))
		{
			istringstream iss(line);
			vector<string> tokens{ istream_iterator<string>{iss}, istream_iterator<string>{} };
			int source = stoi(tokens[0]) - 1; // - 1 because vertices in the source file are labelled 1 to N but we need 0 to N - 1
			int dest = stoi(tokens[1]) - 1;
			double weight = stod(tokens[2]);

			if (lower == false && source < middle && dest < middle) {

				edges.push_back({ source, make_pair(dest, weight) });

				// find the minimum edge weight
				if (weight < minWeight)
					minWeight = weight;
				// find the maximum edge weight
				if (weight > maxWeight)
					maxWeight = weight;

				// find the max vertex label
				if (source + 1 > N)
					N = source + 1;
				if (dest + 1 > N)
					N = dest + 1;
			}
			else if (lower == true && source > middle && dest > middle) {

				edges.push_back({ source - middle, make_pair(dest - middle, weight) });

				// find the minimum edge weight
				if (weight < minWeight)
					minWeight = weight;
				// find the maximum edge weight
				if (weight > maxWeight)
					maxWeight = weight;

				// find the max vertex label
				if (source - middle + 1 > N)
					N = source - middle + 1;
				if (dest - middle + 1 > N)
					N = dest - middle + 1;
			}
		}
		file.close();
	}
	else
		std::cout << " Unable to open file";

	// construct graph
	StandardGraph<double> graph(edges, N);

	std::cout << "\n Graph loaded:";
	std::cout << "\n  - Number of vertices: " << N;
	std::cout << "\n  - Number of edges: " << edges.size();
	std::cout << "\n  - Minimum edge weight: " << minWeight;
	std::cout << "\n  - Maximum edge weight: " << maxWeight << "\n";

	return graph;
}

// object to build random graphs
class GraphGenerator {
private:
	std::default_random_engine r_eng;

public:
	// return a random into between min and max inclusive
	int getRandomInt(int min, int max) {
		std::uniform_int_distribution<int> uniform_dist(min, max);
		return uniform_dist(r_eng);
	}

	// return a random double between min and max inclusive
	double getRandomDouble(double min, double max) {
		std::uniform_real_distribution<double> uniform_dist(min, max);
		return uniform_dist(r_eng);
	}

	// return a random double from the gaussian distribution defined by the mean and standard deviation
	double getGaussianDouble(double mean, double stand_dev) {
		std::normal_distribution<double> distribution(mean, stand_dev);
		return distribution(r_eng);
	}

	// return a directed integer weighted graph with the input parameters
	// size must be an integer greater than 0, density must be a probablity value [0, 1] and min must be less than max
	StandardGraph<int> generateERGraphInt(int size, double density, int minWeight, int maxWeight)
	{
		vector<Edge<int>> edges;

		// consider every possible edge in the graph
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				// if edge is not to itself and independent probability is less than p, add edge
				if ((i != j) && (getRandomDouble(0, 1) < density)) {
					int weight = getRandomInt(minWeight, maxWeight);
					edges.push_back({ i, make_pair(j, weight) });
				}
			}
		}

		StandardGraph<int> graph(edges, size);

		std::cout << "  - Graph created \n";

		return graph;
	}

	// return a directed float weighted graph with the input parameters
	// size must be an integer greater than 0, density must be a probablity value [0, 1] and min must be less than max
	StandardGraph<float> generateERGraphFloat(int size, double density, double minWeight, double maxWeight)
	{
		vector<Edge<float>> edges;

		// consider every possible edge in the graph
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				// if edge is not to itself and independent probability is less than p, add edge
				if ((i != j) && (getRandomDouble(0, 1) < density)) {
					float weight = getRandomDouble(minWeight, maxWeight);
					edges.push_back({ i, make_pair(j, weight) });
				}
			}
		}

		StandardGraph<float> graph(edges, size);

		std::cout << "  - Graph created \n";

		return graph;
	}

	// return a directed double weighted graph with the input parameters
	// size must be an integer greater than 0, density must be a probablity value [0, 1] and min must be less than max
	StandardGraph<double> generateERGraphDouble(int size, double density, double minWeight, double maxWeight)
	{
		vector<Edge<double>> edges;

		// consider every possible edge in the graph
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				// if edge is not to itself and independent probability is less than p, add edge
				if ((i != j) && (getRandomDouble(0, 1) < density)) {
					double weight = getRandomDouble(minWeight, maxWeight);
					edges.push_back({ i, make_pair(j, weight) });
				}
			}
		}

		StandardGraph<double> graph(edges, size);

		std::cout << "  - Graph created \n";

		return graph;
	}

	// return a directed double weighted graph with the input parameters
	// size must be an integer greater than 0, density must be a probablity value [0, 1]
	StandardGraph<double> generateERGraphGaussian(int size, double density, double mean, double stand_dev)
	{
		vector<Edge<double>> edges;

		// consider every possible edge in the graph
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				// if edge is not to itself and independent probability is less than p, add edge
				if ((i != j) && (getRandomDouble(0, 1) < density)) {
					double weight = getGaussianDouble(mean, stand_dev);
					edges.push_back({ i, make_pair(j, weight) });
				}
			}
		}

		StandardGraph<double> graph(edges, size);

		std::cout << "  - Graph created \n";

		return graph;
	}
};




// ALGORITHMS - APSP algorithms

// comparator for determining priority for priority queue (shortest edge comes first) for dijkstra based algorithms
template <class D>
struct prioritize {
public:
	bool operator()(const pair<int, D> &p1, const pair<int, D> &p2)
	{
		return p1.second > p2.second;
	}
};

// object to perform binary Dijkstra's and Fibonacci Dijkstra's
template <class Dtype>
struct Dijkstra {
private:
	priority_queue<pair<int, Dtype>, vector<pair<int, Dtype>>, prioritize<Dtype>> queue; // priority queue for binary Dijkstra's
	FibQueue<pair<int,Dtype>> fib_queue; // fibonacci queue for Fibonacci Dijkstra's

public:

	// BINARY DIJKSTRA FUNCTIONS

	// takes distance vector and fills vector with distances from source - distance list overload
	void SSSP_dijkstra(const StandardGraph<Dtype> &graph, int source, vector<Dtype> &distances)
	{
		int N = graph.size;

		vector<bool> visited(N, false); // determines whether the vertex has been visited or not 

		distances[source] = 0; // set the distance of the source to itself to 0
		queue.push(make_pair(source, distances[source])); // push the source to the queue

		while (!queue.empty())
		{
			pair<int, Dtype> current_pair = queue.top(); // get the next vertex
			queue.pop(); // remove the vertex from the queue

			int current_vertex = current_pair.first;
			Dtype current_distance = current_pair.second;

			if (visited[current_vertex]) // if the vertex is already visited, go to the next vertex
				continue;

			visited[current_vertex] = true; // set the the vertex to visited

			// iterate through all adjacent vertices 
			for (int i = 0; i < graph.list[current_vertex].size(); i++) {

				int adjacent_vertex = graph.list[current_vertex][i].first;

				// if the adjacent vertex is not visited then continue
				if (!visited[adjacent_vertex]) {

					Dtype adjacent_weight = graph.list[current_vertex][i].second;

					// if the current vertex distance + distance from the current vertex to the adjacent vertex
					// is shorter than the shortest known distance to the adjacent vertex then continue
					if (current_distance + adjacent_weight < distances[adjacent_vertex])
					{
						distances[adjacent_vertex] = adjacent_weight + current_distance; // replace distance to the adjacent vertex

						// create a pair that is the adjacent vertex and the new distance to that vertex
						pair<int, Dtype> pair = make_pair(adjacent_vertex, (distances[adjacent_vertex]));
						queue.push(pair); // add the pair to the queue
					}
				}
			}
		}
		
		// ensure the queue is empty for reuse - create new queue because priority queues can't be cleared for some reason
		queue = priority_queue<pair<int, Dtype>, vector<pair<int, Dtype>>, prioritize<Dtype>>();
	}

	// takes distance vector and fills vector with distances from source - distance matrix overload
	void SSSP_dijkstra(const DistanceMatrix<Dtype> &graph, int source, vector<Dtype> &distances)
	{
		int N = graph.size;

		vector<bool> visited(N, false); // determines whether the vertex has been visited or not 

		distances[source] = 0; // set the distance of the source to itself to 0
		queue.push(make_pair(source, distances[source])); // push the source to the queue

		while (!queue.empty())
		{
			pair<int, Dtype> current_pair = queue.top(); // get the next vertex
			queue.pop(); // remove the vertex from the queue

			int current_vertex = current_pair.first;
			Dtype current_distance = current_pair.second;

			if (visited[current_vertex]) // if the vertex is already visited, go to the next vertex
				continue;

			visited[current_vertex] = true; // set the the vertex to visited

			// iterate through all other vertices 
			for (int adjacent_vertex = 0; adjacent_vertex < N; adjacent_vertex++) {

				if (visited[adjacent_vertex]) // if the adjacent vertex is already visited go to the next adjacent vertex
					continue;

				Dtype adjacent_weight = graph.getEdge(current_vertex, adjacent_vertex); // get the weight of the possible edge to the other vertex

				if (adjacent_weight == std::numeric_limits<Dtype>::max()) // if the weight is infinity the edge does not exist
					continue;
				
				// if the current vertex distance + distance from the current vertex to the adjacent vertex
				// is shorter than the shortest known distance to the adjacent vertex then continue
				if (current_distance + adjacent_weight < distances[adjacent_vertex])
				{
					distances[adjacent_vertex] = adjacent_weight + current_distance; // replace that distance

					// create a pair that is the adjacent vertex and the new distance to that vertex
					pair<int, Dtype> pair = make_pair(adjacent_vertex, (distances[adjacent_vertex]));
					queue.push(pair); // add the pair to the queue
				}
			}
		}

		// ensure the queue is empty for reuse - create new queue because priority queues can't be cleared for some reason
		queue = priority_queue<pair<int, Dtype>, vector<pair<int, Dtype>>, prioritize<Dtype>>();
	}

	// returns APSP distance matrix solution of the graph - distance list overload
	DistanceMatrix<Dtype> APSP_dijkstra(const StandardGraph<Dtype> &graph)
	{
		int N = graph.size;

		DistanceMatrix<Dtype> output(N); // container for the APSP solution

		for (int i = 0; i < N; i++) {

			vector<Dtype> distances(N, std::numeric_limits<Dtype>::max()); // stores shortest distances, initial distances of infinity
			SSSP_dijkstra(graph, i, distances); // get SSSP solution for vertex i

			for (int j = 0; j < N; j++) { // fill SSSP solution of vertex i into the APSP solution
				if (distances[j] != std::numeric_limits<Dtype>::max())
					output.setEdge(i, j, distances[j]);
			}
		}

		return output;
	}

	// returns APSP distance matrix solution of the graph - distance matrix overload
	DistanceMatrix<Dtype> APSP_dijkstra(const DistanceMatrix<Dtype> &graph)
	{
		int N = graph.size;

		DistanceMatrix<Dtype> output(N); // container for the APSP solution

		for (int i = 0; i < N; i++) {

			vector<Dtype> distances(N, std::numeric_limits<Dtype>::max()); // stores shortest distances, initial distances of infinity
			SSSP_dijkstra(graph, i, distances); // get SSSP solution for vertex i

			for (int j = 0; j < N; j++) { // fill SSSP solution of vertex i into the APSP solution
				output.setEdge(i, j, distances[j]);
			}
		}

		return output;
	}



	// FIBONACCI DIJKSTRA FUNCTIONS

	// takes distance vector and fills vector with distances from source - distance list overload
	void SSSP_dijkstra_fibqueue(const StandardGraph<Dtype> &graph, int source, vector<Dtype> &distances)
	{
		int N = graph.size;

		vector<bool> visited(N, false); // determines whether the vertex has been visited or not 

		distances[source] = 0; // set the distance of the source to itself to 0
		fib_queue.push(make_pair(source, distances[source])); // push the source to the queue

		while (!fib_queue.empty())
		{
			pair<int, Dtype> current_pair = fib_queue.top(); // get the next vertex
			fib_queue.pop(); // remove the vertex from the queue

			int current_vertex = current_pair.first;
			Dtype current_distance = current_pair.second;

			if (visited[current_vertex]) // if the vertex is already visited, go to the next vertex
				continue;
			visited[current_vertex] = true; // set the the vertex to visited

			// iterate through all adjacent vertices 
			for (int i = 0; i < graph.list[current_vertex].size(); i++) {

				int adjacent_vertex = graph.list[current_vertex][i].first;

				// if the adjacent vertex is not visited, go to the next adjacent vertex
				if (!visited[adjacent_vertex]) {

					Dtype adjacent_weight = graph.list[current_vertex][i].second;

					// if the current vertex distance + distance from the current vertex to the adjacent vertex
					// is shorter than the shortest known distance to the adjacent vertex then continue
					if (current_distance + adjacent_weight < distances[adjacent_vertex])
					{
						distances[adjacent_vertex] = adjacent_weight + current_distance; // replace distance to the adjacent vertex

						// create a pair that is the adjacent vertex and the new distance to that vertex
						pair<int, Dtype> pair = make_pair(adjacent_vertex, (distances[adjacent_vertex]));
						fib_queue.push(pair); // add the pair to the queue
					}
				}
			}
		}

		// ensure the queue is empty for reuse
		fib_queue.clear();
	}

	// takes distance vector and fills vector with distances from source - distance matrix overload
	void SSSP_dijkstra_fibqueue(const DistanceMatrix<Dtype> &graph, int source, vector<Dtype> &distances) 
	{
		int N = graph.size;

		vector<bool> visited(N, false); // determines whether the vertex has been visited or not 

		distances[source] = 0; // set the distance of the source to itself to 0
		fib_queue.push(make_pair(source, distances[source])); // push the source to the queue

		while (!fib_queue.empty())
		{
			pair<int, Dtype> current_pair = fib_queue.top(); // get the next vertex
			fib_queue.pop(); // remove the vertex from the queue

			int current_vertex = current_pair.first;
			Dtype current_distance = current_pair.second;

			if (visited[current_vertex]) // if the vertex is already visited, go to the next vertex
				continue;

			visited[current_vertex] = true; // set the the vertex to visited

			// iterate through all other vertices 
			for (int adjacent_vertex = 0; adjacent_vertex < N; adjacent_vertex++) {

				if (visited[adjacent_vertex]) // if the adjacent vertex is already visited go to the next adjacent vertex
					continue;

				Dtype adjacent_weight = graph.getEdge(current_vertex, adjacent_vertex); // get the weight of the possible edge to the other vertex

				if (adjacent_weight == std::numeric_limits<Dtype>::max()) // if the weight is infinity the edge does not exist
					continue;

				// if the current vertex distance + distance from the current vertex to the adjacent vertex
				// is shorter than the shortest known distance to the adjacent vertex then continue
				if (current_distance + adjacent_weight < distances[adjacent_vertex])
				{
					distances[adjacent_vertex] = adjacent_weight + current_distance; // replace that distance

					// create a pair that is the adjacent vertex and the new distance to that vertex
					pair<int, Dtype> pair = make_pair(adjacent_vertex, (distances[adjacent_vertex]));
					fib_queue.push(pair); // add the pair to the queue
				}
			}
		}

		// ensure the queue is empty for reuse
		fib_queue.clear();
	}

	// returns APSP distance matrix solution of the graph - distance list overload
	DistanceMatrix<Dtype> APSP_dijkstra_fibqueue(const StandardGraph<Dtype> &graph)
	{
		int N = graph.size;

		DistanceMatrix<Dtype> output(N); // container for the APSP solution

		for (int i = 0; i < N; i++) {

			vector<Dtype> distances(N, std::numeric_limits<Dtype>::max()); // stores shortest distances, initial distances of infinity
			SSSP_dijkstra_fibqueue(graph, i, distances); // get SSSP solution for vertex i

			for (int j = 0; j < N; j++) { // fill SSSP solution of vertex i into the APSP solution
				output.setEdge(i, j, distances[j]);
			}
		}

		return output;
	}

	// returns APSP distance matrix solution of the graph - distance matrix overload
	DistanceMatrix<Dtype> APSP_dijkstra_fibqueue(const DistanceMatrix<Dtype> &graph)
	{
		int N = graph.size;

		DistanceMatrix<Dtype> output(N); // container for the APSP solution

		for (int i = 0; i < N; i++) {

			vector<Dtype> distances(N, std::numeric_limits<Dtype>::max()); // stores shortest distances, initial distances of infinity
			SSSP_dijkstra_fibqueue(graph, i, distances); // get SSSP solution for vertex i

			for (int j = 0; j < N; j++) { // fill SSSP solution of vertex i into the APSP solution
				output.setEdge(i, j, distances[j]);
			}
		}

		return output;
	}
};

// object to perform Breadth-First search
struct BreadthFirst {
private:
	queue<int> queue;

public:

	// takes distance vector and fills vector with distances from source
	void SSSP_breadth_first(const StandardGraph<int> &graph, int source, vector<int> &distances)
	{
		int N = graph.size;

		vector<bool> visited(N, false); // determines whether the vertex has been visited or not 

		visited[source] = true; // mark the source as visited
		distances[source] = 0; // set the distance of the source to itself to 0
		queue.push(source); // push the source to the queue

		while (!queue.empty())
		{
			int current_vertex = queue.front(); // get the next vertex
			queue.pop(); // remove the vertex from the queue

			// iterate through all adjacent vertices 
			for (int i = 0; i < graph.list[current_vertex].size(); i++) {

				int adjacent_vertex = graph.list[current_vertex][i].first;

				// if the adjacent vertex is not visited then continue
				if (!visited[adjacent_vertex]) {
					visited[adjacent_vertex] = true; // mark the adjacent vertex as visited
					distances[adjacent_vertex] = distances[current_vertex] + 1; // set the distance to the adjacent vertex

					queue.push(adjacent_vertex); // add the adjacent vertex to the queue
				}
			}
		}
	}

	// returns APSP distance matrix solution of the graph - distance list overload
	DistanceMatrix<int> APSP_breadth_first(const StandardGraph<int> &graph)
	{
		int N = graph.size;

		DistanceMatrix<int> output(N); // container for the APSP solution

		for (int i = 0; i < N; i++) {

			vector<int> distances(N, std::numeric_limits<int>::max()); // stores shortest distances, initial distances of infinity
			SSSP_breadth_first(graph, i, distances); // get SSSP solution for vertex i

			for (int j = 0; j < N; j++) { // fill SSSP solution of vertex i into the APSP solution
				output.setEdge(i, j, distances[j]);
			}
		}

		return output;
	}
};

// object to perform BellmanFord
template <class Dtype>
struct BellmanFord {
public:

	// takes distance vector and fills vector with distances from source - distance list overload
	void SSSP_bellman_ford(const StandardGraph<Dtype> &graph, int source, vector<Dtype> &distances)
	{
		int N = graph.size;

		distances[source] = 0; // set the distance of the source to itself as 0

		for (int i = 1; i <= N - 1; i++) { // repeat N - 1 times
			for (int j = 0; j < N; j++) { // iterate over every vertex

				// if tentative distance to j is not known, go to next vertex
				if (distances[j] == std::numeric_limits<Dtype>::max())
					continue;

				// relax every edge where j is the tail
				for (const pair<int, Dtype> &edge : graph.list[j]) { 
					int dest = edge.first;
					Dtype weight = edge.second;

					// if distance to j plus distance from j to the adjacent vertex
					// is less than tentative distance to adjacent vertex then update that distance
					if (distances[j] + weight < distances[dest]) {
						distances[dest] = distances[j] + weight;
					}
				}
			}
		}
	}

	// takes distance vector and fills vector with distances from source - distance matrix overload
	void SSSP_bellman_ford(const DistanceMatrix<Dtype> &graph, int source, vector<Dtype> &distances)
	{
		int N = graph.size;

		distances[source] = 0; // set the distance of the source to itself as 0

		for (int i = 1; i <= N - 1; i++) { // repeat N - 1 times
			for (int j = 0; j < N; j++) { // iterate over every vertex

				// if tentative distance to j is not known, go to next vertex
				if (distances[j] == std::numeric_limits<Dtype>::max())
					continue;

				for (int dest = 0; dest < N; dest++) { // iterate over every possible edge leaving j

					Dtype weight = graph.getEdge(j, dest);

					// if distance from j to dest is infinity that edge does not exist, so go to next edge
					if (weight == std::numeric_limits<Dtype>::max())
						continue;

					// if distance to j plus distance from j to the adjacent vertex
					// is less than tentative distance to adjacent vertex then update that distance
					if (distances[j] + weight < distances[dest]) {
						distances[dest] = distances[j] + weight;
					}
				}
			}
		}
	}

	// returns APSP distance matrix solution of the graph - distance list overload
	DistanceMatrix<Dtype> APSP_bellman_ford(const StandardGraph<Dtype> &graph)
	{
		int N = graph.size;

		DistanceMatrix<Dtype> output(N); // container for the APSP solution

		for (int i = 0; i < N; i++) {

			vector<Dtype> distances(N, std::numeric_limits<Dtype>::max()); // stores shortest distances, initial distances of infinity
			SSSP_bellman_ford(graph, i, distances); // get SSSP solution for vertex i

			for (int j = 0; j < N; j++) { // fill SSSP solution of vertex i into the APSP solution
				output.setEdge(i, j, distances[j]);
			}
		}

		return output;
	}

	// returns APSP distance matrix solution of the graph - distance matrix overload
	DistanceMatrix<Dtype> APSP_bellman_ford(const DistanceMatrix<Dtype> &graph)
	{
		int N = graph.size;

		DistanceMatrix<Dtype> output(N); // container for the APSP solution

		for (int i = 0; i < N; i++) {

			vector<Dtype> distances(N, std::numeric_limits<Dtype>::max()); // stores shortest distances, initial distances of infinity
			SSSP_bellman_ford(graph, i, distances); // get SSSP solution for vertex i

			for (int j = 0; j < N; j++) { // fill SSSP solution of vertex i into the APSP solution
				output.setEdge(i, j, distances[j]);
			}
		}

		return output;
	}
};

// object to perform Johnson's Algorithm
template <class Dtype>
struct JohnsonDijkstra {
public:
	// returns APSP distance matrix solution of the graph - distance list overload
	DistanceMatrix<Dtype> johnson_dijkstra(StandardGraph<Dtype> &graph)
	{
		// copy the graph. does not need to be done if we are not using the same graph for other algorithms
		StandardGraph<Dtype> copy = graph.deepcopy();

		// add new vertex q (vertex N - 1 in the new graph)
		copy.addVertex();
		int q = copy.size - 1;
		for (int i = 0; i < q; i++) { // add an edge with weight 0 from q to every other vertex
			copy.addEdge(q, i, 0);
		}

		// run SSSP bellman ford from q to get distance array h
		BellmanFord<Dtype> bellman_ford;
		vector<Dtype> h(copy.size, std::numeric_limits<Dtype>::max());
		bellman_ford.SSSP_bellman_ford(copy, q, h);

		// remove q from the graph
		copy.removeVertex(q);

		// reweight every edge
		for (int i = 0; i < copy.size; i++) { // iterate over every vertex
			for (const pair<int, Dtype> &edge : copy.list[i]) { // and every edge leaving that vertex
				copy.reweightEdge(i, edge.first, edge.second + h[i] - h[edge.first]);
			}
		}

		// run APSP dijkstra
		Dijkstra<Dtype> dijkstra;
		DistanceMatrix<Dtype> output = dijkstra.APSP_dijkstra(copy);

		// reverse the reweighting to get the distances in the original graph
		for (int i = 0; i < output.size; i++) { // iterate over every possible edge
			for (int j = 0; j < output.size; j++) { 
				// if that edge exists, reweight it
				if (output.getEdge(i, j) != std::numeric_limits<Dtype>::max()) { 
					output.setEdge(i, j, output.getEdge(i, j) - h[i] + h[j]);
				}
			}
		}

		return output;
	}

	// returns APSP distance matrix solution of the graph - distance matrix overload
	DistanceMatrix<Dtype> johnson_dijkstra(DistanceMatrix<Dtype> &graph)
	{
		// create a copy of the same size of graph with an extra vertex q
		DistanceMatrix<Dtype> copy(graph.size + 1);
		int q = copy.size - 1; // q is the last row and last column in the new distance matrix

		// fill the copy with the edges in the graph
		for (int i = 0; i < graph.size; i++) {
			for (int j = 0; j < graph.size; j++) {
				// if that edge exists, copy it
				if (graph.getEdge(i, j) != std::numeric_limits<Dtype>::max()) {
					copy.setEdge(i, j, graph.getEdge(i, j));
				}
			}
		}

		// add an edge with weight 0 from q to every other vertex
		for (int i = 0; i < q; i++) { 
			copy.setEdge(q, i, 0);
		}

		// run SSSP bellman ford from q to get distance array h
		BellmanFord<Dtype> bellman_ford;
		vector<Dtype> h(copy.size, std::numeric_limits<Dtype>::max());
		bellman_ford.SSSP_bellman_ford(copy, q, h);

		// remove q from the graph - we need to make another copy of the graph without q
		DistanceMatrix<Dtype> copy2(copy.size - 1);

		// fill the second copy with the edges in the copy
		for (int i = 0; i < copy2.size; i++) {
			for (int j = 0; j < copy2.size; j++) {
				copy2.setEdge(i, j, copy.getEdge(i, j));
			}
		}

		// reweight every edge in the second copy
		for (int i = 0; i < copy2.size; i++) { // iterate over every possible edge
			for (int j = 0; j < copy2.size; j++) {
				// if that edge exists, reweight it
				if (copy2.getEdge(i, j) != std::numeric_limits<Dtype>::max()) {
					copy2.setEdge(i, j, copy2.getEdge(i, j) + h[i] - h[j]);
				}
			}
		}

		// run APSP dijkstra
		Dijkstra<Dtype> dijkstra;
		DistanceMatrix<Dtype> output = dijkstra.APSP_dijkstra(copy2);

		// reverse the reweighting to get the distances in the original graph
		for (int i = 0; i < output.size; i++) { // iterate over every possible edge
			for (int j = 0; j < output.size; j++) {
				// if that edge exists, reweight it
				if (output.getEdge(i, j) != std::numeric_limits<Dtype>::max()) {
					output.setEdge(i, j, output.getEdge(i, j) - h[i] + h[j]);
				}
			}
		}

		return output;
	}
};

// object to perform Floyd-Warshall algorithm
template <class Dtype>
struct FloydWarshall {
public:
	// returns APSP distance matrix solution of the graph - distance list overload
	DistanceMatrix<Dtype> floyd_warshall(const StandardGraph<Dtype> &graph)
	{
		int N = graph.size;

		// create empty distance matrix
		DistanceMatrix<Dtype> output(N);

		// copy input graph into distance matrix
		output.fill_with_graph(graph);

		// also set the distance of every vertex to itself as 0
		for (int i = 0; i < N; i++) {
			output.setEdge(i, i, 0);
		}

		for (int k = 0; k < N; k++) {
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					// get the length of the path from i to j through k
					Dtype summed = output.getEdge(i, k) + output.getEdge(k, j);

					// catch additon overflow
					if (output.getEdge(k, j) > 0 && output.getEdge(i, k) > std::numeric_limits<Dtype>::max() - output.getEdge(k, j)) {
						summed = std::numeric_limits<Dtype>::max();
					}

					// if i to j through k is less than the tentative path from i to j, update that distance
					if (summed < output.getEdge(i, j)) {
						output.setEdge(i, j, summed);
					}
				}
			}
		}

		return output;
	}

	// returns APSP distance matrix solution of the graph - distance matrix overload
	DistanceMatrix<Dtype> floyd_warshall(const DistanceMatrix<Dtype> &graph)
	{
		int N = graph.size;

		// create empty distance matrix
		DistanceMatrix<Dtype> output(N);

		// copy distance matrix into output. copying does not need to be done if the orignal graph is not going to be reused
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				// also set the distance of every vertex to itself as 0
				if (i == j) {
					output.setEdge(i, j, 0);
				}

				// if the edge exists, copy it
				if (graph.getEdge(i, j) != std::numeric_limits<Dtype>::max()) {
					output.setEdge(i, j, graph.getEdge(i, j));
				}
			}
		}

		for (int k = 0; k < N; k++) {
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					// get the length of the path from i to j through k
					Dtype summed = output.getEdge(i, k) + output.getEdge(k, j);

					// catch additon overflow
					if (output.getEdge(k, j) > 0 && output.getEdge(i, k) > std::numeric_limits<Dtype>::max() - output.getEdge(k, j)) {
						summed = std::numeric_limits<Dtype>::max();
					}

					// if i to j through k is less than the tentative path from i to j, update that distance
					if (summed < output.getEdge(i, j)) {
						output.setEdge(i, j, summed);
					}
				}
			}
		}

		return output;
	}
};





// RUNNING THE PROGRAM - get input from user, verify algorithm correctness, verify graph generation and run batches

// fill variables for random graph generation with input from user - weights of uniform range between min and max
void get_parameters(int &size, double &density, double &minweight, double &maxweight, int &batch_size) {
	std::cout << "\n Number of vertices (positive integer): ";
	cin >> size;

	std::cout << "\n Density (0 to 1): ";
	cin >> density;

	std::cout << "\n Min weight: ";
	cin >> minweight;

	std::cout << "\n Max weight: ";
	cin >> maxweight;

	std::cout << "\n Batch size (positive integer): ";
	cin >> batch_size;
}

// fill variables for random graph generation with input from user - weights of gaussian distribution centered on mean
void get_parameters_gaussian(int &size, double &density, double &mean, double &sd, int &batch_size) {
	std::cout << "\n Number of vertices (positive integer): ";
	cin >> size;

	std::cout << "\n Density (0 to 1): ";
	cin >> density;

	std::cout << "\n Mean: ";
	cin >> mean;

	std::cout << "\n Standard Deviation: ";
	cin >> sd;

	std::cout << "\n Batch size (positive integer): ";
	cin >> batch_size;
}

// fill variables for random graph generation with input from user - no weight distribution
void get_parameters_no_weights(int &size, double &density, int &batch_size) {
	std::cout << "\n Number of vertices (positive integer): ";
	cin >> size;

	std::cout << "\n Density (0 to 1): ";
	cin >> density;

	std::cout << "\n Batch size (positive integer): ";
	cin >> batch_size;
}


// print verification that the algorithms produce the correct result for a known integer weighted graph
void verify_algorithms()
{
	// declare algorithm structures
	Dijkstra<int> dijkstra_int;
	BellmanFord<int> bellman_ford_int;
	JohnsonDijkstra<int> johnson_dijkstra_int;
	FloydWarshall<int> floyd_warshall_int;

	// load graphs
	StandardGraph<int> graph = loadKnownGraph();
	StandardGraph<int> solution = loadKnownSolution();;

	std::cout << "\n VERIFYING ALGORITHMS ON GRAPHS IN DISTANCE LIST FORMAT\n";

	std::cout << " Graph in distance list form: \n";
	graph.printGraph();

	std::cout << " Known solution of graph: \n";
	DistanceMatrix<int> solution_dist(solution.size);
	solution_dist.fill_with_graph(solution);
	solution_dist.printMatrix();

	std::cout << " APSP Dijkstra binary queue solution: \n";
	DistanceMatrix<int> dijkstra_output = dijkstra_int.APSP_dijkstra(graph);
	dijkstra_output.printMatrix();

	std::cout << " APSP Dijkstra fibonaci queue solution: \n";
	DistanceMatrix<int> dijkstra_output_fib = dijkstra_int.APSP_dijkstra_fibqueue(graph);
	dijkstra_output_fib.printMatrix();

	std::cout << " APSP Bellman-Ford solution: \n";
	DistanceMatrix<int> bellman_output = bellman_ford_int.APSP_bellman_ford(graph);
	bellman_output.printMatrix();

	std::cout << " Johnson-Dijkstra solution: \n";
	DistanceMatrix<int> johnson_output = johnson_dijkstra_int.johnson_dijkstra(graph);
	johnson_output.printMatrix();

	std::cout << " Floyd-Warshall solution: \n";
	DistanceMatrix<int> floyd_output = floyd_warshall_int.floyd_warshall(graph);
	floyd_output.printMatrix();
}

// print verification that the algorithms produce the correct result for a known integer weighted graph while using a distance matrix to store the graph
void verify_algorithms_distance_matrix()
{
	// declare algorithm structures
	Dijkstra<int> dijkstra_int;
	BellmanFord<int> bellman_ford_int;
	JohnsonDijkstra<int> johnson_dijkstra_int;
	FloydWarshall<int> floyd_warshall_int;

	// load graphs
	StandardGraph<int> distance_list_input = loadKnownGraph();
	StandardGraph<int> distance_list_solution = loadKnownSolution();;

	// copy loaded graphs to distance matrix format
	DistanceMatrix<int> graph(distance_list_input.size);
	DistanceMatrix<int> solution(distance_list_solution.size);
	graph.fill_with_graph(distance_list_input);
	solution.fill_with_graph(distance_list_solution);

	std::cout << "\n VERIFYING ALGORITHMS ON GRAPHS IN DISTANCE MATRIX FORMAT\n";

	std::cout << " Graph in distance matrix form: \n";
	graph.printMatrix();

	std::cout << " Known solution of graph: \n";
	solution.printMatrix();

	std::cout << " APSP Dijkstra binary queue solution: \n";
	DistanceMatrix<int> dijkstra_output = dijkstra_int.APSP_dijkstra(graph);
	dijkstra_output.printMatrix();

	std::cout << " APSP Dijkstra fibonaci queue solution: \n";
	DistanceMatrix<int> dijkstra_output_fib = dijkstra_int.APSP_dijkstra_fibqueue(graph);
	dijkstra_output_fib.printMatrix();

	std::cout << " APSP Bellman-Ford solution: \n";
	DistanceMatrix<int> bellman_output = bellman_ford_int.APSP_bellman_ford(graph);
	bellman_output.printMatrix();

	std::cout << " Johnson-Dijkstra solution: \n";
	DistanceMatrix<int> johnson_output = johnson_dijkstra_int.johnson_dijkstra(graph);
	johnson_output.printMatrix();

	std::cout << " Floyd-Warshall solution: \n";
	DistanceMatrix<int> floyd_output = floyd_warshall_int.floyd_warshall(graph);
	floyd_output.printMatrix();
}

// print verification that the algorithms produce the correct result for a known integer weighted graph with negative weights
void verify_algorithms_negatives()
{
	// declare algorithm structures
	BellmanFord<int> bellman_ford_int;
	JohnsonDijkstra<int> johnson_dijkstra_int;
	FloydWarshall<int> floyd_warshall_int;

	// load graphs
	StandardGraph<int> graph = loadKnownGraphNegatives();
	StandardGraph<int> solution = loadKnownSolutionNegatives();

	std::cout << "\n VERIFYING ALGORITHMS ON GRAPHS WITH NEGATIVE WEIGHTS\n";

	std::cout << " Graph in distance list form: \n";
	graph.printGraph();
	std::cout << "\n";

	std::cout << " Known solution of graph: \n";
	DistanceMatrix<int> solution_dist(solution.size);
	solution_dist.fill_with_graph(solution);
	solution_dist.printMatrix();


	std::cout << " Floyd-Warshall solution: \n";
	DistanceMatrix<int> floyd_output = floyd_warshall_int.floyd_warshall(graph);
	floyd_output.printMatrix();

	std::cout << " APSP Bellman-Ford solution: \n";
	DistanceMatrix<int> bellman_output = bellman_ford_int.APSP_bellman_ford(graph);
	bellman_output.printMatrix();

	std::cout << " Johnson-Dijkstra solution: \n";
	DistanceMatrix<int> johnson_output = johnson_dijkstra_int.johnson_dijkstra(graph);
	johnson_output.printMatrix();
}

// print verification that the algorithms produce the correct result for known unweighted graphs - unused
void verify_algorithms_unweighted()
{
	// declare algorithm structures
	Dijkstra<int> dijkstra_int;
	BellmanFord<int> bellman_ford_int;
	JohnsonDijkstra<int> johnson_dijkstra_int;
	FloydWarshall<int> floyd_warshall_int;
	BreadthFirst breadth_first;

	// load graphs
	StandardGraph<int> graph = loadKnownGraphUnweighted();
	StandardGraph<int> solution = loadKnownSolutionUnweighted();;

	std::cout << " VERIFYING ALGORITHMS ON UNWEIGHTED GRAPHS\n";

	std::cout << " Graph in distance list form: \n";
	graph.printGraph();
	std::cout << "\n";

	std::cout << " Known solution of graph: \n";
	DistanceMatrix<int> solution_dist(solution.size);
	solution_dist.fill_with_graph(solution);
	solution_dist.printMatrix();


	std::cout << "APSP Dijkstra binary queue solution: \n";
	DistanceMatrix<int> dijkstra_output = dijkstra_int.APSP_dijkstra(graph);
	dijkstra_output.printMatrix();

	std::cout << "APSP Dijkstra fibonaci queue solution: \n";
	DistanceMatrix<int> dijkstra_output_fib = dijkstra_int.APSP_dijkstra_fibqueue(graph);
	dijkstra_output_fib.printMatrix();

	std::cout << "APSP Bellman-Ford solution: \n";
	DistanceMatrix<int> bellman_output = bellman_ford_int.APSP_bellman_ford(graph);
	bellman_output.printMatrix();

	std::cout << "Johnson-Dijkstra solution: \n";
	DistanceMatrix<int> johnson_output = johnson_dijkstra_int.johnson_dijkstra(graph);
	johnson_output.printMatrix();

	std::cout << "Floyd-Warshall solution: \n";
	DistanceMatrix<int> floyd_output = floyd_warshall_int.floyd_warshall(graph);
	floyd_output.printMatrix();

	std::cout << "Breadth-First solution: \n";
	DistanceMatrix<int> breadth_output = breadth_first.APSP_breadth_first(graph);
	breadth_output.printMatrix();
}

// print verification that the algorithms produce the correct result for randomly generated graphs with different parameters
void verify_graph_generation()
{
	std::cout << "\n VISUALISING RANDOMLY GENERATED GRAPHS\n";

	GraphGenerator graph_gen;

	std::cout << " Sparse int ER graph with uniform weights from 0 to 10\n";
	StandardGraph<int> graph = graph_gen.generateERGraphInt(8, 0.2, 0, 10);
	graph.printGraph();

	std::cout << " Dense int ER graph with uniform weights from 0 to 10\n";
	graph = graph_gen.generateERGraphInt(8, 0.95, 0, 10);
	graph.printGraph();

	std::cout << " Medium density double ER graph (uniform weights 0 to 1)\n";
	StandardGraph<double> graph_double = graph_gen.generateERGraphDouble(8, 0.5, 0, 1);
	graph_double.printGraph();

	std::cout << " Medium density double ER graph with gaussian weights (mean = 10, sd = 2)\n";
	graph_double = graph_gen.generateERGraphGaussian(8, 0.5, 10, 2);
	graph_double.printGraph();
}


// performance testing on graphs with int weights
void run_batch_int()
{
	GraphGenerator graph_gen;

	Dijkstra<int> dijkstra;
	JohnsonDijkstra<int> johnson_dijkstra;
	FloydWarshall<int> floyd_warshall;

	auto t_start = Clock::now();
	auto t_end = Clock::now();

	double total_times[10] = { 0 };

	int size;
	double density;
	double minweight;
	double maxweight;
	int batch_size;

	get_parameters(size, density, minweight, maxweight, batch_size);

	for (int i = 0; i < batch_size; i++)
	{
		std::cout << "\n Iteration " << i + 1 << " started:\n";

		StandardGraph<int> graph = graph_gen.generateERGraphInt(size, density, minweight, maxweight);

		// Binary Dijkstra
		std::cout << "\n  - Binary Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra(graph);
		t_end = Clock::now();

		total_times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Fibonacci Dijkstra
		std::cout << "\n  - Fibonacci Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra_fibqueue(graph);
		t_end = Clock::now();

		total_times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Johnson Dijkstra
		std::cout << "\n  - Johnson Dijkstra started\n";
		t_start = Clock::now();
		johnson_dijkstra.johnson_dijkstra(graph);
		t_end = Clock::now();

		total_times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// FloydWarshall
		std::cout << "\n  - Floyd Warshall started\n";
		t_start = Clock::now();
		floyd_warshall.floyd_warshall(graph);
		t_end = Clock::now();

		total_times[3] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		std::cout << "\n  - Iteration " << i + 1 << " ended\n";
	}

	std::cout << "\n Results:\n";

	std::cout << "  - Time to solve Binary Heap Dijkstra: " << int(total_times[0] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Fibonacci Heap: " << int(total_times[1] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Johnson Dijkstra: " << int(total_times[2] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Floyd Warshall: " << int(total_times[3] / batch_size) << " seconds \n";

	std::cout << "\n";
}

// performance testing on graphs with float weights
void run_batch_float()
{
	GraphGenerator graph_gen;

	Dijkstra<float> dijkstra;
	JohnsonDijkstra<float> johnson_dijkstra;
	FloydWarshall<float> floyd_warshall;

	auto t_start = Clock::now();
	auto t_end = Clock::now();

	double total_times[10] = { 0 };

	int size;
	double density;
	double minweight;
	double maxweight;
	int batch_size;

	get_parameters(size, density, minweight, maxweight, batch_size);

	for (int i = 0; i < batch_size; i++)
	{
		std::cout << "\n Iteration " << i + 1 << " started:\n";

		StandardGraph<float> graph = graph_gen.generateERGraphFloat(size, density, minweight, maxweight);

		// Binary Dijkstra
		std::cout << "\n  - Binary Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra(graph);
		t_end = Clock::now();

		total_times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Fibonacci Dijkstra
		std::cout << "\n  - Fibonacci Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra_fibqueue(graph);
		t_end = Clock::now();

		total_times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Johnson Dijkstra
		std::cout << "\n  - Johnson Dijkstra started\n";
		t_start = Clock::now();
		johnson_dijkstra.johnson_dijkstra(graph);
		t_end = Clock::now();

		total_times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// FloydWarshall
		std::cout << "\n  - Floyd Warshall started\n";
		t_start = Clock::now();
		floyd_warshall.floyd_warshall(graph);
		t_end = Clock::now();

		total_times[3] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		std::cout << "\n  - Iteration " << i + 1 << " ended\n";
	}

	std::cout << "\n Results:\n";

	std::cout << "  - Time to solve Binary Heap Dijkstra: " << int(total_times[0] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Fibonacci Heap: " << int(total_times[1] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Johnson Dijkstra: " << int(total_times[2] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Floyd Warshall: " << int(total_times[3] / batch_size) << " seconds \n";

	std::cout << "\n";
}

// performance testing on graphs with double weights
void run_batch_double()
{
	GraphGenerator graph_gen;

	Dijkstra<double> dijkstra;
	JohnsonDijkstra<double> johnson_dijkstra;
	FloydWarshall<double> floyd_warshall;

	auto t_start = Clock::now();
	auto t_end = Clock::now();

	double total_times[10] = { 0 };

	int size;
	double density;
	double minweight;
	double maxweight;
	int batch_size;

	get_parameters(size, density, minweight, maxweight, batch_size);

	for (int i = 0; i < batch_size; i++)
	{
		std::cout << "\n Iteration " << i + 1 << " started:\n";

		StandardGraph<double> graph = graph_gen.generateERGraphDouble(size, density, minweight, maxweight);

		// Binary Dijkstra
		std::cout << "\n  - Binary Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra(graph);
		t_end = Clock::now();

		total_times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Fibonacci Dijkstra
		std::cout << "\n  - Fibonacci Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra_fibqueue(graph);
		t_end = Clock::now();

		total_times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Johnson Dijkstra
		std::cout << "\n  - Johnson Dijkstra started\n";
		t_start = Clock::now();
		johnson_dijkstra.johnson_dijkstra(graph);
		t_end = Clock::now();

		total_times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// FloydWarshall
		std::cout << "\n  - Floyd Warshall started\n";
		t_start = Clock::now();
		floyd_warshall.floyd_warshall(graph);
		t_end = Clock::now();

		total_times[3] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		std::cout << "\n  - Iteration " << i + 1 << " ended\n";
	}

	std::cout << "\n Results:\n";

	std::cout << "  - Time to solve Binary Heap Dijkstra: " << int(total_times[0] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Fibonacci Heap: " << int(total_times[1] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Johnson Dijkstra: " << int(total_times[2] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Floyd Warshall: " << int(total_times[3] / batch_size) << " seconds \n";

	std::cout << "\n";
}

// performance testing on graphs with double weights - inluding the Bellman-Ford algorithm
void run_batch_double_with_bellman()
{
	GraphGenerator graph_gen;

	Dijkstra<double> dijkstra;
	BellmanFord<double> bellman_ford;
	JohnsonDijkstra<double> johnson_dijkstra;
	FloydWarshall<double> floyd_warshall;

	auto t_start = Clock::now();
	auto t_end = Clock::now();

	double total_times[10] = { 0 };

	int size;
	double density;
	double minweight;
	double maxweight;
	int batch_size;

	get_parameters(size, density, minweight, maxweight, batch_size);

	for (int i = 0; i < batch_size; i++)
	{
		std::cout << "\n Iteration " << i + 1 << " started:\n";

		StandardGraph<double> graph = graph_gen.generateERGraphDouble(size, density, minweight, maxweight);

		// Binary Dijkstra
		std::cout << "\n  - Binary Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra(graph);
		t_end = Clock::now();

		total_times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Fibonacci Dijkstra
		std::cout << "\n  - Fibonacci Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra_fibqueue(graph);
		t_end = Clock::now();

		total_times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Bellman Ford
		std::cout << "\n  - Bellman Ford started\n";
		t_start = Clock::now();
		bellman_ford.APSP_bellman_ford(graph);
		t_end = Clock::now();

		total_times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Johnson Dijkstra
		std::cout << "\n  - Johnson Dijkstra started\n";
		t_start = Clock::now();
		johnson_dijkstra.johnson_dijkstra(graph);
		t_end = Clock::now();

		total_times[3] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// FloydWarshall
		std::cout << "\n  - Floyd Warshall started\n";
		t_start = Clock::now();
		floyd_warshall.floyd_warshall(graph);
		t_end = Clock::now();

		total_times[4] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		std::cout << "\n  - Iteration " << i + 1 << " ended\n";
	}

	std::cout << "\n Results:\n";

	std::cout << "  - Time to solve Binary Heap Dijkstra: " << int(total_times[0] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Fibonacci Heap: " << int(total_times[1] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Bellman Ford: " << int(total_times[2] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Johnson Dijkstra: " << int(total_times[3] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Floyd Warshall: " << int(total_times[4] / batch_size) << " seconds \n";

	std::cout << "\n";
}

// performance testing on graphs with gaussian distributed double weights
void run_batch_gaussian()
{
	GraphGenerator graph_gen;

	Dijkstra<double> dijkstra;
	JohnsonDijkstra<double> johnson_dijkstra;
	FloydWarshall<double> floyd_warshall;

	auto t_start = Clock::now();
	auto t_end = Clock::now();

	double total_times[10] = { 0 };

	int size;
	double density;
	double mean;
	double sd;
	int batch_size;

	get_parameters_gaussian(size, density, mean, sd, batch_size);

	for (int i = 0; i < batch_size; i++)
	{
		std::cout << "\n Iteration " << i + 1 << " started:\n";

		StandardGraph<double> graph = graph_gen.generateERGraphGaussian(size, density, mean, sd);

		// Binary Dijkstra
		std::cout << "\n  - Binary Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra(graph);
		t_end = Clock::now();

		total_times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Fibonacci Dijkstra
		std::cout << "\n  - Fibonacci Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra_fibqueue(graph);
		t_end = Clock::now();

		total_times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Johnson Dijkstra
		std::cout << "\n  - Johnson Dijkstra started\n";
		t_start = Clock::now();
		johnson_dijkstra.johnson_dijkstra(graph);
		t_end = Clock::now();

		total_times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// FloydWarshall
		std::cout << "\n  - Floyd Warshall started\n";
		t_start = Clock::now();
		floyd_warshall.floyd_warshall(graph);
		t_end = Clock::now();

		total_times[3] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		std::cout << "\n  - Iteration " << i + 1 << " ended\n";
	}

	std::cout << "\n Results:\n";

	std::cout << "  - Time to solve Binary Heap Dijkstra: " << int(total_times[0] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Fibonacci Heap: " << int(total_times[1] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Johnson Dijkstra: " << int(total_times[2] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Floyd Warshall: " << int(total_times[3] / batch_size) << " seconds \n";

	std::cout << "\n";
}

// performance testing on constant 1 weighted graphs, analagous to unweighted, to compare breadth-first search to other algorithms
void run_batch_unweighted()
{
	GraphGenerator graph_gen;

	Dijkstra<int> dijkstra;
	JohnsonDijkstra<int> johnson_dijkstra;
	FloydWarshall<int> floyd_warshall;
	BreadthFirst breadth_first;

	auto t_start = Clock::now();
	auto t_end = Clock::now();

	double total_times[10] = { 0 };

	int size;
	double density;
	int batch_size;

	get_parameters_no_weights(size, density, batch_size);

	for (int i = 0; i < batch_size; i++)
	{
		std::cout << "\n Iteration " << i + 1 << " started:\n";

		StandardGraph<int> graph = graph_gen.generateERGraphInt(size, density, 1, 1);

		// Binary Dijkstra
		std::cout << "\n  - Binary Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra(graph);
		t_end = Clock::now();

		total_times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Fibonacci Dijkstra
		std::cout << "\n  - Fibonacci Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra_fibqueue(graph);
		t_end = Clock::now();

		total_times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Johnson Dijkstra
		std::cout << "\n  - Johnson Dijkstra started\n";
		t_start = Clock::now();
		johnson_dijkstra.johnson_dijkstra(graph);
		t_end = Clock::now();

		total_times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// FloydWarshall
		std::cout << "\n  - Floyd Warshall started\n";
		t_start = Clock::now();
		floyd_warshall.floyd_warshall(graph);
		t_end = Clock::now();

		total_times[3] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// BreadthFirst
		std::cout << "\n  - Breadth First started\n";
		t_start = Clock::now();
		breadth_first.APSP_breadth_first(graph);
		t_end = Clock::now();

		total_times[4] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		std::cout << "\n  - Iteration " << i + 1 << " ended\n";
	}

	std::cout << "\n Results:\n";

	std::cout << "  - Time to solve Binary Heap Dijkstra: " << int(total_times[0] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Fibonacci Heap: " << int(total_times[1] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Johnson Dijkstra: " << int(total_times[2] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Floyd Warshall: " << int(total_times[3] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Breadth First: " << int(total_times[4] / batch_size) << " seconds \n";

	std::cout << "\n";
}

// performance testing on graphs with double weights using ditsance matrices
void run_batch_distance_matrix()
{
	GraphGenerator graph_gen;

	Dijkstra<double> dijkstra;
	JohnsonDijkstra<double> johnson_dijkstra;
	FloydWarshall<double> floyd_warshall;

	auto t_start = Clock::now();
	auto t_end = Clock::now();

	double total_times[10] = { 0 };

	int size;
	double density;
	double minweight;
	double maxweight;
	int batch_size;

	get_parameters(size, density, minweight, maxweight, batch_size);

	for (int i = 0; i < batch_size; i++)
	{
		std::cout << "\n Iteration " << i + 1 << " started:\n";

		StandardGraph<double> distance_list = graph_gen.generateERGraphDouble(size, density, minweight, maxweight);
		DistanceMatrix<double> graph(distance_list.size);
		graph.fill_with_graph(distance_list);

		// Binary Dijkstra
		std::cout << "\n  - Binary Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra(graph);
		t_end = Clock::now();

		total_times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Fibonacci Dijkstra
		std::cout << "\n  - Fibonacci Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra_fibqueue(graph);
		t_end = Clock::now();

		total_times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Johnson Dijkstra
		std::cout << "\n  - Johnson Dijkstra started\n";
		t_start = Clock::now();
		johnson_dijkstra.johnson_dijkstra(graph);
		t_end = Clock::now();

		total_times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// FloydWarshall
		std::cout << "\n  - Floyd Warshall started\n";
		t_start = Clock::now();
		floyd_warshall.floyd_warshall(graph);
		t_end = Clock::now();

		total_times[3] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		std::cout << "\n  - Iteration " << i + 1 << " ended\n";
	}

	std::cout << "\n Results:\n";

	std::cout << "  - Time to solve Binary Heap Dijkstra: " << int(total_times[0] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Fibonacci Heap: " << int(total_times[1] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Johnson Dijkstra: " << int(total_times[2] / batch_size) << " seconds \n";

	std::cout << "  - Time to solve Floyd Warshall: " << int(total_times[3] / batch_size) << " seconds \n";

	std::cout << "\n";
}

// performance testing on graphs of real data loaded from files
void run_on_real_graph()
{
	Dijkstra<double> dijkstra;
	JohnsonDijkstra<double> johnson_dijkstra;
	FloydWarshall<double> floyd_warshall;

	auto t_start = Clock::now();
	auto t_end = Clock::now();

	double total_times[10] = { 0 };

	string filename;
	bool lower;
	cout << "\n Enter a file path: ";
	cin >> filename;
	StandardGraph<double> graph = loadRealGraph(filename);

	// Binary Dijkstra
	std::cout << "\n  - Binary Heap Dijkstra started\n";
	t_start = Clock::now();
	dijkstra.APSP_dijkstra(graph);
	t_end = Clock::now();

	total_times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

	// Fibonacci Dijkstra
	std::cout << "\n  - Fibonacci Heap Dijkstra started\n";
	t_start = Clock::now();
	dijkstra.APSP_dijkstra_fibqueue(graph);
	t_end = Clock::now();

	total_times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();


	// Johnson Dijkstra
	std::cout << "\n  - Johnson Dijkstra started\n";
	t_start = Clock::now();
	johnson_dijkstra.johnson_dijkstra(graph);
	t_end = Clock::now();

	total_times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

	// FloydWarshall
	std::cout << "\n  - Floyd Warshall started\n";
	t_start = Clock::now();
	floyd_warshall.floyd_warshall(graph);
	t_end = Clock::now();

	total_times[3] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();


	std::cout << "\n Results:\n";

	std::cout << "  - Time to solve Binary Heap Dijkstra: " << int(total_times[0]) << " seconds \n";

	std::cout << "  - Time to solve Fibonacci Heap: " << int(total_times[1]) << " seconds \n";

	std::cout << "  - Time to solve Johnson Dijkstra: " << int(total_times[2]) << " seconds \n";

	std::cout << "  - Time to solve Floyd Warshall: " << int(total_times[3]) << " seconds \n";

	std::cout << "\n";
}

// performance testing on graphs of real data loaded from files - use when the graph is too big for the memory
void run_on_real_graph_split()
{
	Dijkstra<double> dijkstra;
	JohnsonDijkstra<double> johnson_dijkstra;
	FloydWarshall<double> floyd_warshall;

	auto t_start = Clock::now();
	auto t_end = Clock::now();

	double total_times[10] = { 0 };

	string filename;
	bool lower;
	cout << "\n Enter a file path: ";
	cin >> filename;
	cout << "\n Upper or lower half? (0 or 1): ";
	cin >> lower;
	StandardGraph<double> graph = loadRealGraphSplit(filename, 16726, lower);

	// Binary Dijkstra
	std::cout << "\n  - Binary Heap Dijkstra started\n";
	t_start = Clock::now();
	dijkstra.APSP_dijkstra(graph);
	t_end = Clock::now();

	total_times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

	// Fibonacci Dijkstra
	std::cout << "\n  - Fibonacci Heap Dijkstra started\n";
	t_start = Clock::now();
	dijkstra.APSP_dijkstra_fibqueue(graph);
	t_end = Clock::now();

	total_times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();


	// Johnson Dijkstra
	std::cout << "\n  - Johnson Dijkstra started\n";
	t_start = Clock::now();
	johnson_dijkstra.johnson_dijkstra(graph);
	t_end = Clock::now();

	total_times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

	// FloydWarshall
	std::cout << "\n  - Floyd Warshall started\n";
	t_start = Clock::now();
	floyd_warshall.floyd_warshall(graph);
	t_end = Clock::now();

	total_times[3] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();


	std::cout << "\n Results:\n";

	std::cout << "  - Time to solve Binary Heap Dijkstra: " << int(total_times[0]) << " seconds \n";

	std::cout << "  - Time to solve Fibonacci Heap: " << int(total_times[1]) << " seconds \n";

	std::cout << "  - Time to solve Johnson Dijkstra: " << int(total_times[2]) << " seconds \n";

	std::cout << "  - Time to solve Floyd Warshall: " << int(total_times[3]) << " seconds \n";

	std::cout << "\n";
}

// choose a batch operation
void run_batch() {

	cout << "\n CHOOSE A BATCH OPERATION FROM THE LIST \n";

	cout << " 1) run_batch_double()\n";
	cout << " 2) run_batch_float()\n";
	cout << " 3) run_batch_int()\n";
	cout << " 4) run_batch_double_with_bellman()\n";
	cout << " 5) run_batch_distance_matrix()\n";
	cout << " 6) run_batch_gaussian()\n";
	cout << " 7) run_batch_unweighted()\n";
	cout << " 8) run_on_real_graph()\n";

	cout << "\n Choice: ";

	int choice;
	while (true)
	{
		cin >> choice;

		if (choice == 1) {
			run_batch_double();
			break;
		} 
		else if (choice == 2) {
			run_batch_float();
			break;
		}
		else if (choice == 3) {
			run_batch_int();
			break;
		}
		else if (choice == 4) {
			run_batch_double_with_bellman();
			break;
		}
		else if (choice == 5) {
			run_batch_distance_matrix();
			break;
		}
		else if (choice == 6) {
			run_batch_gaussian();
			break;
		}
		else if (choice == 7) {
			run_batch_unweighted();
			break;
		}
		else if (choice == 8) {
			run_on_real_graph();
			break;
		}
		else {
			cout << "\n Please enter a valid number: ";
		}
	}
}


int main()
{
	cout << "\n --- for information on how to use this program see the \"Implementation\" section in the writeup --- \n";

	cout << "\n CHOOSE A FUNCTION \n";

	cout << " 1) verify implementation of APSP algorithms \n";
	cout << " 2) verify random graph generation \n";
	cout << " 3) do performance testing of APSP algorithms \n";

	cout << "\n Choice: ";

	int choice;
	while (true)
	{
		cin >> choice;

		if (choice == 1) {
			verify_algorithms();
			verify_algorithms_distance_matrix();
			verify_algorithms_negatives();
			verify_algorithms_unweighted();
			break;
		}
		else if (choice == 2) {
			verify_graph_generation();
			break;
		}
		else if (choice == 3) {
			run_batch();
			break;
		}
		else {
			cout << "\n Please enter a valid number: ";
		}
	}

	int any;
	std::cout << "Close with any key: ";
	cin >> any;

	return 0;
}
