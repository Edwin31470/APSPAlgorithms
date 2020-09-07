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

// structures to contain edge information
template <class T>
struct Edge {
	int source;
	pair<int, T> dest; // first is destination node, second is edge weight
};

// standard graph object - weighted directed graph with optional directed flag
template <class T>
class StandardGraph
{
public:
	// construct a vector of vectors to represent an adjacency list
	vector<vector<pair<int, T>>> list;
	int size;
	int num_edges;
	bool directed;

	// construct the graph from edges - directed unless otherwise specified
	StandardGraph(vector<Edge<T>> const &edges, int N, bool dir = true)
	{
		size = N;
		num_edges = edges.size();
		directed = dir;

		// resize the vector to N elements of type vector<int>
		list.resize(N);

		// add edges to the directed graph
		for (auto &edge : edges)
		{
			// insert at the end
			list[edge.source].push_back(edge.dest);

			// if undirected, add reverse edge
			if (!directed)
				list[edge.dest.first].push_back(make_pair(edge.source, edge.dest.second));
		}
	}

	// print the contents of this graph to console
	void printGraph()
	{
		for (int i = 0; i < size; i++)
		{
			// print current vertex number
			std::cout << i << " --> ";

			// print all neighboring vertices of vertex i
			for (const pair<int, T> &dest : list[i])
				std::cout << dest.first << "/" << dest.second << " ";
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
		vector<pair<int, T>> node;
		list.push_back(node);

		size = size + 1;
	}

	// removes all edges associated with a vertex. does not actually remove the vertex label unless the vertex is the last vertex
	void removeVertex(int label)
	{
		for (int i = 0; i < size; i++) { // iterate over every node
			vector<pair<int, T>> node = list[i]; // list[i] must be copied to avoid changing container size as edges are deleted
			for (pair<int, T> edge : node) { // and every edge leaving that node

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

	// add an edge from source to dest with a weight - directed edge unless graph is undirected
	void addEdge(int source, int dest, T weight)
	{
		list[source].push_back(make_pair(dest, weight));

		if (!directed)
			list[dest].push_back(make_pair(source, weight));
	}

	// remove edge from source to dest - if directed remove edge from dest to source too
	void removeEdge(int source, int dest)
	{
		list[source].erase(remove_if(
			list[source].begin(),
			list[source].end(),
			[dest](pair<int, T> edge) {return edge.first == dest; }
		));

		if (!directed) {
			list[dest].erase(remove_if(
				list[dest].begin(),
				list[dest].end(),
				[source](pair<int, T> edge) {return edge.first == source; }
			));
		}
		
	}

	// reweight edge from source to dest - if directed reweight edge from dest to source too
	void reweightEdge(int source, int dest, T weight)
	{
		for (pair<int, T> &edge : list[source]) {
			if (edge.first == dest) {
				edge.second = weight;
			}
		}

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
	vector<T> matrix;
	int size;
	
	DistanceMatrix(int N)
	{
		size = N;

		matrix.resize(N*N);

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				// set trace to 0
				if (i == j) {
					setEdge(i, j, 0);
				}
				// set all other values to positive maximum
				else {
					setEdge(i, j, std::numeric_limits<T>::max());
				}
			}
		}
	}

	DistanceMatrix(vector<Edge<T>> const &edges, int N)
	{
		size = N;

		matrix.resize(N*N);

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				// set trace to 0
				if (i == j) {
					setEdge(i, j, 0);
				}
				// set all other values to positive maximum
				else {
					setEdge(i, j, std::numeric_limits<T>::max());
				}
			}
		}

		// add edges to the matrix
		for (const Edge<T> &edge : edges)
		{
			setEdge(edge.source, edge.dest.first, edge.dest.second);
		}
	}

	void setEdge(int source, int dest, T weight) {
		matrix[(source * size) + dest] = weight;
	}

	T getEdge(int source, int dest) {
		return matrix[(source * size) + dest];
	}

	// put adjacency list graph into distance matrix
	void fill_with_graph(const StandardGraph<T> &graph)
	{
		for (int i = 0; i < size; i++) {
			for (const pair<int, T> &dest : graph.list[i]) {
				setEdge(i, dest.first, dest.second);
			}
		}
	}

	void printMatrix()
	{
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				if (getEdge(i, j) == std::numeric_limits<T>::max()) {
					std::cout << "-  ";
				}
				else {
					std::cout << getEdge(i, j) << " ";
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



// UTILITY FUNCTIONS

// returns an int graph of specific values
StandardGraph<int> loadKnownGraph()
{
	// Starting graph:
	// 0 7 3 - - - - -
	// 3 0 7 - - - - -
	// 9 - 0 - - 4 4 -
	// - - - 0 - - - -
	// - - - - 0 - - 8
	// - 6 - 2 - 0 5 -
	// - - - - - - 0 7
	// - - - - - - - 0

	vector<Edge<int>> edges =
	{
		{ 0, make_pair(1, 7) }, { 0, make_pair(2, 3) },
		{ 1, make_pair(0, 3) }, { 1, make_pair(2, 7) },
		{ 2, make_pair(0, 9) }, { 2, make_pair(5, 4) }, { 2, make_pair(6, 4) },

		{ 4, make_pair(7, 8) },
		{ 5, make_pair(1, 6) }, { 5, make_pair(3, 2) }, { 5, make_pair(6, 5) },
		{ 6, make_pair(7, 7) }
	};

	// Number of nodes in the graph
	int N = 8;

	// construct graph
	StandardGraph<int> graph(edges, N);

	return graph;
}

// returns a double graph of specific
StandardGraph<double> loadKnownDoubleGraph()
{
	// Starting graph:
	// 0 7 3 - - - - -
	// 3 0 7 - - - - -
	// 9 - 0 - - 4 4 -
	// - - - 0 - - - -
	// - - - - 0 - - 8
	// - 6 - 2 - 0 5 -
	// - - - - - - 0 7
	// - - - - - - - 0

	vector<Edge<double>> edges =
	{
		{ 0, make_pair(1, 7.1) }, { 0, make_pair(2, 3.1) },
		{ 1, make_pair(0, 3.1) }, { 1, make_pair(2, 7.1) },
		{ 2, make_pair(0, 9.1) }, { 2, make_pair(5, 4.1) }, { 2, make_pair(6, 4.1) },

		{ 4, make_pair(7, 8.1) },
		{ 5, make_pair(1, 6.1) }, { 5, make_pair(3, 2.1) }, { 5, make_pair(6, 5.1) },
		{ 6, make_pair(7, 7.1) }
	};

	// Number of nodes in the graph
	int N = 8;

	// construct graph
	StandardGraph<double> graph(edges, N);

	return graph;
}

// returns an int graph of specific values including negatives
StandardGraph<int> loadKnownGraphNegatives()
{
	// Starting graph:
	// 0  -  -  2
	// 6  0  3  -
	// 4  -  0  5 
	// -  -7 -3 0 

	vector<Edge<int>> edges =
	{
		{ 0, make_pair(3, 2) },
		{ 1, make_pair(0, 6) }, { 1, make_pair(2, 3) },
		{ 2, make_pair(0, 4) }, { 2, make_pair(3, 5) },
		{ 3, make_pair(1, -7) }, { 3, make_pair(2, -3) }
	};

	// Number of nodes in the graph
	int N = 4;

	// construct graph
	StandardGraph<int> graph(edges, N);

	return graph;
}

// returns an undirected int graph of specific values
StandardGraph<int> loadKnownUnidrectedGraph()
{
	/*
	// Starting graph:
	// 0 - - - - - - -
	// 4 0 - - - - - -
	// - 8 0 - - - - -
	// - - 7 0 - - - -
	// - - 2 - 0 - - -
	// 8 9 - - 7 0 - -
	// - - - - 5 1 0 -
	// - - 4 - - - 2 0

	vector<Edge> edges =
	{
		{ 1, make_pair(0, 4) },
		{ 2, make_pair(1, 8) },
		{ 3, make_pair(2, 7) },
		{ 4, make_pair(2, 2) },
		{ 5, make_pair(0, 8) }, { 5, make_pair(1, 9) }, { 5, make_pair(4, 7) },
		{ 6, make_pair(4, 5) }, { 6, make_pair(5, 1) },
		{ 7, make_pair(2, 4) }, { 7, make_pair(6, 2) }
	};

	// Number of nodes in the graph
	int N = 8;

	*/

	vector<Edge<int>> edges =
	{
		{ 0, make_pair(1, 78) }, { 0, make_pair(2, 21) }, { 0, make_pair(10, 5001) },
		{ 1, make_pair(2, 3) }, { 1, make_pair(6, 89) }, { 1, make_pair(7, 829) }, { 1, make_pair(9, 5) }, { 1, make_pair(10, 290) },
		
		{ 3, make_pair(4, 157) }, { 3, make_pair(6, 41) },
		{ 4, make_pair(5, 1688) }, { 4, make_pair(6, 7) },
		{ 5, make_pair(6, 3622) },
		{ 6, make_pair(7, 36025) },
		{ 7, make_pair(8, 11) },
		{ 8, make_pair(9, 19) }
	};

	// Number of nodes in the graph
	int N = 11;

	// construct graph
	StandardGraph<int> graph(edges, N, false);

	return graph;
}

// returns the known APSP solution of graph from loadKnownGraph()
StandardGraph<int> loadKnownSolution()
{

	// Solution:
	// 0  7  3  9  -  7  7  14
	// 3  0  6  12 -  10 10 17
	// 9  10 0  6  -  4  4  11
	// -  -  -  0  -  -  -  -
	// -  -  -  -  0  -  -  8
	// 9  6  12 2  -  0  5  12
	// -  -  -  -  -  -  0  7
	// -  -  -  -  -  -  -  0

	vector<Edge<int>> edges =
	{
		{ 0, make_pair(1, 7) }, { 0, make_pair(2, 3) }, { 0, make_pair(3, 9) },
		{ 0, make_pair(5, 7) }, { 0, make_pair(6, 7) }, { 0, make_pair(7, 14) },

		{ 1, make_pair(0, 3) }, { 1, make_pair(2, 6) }, { 1, make_pair(3, 12) },
		{ 1, make_pair(5, 10) }, { 1, make_pair(6, 10) }, { 1, make_pair(7, 17) },

		{ 2, make_pair(0, 9) }, { 2, make_pair(1, 10) }, { 2, make_pair(3, 6) },
		{ 2, make_pair(5, 4) }, { 2, make_pair(6, 4) }, { 2, make_pair(7, 11) },

		{ 4, make_pair(7, 8) },

		{ 5, make_pair(0, 9) }, { 5, make_pair(1, 6) }, { 5, make_pair(2, 12) },
		{ 5, make_pair(3, 2) }, { 5, make_pair(6, 5) }, { 5, make_pair(7, 12) },

		{ 6, make_pair(7, 7) },
	};

	// Number of nodes in the graph
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
		{ 0, make_pair(1, -5) }, { 0, make_pair(2, -2) }, { 0, make_pair(3, 2) },
		{ 1, make_pair(0, 6) }, { 1, make_pair(2, 3) }, { 1, make_pair(3, 8) },
		{ 2, make_pair(0, 4) }, { 2, make_pair(1, -2) }, { 2, make_pair(3, 5) },
		{ 3, make_pair(0, -1) }, { 3, make_pair(1, -7) }, { 3, make_pair(2, -4) },
	};

	// Number of nodes in the graph
	int N = 4;

	// construct graph
	StandardGraph<int> graph(edges, N);

	return graph;
}

// builds and returns a graph contained in a file - takes a file path - for file format see writeup
StandardGraph<double> loadRealGraph(string filename)
{
	vector<Edge<double>> edges;
	// Number of nodes in the graph
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

			// find the max vertex label - dest is needed as it is possible for no edges to leave the max vertex
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
		std::cout << "Unable to open file";

	// construct graph
	StandardGraph<double> graph(edges, N);

	std::cout << "\nGraph loaded:";
	std::cout << "\n - Number of vertices: " << N;
	std::cout << "\n - Number of edges: " << edges.size();
	std::cout << "\n - Minimum edge weight: " << minWeight;
	std::cout << "\n - Maximum edge weight: " << maxWeight << "\n";

	return graph;
}

// object to build random graphs
class GraphGenerator {
private:
	std::default_random_engine r_eng;

public:
	int getRandomInt(int min, int max) {
		std::uniform_int_distribution<int> uniform_dist(min, max);
		return uniform_dist(r_eng);
	}

	double getRandomDouble(double min, double max) {
		std::uniform_real_distribution<double> uniform_dist(min, max);
		return uniform_dist(r_eng);
	}

	double getGaussianDouble(double mean, double stand_dev) {
		std::normal_distribution<double> distribution(mean, stand_dev);
		return distribution(r_eng);
	}


	// size must be an integer greater than 0 and density must be a probablity value [0, 1]
	StandardGraph<int> generateERGraphInt(int size, double density, int minWeight, int maxWeight)
	{
		vector<Edge<int>> edges;

		// add every edge with probability of density
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				// if edge is not to itself and independent probability is less than p, add edge
				if ((i != j) && (getRandomDouble(0, 1) < density)) {
					int weight = getRandomInt(minWeight, maxWeight);
					//std::cout << "Rand weight: " << weight << "\n";
					edges.push_back({ i, make_pair(j, weight) });
				}
			}
		}

		StandardGraph<int> graph(edges, size);

		std::cout << " - Graph created \n";

		return graph;
	}

	// size must be an integer greater than 0 and density must be a probablity value [0, 1]
	StandardGraph<float> generateERGraphFloat(int size, double density, double minWeight, double maxWeight)
	{
		vector<Edge<float>> edges;

		// add every edge with probability of density
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

		std::cout << " - Graph created \n";

		return graph;
	}

	// size must be an integer greater than 0 and density must be a probablity value [0, 1]
	StandardGraph<double> generateERGraphDouble(int size, double density, double minWeight, double maxWeight)
	{
		vector<Edge<double>> edges;

		// add every edge with probability of density
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				// if edge is not to itself and independent probability is less than p, add edge
				if ((i != j) && (getRandomDouble(0, 1) < density)) {
					double weight = getRandomDouble(minWeight, maxWeight);
					//std::cout << "Rand weight: " << weight << "\n";
					edges.push_back({ i, make_pair(j, weight) });
				}
			}
		}

		StandardGraph<double> graph(edges, size);

		std::cout << " - Graph created \n";

		return graph;
	}

	// size must be an integer greater than 0 and density must be a probablity value [0, 1]
	StandardGraph<double> generateERGraphGaussian(int size, double density, double mean, double stand_dev)
	{
		vector<Edge<double>> edges;

		// add every edge with probability of density
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				// if edge is not to itself and independent probability is less than p, add edge
				if ((i != j) && (getRandomDouble(0, 1) < density)) {
					double weight = getGaussianDouble(mean, stand_dev);
					if (weight < 0) {
						std::cout << " - warning: negative edge detected\n";
						weight = 0;
					}
					edges.push_back({ i, make_pair(j, weight) });
				}
			}
		}

		StandardGraph<double> graph(edges, size);

		std::cout << " - Graph created \n";

		return graph;
	}

	// size must be an integer greater than 0 and density must be a probablity value [0, 1]
	void generateBAGraphInt(int size, int m, int minWeight, int maxWeight)
	{

	}
};




// ALGORITHMS

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
	priority_queue<pair<int, Dtype>, vector<pair<int, Dtype>>, prioritize<Dtype>> queue; // priority queue to store vertex weight pairs 
	FibQueue<pair<int,Dtype>> fib_queue; // fibonacci queue to store vertex weight pairs

public:
	// takes distance vector and fills vector with distances from source
	void SSSP_dijkstra(const StandardGraph<Dtype> &graph, int source, vector<Dtype> &distances) //Algorithm for SSSP  
	{
		int N = graph.size;

		vector<bool> visited(N, false); // determines whether the node has been visited or not 

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

				// if the adjacent vertex is not visited, go to the next adjacent vertex
				if (!visited[adjacent_vertex]) {

					Dtype adjacent_weight = graph.list[current_vertex][i].second;

					// if the current vertex distance + distance from the current vertex to the adjacent vertex
					// is shorter than the shortest known distance to the adjacent vertex, replace that distance
					if (current_distance + adjacent_weight < distances[adjacent_vertex])
					{
						distances[adjacent_vertex] = adjacent_weight + current_distance; // set the new distance

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

	// takes distance vector and fills vector with distances from source
	void SSSP_dijkstra_fibqueue(const StandardGraph<Dtype> &graph, int source, vector<Dtype> &distances) //Algorithm for SSSP  
	{
		int N = graph.size;

		vector<bool> visited(N, false); // determines whether the node has been visited or not 

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
					// is shorter than the shortest known distance to the adjacent vertex, replace that distance
					if (current_distance + adjacent_weight < distances[adjacent_vertex])
					{
						distances[adjacent_vertex] = adjacent_weight + current_distance; // set the new distance

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

	// returns APSP distance matrix of given graph
	DistanceMatrix<Dtype> APSP_dijkstra(const StandardGraph<Dtype> &graph)
	{
		int N = graph.size;

		DistanceMatrix<Dtype> output(N);

		for (int i = 0; i < N; i++) {

			vector<Dtype> distances(N, std::numeric_limits<Dtype>::max()); //Stores shortest distance, initial distances of infinity
			SSSP_dijkstra(graph, i, distances);

			for (int j = 0; j < N; j++) {
				if (distances[j] != std::numeric_limits<int>::max())
					output.setEdge(i, j, distances[j]);
			}
		}

		return output;
	}

	DistanceMatrix<Dtype> APSP_dijkstra_fibqueue(const StandardGraph<Dtype> &graph)
	{
		int N = graph.size;

		DistanceMatrix<Dtype> output(N);

		for (int i = 0; i < N; i++) {

			vector<Dtype> distances(N, std::numeric_limits<Dtype>::max()); //Stores shortest distance, initial distances of infinity
			SSSP_dijkstra_fibqueue(graph, i, distances);

			for (int j = 0; j < N; j++) {
				if (distances[j] != std::numeric_limits<int>::max())
					output.setEdge(i, j, distances[j]);
			}
		}

		return output;
	}
};

// object to perform BreadFirst search - needs to be updated
struct BreadthFirst {
private:
	//priority_queue<int, vector<int>, std::greater<int>> queue;
	queue<int> queue;

public:

	void SSSP_breadth_first(const StandardGraph<int> &graph, int source, vector<int> &distances) //Algorithm for SSSP  
	{
		int N = graph.size;

		vector<bool> visited(N, false); // determines whether the node has been visited or not 

		visited[source] = true;
		distances[source] = 0;
		queue.push(source);

		while (!queue.empty())
		{
			int current_vertex = queue.front();
			queue.pop();

			// iterate through all adjacent vertices 
			for (int i = 0; i < graph.list[current_vertex].size(); i++) {

				int adjacent_vertex = graph.list[current_vertex][i].first;

				// if the adjacent vertex is not visited, go to the next adjacent vertex
				if (!visited[adjacent_vertex]) {
					visited[adjacent_vertex] = true;
					distances[adjacent_vertex] = distances[current_vertex] + 1; // set the new distance

					queue.push(adjacent_vertex); // add the adjacent vertex to the queue
				}
			}
		}
	}

	DistanceMatrix<int> APSP_breadth_first(StandardGraph<int> graph)
	{
		int N = graph.size;

		DistanceMatrix<int> output(N);

		for (int i = 0; i < N; i++) {

			vector<int> distances(N, std::numeric_limits<int>::max()); //Stores shortest distance, initial distances of infinity
			SSSP_breadth_first(graph, i, distances);

			for (int j = 0; j < N; j++) {
				if (distances[j] != std::numeric_limits<int>::max())
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
	void SSSP_bellman_ford(const StandardGraph<Dtype> &graph, int source, vector<Dtype> &distances)
	{
		int N = graph.size;

		distances[source] = 0;

		for (int i = 1; i <= N - 1; i++) { // repeat N - 1 times
			for (int j = 0; j < N; j++) { // iterate over every node
				for (const pair<int, Dtype> &edge : graph.list[j]) { // relax every edge where that node is the source
					int dest = edge.first;
					Dtype weight = edge.second;

					if (distances[j] != std::numeric_limits<Dtype>::max() && distances[j] + weight < distances[dest]) {
						distances[dest] = distances[j] + weight;
					}
				}
			}
		}
	}

	DistanceMatrix<Dtype> APSP_bellman_ford(const StandardGraph<Dtype> &graph)
	{
		int N = graph.size;

		DistanceMatrix<Dtype> output(N);

		for (int i = 0; i < N; i++) {

			vector<Dtype> distances(N, std::numeric_limits<Dtype>::max());
			SSSP_bellman_ford(graph, i, distances);

			for (int j = 0; j < N; j++) {
				if (distances[j] != std::numeric_limits<int>::max())
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
	DistanceMatrix<Dtype> johnson_dijkstra(StandardGraph<Dtype> &graph)
	{
		// copy the graph. does not need to be done if we are not using the same graph for other algorithms
		StandardGraph<Dtype> copy = graph.deepcopy();

		// add new vertex q (equivalent to the last node in the new graph)
		copy.addVertex();
		int q = copy.size - 1;
		for (int i = 0; i < q; i++) // add an edge with weight 0 from q to every other node
		{
			copy.addEdge(q, i, 0);
		}

		//cout << "\nAdded q: \n";
		//copy.printGraph();

		// run SSSP bellman ford from q to get distance array h
		BellmanFord<Dtype> bellman_ford;
		vector<Dtype> h(copy.size, std::numeric_limits<Dtype>::max());
		bellman_ford.SSSP_bellman_ford(copy, q, h);
		
		//cout << "\nDistance Array: \n";
		//for (Dtype i : h) { cout << i << endl; }

		// remove q from the graph
		copy.removeVertex(q);

		//cout << "\nRemoved q: \n";
		//copy.printGraph();

		// reweight every edge
		for (int i = 0; i < copy.size; i++) { // iterate over every node
			for (const pair<int, Dtype> &edge : copy.list[i]) { // and every edge leaving that node
				copy.reweightEdge(i, edge.first, edge.second + h[i] - h[edge.first]);
				// cout << edge.second << " + " << h[i] << " - " << h[edge.first] << " = " << edge.second + h[i] - h[edge.first] << endl;
			}
		}

		//cout << "\nReweighted: \n";
		//copy.printGraph();

		// run APSP dijkstra
		Dijkstra<Dtype> dijkstra;
		DistanceMatrix<Dtype> output = dijkstra.APSP_dijkstra(copy);

		// reverse the reweighting to get original edge lengths
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
	DistanceMatrix<Dtype> floyd_warshall(const StandardGraph<Dtype> &graph)
	{
		int N = graph.size;

		// create empty distance matrix
		DistanceMatrix<Dtype> output(N);

		// copy input graph into distance matrix
		output.fill_with_graph(graph);

		for (int k = 0; k < N; k++) {
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					Dtype summed = output.getEdge(i, k) + output.getEdge(k, j);

					// catch additon overflow
					if (output.getEdge(k, j) > 0 && output.getEdge(i, k) > std::numeric_limits<Dtype>::max() - output.getEdge(k, j)) {
						summed = std::numeric_limits<Dtype>::max();
					}

					if (output.getEdge(i, j) > summed) {
						output.setEdge(i, j, summed);
					}
				}
			}
		}

		return output;
	}
};





// RUNNING THE PROGRAM

void get_parameters(int &size, double &density, double &minweight, double &maxweight, int &batch_size) {
	std::cout << "\nNumber of vertices (positive integer): ";
	cin >> size;

	std::cout << "\nDensity (0 to 1): ";
	cin >> density;

	std::cout << "\nMin weight: ";
	cin >> minweight;

	std::cout << "\nMax weight: ";
	cin >> maxweight;

	std::cout << "\nBatch size (positive integer): ";
	cin >> batch_size;
}

void get_parameters_gaussian(int &size, double &density, double &mean, double &sd, int &batch_size) {
	std::cout << "\nNumber of vertices (positive integer): ";
	cin >> size;

	std::cout << "\nDensity (0 to 1): ";
	cin >> density;

	std::cout << "\nMean: ";
	cin >> mean;

	std::cout << "\nStandard Deviation: ";
	cin >> sd;

	std::cout << "\nBatch size (positive integer): ";
	cin >> batch_size;
}

void get_parameters_no_weights(int &size, double &density, int &batch_size) {
	std::cout << "\nNumber of vertices (0 to 5000): ";
	cin >> size;

	std::cout << "\nDensity (0 to 1): ";
	cin >> density;

	std::cout << "\nBatch size (0 to x): ";
	cin >> batch_size;
}

// performance testing on graphs with double weigts - inluding the Bellman-Ford algorithm
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
		std::cout << "\nIteration " << i + 1 << " started:\n";

		StandardGraph<double> graph = graph_gen.generateERGraphDouble(size, density, minweight, maxweight);

		// Binary Dijkstra
		std::cout << "\n - Binary Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra(graph);
		t_end = Clock::now();

		total_times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Fibonacci Dijkstra
		std::cout << "\n - Fibonacci Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra_fibqueue(graph);
		t_end = Clock::now();

		total_times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Bellman Ford
		std::cout << "\n - Bellman Ford started\n";
		t_start = Clock::now();
		bellman_ford.APSP_bellman_ford(graph);
		t_end = Clock::now();

		total_times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Johnson Dijkstra
		std::cout << "\n - Johnson Dijkstra started\n";
		t_start = Clock::now();
		johnson_dijkstra.johnson_dijkstra(graph);
		t_end = Clock::now();

		total_times[3] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
		
		// FloydWarshall
		std::cout << "\n - Floyd Warshall started\n";
		t_start = Clock::now();
		floyd_warshall.floyd_warshall(graph);
		t_end = Clock::now();
		
		total_times[4] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		std::cout << "\n - Iteration " << i + 1 << " ended\n";
	}

	std::cout << "\nResults:\n";

	std::cout << " - Time to solve Binary Heap Dijkstra: " << int(total_times[0]/batch_size) << " seconds \n";

	std::cout << " - Time to solve Fibonacci Heap: " << int(total_times[1] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Bellman Ford: " << int(total_times[2] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Johnson Dijkstra: " << int(total_times[3] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Floyd Warshall: " << int(total_times[4] / batch_size) << " seconds \n";

	std::cout << "\n";
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
		std::cout << "\nIteration " << i + 1 << " started:\n";

		StandardGraph<int> graph = graph_gen.generateERGraphInt(size, density, minweight, maxweight);

		// Binary Dijkstra
		std::cout << "\n - Binary Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra(graph);
		t_end = Clock::now();

		total_times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Fibonacci Dijkstra
		std::cout << "\n - Fibonacci Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra_fibqueue(graph);
		t_end = Clock::now();

		total_times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Johnson Dijkstra
		std::cout << "\n - Johnson Dijkstra started\n";
		t_start = Clock::now();
		johnson_dijkstra.johnson_dijkstra(graph);
		t_end = Clock::now();

		total_times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// FloydWarshall
		std::cout << "\n - Floyd Warshall started\n";
		t_start = Clock::now();
		floyd_warshall.floyd_warshall(graph);
		t_end = Clock::now();

		total_times[3] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		std::cout << "\n - Iteration " << i + 1 << " ended\n";
	}

	std::cout << "\nResults:\n";

	std::cout << " - Time to solve Binary Heap Dijkstra: " << int(total_times[0] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Fibonacci Heap: " << int(total_times[1] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Johnson Dijkstra: " << int(total_times[2] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Floyd Warshall: " << int(total_times[3] / batch_size) << " seconds \n";

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
		std::cout << "\nIteration " << i + 1 << " started:\n";

		StandardGraph<float> graph = graph_gen.generateERGraphFloat(size, density, minweight, maxweight);

		// Binary Dijkstra
		std::cout << "\n - Binary Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra(graph);
		t_end = Clock::now();

		total_times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Fibonacci Dijkstra
		std::cout << "\n - Fibonacci Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra_fibqueue(graph);
		t_end = Clock::now();

		total_times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Johnson Dijkstra
		std::cout << "\n - Johnson Dijkstra started\n";
		t_start = Clock::now();
		johnson_dijkstra.johnson_dijkstra(graph);
		t_end = Clock::now();

		total_times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// FloydWarshall
		std::cout << "\n - Floyd Warshall started\n";
		t_start = Clock::now();
		floyd_warshall.floyd_warshall(graph);
		t_end = Clock::now();

		total_times[3] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		std::cout << "\n - Iteration " << i + 1 << " ended\n";
	}

	std::cout << "\nResults:\n";

	std::cout << " - Time to solve Binary Heap Dijkstra: " << int(total_times[0] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Fibonacci Heap: " << int(total_times[1] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Johnson Dijkstra: " << int(total_times[2] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Floyd Warshall: " << int(total_times[3] / batch_size) << " seconds \n";

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
		std::cout << "\nIteration " << i + 1 << " started:\n";

		StandardGraph<double> graph = graph_gen.generateERGraphDouble(size, density, minweight, maxweight);

		// Binary Dijkstra
		std::cout << "\n - Binary Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra(graph);
		t_end = Clock::now();

		total_times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Fibonacci Dijkstra
		std::cout << "\n - Fibonacci Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra_fibqueue(graph);
		t_end = Clock::now();

		total_times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Johnson Dijkstra
		std::cout << "\n - Johnson Dijkstra started\n";
		t_start = Clock::now();
		johnson_dijkstra.johnson_dijkstra(graph);
		t_end = Clock::now();

		total_times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// FloydWarshall
		std::cout << "\n - Floyd Warshall started\n";
		t_start = Clock::now();
		floyd_warshall.floyd_warshall(graph);
		t_end = Clock::now();

		total_times[3] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		std::cout << "\n - Iteration " << i + 1 << " ended\n";
	}

	std::cout << "\nResults:\n";

	std::cout << " - Time to solve Binary Heap Dijkstra: " << int(total_times[0] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Fibonacci Heap: " << int(total_times[1] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Johnson Dijkstra: " << int(total_times[2] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Floyd Warshall: " << int(total_times[3] / batch_size) << " seconds \n";

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
		std::cout << "\nIteration " << i + 1 << " started:\n";

		StandardGraph<double> graph = graph_gen.generateERGraphGaussian(size, density, mean, sd);

		// Binary Dijkstra
		std::cout << "\n - Binary Heap Dijkstra started\n";
		t_start = Clock::now();
		//dijkstra.APSP_dijkstra(graph);
		t_end = Clock::now();

		total_times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Fibonacci Dijkstra
		std::cout << "\n - Fibonacci Heap Dijkstra started\n";
		t_start = Clock::now();
		//dijkstra.APSP_dijkstra_fibqueue(graph);
		t_end = Clock::now();

		total_times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Johnson Dijkstra
		std::cout << "\n - Johnson Dijkstra started\n";
		t_start = Clock::now();
		johnson_dijkstra.johnson_dijkstra(graph);
		t_end = Clock::now();

		total_times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// FloydWarshall
		std::cout << "\n - Floyd Warshall started\n";
		t_start = Clock::now();
		//floyd_warshall.floyd_warshall(graph);
		t_end = Clock::now();

		total_times[3] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		std::cout << "\n - Iteration " << i + 1 << " ended\n";
	}

	std::cout << "\nResults:\n";

	std::cout << " - Time to solve Binary Heap Dijkstra: " << int(total_times[0] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Fibonacci Heap: " << int(total_times[1] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Johnson Dijkstra: " << int(total_times[2] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Floyd Warshall: " << int(total_times[3] / batch_size) << " seconds \n";

	std::cout << "\n";
}

// performance testing on constant weighted graphs, analagous to unweighted, to compare breadth-first search to other algorithms
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
		std::cout << "\nIteration " << i + 1 << " started:\n";

		StandardGraph<int> graph = graph_gen.generateERGraphInt(size, density, 1, 1);

		// Binary Dijkstra
		std::cout << "\n - Binary Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra(graph);
		t_end = Clock::now();

		total_times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Fibonacci Dijkstra
		std::cout << "\n - Fibonacci Heap Dijkstra started\n";
		t_start = Clock::now();
		dijkstra.APSP_dijkstra_fibqueue(graph);
		t_end = Clock::now();

		total_times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// Johnson Dijkstra
		std::cout << "\n - Johnson Dijkstra started\n";
		t_start = Clock::now();
		johnson_dijkstra.johnson_dijkstra(graph);
		t_end = Clock::now();

		total_times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// FloydWarshall
		std::cout << "\n - Floyd Warshall started\n";
		t_start = Clock::now();
		floyd_warshall.floyd_warshall(graph);
		t_end = Clock::now();

		total_times[3] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		// BreadthFirst
		std::cout << "\n - Breadth First started\n";
		t_start = Clock::now();
		breadth_first.APSP_breadth_first(graph);
		t_end = Clock::now();

		total_times[4] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

		std::cout << "\n - Iteration " << i + 1 << " ended\n";
	}

	std::cout << "\nResults:\n";

	std::cout << " - Time to solve Binary Heap Dijkstra: " << int(total_times[0] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Fibonacci Heap: " << int(total_times[1] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Johnson Dijkstra: " << int(total_times[2] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Floyd Warshall: " << int(total_times[3] / batch_size) << " seconds \n";

	std::cout << " - Time to solve Breadth First: " << int(total_times[4] / batch_size) << " seconds \n";

	std::cout << "\n";
}

// performance testing on graphs of real data loaded from files
void run_on_real_graph(string filename)
{
	Dijkstra<double> dijkstra;
	JohnsonDijkstra<double> johnson_dijkstra;
	FloydWarshall<double> floyd_warshall;

	auto t_start = Clock::now();
	auto t_end = Clock::now();

	double total_times[10] = { 0 };

	StandardGraph<double> graph = loadRealGraph(filename);

	// Binary Dijkstra
	std::cout << "\n - Binary Heap Dijkstra started\n";
	t_start = Clock::now();
	dijkstra.APSP_dijkstra(graph);
	t_end = Clock::now();

	total_times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

	// Fibonacci Dijkstra
	std::cout << "\n - Fibonacci Heap Dijkstra started\n";
	t_start = Clock::now();
	dijkstra.APSP_dijkstra_fibqueue(graph);
	t_end = Clock::now();

	total_times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();


	// Johnson Dijkstra
	std::cout << "\n - Johnson Dijkstra started\n";
	t_start = Clock::now();
	johnson_dijkstra.johnson_dijkstra(graph);
	t_end = Clock::now();

	total_times[2] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

	// FloydWarshall
	std::cout << "\n - Floyd Warshall started\n";
	t_start = Clock::now();
	//floyd_warshall.floyd_warshall(graph);
	t_end = Clock::now();

	total_times[3] += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();


	std::cout << "\nResults:\n";

	std::cout << " - Time to solve Binary Heap Dijkstra: " << int(total_times[0]) << " seconds \n";

	std::cout << " - Time to solve Fibonacci Heap: " << int(total_times[1]) << " seconds \n";

	std::cout << " - Time to solve Johnson Dijkstra: " << int(total_times[2]) << " seconds \n";

	std::cout << " - Time to solve Floyd Warshall: " << int(total_times[3]) << " seconds \n";

	std::cout << "\n";
}


// print verification that the algorithms produce the correct result for known integer weighted graphs
void verify_algorithms_int()
{
	// declare algorithm structures
	Dijkstra<int> dijkstra_int;
	BellmanFord<int> bellman_ford_int;
	JohnsonDijkstra<int> johnson_dijkstra_int;
	FloydWarshall<int> floyd_warshall_int;

	// load graphs
	StandardGraph<int> graph = loadKnownGraph();
	StandardGraph<int> solution = loadKnownSolution();;

	std::cout << "VERIFYING ALGORITHMS\n";

	// print adjacency list representation of graph
	std::cout << "Adjacency List of graph: \n";
	graph.printGraph();
	std::cout << "\n";

	// print distance matrix representation of graph
	std::cout << "Distance Matrix of graph: \n";
	DistanceMatrix<int> dist(graph.size);
	dist.fill_with_graph(graph);
	dist.printMatrix();

	std::cout << "Distance Matrix of solution: \n";
	DistanceMatrix<int> solution_dist(solution.size);
	solution_dist.fill_with_graph(solution);
	solution_dist.printMatrix();


	std::cout << "Floyd-Warshall solution: \n";
	DistanceMatrix<int> floyd_output = floyd_warshall_int.floyd_warshall(graph);
	floyd_output.printMatrix();


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
}

// print verification that the algorithms produce the correct result for known integer weighted graphs with negative weights
void verify_algorithms_negatives()
{
	// declare algorithm structures
	BellmanFord<int> bellman_ford_int;
	JohnsonDijkstra<int> johnson_dijkstra_int;
	FloydWarshall<int> floyd_warshall_int;

	// load graphs
	StandardGraph<int> graph = loadKnownGraphNegatives();
	StandardGraph<int> solution = loadKnownSolutionNegatives();

	std::cout << "VERIFYING ALGORITHMS ON GRAPHS WITH NEGATIVE WEIGHTS\n";

	// print adjacency list representation of graph
	std::cout << "Adjacency List of graph: \n";
	graph.printGraph();
	std::cout << "\n";

	// print distance matrix representation of graph
	std::cout << "Distance Matrix of graph: \n";
	DistanceMatrix<int> dist(graph.size);
	dist.fill_with_graph(graph);
	dist.printMatrix();

	std::cout << "Distance Matrix of solution: \n";
	DistanceMatrix<int> solution_dist(solution.size);
	solution_dist.fill_with_graph(solution);
	solution_dist.printMatrix();


	std::cout << "Floyd-Warshall solution: \n";
	DistanceMatrix<int> floyd_output = floyd_warshall_int.floyd_warshall(graph);
	floyd_output.printMatrix();

	std::cout << "APSP Bellman-Ford solution: \n";
	DistanceMatrix<int> bellman_output = bellman_ford_int.APSP_bellman_ford(graph);
	bellman_output.printMatrix();

	std::cout << "Johnson-Dijkstra solution: \n";
	DistanceMatrix<int> johnson_output = johnson_dijkstra_int.johnson_dijkstra(graph);
	johnson_output.printMatrix();
}

// print verification that the algorithms produce the correct result for known double weighted graphs
void verify_algorithms_doubles()
{
	// declare algorithm structures
	Dijkstra<double> dijkstra_double;
	BellmanFord<double> bellman_ford_double;
	JohnsonDijkstra<double> johnson_dijkstra_double;
	FloydWarshall<double> floyd_warshall_double;

	// load graphs
	StandardGraph<double> double_graph = loadKnownDoubleGraph();
	StandardGraph<int> solution = loadKnownSolution();

	std::cout << "VERIFYING ALGORITHMS ON GRAPHS WITH DOUBLE WEIGHTS\n";

	// print adjacency list representation of graph
	std::cout << "Adjacency List of double graph: \n";
	double_graph.printGraph();
	std::cout << "\n";

	// print distance matrix representation of graph
	std::cout << "Distance Matrix of double graph: \n";
	DistanceMatrix<double> dist(double_graph.size);
	dist.fill_with_graph(double_graph);
	dist.printMatrix();

	std::cout << "Distance Matrix of solution: \n";
	DistanceMatrix<int> solution_dist(solution.size);
	solution_dist.fill_with_graph(solution);
	solution_dist.printMatrix();

	std::cout << "Floyd-Warshall solution: \n";
	DistanceMatrix<double> floyd_double_output = floyd_warshall_double.floyd_warshall(double_graph);
	floyd_double_output.printMatrix();

	std::cout << "APSP Dijkstra solution binary queue: \n";
	DistanceMatrix<double> dijkstra_double_output = dijkstra_double.APSP_dijkstra(double_graph);
	dijkstra_double_output.printMatrix();

	std::cout << "APSP Dijkstra fibonaci queue solution: \n";
	DistanceMatrix<double> dijkstra_output_fib_double = dijkstra_double.APSP_dijkstra_fibqueue(double_graph);
	dijkstra_output_fib_double.printMatrix();

	std::cout << "APSP Bellman-Ford solution: \n";
	DistanceMatrix<double> bellman_output = bellman_ford_double.APSP_bellman_ford(double_graph);
	bellman_output.printMatrix();

	std::cout << "Johnson-Dijkstra solution: \n";
	DistanceMatrix<double> johnson_output = johnson_dijkstra_double.johnson_dijkstra(double_graph);
	johnson_output.printMatrix();
}

// print verification that the algorithms produce the correct result for known integer weighted graphs
void verify_algorithms_unweighted()
{
	// declare algorithm structures
	Dijkstra<int> dijkstra_int;
	BellmanFord<int> bellman_ford_int;
	JohnsonDijkstra<int> johnson_dijkstra_int;
	FloydWarshall<int> floyd_warshall_int;
	BreadthFirst breadth_first;

	// load graphs
	StandardGraph<int> graph = loadKnownGraph();
	StandardGraph<int> solution = loadKnownSolution();;

	std::cout << "VERIFYING ALGORITHMS\n";

	// print adjacency list representation of graph
	std::cout << "Adjacency List of graph: \n";
	graph.printGraph();
	std::cout << "\n";

	// print distance matrix representation of graph
	std::cout << "Distance Matrix of graph: \n";
	DistanceMatrix<int> dist(graph.size);
	dist.fill_with_graph(graph);
	dist.printMatrix();

	std::cout << "Distance Matrix of solution: \n";
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

// print verification that the algorithms produce the correct result for randomly generated graphs with integer weights
void verify_graph_generation_int()
{
	std::cout << "\nVERIFY GENERATED GRAPHS ARE CORRECT\n";

	GraphGenerator graph_gen;

	std::cout << "sparse int ER graph (uniform weights 0 to 10)\n";
	StandardGraph<int> graph = graph_gen.generateERGraphInt(8, 0.2, 0, 10);
	graph.printGraph();

	std::cout << "dense int ER graph (uniforms weights 0 to 10)\n";
	graph = graph_gen.generateERGraphInt(8, 0.95, 0, 10);
	graph.printGraph();

	std::cout << "sparse double ER graph (uniform weights 0 to 1)\n";
	StandardGraph<double> graph_double = graph_gen.generateERGraphDouble(8, 0.2, 0, 1);
	graph_double.printGraph();

	std::cout << "dense double ER graph (uniform weights 0 to 1)\n";
	graph_double = graph_gen.generateERGraphDouble(8, 1, 0, 1);
	graph_double.printGraph();

	std::cout << "gaussian ER graph (mean = 10, sd = 1)\n";
	graph_double = graph_gen.generateERGraphGaussian(8, 0.5, 10, 1);
	graph_double.printGraph();
}

// print verification that the algorithms produce the correct result for randomly generated graphs with double weights
void verify_graph_generation_doubles()
{
	std::cout << "\nVERIFY GENERATED GRAPHS ARE CORRECT\n";

	GraphGenerator graph_gen;

	std::cout << "sparse double ER graph (uniform weights 0 to 100)\n";
	StandardGraph<double> graph_double = graph_gen.generateERGraphDouble(10, 0.2, 0, 100);
	graph_double.printGraph();

	std::cout << "dense double ER graph (uniform weights 0 to 100)\n";
	graph_double = graph_gen.generateERGraphDouble(10, 0.95, 0, 100);
	graph_double.printGraph();

	std::cout << "gaussian ER graph (mean = 10, sd = 1)\n";
	graph_double = graph_gen.generateERGraphGaussian(10, 0.95, 100, 10);
	graph_double.printGraph();
}



int main()
{
	// TODO - function to choose a function

	//verify_algorithms_int();
	//verify_algorithms_doubles();
	//verify_algorithms_negatives();
	verify_algorithms_unweighted();
	//verify_graph_generation();
	//verify_graph_generation_doubles();


	std::cout << "\n\n" << "TESTING ON RANDOM GRAPHS\n";
	//run_batch_double();
	//run_batch_float();
	//run_batch_int();
	//run_batch_gaussian();
	//run_batch_with_bellman();
	//run_batch_unweighted();
	//run_on_real_graph("openflights.txt");

	int any;
	std::cout << "Close with any key: ";
	cin >> any;

	return 0;
}


// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
