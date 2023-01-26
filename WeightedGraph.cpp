#include <bits/stdc++.h>
using namespace std;

class WeightedGraph {
    int n, m;

    // pairs have the following form : (cost, neighbour)
    vector<vector<pair<int, int>>> adj;
    
    //(cost, firstNode, secondNode)
    vector<tuple<int, int, int>> edges;
    vector<bool> viz;
    
    // the representant of each vertex during kruskal
    int rep(int x, vector<int>& parent);

    // union of disjoint sets during kruskal
    void setUnion(int x, int y, vector<int>& parent, vector<int>& height);
    
    // dfs used in topSort()
    void topSortDfs(int node, stack<int>& topSortStack);

public:
    static const int inf;
    enum type{directed, undirected};   
    
    // returns the number of nodes in the same connected component as source, excluding source
    int countNodesAccessibleFrom(int source);

    // returns the distance matrix from the adjecency lists 
    vector<vector<int>> toDistanceMatrix();

    // edges must be in the form {cost, firstNode, secondNode}
    WeightedGraph(int _n, int _m, vector<tuple<int, int, int>>& edges, WeightedGraph::type type);

    // returns (total_cost, vector_of_edges)
    pair<int, vector<pair<int, int>>> kruskal();
    
    // returns (total_cost, vector_of_edges)
    pair<int, vector<pair<int, int>>> prim();

    // returns (dist_vector, parent_vector)
    pair<vector<int>, vector<int>> dijkstra(int source);

    // returns (dist_from_src_to_dest, parent_vector)
    pair<int, vector<int>> dijkstra2(int source, int dest);

    // returns (dist_vector, parent_vector)
    pair<vector<int>, vector<int>> distanceDAG(int source);

    // returns dist_vector, parent_vector) -------> if negative cycle returns a pair of empty vectors
    pair<vector<int>, vector<int>> bellmanFord(int source);

    // returns (dist_matrix, pred_matrix) --------> if negative cycle returns a pair of empty matrices
    pair<vector<vector<int>>, vector<vector<int>>> royFloyd();

    // returns a topological sort of the graph
    vector<int> topSort();

    // calculates the minimum cost to connect all nodes to a source (directly or indirectly), given a set of sources
    // returns (source_with_most_nodes, edges of the partial graph)
    pair<int, vector<pair<int, int>>> multiPrim(unordered_set<int>& sources);
};

const int WeightedGraph::inf = INT_MAX;

WeightedGraph::WeightedGraph(int _n, int _m, vector<tuple<int, int, int>>& _edges, WeightedGraph::type type): n(_n), m(_m), edges(_edges) {
    adj.resize(n + 1);
    viz.resize(n + 1);
    if(type == directed)
        for(auto edge : _edges) 
            adj[get<1>(edge)].push_back({get<0>(edge), get<2>(edge)});

    else
        for(auto edge : _edges) {
            adj[get<1>(edge)].push_back({get<0>(edge), get<2>(edge)});
            adj[get<2>(edge)].push_back({get<0>(edge), get<1>(edge)});
        }
}

void WeightedGraph::topSortDfs(int node, stack<int>& topSortStack) {
    viz[node] = true;
    for(auto nextPair : adj[node]) 
        if(!viz[nextPair.second])
            topSortDfs(nextPair.second, topSortStack);

    topSortStack.push(node);
}

vector<int> WeightedGraph::topSort() {
    stack<int> topSortStack;
    fill(viz.begin(), viz.end(), false);
    for(int i = 1; i <= n; ++ i) {
        if(!viz[i])
            topSortDfs(i, topSortStack);
    }

    vector<int> result;
    while(!topSortStack.empty()) {
        result.push_back(topSortStack.top());
        topSortStack.pop();
    }

    return result;
}

int WeightedGraph::rep(int x, vector<int>& parent) {
    if (!parent[x])
        return x;
    
    int representative = rep(parent[x], parent);
    parent[x] = representative;
    return representative;
}

void WeightedGraph::setUnion(int x, int y, vector<int>& parent, vector<int>& height) {
    int repx = rep(x, parent), repy = rep(y, parent);

    if (height[repx] < height[repy])
        parent[repx] = repy;

    else if (height[repy] < height[repx])
        parent[repy] = repx;

    else {
        parent[repx] = repy;
        ++ height[repy];
    }
}

pair<int, vector<pair<int, int>>> WeightedGraph::kruskal() {
    vector<int> parent(n + 1, 0), height(n + 1, 0);
    sort(edges.begin(), edges.end(), [](auto a, auto b) {return get<0>(a) < get<0>(b);});

    int totalCost = 0;
    vector<pair<int, int>> sol;
    for(auto edge : edges) {
        int node1 = get<1>(edge), node2 = get<2>(edge), cost = get<0>(edge);
        if(rep(node1, parent) == rep(node2, parent))
            continue;
        
        setUnion(node1, node2, parent, height);
        totalCost += cost;
        sol.push_back({node1, node2});
    }

    return make_pair(totalCost, sol);
}

pair<int, vector<pair<int, int>>> WeightedGraph::prim() {
    vector<bool> inTree(n + 1, false);
    priority_queue<tuple<int, int, int>, vector<tuple<int, int, int>>, greater<tuple<int, int, int>>> q;
    vector<pair<int, int>> sol;
    int totalCost = 0;

     q.push({0, -1, 1});
    while(!q.empty() && sol.size() < n - 1) {
        int currCost = get<0>(q.top()), prevNode = get<1>(q.top()), currNode = get<2>(q.top());
        q.pop();
        
        if(inTree[currNode])
            continue;

        totalCost += currCost;
        if(prevNode != -1)
            sol.push_back({prevNode, currNode});

        inTree[currNode] = true;
        for(auto nextPair : adj[currNode]) {
            int nextNode = nextPair.second, nextCost = nextPair.first;
            q.push({nextCost, currNode, nextNode});
        }
    }

    return make_pair(totalCost, sol);
}

// first vector is distance, the second is previous node of each node
pair<vector<int>, vector<int>> WeightedGraph::dijkstra(int source) {
    vector<int> dist(n + 1, inf), prec(n + 1, 0);
    fill(viz.begin(), viz.end(), false);
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> q;
    dist[source] = 0;
    q.push({0, source});
    while(!q.empty()) {
        int node = q.top().second, currDist = q.top().first;
        q.pop();

        if(viz[node])
            continue;
        
        viz[node] = true;
        for(auto& nextPair : adj[node]) {
            if(dist[nextPair.second] > dist[node] + nextPair.first) {
                dist[nextPair.second] = dist[node] + nextPair.first;
                prec[nextPair.second] = node;
                q.push({dist[nextPair.second], nextPair.second});
            }
        }
    }

    return make_pair(dist, prec);
}

pair<int, vector<int>> WeightedGraph::dijkstra2(int source, int dest) {
    vector<int> dist(n + 1, inf), prec(n + 1, 0);
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> q;
    dist[source] = 0;
    q.push({0, source});
    while(!q.empty()) {
        int node = q.top().second, currDist = q.top().first;
        q.pop();

        if(node == dest) {
            return make_pair(currDist, prec);
        }

        for(auto& nextPair : adj[node]) {
            if(dist[nextPair.second] > dist[node] + nextPair.first) {
                dist[nextPair.second] = dist[node] + nextPair.first;
                prec[nextPair.second] = node;
                q.push({dist[nextPair.second], nextPair.second});
            }
        }
    }

    return make_pair(-1, vector<int>());
}

pair<vector<int>, vector<int>> WeightedGraph::distanceDAG(int source) {
    vector<int> dist(n + 1, inf), prec(n + 1, 0);
    dist[source] = 0;
    vector<int> sorted = topSort();

    for(auto node : sorted) {
        if(dist[node] == inf)
            continue;
        for(auto nextPair : adj[node]) {
            int nextNode = nextPair.second, nextDist = nextPair.first;
            if(dist[nextNode] > dist[node] + nextDist) {
                dist[nextNode] = dist[node] + nextDist;
                prec[nextNode] = node;
            }
        }
    }

    return make_pair(dist, prec);
}

// first vector is distance vector; second vector is previous node vector
pair<vector<int>, vector<int>> WeightedGraph::bellmanFord(int source) {
    queue<int> q;
    vector<int> dist(n + 1, inf), prec(n + 1, 0), count(n + 1, 0);
    dist[source] = 0;
    q.push(source);

    while(!q.empty()) {
        int currNode = q.front();
        q.pop();

        ++ count[currNode];
        if(count[currNode] == n)
            return make_pair(vector<int>(), vector<int>());     // negative cycle

        for(auto nextPair : adj[currNode]) {
            int nextNode = nextPair.second, nextDist = nextPair.first;
            if(dist[nextNode] > dist[currNode] + nextDist) {
                dist[nextNode] = dist[currNode] + nextDist;
                prec[nextNode] = currNode;
                q.push(nextNode);
            }
        }
    }

    return make_pair(dist, prec);
}

vector<vector<int>> WeightedGraph::toDistanceMatrix() {
    vector<vector<int>> matrix(n + 1, vector<int>(n + 1, inf));
    for(int i = 1; i <= n; ++ i)
        matrix[i][i] = 0;
    
    for(int i = 1; i <= n; ++ i)
        for(auto Pair : adj[i]) 
            matrix[i][Pair.second] = Pair.first;
    
    return matrix;
}

pair<vector<vector<int>>, vector<vector<int>>> WeightedGraph::royFloyd() {
    auto w = toDistanceMatrix();
    vector<vector<int>> prec(n + 1, vector<int>(n + 1)), dist(w);
    
    for(int i = 1; i <= n; ++ i) 
        for(int j = 1; j <= n; ++ j) {
            if(i == j)
                continue;

            if(dist[i][j] == inf)
                prec[i][j] = 0;
            else
                prec[i][j] = i;
        }
            
    

    for(int k = 1; k <= n; ++ k) 
        for(int i = 1; i <= n; ++ i)
            for(int j = 1; j <= n; ++ j)
                if(dist[i][k] != inf && dist[k][j] != inf && dist[i][j] > dist[i][k] + dist[k][j]) {
                    dist[i][j] = dist[i][k] + dist[k][j];
                    prec[i][j] = prec[k][j];
                }

    for(int i = 1; i <= n; ++ i)
        if(dist[i][i] < 0)
            return make_pair(vector<vector<int>>(), vector<vector<int>>());

    return make_pair(dist, prec);
}

pair<int, vector<pair<int, int>>> WeightedGraph::multiPrim(unordered_set<int>& sources) {
    vector<bool> inTree(n + 1, false);
    unordered_map<int, int> sourceNodes;  // tine minte sursa din care pleci cumva si ea o sa se transmita din nod in nod ish :()
    priority_queue<tuple<int, int, int, int>, vector<tuple<int, int, int, int>>, greater<tuple<int, int, int, int>>> q;
    vector<pair<int, int>> sol;
    int totalCost = 0;
    for(auto source : sources) {
        sourceNodes[source] = 0;
        q.push({0, -1, source, source});
    }

    while(!q.empty() && sol.size() < n - (int) sources.size()) {
        int currCost = get<0>(q.top()), prevNode = get<1>(q.top()), currNode = get<2>(q.top()), ps = get<3>(q.top());
        q.pop();
        
        if(inTree[currNode])
            continue;

        totalCost += currCost;
        if(prevNode != -1)
            sol.push_back({prevNode, currNode});

        sourceNodes[ps] ++;
        inTree[currNode] = true;
        for(auto nextPair : adj[currNode]) {
            int nextNode = nextPair.second, nextCost = nextPair.first;
            if(sources.find(currNode) != sources.end() && sources.find(nextNode) != sources.end())
                continue;
            q.push({nextCost, currNode, nextNode, ps});
        }
    }

    int maxPowerSource, maxPowerSourceNodes = 0;
    for(auto it : sourceNodes)
        if(maxPowerSourceNodes < it.second) {
            maxPowerSourceNodes = it.second;
            maxPowerSource = it.first;
        }

    return make_pair(maxPowerSource, sol);
}

int WeightedGraph::countNodesAccessibleFrom(int source) {
    fill(viz.begin(), viz.end(), false);
    queue<int> q;
    viz[source] = true;
    int noVisited = 0;

    while(!q.empty()) {
        int currNode = q.front();
        q.pop();

        ++ noVisited;

        for(auto nextPair : adj[currNode]) {
            if(!viz[nextPair.second]) {
                viz[nextPair.second] = true;
                q.push(nextPair.second);
            }
        }
    }   

    return noVisited - 1;
}

int main() {

    return 0;
}