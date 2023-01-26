#include <bits/stdc++.h>
using namespace std;

class UnweightedGraph {
    int n, m;
    vector<vector<int>> adj;
    vector<bool> viz;
    bool cycle;

    // dfs functions useful for public functions
    void hasCycleDfs(int node, int parent, vector<bool>& added);
    void topSortDfs(int node, stack<int>& topSortStack);
    void criticalNodesDfs(int node, int parent, vector<int>& result, vector<int>& lvl, vector<int>& minLvl, int root);
    void criticalEdgesDfs(int node, int parent, vector<pair<int, int>>& result, vector<int>& lvl, vector<int>& minLvl);
    void CCDFS(int node, vector<int>& component);

    // computes an eulerian cycle/path on an undirected graph -> used in getEulerianCycle()
    deque<int> getEulerianUndirected(bool requiresCycle);

    // computes an eulerian cycle/path on an undirected graph -> used in getEulerianCycle()
    deque<int> getEulerianDirected(bool requiresCycle);

    // useful in getEulerianUndirected; transforms the input in a string like "first - second"
    string edgeToString(int first, int second) const;       
public:
    enum type{directed, undirected};
private: 
    type Type;
public:

    // adds an edge to the graph (x - y if undirected; x -> y if directed)
    void addEdge(int x, int y);

    // returns a vector that contains the connected components of the graph
    vector<vector<int>> getConnectedComponents();

    // n is the number of nodes in the graph
    // m is the number of edges
    // each pair in edges must have the following form: (start_node, end_node)
    UnweightedGraph(int _n, int _m, vector<pair<int, int>>& edges, UnweightedGraph::type type);

    // simple implementation of dfs (completely useless on it's own)
    void dfs(int node);

    // returns the maximum distance between the source and any other node in the graph
    int bfs(int startNode);

    // returns (distance, path_from_src_to_dest) ----> if there is no path between them returns (-1, empty vector) 
    pair<int, vector<int>> bfs2(int src, int dest);

    // returns the adjacency list of the transpose graph
    vector<vector<int>> transpose();

    bool isBipartite();

    // returns a topological sort of the graph
    vector<int> topSort();

    bool hasCycle();

    // returns the critical nodes of the graph
    vector<int> criticalNodes();

    // the elements of a pair represents the ends of a critical edge
    vector<pair<int, int>> criticalEdges();

    // returns a vector containing Strongly Connected Components
    vector<vector<int>> getSCComponents();

    // the dfs function that accumulates the nodes of a Strongly Connected Component in component parameter
    void getSCC(int node, vector<vector<int>>& adjt, vector<int>& component);

    // returns a vector containing the Biconnected Components (each biconnected component is defined by it's edges)
    vector<vector<pair<int, int>>> getBCComponents();

    // the dfs function that accumulates the edges of a Biconnected Component
    void getBCC(int node, int parent, vector<vector<pair<int, int>>>& result, vector<int>& lvl, vector<int>& minLvl, stack<pair<int, int>>& st, int root);
    
    // returns the euleria cycle (nodes) -------> if no eulerian cycle returns empty vector
    deque<int> getEulerianCycle();

    // returns the euleria paht (nodes) -------> if no eulerian paht returns empty vector
    deque<int> getEulerianPath();
};

UnweightedGraph::UnweightedGraph(int _n, int _m, vector<pair<int, int>>& edges, UnweightedGraph::type type): n(_n), m(_m), Type(type) {
    adj.resize(n + 1);
    viz.resize(n + 1, false);
    
    if(type == directed)
        for(auto edge : edges) {
            adj[edge.first].push_back(edge.second);
        }
    
    else
        for(auto edge : edges) {
            adj[edge.first].push_back(edge.second);
            adj[edge.second].push_back(edge.first);
        }
}

string UnweightedGraph::edgeToString(int first, int second) const {
    if(first > second)
        return to_string(second) + "-" + to_string(first);
    
    return to_string(first) + "-" + to_string(second);
}

void UnweightedGraph::dfs(int node)  {
        viz[node] = true;
        for(auto nextNode : adj[node])
            if(!viz[nextNode])
                dfs(nextNode);
}

int UnweightedGraph::bfs(int startNode) {
    fill(viz.begin(), viz.end(), false);
    vector<int> dist(n + 1, 0);
    queue<int> q;
    q.push(startNode);
    viz[startNode] = true;
    dist[startNode] = 0;
    int maxDist = 0;

    while(!q.empty()) {
        int currNode = q.front();
        q.pop();
        maxDist = max(maxDist, dist[currNode]);

        for(auto nextNode : adj[currNode]) {
            if(!viz[nextNode]) {
                viz[nextNode] = true;
                dist[nextNode] = dist[currNode] + 1;
                q.push(nextNode);
            }
        }
    }

    return maxDist;
}

vector<vector<int>> UnweightedGraph::transpose() {
    vector<vector<int>> adjt(n + 1);
    for(int i = 1; i <= n; ++ i){
        for(auto neighbour : adj[i])
            adjt[neighbour].push_back(i);
    }

    return adjt;
}

bool UnweightedGraph::isBipartite() {
    queue<int> q;
    vector<int> color(n + 1, 0);

    for(int i = 1; i <= n; ++ i) {
        if(color[i])
            continue;
        
        q.push(i);
        color[i] = 1;
        while(!q.empty()) {
            int currNode = q.front();
            q.pop();

            for(auto nextNode : adj[currNode])
                if(color[nextNode] == color[currNode])
                    return false;
                else if(!color[nextNode]) {
                    color[nextNode] = 3 - color[currNode];
                    q.push(nextNode);
                }
        }
    }

    return true;
}

void UnweightedGraph::topSortDfs(int node, stack<int>& topSortStack) {
    viz[node] = true;
    for(auto nextNode : adj[node]) 
        if(!viz[nextNode])
            topSortDfs(nextNode, topSortStack);

    topSortStack.push(node);
}

vector<int> UnweightedGraph::topSort() {
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

void UnweightedGraph::hasCycleDfs(int node, int parent, vector<bool>& added) {
    viz[node] = true;
    for(auto nextNode : adj[node]) {
        if(nextNode == parent)
            continue;

        if(!viz[nextNode]) {
            added[nextNode] = true;
            hasCycleDfs(nextNode, node, added);
            added[nextNode] = false;
        }

        if(added[nextNode]) {
            cycle = true;
            return;
        }
    }
}

bool UnweightedGraph::hasCycle() {
    vector<bool> added(n + 1, false);
    fill(viz.begin(), viz.end(), false);
    cycle = false;

    for(int i = 1; i <= n; ++ i) {
        if(!viz[i]) {
            added[i] = true;
            hasCycleDfs(i, 0, added);
            added[i] = false;
        }
        if(cycle)
            break;
    }

    return cycle;
}

void UnweightedGraph::criticalNodesDfs(int node, int parent, vector<int>& result, vector<int>& lvl, vector<int>& minLvl, int root) {
    viz[node] = true;
    lvl[node] = lvl[parent] + 1;
    minLvl[node] = lvl[node];
    bool critical = false;
    int noChildren = 0;   

    for(auto nextNode : adj[node]) {
        if(nextNode == parent)
            continue; 

        if(!viz[nextNode]) {
            ++ noChildren;
            criticalNodesDfs(nextNode, node, result, lvl, minLvl, root);  
            minLvl[node] = min(minLvl[node], minLvl[nextNode]);
            if(node != root && minLvl[nextNode] >= lvl[node])
                critical = true;
        }

        else
            minLvl[node] = min(minLvl[node], lvl[nextNode]);
    }

    if(node == root && noChildren > 1)
        critical = true;

    if(critical)
        result.push_back(node);
}

vector<int> UnweightedGraph::criticalNodes() {
    vector<int> result, lvl(n + 1, 0), minLvl(n + 1, 0);
    fill(viz.begin(), viz.end(), false);

    criticalNodesDfs(1, 0, result, lvl, minLvl, 1);

    return result;
}

void UnweightedGraph::criticalEdgesDfs(int node, int parent, vector<pair<int, int>>& result, vector<int>& lvl, vector<int>& minLvl) {
    viz[node] = true;
    lvl[node] = lvl[parent] + 1;
    minLvl[node] = lvl[node];
    int noChildren = 0;    

    for(auto nextNode : adj[node]) {
        if(nextNode == parent)
            continue; 

        if(!viz[nextNode]) {
            ++ noChildren;
            criticalEdgesDfs(nextNode, node, result, lvl, minLvl);
            minLvl[node] = min(minLvl[node], minLvl[nextNode]);
            if(minLvl[nextNode] > lvl[node])
                result.push_back({node, nextNode});
        }

        else
            minLvl[node] = min(minLvl[node], lvl[nextNode]);
    } 
}

vector<pair<int, int>> UnweightedGraph::criticalEdges() {
    vector<pair<int, int>> result;
    vector<int> lvl(n + 1, 0), minLvl(n + 1, 0);
    fill(viz.begin(), viz.end(), false);

    criticalEdgesDfs(1, 0, result, lvl, minLvl);

    return result;
}

void UnweightedGraph::getSCC(int node, vector<vector<int>>& adjt, vector<int>& component) {
    viz[node] = false;
    component.push_back(node);
    for(auto nextNode : adjt[node])
        if(viz[nextNode])
            getSCC(nextNode, adjt, component);
}

vector<vector<int>> UnweightedGraph::getSCComponents() {
    vector<vector<int>> components;
    vector<int> sorted = topSort();
    vector<vector<int>> adjt = transpose();
    vector<int> currComponent;

    for(auto node : sorted) 
        if(viz[node]) {
            getSCC(node, adjt, currComponent);
            components.emplace_back(currComponent);
            currComponent = vector<int>();
        }

    return components;
}

void UnweightedGraph::getBCC(int node, int parent, vector<vector<pair<int, int>>>& result, vector<int>& lvl, vector<int>& minLvl, stack<pair<int, int>>& st, int root) {
    viz[node] = true;
    lvl[node] = lvl[parent] + 1;
    minLvl[node] = lvl[node];
    bool critical = false;
    int noChildren = 0;    

    for(auto nextNode : adj[node]) {
        if(nextNode == parent)
            continue; 

        if(!viz[nextNode]) {
            st.push({node, nextNode});
            ++ noChildren;
            getBCC(nextNode, node, result, lvl, minLvl, st, root);
            minLvl[node] = min(minLvl[node], minLvl[nextNode]);
            critical = (node != root && minLvl[nextNode] >= lvl[node] && !st.empty()) || (node == root && noChildren > 1 && !st.empty());

            if(critical) {
                vector<pair<int, int>> partialResult;
                while(!st.empty()) {
                    auto topElem = st.top();
                    partialResult.push_back(topElem);
                    st.pop();
                    if(topElem == make_pair(node, nextNode))
                        break;
                }

                result.push_back(partialResult);
            }
        }
        
        else {  
            minLvl[node] = min(minLvl[node], lvl[nextNode]);
            if(lvl[nextNode] < lvl[node])   
                st.push({node, nextNode});
        }
           
    }
}

vector<vector<pair<int, int>>> UnweightedGraph::getBCComponents() {
    vector<vector<pair<int, int>>> result;
    stack<pair<int, int>> stackBCC;
    fill(viz.begin(), viz.end(), false);
    vector<int> lvl(n + 1, 0), minLvl(n + 1, 0);
    
    for(int i = 1; i <= n; ++ i) {
        if(!viz[i]) {
            stackBCC.push({0, i});
            getBCC(i, 0, result, lvl, minLvl, stackBCC, i);
        }
        
        if(!stackBCC.empty() && stackBCC.top() != make_pair(0, i)) {
            vector<pair<int, int>> partialResult;
            while(!stackBCC.empty() && stackBCC.top() != make_pair(0, i)) {
                partialResult.push_back(stackBCC.top());
                stackBCC.pop();
            }
            stackBCC.pop();
            result.push_back(partialResult);
        }
    }
        
    return result;
}

deque<int> UnweightedGraph::getEulerianUndirected(bool requiresCycle) {
    unordered_set<string> usedEdges;
    vector<int> deg(n + 1, 0);
    vector<vector<int>> adjTemp = adj;

    int odd = 0;
    int startNode = 1;
    for(int i = 1; i <= n; ++ i) {
        deg[i] = adj[i].size();
        if(deg[i] % 2 != 0) {
            startNode = i;
            ++ odd;
        }
    }

    if((odd && requiresCycle) || odd > 2)
        return deque<int>();
    
    deque<int> result;
    stack<int> Stack;
    Stack.push(startNode);
    int currNode = startNode, nextNode;

    while(!Stack.empty()){
        nextNode = adjTemp[currNode].back();
        while(adjTemp[currNode].size() && usedEdges.find(edgeToString(currNode, nextNode)) != usedEdges.end()) {
            adjTemp[currNode].pop_back();
            nextNode = adjTemp[currNode].back();
        }

        if(!adjTemp[currNode].size()) {
            result.push_back(currNode);
            currNode = Stack.top();
            Stack.pop();
        }

        else {
            Stack.push(currNode);                
            if(deg[nextNode] % 2 != 0 && (int)adjTemp[currNode].size() >= 2) {
                int lastPos = adjTemp[currNode].size() - 1;
                swap(adjTemp[currNode][lastPos], adjTemp[currNode][lastPos - 1]);
                nextNode = adjTemp[currNode].back();
            }

            adjTemp[currNode].pop_back();
            usedEdges.insert(edgeToString(currNode, nextNode));
            currNode = nextNode;
        }
    }

    return result;
}

deque<int> UnweightedGraph::getEulerianDirected(bool requiresCycle) {
    vector<int> degIn(n + 1, 0), degOut(n + 1, 0);
    vector<vector<int>> adjTemp = adj;

    for(int i = 1; i <= n; ++ i) {
        degOut[i] += adj[i].size();
        for(auto nextNode : adj[i])
            degIn[nextNode] ++;
    }

    int startNode = 1;
    int unbalanced = 0;
    for(int i = 1; i <= n; ++ i) 
        if(degIn[i] != degOut[i]) {
            unbalanced ++;
            if(degOut[i] > degIn[i])
                startNode = i;
        }

    if((unbalanced && requiresCycle) || unbalanced > 2)
        return deque<int>();
    
    stack<int> Stack;
    deque<int> result;
    Stack.push(startNode);
    int currNode = startNode, nextNode;
    while(!Stack.empty()) {
        if(!adjTemp[currNode].size()) {
            result.push_front(currNode);
            currNode = Stack.top();
            Stack.pop();
        }

        else {
            Stack.push(currNode);
            nextNode = adjTemp[currNode].back();
            if(degIn[nextNode] != degOut[nextNode] && (int) adjTemp[currNode].size() >= 2) {
                int lastPos = adjTemp[currNode].size() - 1;
                swap(adjTemp[currNode][lastPos], adjTemp[currNode][lastPos - 1]);
                nextNode = adjTemp[currNode].back();
            }
            adjTemp[currNode].pop_back();
            currNode = nextNode;
        }
    }

    return result;
}

deque<int> UnweightedGraph::getEulerianCycle() {
    if(Type == undirected)
        return getEulerianUndirected(true);
    return getEulerianDirected(true);
}

deque<int> UnweightedGraph::getEulerianPath() {
    if(Type == undirected)
        return getEulerianUndirected(false);
    return getEulerianDirected(false);
}

void UnweightedGraph::CCDFS(int node, vector<int>& component) {
    component.push_back(node);
    viz[node] = true;
    for(auto nextNode : adj[node])
        if(!viz[nextNode])
            CCDFS(nextNode, component);
}

vector<vector<int>> UnweightedGraph::getConnectedComponents() {
    fill(viz.begin(), viz.end(), false);
    vector<vector<int>> result;
    for(int i = 1; i <= n; ++ i)
        if(!viz[i]) {
            result.push_back(vector<int>());
            CCDFS(i, result.back());
        }
            
    return result;
}

void UnweightedGraph::addEdge(int x, int y) {
    if(Type == undirected) {
        adj[x].push_back(y);
        adj[y].push_back(x);
    }

    else
        adj[x].push_back(y);
    
    m ++;
}

pair<int, vector<int>> UnweightedGraph::bfs2(int src, int dest) {
    if(src == dest)
        return make_pair(0, vector<int>({src}));

    vector<int> prevNode(n + 1, 0);
    vector<int> dist(n + 1);
    fill(viz.begin(), viz.end(), false);
    queue<int> q;
    q.push(src);
    viz[src] = true;
    dist[src] = 0;

    while(!q.empty()) {
        int currNode = q.front();
        q.pop();

        if(currNode == dest)
            break;

        for(auto nextNode : adj[currNode]) {
            if(!viz[nextNode]) {
                viz[nextNode] = true;
                dist[nextNode] = dist[currNode] + 1;
                prevNode[nextNode] = currNode;
                q.push(nextNode);
            }
        }
    }

    if(!viz[dest])
        return make_pair(-1, vector<int>());

    vector<int> result;
    int node = dest;
    do {
        result.push_back(node);
        node = prevNode[node];
    }while(node);

    reverse(result.begin(), result.end());
    return make_pair(dist[dest], result);
}

string coordsToString(int i, int j) {
    return to_string(i) + "-" + to_string(j);
}

bool coordsOk(int i, int j, int lowi, int highi, int lowj, int highj) {
    return (lowi <= i && i <= highi && lowj <= j && j <= highj);
}

int main() {

    return 0;
}