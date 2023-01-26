#include <bits/stdc++.h>
using namespace std;

class FlowNetwork {
    static const int inf;
    int n, m, s, t, maxFlow = 0;
    struct edge {
        int src, dst, capacity, flow;
        edge* back;
    };
    vector<vector<edge*>> adj;
    vector<edge*> parentEdge;
    void appendEdge(int x, int y, int capacity) {
        edge* e1 = new edge;
        edge* e2 = new edge;

        e1->back = e2;
        e2->back = e1;
    
        e1->src = x;
        e1->dst = y;
        e1->capacity = capacity;
        e2->flow = e1->flow = 0;
        e2->src = y;
        e2->dst = x;
        e2->capacity = 0;

        adj[x].push_back(e1);
        adj[y].push_back(e2);
    }

    // computes the maxFlow of the FlowNetwork
    // once computed, the value of the max flow is stored in the maxFlow variable
    void computeMaxFlow();

public:
    // If you want bipartite matching, you must give the capacity of 1 to all edges
    // an edge from edges vector has the following form : (src, dst, capacity)
    // n is the number of nodes in the network, including source and terminal
    // m is the number of edges in the network
    FlowNetwork(int _n, int _m, int source, int term, vector<tuple<int, int, int>>& edges);
    void changeSource(int _s) { s = _s; }
    void changeTerminal(int _t) { t = _t; }

    // builds paths from src to term and checks if term is reachable from src
    bool bfs(int src, int term);

    // returns the maxFlow; if it's not computed, it computes it and then returns it
    int getMaxFlow();

    // returns the direct edges of the mincut
    vector<pair<int, int>> getMinCut();

    // returns the edges chosen in the maximum Matching algorithm
    vector<pair<int, int>> findMatching();
};

const int FlowNetwork::inf = INT_MAX;

FlowNetwork::FlowNetwork(int _n, int _m, int source, int term, vector<tuple<int, int, int>>& edges): n(_n), m(_m), s(source), t(term) {
    adj.resize(n + 1);
    parentEdge.resize(n + 1);
    maxFlow = 0;

    int x, y, c;
    for(auto edge : edges) {
        x = get<0>(edge);
        y = get<1>(edge);
        c = get<2>(edge);
        appendEdge(x, y, c);
    }
}

bool FlowNetwork::bfs(int src, int term) {
    fill(parentEdge.begin(), parentEdge.end(), nullptr);
    queue<int> q;
    q.push(src);

    int currNode;
    while(!q.empty()) {
        currNode = q.front();
        q.pop();        

        for(auto nextEdge : adj[currNode])  {
            int nextNode = nextEdge->dst;
            if(!parentEdge[nextNode] && nextEdge->capacity - nextEdge->flow > 0) {
                parentEdge[nextNode] = nextEdge;
                if(nextNode != term) {
                    q.push(nextNode);
                }
            }
        }
    }   

    return parentEdge[term] != nullptr;
}

void FlowNetwork::computeMaxFlow() {
    while(bfs(s, t)) {
        for(auto terminalEdge : adj[t]) {
            if(terminalEdge->back->capacity - terminalEdge->back->flow <= 0 || !parentEdge[terminalEdge->dst])        
                continue;
            
            parentEdge[t] = terminalEdge->back;
            int minValOnPath = inf;
            for(auto current = t; current != s; current = parentEdge[current]->src)
                minValOnPath = min(minValOnPath, parentEdge[current]->capacity - parentEdge[current]->flow);
            
            if(!minValOnPath)    
                continue;
            
            maxFlow += minValOnPath;
            
            for(auto current = t; current != s; current = parentEdge[current]->src) {
                parentEdge[current]->flow += minValOnPath;
                parentEdge[current]->back->flow -= minValOnPath;
            }
        }
    }
}

int FlowNetwork::getMaxFlow(){
    if(!maxFlow)
        computeMaxFlow();
    return maxFlow;
}

vector<pair<int, int>> FlowNetwork::getMinCut() {
    if(!maxFlow)
        computeMaxFlow();
    
    vector<bool> viz(n + 1);
    queue<int> q;
    q.push(s);
    viz[s] = true;

    vector<pair<int, int>> result;
    int currNode;
    while(!q.empty()){
        currNode = q.front();
        q.pop();

        for(auto nextEdge : adj[currNode]) {
            if(nextEdge->flow && nextEdge->capacity - nextEdge->flow == 0)
                result.push_back({nextEdge->src, nextEdge->dst});
            
            else if(!viz[nextEdge->dst] && nextEdge->capacity != 0) {
                q.push(nextEdge->dst);
                viz[nextEdge->dst] = true;
            }
        }
    }

    return result;
}

vector<pair<int, int>> FlowNetwork::findMatching() {
    if(!maxFlow)
        computeMaxFlow();

    vector<pair<int, int>> result;
    
    int nrLeft = (int) adj[s].size();

    for(auto sourceEdge : adj[s]) {
        int leftNode = sourceEdge->dst;
        for(auto edge : adj[leftNode]) {
            if(edge->dst == s)
                continue;
            if(edge->capacity - edge->flow == 0)
                result.push_back({leftNode, edge->dst});
        }
    }

    return result;
}

int main() {

    return 0;
}