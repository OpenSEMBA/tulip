#include "Graph.h"

#include <algorithm>
#include <functional>
#include <map>
#include <queue>
#include <set>
#include <sstream>

namespace step2gmsh {

std::vector<Graph::Node> Graph::roots() const {
    std::vector<Node> r;
    for (const auto& node : nodes_) {
        bool isChild = false;
        for (const auto& edge : edges_) {
            if (edge.second == node) {
                isChild = true;
                break;
            }
        }
        if (!isChild) {
            r.push_back(node);
        }
    }
    return r;
}

const std::vector<Graph::Node>& Graph::nodes() const { return nodes_; }
const std::vector<Graph::Edge>& Graph::edges() const { return edges_; }

void Graph::setNodes(std::vector<Node> nodes) { nodes_ = std::move(nodes); }
void Graph::setEdges(std::vector<Edge> edges) { edges_ = std::move(edges); }

void Graph::addNode(const Node& node) {
    if (std::find(nodes_.begin(), nodes_.end(), node) == nodes_.end()) {
        nodes_.push_back(node);
    }
}

void Graph::addEdge(const Node& source, const Node& dest) {
    addNode(source);
    addNode(dest);
    Edge e{source, dest};
    if (std::find(edges_.begin(), edges_.end(), e) == edges_.end()) {
        edges_.push_back(e);
    }
}

std::map<Graph::Node, std::vector<Graph::Node>> Graph::getConnections() const {
    std::map<Node, std::vector<Node>> connections;
    for (const auto& node : nodes_) {
        connections[node] = {};
    }
    for (const auto& edge : edges_) {
        connections[edge.first].push_back(edge.second);
    }
    return connections;
}

std::vector<Graph::Node> Graph::getParentNodes() const {
    std::vector<Node> parents;
    for (const auto& edge : edges_) {
        parents.push_back(edge.first);
    }
    return parents;
}

std::vector<Graph::Node> Graph::getChildNodes() const {
    std::vector<Node> children;
    for (const auto& edge : edges_) {
        children.push_back(edge.second);
    }
    return children;
}

void Graph::pruneToLongestPaths() {
    auto connections = getConnections();
    auto childNodes = getChildNodes();

    std::vector<Node> rootNodes;
    for (const auto& node : nodes_) {
        if (std::find(childNodes.begin(), childNodes.end(), node) == childNodes.end()) {
            rootNodes.push_back(node);
        }
    }

    std::vector<std::vector<Node>> allPaths;
    std::function<void(const Node&, std::vector<Node>)> dfs =
        [&](const Node& node, std::vector<Node> path) {
            path.push_back(node);
            auto it = connections.find(node);
            if (it == connections.end() || it->second.empty()) {
                allPaths.push_back(path);
                return;
            }
            for (const auto& child : it->second) {
                dfs(child, path);
            }
        };

    for (const auto& root : rootNodes) {
        dfs(root, {});
    }

    std::map<Node, std::vector<Node>> leafToPath;
    for (const auto& path : allPaths) {
        const Node& leaf = path.back();
        auto it = leafToPath.find(leaf);
        if (it == leafToPath.end() || path.size() > it->second.size()) {
            leafToPath[leaf] = path;
        }
    }

    std::set<Node> newNodes;
    std::set<Edge> newEdges;
    for (const auto& [leaf, path] : leafToPath) {
        for (const auto& n : path) {
            newNodes.insert(n);
        }
        for (std::size_t i = 0; i + 1 < path.size(); ++i) {
            newEdges.insert({path[i], path[i + 1]});
        }
    }

    nodes_ = std::vector<Node>(newNodes.begin(), newNodes.end());
    edges_ = std::vector<Edge>(newEdges.begin(), newEdges.end());
}

std::map<Graph::Node, std::vector<Graph::Node>> Graph::getAdjacencyTree() const {
    std::map<Node, std::vector<Node>> tree;
    for (const auto& root : roots()) {
        tree[""].push_back(root);
    }
    for (const auto& edge : edges_) {
        tree[edge.first].push_back(edge.second);
    }
    return tree;
}

std::vector<Graph::Node> Graph::getNodesByLevels() const {
    auto tree = getAdjacencyTree();
    std::queue<std::pair<Node, int>> q;
    q.push({"", 0});
    std::vector<Node> nodeList;
    while (!q.empty()) {
        auto [node, level] = q.front();
        q.pop();
        nodeList.push_back(node);
        auto it = tree.find(node);
        if (it != tree.end()) {
            for (const auto& child : it->second) {
                q.push({child, level + 1});
            }
        }
    }
    // Remove the first element (empty-string sentinel root)
    return std::vector<Node>(nodeList.begin() + 1, nodeList.end());
}

void Graph::reorderData() {
    std::sort(edges_.begin(), edges_.end());
    std::sort(nodes_.begin(), nodes_.end());
}

std::string Graph::toString() const {
    std::ostringstream oss;
    oss << "Graph(Nodes: [";
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        if (i > 0) oss << ", ";
        oss << nodes_[i];
    }
    oss << "],\n Edges: [";
    for (std::size_t i = 0; i < edges_.size(); ++i) {
        if (i > 0) oss << ", ";
        oss << "(" << edges_[i].first << ", " << edges_[i].second << ")";
    }
    oss << "])";
    return oss.str();
}

} // namespace step2gmsh
