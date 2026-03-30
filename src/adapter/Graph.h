#pragma once

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace tulip {

class Graph {
public:
    using Node = std::string;
    using Edge = std::pair<Node, Node>;

    Graph() = default;

    std::vector<Node> roots() const;

    const std::vector<Node>& nodes() const;
    const std::vector<Edge>& edges() const;

    void setNodes(std::vector<Node> nodes);
    void setEdges(std::vector<Edge> edges);

    void addNode(const Node& node);
    void addEdge(const Node& source, const Node& dest);

    std::map<Node, std::vector<Node>> getConnections() const;
    std::vector<Node> getParentNodes() const;
    std::vector<Node> getChildNodes() const;

    void pruneToLongestPaths();

    std::map<Node, std::vector<Node>> getAdjacencyTree() const;
    std::vector<Node> getNodesByLevels() const;

    void reorderData();

    std::string toString() const;

private:
    std::vector<Node> nodes_;
    std::vector<Edge> edges_;
};

} // namespace tulip
