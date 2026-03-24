#include <gtest/gtest.h>

#include "Graph.h"

using namespace step2gmsh;

class GraphTest : public ::testing::Test {
protected:
    Graph graph;
};

TEST_F(GraphTest, addNode) {
    graph.addNode("A");
    auto nodes = graph.nodes();
    EXPECT_NE(std::find(nodes.begin(), nodes.end(), "A"), nodes.end());

    graph.addNode("A");
    int count = (int)std::count(nodes.begin(), nodes.end(), "A");
    EXPECT_EQ(count, 1);
}

TEST_F(GraphTest, addEdge) {
    graph.addEdge("A", "B");
    auto edges = graph.edges();
    EXPECT_NE(std::find(edges.begin(), edges.end(), Graph::Edge("A", "B")), edges.end());

    auto nodes = graph.nodes();
    EXPECT_NE(std::find(nodes.begin(), nodes.end(), "A"), nodes.end());
    EXPECT_NE(std::find(nodes.begin(), nodes.end(), "B"), nodes.end());

    graph.addEdge("A", "B");
    int count = (int)std::count(edges.begin(), edges.end(), Graph::Edge("A", "B"));
    EXPECT_EQ(count, 1);
}

TEST_F(GraphTest, settersAndGetters) {
    std::vector<Graph::Node> nodes = {"X", "Y"};
    std::vector<Graph::Edge> edges = {{"X", "Y"}};
    graph.setNodes(nodes);
    graph.setEdges(edges);
    EXPECT_EQ(graph.nodes(), nodes);
    EXPECT_EQ(graph.edges(), edges);
}

TEST_F(GraphTest, getConnections) {
    graph.addEdge("A", "B");
    graph.addEdge("A", "C");
    graph.addNode("D");

    auto connections = graph.getConnections();
    EXPECT_EQ(connections.at("A"), (std::vector<Graph::Node>{"B", "C"}));
    EXPECT_TRUE(connections.at("B").empty());
    EXPECT_TRUE(connections.at("C").empty());
    EXPECT_TRUE(connections.at("D").empty());
}

TEST_F(GraphTest, toString) {
    graph.addEdge("A", "B");
    std::string s = graph.toString();
    EXPECT_NE(s.find("A"), std::string::npos);
    EXPECT_NE(s.find("B"), std::string::npos);
    EXPECT_NE(s.find("Edges"), std::string::npos);
}

TEST_F(GraphTest, pruneToLongestPaths) {
    graph.setNodes({"A", "B", "C", "D", "E", "F", "G"});
    graph.setEdges({
        {"A", "B"}, {"A", "C"}, {"A", "D"}, {"A", "E"},
        {"B", "C"}, {"B", "E"},
        {"C", "E"},
        {"F", "G"}
    });

    std::vector<Graph::Edge> expectedEdges = {
        {"A", "B"}, {"A", "D"},
        {"B", "C"},
        {"C", "E"},
        {"F", "G"}
    };

    graph.pruneToLongestPaths();

    auto resultEdges = graph.edges();
    std::sort(resultEdges.begin(), resultEdges.end());
    std::sort(expectedEdges.begin(), expectedEdges.end());
    EXPECT_EQ(resultEdges, expectedEdges);
}

TEST_F(GraphTest, getNodesByLevel) {
    graph.setNodes({"A", "B", "C", "D", "E", "F", "G"});
    graph.setEdges({
        {"A", "B"}, {"A", "D"},
        {"B", "C"},
        {"C", "E"},
        {"F", "G"}
    });

    std::vector<Graph::Node> expected = {"A", "F", "B", "D", "G", "C", "E"};
    EXPECT_EQ(graph.getNodesByLevels(), expected);
}

TEST_F(GraphTest, getRoots) {
    graph.setNodes({"A", "B", "C", "D", "E", "F", "G"});
    graph.setEdges({
        {"A", "B"}, {"A", "D"},
        {"B", "C"},
        {"C", "E"},
        {"F", "G"}
    });
    EXPECT_EQ(graph.roots(), (std::vector<Graph::Node>{"A", "F"}));
}
