#pragma once

#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/function.hpp>

#include <iostream>
#include <vector>
#include <set>
#include <map>

namespace pulmtln {

using namespace boost;

using VertexId = int;
using Path = std::vector<VertexId>;
using IdSet = std::set<VertexId>;

struct Vertex_t { 
    VertexId id; 
    bool operator==(const Vertex_t& rhs) const { return id == rhs.id; }
    bool operator<(const Vertex_t& rhs) const { return id < rhs.id; }
};

typedef adjacency_list<vecS, vecS, bidirectionalS, Vertex_t, property<edge_index_t, int>> graph_t;

typedef filtered_graph<
    graph_t, 
    function<bool(graph_t::edge_descriptor)>, 
    function<bool(graph_t::vertex_descriptor)> > ComponentGraph;

typedef property_map<graph_t, vertex_index_t>::type IndexMap;


typedef graph_traits<graph_t>::vertex_iterator vertex_iter;
/// @brief Implements an oriented graph class based on <a href="https://www.boost.org/doc/libs/1_78_0/libs/graph/doc/index.html">the BOOST Graph Library</a> 
class DirectedGraph {
public:
    typedef std::pair<VertexId, VertexId> EdgeIds;
    typedef std::map<VertexId, graph_t::vertex_descriptor> VertexMap;
    typedef std::vector<VertexId> Path;
    typedef std::vector<Path> Paths;

    DirectedGraph() = default;
    
    /// @brief Add vertex identified by VertexId id to the graph.
    void addVertex(const VertexId& id);
    /// @brief Add oriented edge from vertex id1 to vertex id2 to the graph. 
    void addEdge(const VertexId& id1, const VertexId& id2);

    void addClosedPath(const std::vector<VertexId>& vIds);
    
    /// @brief Returns the number of vertices in the graph. 
    std::size_t verticesSize() const { return num_vertices(graph_); }
    /// @brief return the number of edges in the graph.
    std::size_t edgesSize() const { return num_edges(graph_); }

    /// @brief Returns the VertexId corresponding to the graph vertices. 
    IdSet getVertices() const;
    
    /// @brief Returns all the Path in the graph that form a cycle.
    Paths findCycles() const;

    /// @brief Return a graph for each disconnected subgraph within the graph. 
    std::vector<DirectedGraph> split() const;

    /// @brief Returns the boundary of the graph, formed by all the external edges. 
    DirectedGraph getBoundaryGraph() const;

    /// @brief Returns edges as a list of pairs.
    std::vector<std::pair<VertexId, VertexId>> getEdgesAsPairs() const;
    
private:
    VertexMap vertexMap_;
    graph_t graph_;
};


}
