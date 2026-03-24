#include "DirectedGraph.h"

#include <algorithm>
#include <stdexcept>
#include <assert.h>
#include <set>

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/tiernan_all_cycles.hpp>

namespace boost { template <typename G> inline void renumber_vertex_indices(G& g) { throw; } }

namespace pulmtln {

IdSet DirectedGraph::getVertices() const 
{
    IdSet res;
    auto verts = vertices(graph_);
    for (auto& vp = verts; vp.first != vp.second; ++vp.first) {
        res.insert(graph_[*vp.first].id);
    }
    return res;
}

void DirectedGraph::addVertex(const VertexId& id) 
{
    if (vertexMap_.find(id) == vertexMap_.end()) {
        vertexMap_.emplace(id, add_vertex(Vertex_t({ id }), graph_));
    }
}

void DirectedGraph::addEdge(const VertexId& id1, const VertexId& id2) 
{
    if (id1 == id2) {
        throw std::runtime_error("Edges starting and finishing in same vertex are not allowed.");
    }
    
    if (vertexMap_.count(id1) == 0) {
        addVertex(id1);
    }
    if (vertexMap_.count(id2) == 0) {
        addVertex(id2);
    }
    
    add_edge(vertexMap_[id1], vertexMap_[id2], graph_);
}

void DirectedGraph::addClosedPath(const std::vector<VertexId>& vIds)
{
    if (vIds.size() == 1) {
        throw std::runtime_error("Paths must contain at least two vertices.");
    }
    
    for (auto i{0}; i < vIds.size(); ++i) {
        auto iNext{ (i + 1) % vIds.size() };
        addEdge(vIds[i], vIds[iNext]);
    }
}

std::vector<DirectedGraph> DirectedGraph::split() const
{
    if (num_vertices(graph_) == 0) {
        return {};
    }

    std::vector<std::size_t> component(num_vertices(graph_));
    std::size_t num = connected_components(graph_, &component[0]);
    if (num == 1) {
        return { *this };
    }

    IndexMap index = get(vertex_index, graph_);
    std::vector<DirectedGraph> res(num);
    for (std::size_t c = 0; c < num; c++) {
        graph_t::edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = edges(graph_); ei != ei_end; ++ei) {
            VertexId id1 = graph_[index[source(*ei, graph_)]].id;
            VertexId id2 = graph_[index[target(*ei, graph_)]].id;
            std::size_t v1 = vertexMap_.find(id1)->second;
            std::size_t v2 = vertexMap_.find(id2)->second;
            if (component[v1] == c && component[v2] == c) {
                res[c].addEdge(id1, id2);
            }
        }
    }
    return res;
}

DirectedGraph DirectedGraph::getBoundaryGraph() const
{
    DirectedGraph res;
    graph_traits<graph_t>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(graph_); ei != ei_end; ++ei) {
        auto src = source(*ei, graph_);
        auto tgt = target(*ei, graph_);
        if (!edge(tgt, src, graph_).second) {
            res.addEdge(graph_[src].id, graph_[tgt].id);
        }
    }
    return res;
}

struct cycle_recorder
{
    template <typename Path, typename graph_t>
    inline void cycle(const Path& p, const graph_t& g)
    { 
        cycles->push_back({});
        for (auto v : p) {
            cycles->back().push_back(g[v].id);
        }
    }
    std::vector<DirectedGraph::Path>* cycles;
};

std::vector<Path> DirectedGraph::findCycles() const
{
    std::vector<Path> res;
    cycle_recorder vis{&res};
    tiernan_all_cycles(graph_, vis);
    return res;
}


std::vector<std::pair<VertexId,VertexId>> DirectedGraph::getEdgesAsPairs() const
{
    std::vector<std::pair<VertexId, VertexId>> res;

    graph_traits<graph_t>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(graph_); ei != ei_end; ++ei) {
        auto src = source(*ei, graph_);
        auto tgt = target(*ei, graph_);

        auto srcId = graph_[src].id;
        auto tgtId = graph_[tgt].id;

        res.push_back(std::make_pair(srcId, tgtId));
    }

    return res;
}


}