#pragma once

#include "DirectedGraph.h"
#include "Materials.h"

namespace tulip {

struct Domain {
    using Id = int;
    using IdToDomain = std::map<Id, Domain>;

    static constexpr int UNDEFINED_GROUND = -1;

    ConductorId ground{ UNDEFINED_GROUND };
    std::set<ConductorId> conductorIds;
};


// DomainTree is a tree graph with a single root having
// a vertex for each domain.
class DomainTree : private DirectedGraph {
public:
    DomainTree() = default;
    DomainTree(const Domain::IdToDomain&);

    std::set<ConductorId> getConductorsInsideConductor(ConductorId) const;
    std::set<ConductorId> getConductorsInDomain(Domain::Id) const;

    using DirectedGraph::getEdgesAsPairs; 
    using DirectedGraph::verticesSize;

private:
    Domain::IdToDomain idToDomain;
};

}