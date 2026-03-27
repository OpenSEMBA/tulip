#pragma once

#include "DirectedGraph.h"
#include "Model.h"

namespace pulmtln {

struct Domain {
    using Id = int;
    using ElementIds = std::set<int>;
    using IdToDomain = std::map<Id, Domain>;

    MaterialId ground{ -1 };
    std::set<MaterialId> conductorIds;
    ElementIds elems;   
    ElementIds bdrElems;

    static IdToDomain buildDomains(const Model&);
    static Model buildModelForDomain(
        mfem::Mesh& globalMesh,
        const Materials& materials,
        const Domain& domain);
};


// DomainTree is a tree graph with a single root having
// a vertex for each domain.
// The root contains Materialid 0.
class DomainTree : private DirectedGraph {
public:
    DomainTree() = default;
    DomainTree(const Domain::IdToDomain&);

    using DirectedGraph::getEdgesAsPairs; 
    using DirectedGraph::verticesSize;
};

}