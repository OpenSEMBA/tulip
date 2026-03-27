#include "Domain.h"

using namespace mfem;

namespace pulmtln {

DirectedGraph buildMeshGraph(const mfem::Mesh& mesh)
{
	DirectedGraph meshGraph;
	for (auto e{ 0 }; e < mesh.GetNE(); ++e) {
		const mfem::Element* elem{ mesh.GetElement(e) };
		DirectedGraph::Path vs(elem->GetNVertices());
		for (auto i{ 0 }; i < vs.size(); ++i) {
			vs[i] = elem->GetVertices()[i];
		}
		meshGraph.addClosedPath(vs);
	}

	return meshGraph;
}

Domain::ElementIds getBdrElemsInDomain(
	const Material* mat, 
	const IdSet& verticesInDomain, 
	const mfem::Mesh& mesh)
{
	Domain::ElementIds res;
	
	for (auto e{ 0 }; e < mesh.GetNBE(); ++e) {
		const mfem::Element* elem{ mesh.GetBdrElement(e) };
		if (elem->GetAttribute() != mat->attribute) {
			continue;
		}
		IdSet verticesWithAttribute;
		verticesWithAttribute.insert(
			elem->GetVertices(),
			elem->GetVertices() + elem->GetNVertices());
		
		IdSet common;
		std::set_intersection(
			verticesInDomain.begin(), verticesInDomain.end(),
			verticesWithAttribute.begin(), verticesWithAttribute.end(),
			std::inserter(common, common.begin()) );

		if (common.size() < 2) {
			continue;
		}

		res.insert(e);
	}

	return res;
}
	
DomainTree::DomainTree(const Domain::IdToDomain& domains)
{
	// Check conductor 0 must be in a single domain.
	bool foundOnce{ false };
	for (const auto& [id, dom] : domains) {
		if (dom.conductorIds.count(0) && !foundOnce) {
			foundOnce = true;
			continue;
		}
		if (dom.conductorIds.count(0) && foundOnce) {
			throw std::runtime_error(
				"Conductor 0 must be only in a single domain."
			);
		}
	}
	if (!foundOnce) {
		throw std::runtime_error("Conductor 0 is not present");
	}
	
	std::multimap<MaterialId, Domain::Id> condIdToDomId;
	for (const auto& [domId, dom] : domains) {
		addVertex(domId);
		for (const auto& condId : dom.conductorIds) {
			condIdToDomId.emplace(condId, domId);
		}
	}

	int prevCondId{ -1 };
	int prevDomId{ -1 };
	for (const auto& [cId, dId] : condIdToDomId) {
		if (cId == prevCondId) {
			addEdge(prevDomId, dId);
		}
		else {
			prevCondId = cId;
			prevDomId = dId;
			addVertex(dId);
		}
	}

	// Post-conditions
	assert(findCycles().size() == 0);
}

Domain::IdToDomain Domain::buildDomains(const Model& model)
{
	if (model.determineOpenness() != Model::Openness::closed) {
		throw std::runtime_error("Domains can only be determined for closed problems.");
	}

	const mfem::Mesh& mesh{ *model.getMesh() };

	Domain::IdToDomain res;
	Domain::Id id{ 0 };
	for (const auto& domainMeshGraph : buildMeshGraph(mesh).split()) {

		const auto vsInDomain{ domainMeshGraph.getVertices() };

		Domain domain;
		// Determine which elements belong to domain.
		for (auto e{ 0 }; e < mesh.GetNE(); ++e) {
			const auto& elem{ *mesh.GetElement(e) };
			if (std::all_of(
				elem.GetVertices(), 
				elem.GetVertices() + elem.GetNVertices(),
				[&vsInDomain](int i) { return vsInDomain.count(i) == 1; }
			)) {
				domain.elems.insert(e);
			}
		}

		// Determine conductors in domain.
		for (const auto& pec : model.getMaterials().pecs) {
			auto bdrElems = getBdrElemsInDomain(&pec, vsInDomain, mesh);
			if (!bdrElems.empty()) {
				domain.conductorIds.insert(Materials::getMaterialIdFromName(pec.name));
				domain.bdrElems.insert(bdrElems.begin(), bdrElems.end());
			}
		}

		//
		res[id++] = domain;
	}

	// Sets grounds.
	for (auto& [domId, dom] : res) {
		if (dom.conductorIds.count(0) == 1) {
			dom.ground = 0;
		}
	}
	for (const auto& edge : DomainTree{ res }.getEdgesAsPairs()) {
		const auto& c1{ res[edge.first].conductorIds };
		const auto& c2{ res[edge.second].conductorIds };
		std::set<MaterialId> common;
		std::set_intersection(
			c1.begin(), c1.end(),
			c2.begin(), c2.end(),
			std::inserter(common, common.begin())
		);
		assert(common.size() == 1);
		res[edge.second].ground = *common.begin();
	}
	
	// Post conditions.
	// All elements belong to a single domain.
	// Conductor 0 is always ground and is in a single domain.

	return res;
}

Model Domain::buildModelForDomain(
	mfem::Mesh& globalMesh, 
	const Materials& materials, 
	const Domain& domain)
{
	// We must modify attributes in the global mesh to identify elements in subdomain.
	std::map<int, int> elemToAttBackup;
	for (auto e{ 0 }; e < globalMesh.GetNE(); ++e) {
		elemToAttBackup.emplace(e, globalMesh.GetElement(e)->GetAttribute());
	}

	// --
	for (auto e{ 0 }; e < globalMesh.GetNE(); ++e) {
		if (domain.elems.count(e)) {
			globalMesh.GetElement(e)->SetAttribute(1);
		}
		else {
			globalMesh.GetElement(e)->SetAttribute(2);
		}
	}
	globalMesh.SetAttributes();
	globalMesh.Finalize();
	// --
	Array<int> subdomainAttrs(1);
	subdomainAttrs[0] = 1;
	auto domainMesh{ SubMesh::CreateFromDomain(globalMesh, subdomainAttrs) };

	// Sets attributes in domain mesh.
	for (auto e{ 0 }; e < domainMesh.GetNE(); ++e) {
		const auto parentAttribute{
			elemToAttBackup.at(domainMesh.GetParentElementIDMap()[e])
		};
		domainMesh.GetElement(e)->SetAttribute(parentAttribute);
	}

	// Sets bdr attributes in domain mesh.
	auto parentFaceIdMap = SubMeshUtils::BuildFaceMap(
		globalMesh, domainMesh, domainMesh.GetParentElementIDMap());
	const auto faceToBdrElem{ globalMesh.GetFaceToBdrElMap() };
	for (auto b{ 0 }; b < domainMesh.GetNBE(); b++) {
		domainMesh.GetBdrElement(b)->SetAttribute(
			globalMesh.GetBdrAttribute(
				faceToBdrElem[
					parentFaceIdMap[
						domainMesh.GetBdrFace(b)]]
			)
		);
	}
	domainMesh.SetAttributes();
	domainMesh.Finalize();

	// Restores original attributes in global mesh.
	for (auto e{ 0 }; e < globalMesh.GetNE(); ++e) {
		globalMesh.GetElement(e)->SetAttribute(elemToAttBackup.at(e));
	}
	globalMesh.SetAttributes();
	globalMesh.Finalize();

	Model res{ std::move(domainMesh), materials };
	res.setGroundConductorId(domain.ground);

	return res;
}


}