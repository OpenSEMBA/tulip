#include "Model.h"

#include <set>
#include <assert.h>

#include "DirectedGraph.h"

namespace pulmtln {

using namespace mfem;

std::vector<const Element*> getElementsWithAttribute(const Mesh& mesh, int attr)
{
	std::vector<const Element*> res;
	for (auto e{ 0 }; e < mesh.GetNBE(); ++e) {
		const auto elemPtr{ mesh.GetBdrElement(e) };
		if (elemPtr->GetAttribute() == attr) {
			res.push_back(elemPtr);
		}
	}
	return res;
}

std::multimap<int, const Element*> determineClosedLoops(const std::vector<const Element*>& elems)
{
	DirectedGraph g;
	std::map<int, const Element*> vToE;
	for (const auto& e : elems) {
		assert(e->GetType() == Element::Type::SEGMENT);
		int v0{ e->GetVertices()[0] };
		int v1{ e->GetVertices()[1] };
		g.addEdge(v0, v1);
		vToE[v0] = e;
	}

	std::multimap<int, const Element*> res;
	int cycleCount{ 0 };
	for (const auto& cycle : g.findCycles()) {
		for (const auto& v : cycle) {
			res.emplace(cycleCount, vToE.at(v));
		}
		cycleCount++;
	}

	return res;
}

Materials filterOutMaterialsNotPresentInMesh(
	const Materials& materials, 
	const Mesh& mesh)
{
	NameToAttrMap matInMesh;
	NameToAttrMap allMats{ materials.buildNameToAttrMap() };
	for (const auto& [name, attr] : allMats) {
		for (auto e{ 0 }; e < mesh.GetNBE(); e++) {
			if (mesh.GetBdrElement(e)->GetAttribute() == attr) {
				matInMesh.emplace(name, attr);
				break;
			}
		}
	}
	for (const auto& [name, attr] : allMats) {
		for (auto e{ 0 }; e < mesh.GetNE(); e++) {
			if (mesh.GetElement(e)->GetAttribute() == attr) {
				matInMesh.emplace(name, attr);
				break;
			}
		}
	}

	Materials res{ materials };
	res.removeMaterialsNotInList(matInMesh);
	return res;
}

Model::Model(
	Mesh&& mesh,  
	const Materials& materials) :
	mesh_{ std::make_unique<mfem::Mesh>(std::move(mesh)) }
{
	materials_ = filterOutMaterialsNotPresentInMesh(materials, *mesh_);
}

bool elementsFormOpenLoops(const std::vector<const Element*>& elems)
{
	return determineClosedLoops(elems).size() != elems.size();
}

std::size_t Model::numberOfConductors() const
{
	return materials_.buildNameToAttrMapFor<PEC>().size();
}

Model::Openness Model::determineOpenness() const
{
	assert(materials_.openBoundaries.size() <= 1);

	if (materials_.openBoundaries.size() == 0) {
		return Openness::closed;
	}

	for (const auto& m : materials_.openBoundaries) {
		if (!elementsFormOpenLoops(getElementsWithAttribute(*mesh_, m.attribute))) {
			return Openness::open;
		}
	}
	return Openness::semiopen;
}

double Model::getAreaOfMaterial(const std::string& materialName) const
{
	auto materials{ getMaterials().buildNameToAttrMap() };
	if (!materials.count(materialName)) {
		throw std::runtime_error("Unable to determine inner region.");
	}
	int innerRegionTag{ materials.at(materialName) };

	if (getMaterials().isDomainMaterial(materialName)) {
		Mesh m{ *getMesh() };
		double area = 0.0;
		for (int i = 0; i < m.GetNE(); ++i) {
			if (innerRegionTag == m.GetAttribute(i)) {
				area += m.GetElementVolume(i);
			}
		}
		return area;
	}
	else {
		return getMaterials().get<PEC>(materialName).area;
	}

}

Box Model::getBoundingBoxOfMaterial(const std::string& materialName) const
{
	Box res{
		{infinity(), infinity()},
		{-infinity(), -infinity()}
	};

	auto materials{ getMaterials().buildNameToAttrMap() };
	auto tag = materials.at(materialName);
	
	for (int e = 0; e < mesh_->GetNE(); ++e) {
		auto el = mesh_->GetElement(e);
		if (el->GetAttribute() != tag) {
			continue;
		}
		for (int v = 0; v < el->GetNVertices(); ++v) {
			auto vId = el->GetVertices()[v];
			for (int x = 0; x < 2; ++x) {
				auto vx = mesh_->GetVertex(vId)[x];
				res.min[x] = std::min(res.min[x], vx);
				res.max[x] = std::max(res.max[x], vx);
			}
		}
	}

	return res;
}


}