#include "Model.h"

#include <cmath>
#include <set>
#include <assert.h>

#include "DirectedGraph.h"

namespace tulip {

using namespace mfem;

bool isAttributePresentInMesh(const Mesh& mesh, int attr)
{
	for (int e = 0; e < mesh.GetNBE(); ++e) {
		if (mesh.GetBdrElement(e)->GetAttribute() == attr) {
			return true;
		}
	}
	for (int e = 0; e < mesh.GetNE(); ++e) {
		if (mesh.GetElement(e)->GetAttribute() == attr) {
			return true;
		}
	}
	return false;
}


Materials filterOutMaterialsNotPresentInMesh(
	const Materials& materials, 
	const Mesh& mesh)
{
	Materials res;
	for (const auto* m : materials.getAll()) {
		if (!isAttributePresentInMesh(mesh, m->getAttribute())) {
			continue;
		}
		if (const auto* c = dynamic_cast<const Conductor*>(m)) {
			res.addConductor(
				c->getAttribute(),
				c->getConductorId(),
				c->getResistancePerMeter());
		} else if (const auto* d = dynamic_cast<const Dielectric*>(m)) {
			res.addDielectric(d->getAttribute(), d->getRelativePermittivity());
		} else if (dynamic_cast<const Open*>(m) != nullptr) {
			res.addOpenBoundary(m->getAttribute());
		}
	}
	return res;
}

Model::Model(
	Mesh&& mesh,  
	const Materials& materials) :
	materials_{},
	mesh_{ std::make_unique<mfem::Mesh>(std::move(mesh)) }
{
	materials_ = filterOutMaterialsNotPresentInMesh(materials, *mesh_);
	domains_ = buildDomains();
}

Model::Openness Model::determineOpenness() const
{
	if (materials_.getOpenBoundaries().empty()) {
		return Openness::closed;
	} else {
		return Openness::open;
	}
}

DirectedGraph Model::buildMeshGraph() const
{
	DirectedGraph meshGraph;
	for (auto e{ 0 }; e < mesh_->GetNE(); ++e) {
		const mfem::Element* elem{ mesh_->GetElement(e) };
		DirectedGraph::Path vs(elem->GetNVertices());
		for (auto i{ 0 }; i < vs.size(); ++i) {
			vs[i] = elem->GetVertices()[i];
		}
		meshGraph.addClosedPath(vs);
	}

	return meshGraph;
}

std::set<int> getBdrElemsInDomain(
	const Material* mat, 
	const IdSet& verticesInDomain, 
	const mfem::Mesh& mesh)
{
	std::set<int> res;
	
	for (auto e{ 0 }; e < mesh.GetNBE(); ++e) {
		const mfem::Element* elem{ mesh.GetBdrElement(e) };
		if (elem->GetAttribute() != mat->getAttribute()) {
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

Domain::IdToDomain Model::buildDomains() const
{
	Domain::IdToDomain res;
	Domain::Id id{ 0 };

	// Determine conductors in domain.
	for (const auto& domainMeshGraph : buildMeshGraph().split()) {
		const auto vsInDomain{ domainMeshGraph.getVertices() };
		Domain domain;
		for (const auto* pec : getMaterials().getConductors()) {
			auto bdrElems = getBdrElemsInDomain(pec, vsInDomain, *mesh_);
			if (!bdrElems.empty()) {
				domain.conductorIds.insert(pec->getConductorId());
			}
		}
		res[id++] = domain;
	}

	// Sets grounds. root domain is left with the default value.
	for (const auto& edge : DomainTree{ res }.getEdgesAsPairs()) {
		const auto& c1{ res[edge.first].conductorIds };
		const auto& c2{ res[edge.second].conductorIds };
		std::set<ConductorId> common;
		std::set_intersection(
			c1.begin(), c1.end(),
			c2.begin(), c2.end(),
			std::inserter(common, common.begin())
		);
		res[edge.second].ground = *common.begin();
	}

	return res;
}

double getBoundaryElementMeasure(const Mesh& mesh, int bdrElemId)
{
	const auto* elem = mesh.GetBdrElement(bdrElemId);
	if (elem->GetNVertices() != 2) {
		return 0.0;
	}
	const int v0 = elem->GetVertices()[0];
	const int v1 = elem->GetVertices()[1];
	const auto* p0 = mesh.GetVertex(v0);
	const auto* p1 = mesh.GetVertex(v1);

	double sumSq = 0.0;
	for (int d = 0; d < mesh.SpaceDimension(); ++d) {
		double dx = p1[d] - p0[d];
		sumSq += dx * dx;
	}
	return std::sqrt(sumSq);
}

double Model::getAreaOfMaterial(const Material* m) const
{
	const int attr = m->getAttribute();

	if (m->isDomainMaterial()) {
		double area = 0.0;
		for (int i = 0; i < mesh_->GetNE(); ++i) {
			if (mesh_->GetElement(i)->GetAttribute() == attr) {
				area += mesh_->GetElementVolume(i);
			}
		}
		return area;
	} else {
		MFEM_VERIFY(mesh_->Dimension() == 2,
			"Boundary-enclosed area is implemented for 2D meshes only.");

		const int maxAttrInMesh =
			mesh_->bdr_attributes.Size() > 0 ? mesh_->bdr_attributes.Max() : 0;
		const int maxAttr = std::max(attr, maxAttrInMesh);
		Array<int> bdrMarker(maxAttr);
		bdrMarker = 0;
		if (attr > 0) {
			bdrMarker[attr - 1] = 1;
		}

		H1_FECollection fec(1, mesh_->Dimension());
		FiniteElementSpace fes(mesh_.get(), &fec);

		auto g_fun = [](const Vector& x, Vector& g) {
			g.SetSize(2);
			g = 0.0;
			g(0) = x(0);
		};
		VectorFunctionCoefficient gcoef(2, g_fun);

		LinearForm lf(&fes);
		lf.AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(gcoef), bdrMarker);
		lf.Assemble();

		GridFunction one(&fes);
		one = 1.0;

		const double signedArea = lf * one;
		return std::abs(signedArea);
	}
}

bool Model::isOuterRegion(const Material* material) const
{
	if (dynamic_cast<const Open*>(material) != nullptr) {
		return true;
	}

	const auto* dielectric = dynamic_cast<const Dielectric*>(material);
	if (dielectric == nullptr) {
		return false;
	}

	std::set<int> openBoundaryAttributes;
	for (const auto* openBoundary : materials_.getOpenBoundaries()) {
		openBoundaryAttributes.insert(openBoundary->getAttribute());
	}
	if (openBoundaryAttributes.empty()) {
		return false;
	}

	std::set<int> verticesOnOpenBoundary;
	for (int b = 0; b < mesh_->GetNBE(); ++b) {
		const auto* bdrElem = mesh_->GetBdrElement(b);
		if (!openBoundaryAttributes.count(bdrElem->GetAttribute())) {
			continue;
		}
		for (int v = 0; v < bdrElem->GetNVertices(); ++v) {
			verticesOnOpenBoundary.insert(bdrElem->GetVertices()[v]);
		}
	}

	for (int e = 0; e < mesh_->GetNE(); ++e) {
		const auto* elem = mesh_->GetElement(e);
		if (elem->GetAttribute() != dielectric->getAttribute()) {
			continue;
		}
		for (int v = 0; v < elem->GetNVertices(); ++v) {
			if (verticesOnOpenBoundary.count(elem->GetVertices()[v])) {
				return true;
			}
		}
	}

	return false;
}

Box Model::getInnerRegionBoundingBox() const
{
	Box res{
		{infinity(), infinity()},
		{-infinity(), -infinity()}
	};

	for (auto m: materials_.getAll()) {
		if (isOuterRegion(m)) {
			continue;
		}
		
		auto tag = m->getAttribute();
		
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
		
	}
	
	return res;
}


}