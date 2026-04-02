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
			res.addConductor(c->getAttribute(), c->getConductorId(), c->isGround());
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
}

Model::Openness Model::determineOpenness() const
{
	if (materials_.getOpenBoundaries().empty()) {
		return Openness::closed;
	} else {
		return Openness::open;
	}
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
	}

	double length = 0.0;
	for (int i = 0; i < mesh_->GetNBE(); ++i) {
		if (mesh_->GetBdrElement(i)->GetAttribute() == attr) {
			length += getBoundaryElementMeasure(*mesh_, i);
		}
	}
	return length;

}

Box Model::getInnerRegionBoundingBox() const
{
	Box res{
		{infinity(), infinity()},
		{-infinity(), -infinity()}
	};

	for (auto m: materials_.getAll()) {
		if (m->isOuterRegion()) {
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