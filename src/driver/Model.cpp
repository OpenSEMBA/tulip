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
			res.addConductor(c->getAttribute(), c->getConductorId());
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