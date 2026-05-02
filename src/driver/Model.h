#pragma once

#include "AttrToValueMap.h"
#include "Materials.h"
#include "Domain.h"
#include "mfem.hpp"

#include <array>
#include <memory>

namespace tulip {

struct Box {
	std::array<double,2> min, max;

	double area() const 
	{
		return (max[0] - min[0]) * (max[1] - min[1]);
	}

	bool isWithinBox(const mfem::Vector& point) const 
	{
		return (point[0] >= min[0] && point[0] <= max[0]) &&
			(point[1] >= min[1] && point[1] <= max[1]);
	}

	void displace(const std::array<double,2>& point) {
		for (int i = 0; i < 2; ++i) {
			min[i] += point[i];
			max[i] += point[i];
		}
	}

	bool operator==(const Box& rhs) const
	{
		return min == rhs.min && max == rhs.max;
	}
};

class Model {
public:
	enum class Openness {
		open,     // The most external boundary is open.
		closed    // The most external boundary is a conductor.
	};

	Model() = default;
	Model(
		mfem::Mesh&& mesh,         // Model gets ownership of mesh.
		const Materials& materials // Stores only materials present in mesh.
	);

	mfem::Mesh* getMesh() { return mesh_.get(); }
	const mfem::Mesh* getMesh() const { return mesh_.get(); }

	const Materials& getMaterials() const { return materials_; }
	std::size_t numberOfConductors() const;
	
	void setGroundConductorId(ConductorId id) { groundConductorId_ = id; }
	ConductorId getGroundConductorId() const { return groundConductorId_; }

	const Domain::IdToDomain& getDomains() const { return domains_; }

	Openness determineOpenness() const;
	
	double getAreaOfMaterial(const Material* m) const;
	Box getInnerRegionBoundingBox() const;
	
    bool isOuterRegion(const Material*) const;

private:
	Materials materials_;
	std::unique_ptr<mfem::Mesh> mesh_;
	Domain::IdToDomain domains_;
	ConductorId groundConductorId_{ 0 };
	
	Domain::IdToDomain buildDomains() const;
	DirectedGraph buildMeshGraph() const;
};

}