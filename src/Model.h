#pragma once

#include "AttrToValueMap.h"
#include "Materials.h"

namespace pulmtln {

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

	bool operator==(const Box& rhs) const
	{
		return min == rhs.min && max == rhs.max;
	}
};

class Model {
public:
	enum class Openness {
		open,     // The most external boundary is open.
		semiopen, // The most external boundaries are open and a conductor.
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
	
	void setGroundConductorId(MaterialId id) { groundConductorId_ = id; }
	MaterialId getGroundConductorId() const { return groundConductorId_; }

	Openness determineOpenness() const;
	
	double getAreaOfMaterial(const std::string& materialName) const;
	Box getBoundingBoxOfMaterial(const std::string& materialName) const;

private:
	Materials materials_;
	std::unique_ptr<mfem::Mesh> mesh_;
	MaterialId groundConductorId_{ 0 };
};

}