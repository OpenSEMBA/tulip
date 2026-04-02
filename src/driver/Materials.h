#pragma once

#include "AttrToValueMap.h"

#include <string>
#include <list>
#include <memory>
#include <vector>

namespace tulip {

using ConductorId = int;
using Attribute = int;

// This are materials used by the solver. 
// There are only three options: 
// - Conductor (treated as PEC), 
// - Dielectric.
// - Open.
class Material {
public:
	virtual ~Material() = default;
	Attribute getAttribute() const { return attribute; }
	void setAsOuterRegion() { isOuterRegion_ = true; }
	bool isOuterRegion() const { return isOuterRegion_; }
	virtual bool isDomainMaterial() const = 0;

protected:
	explicit Material(Attribute attr) : attribute(attr) {}
	Material(const Material&) = default;
	Material& operator=(const Material&) = default;
	Material(Material&&) = default;
	Material& operator=(Material&&) = default;

private:
	Attribute attribute = -1;
	bool isOuterRegion_ = false;
};

class Conductor : public Material {
public:
	explicit Conductor(Attribute attr, ConductorId id, bool isGround = false)
		: Material(attr), conductorId(id), isGround_(isGround)
	{}

	ConductorId getConductorId() const {
		return conductorId;
	}
	bool isGround() const { return isGround_; }
	void setAsGround() { isGround_ = true; }
	bool isDomainMaterial() const { return false; }
private:	
	ConductorId conductorId = -1;
	bool isGround_ = false;
};

class Dielectric : public Material {
public:
	explicit Dielectric(Attribute attr, double epsr = 1.0)
		: Material(attr), relativePermittivity(epsr)
	{}

	double getRelativePermittivity() const {
		return relativePermittivity;
	}
	bool isDomainMaterial() const { return true; }
private:
	double relativePermittivity = 1.0;
};

class Open : public Material {
	public:
	explicit Open(Attribute attr) : Material(attr) {}
	bool isDomainMaterial() const { return false; }
};


class Materials {
public: 
	void addConductor(Attribute attribute, ConductorId id, bool isGround = false);
	void addDielectric(Attribute attribute, double relativePermittivity);
	void addOpenBoundary(Attribute attribute);

	void removeMaterialsNotInList(const NameToAttrMap allowedMaterials);
	
	// List of conductors sorted by their conductorId.
	std::list<const Conductor*> getConductors() const;

	std::list<const Dielectric*> getDielectrics() const;
	
	std::list<const Open*> getOpenBoundaries() const;

	std::list<const Material*> getAll() const;
	
	bool hasDielectrics() const;
	
private:
	std::list<std::unique_ptr<Material>> materials_;
};

}