#pragma once

#include "AttrToValueMap.h"

#include <string>
#include <list>
#include <memory>

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
	Attribute getAttribute() const {
		return attribute;
	}
private:
	Attribute attribute = -1;
};

class Conductor : public Material {
public:
	ConductorId getConductorId() const {
		return conductorId;
	}
private:	
	ConductorId conductorId = -1;
};

class Dielectric : public Material {
public:
	double getRelativePermittivity() const {
		return relativePermittivity;
	}
	void setAsOuterRegion() { isOuterRegion = true; }
	bool isOuterRegion() const { return isOuterRegion; }
private:
	double relativePermittivity = 1.0;
	bool isOuterRegion = false;
};

class Open : public Material {
};


class Materials {
public: 
	
	// List of conductors sorted by their conductorId.
	std::list<const Conductor*> getConductors() const;

	std::list<const Dielectric*> getDielectrics() const;
	
	std::list<const Open*> getOpenBoundaries() const;
	
	bool hasDielectrics() const;
	
private:
	std::list<std::unique_ptr<Material>> materials_;
};

}