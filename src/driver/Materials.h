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
	Attribute getAttribute() const { return attribute; }
	void setAsOuterRegion() { isOuterRegion_ = true; }
	bool isOuterRegion() const { return isOuterRegion_; }
private:
	Attribute attribute = -1;
	bool isOuterRegion_ = false;
};

class Conductor : public Material {
public:
	ConductorId getConductorId() const {
		return conductorId;
	}
	bool isGround() const;
	void setAsGround() { isGround_ = true; }
private:	
	ConductorId conductorId = -1;
	bool isGround_ = false;
};

class Dielectric : public Material {
public:
	double getRelativePermittivity() const {
		return relativePermittivity;
	}
private:
	double relativePermittivity = 1.0;
};

class Open : public Material {
};


class Materials {
public: 
	
	// List of conductors sorted by their conductorId.
	std::list<const Conductor*> getConductors() const; // TODO

	std::list<const Dielectric*> getDielectrics() const; // TODO
	
	std::list<const Open*> getOpenBoundaries() const; // TODO

	std::list<const Material*> getAll() const; // TODO
	
	bool hasDielectrics() const;
	
private:
	std::list<std::unique_ptr<Material>> materials_;
};

}