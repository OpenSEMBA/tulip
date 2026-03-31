#pragma once

#include "AttrToValueMap.h"

#include <string>
#include <list>
#include <memory>

namespace tulip {

using LayerId = int;
using Attribute = int;
using IdToAttrMap = std::map<LayerId, Attribute>;

// This are materials used by the solver. 
// There are only three options: 
// - Conductor (treated as PEC), 
// - Dielectric.
// - Open.
class Material {
public:
	virtual bool isDomainMaterial() const = 0;

private:
	LayerId id = -1;
	Attribute attribute = -1;
};

class Conductor : public Material {
public:	
	bool isDomainMaterial() const { return false; }
};

class Dielectric : public Material {
public:
	bool isDomainMaterial() const { return true; }

private:
	double relativePermittivity = 1.0;
};

class Open : public Material {
	bool isDomainMaterial() const { return false; }
};


class Materials {
public: 
	template <class T>
	NameToAttrMap buildNameToAttrMapFor() const
	{
		NameToAttrMap res;
		// TODO 
		return res;
	}

	template <class T>
	IdToAttrMap buildIdToAttrMapFor() const
	{
		IdToAttrMap res;
		// TODO
		return res;
	}

	template <class T>
	const T& getByLayerId(const LayerId id) const
	{
		// TODO
	}

	void removeMaterialsNotInList(const NameToAttrMap allowedMaterials);
	bool hasDielectrics() const;
	
private:
	std::list<std::unique_ptr<Material>> materials_;
};

}