#pragma once

#include "AttrToValueMap.h"

#include <string>
#include <vector>
#include <memory>

namespace tulip {

using MaterialId = int;
using Attribute = int;
using IdToAttrMap = std::map<MaterialId, Attribute>;

class Material {
public:
	virtual bool isDomainMaterial() const = 0;

private:
	std::string name;
	Attribute attribute = -1;
};

class Conductor : public Material {
public:	
	bool isDomainMaterial() const { return false; }

private:
	double resistancePerMeter = 0.0;
	
	// Area enclosed by the outermost curve. 
	// It is different from area if the layer contains holes.
	double enclosedArea = -1.0; 
};

class Shield : public Conductor {
public: 
	bool isDomainMaterial() const { return false; }

private:
	enum class Direction { 
		both, 
		inwards, 
		outwards
	};
	double inductancePerMeter = 0.0;
	Direction direction = Direction::both;

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

	Material&
	const T& getByName(const std::string& name) const
	{
		// TODO
	}

	template <class T>
	const T& getByMaterialId(const MaterialId id) const
	{
		// TODO
	}

	NameToAttrMap buildNameToAttrMap() const;
	void removeMaterialsNotInList(const NameToAttrMap allowedMaterials);
	bool hasDielectrics() const;
	
private:
	std::map<MaterialId, std::unique_ptr<Material>> materials_;
};

}