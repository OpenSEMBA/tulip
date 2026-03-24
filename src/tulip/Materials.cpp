#include "Materials.h"
#include "constants.h"

namespace pulmtln {


NameToAttrMap Materials::buildNameToAttrMap() const
{
	NameToAttrMap res{ buildNameToAttrMapFor<PEC>() };
	{
		auto aux{ buildNameToAttrMapFor<OpenBoundary>() };
		res.insert(aux.begin(), aux.end());
	}
	{
		auto aux{ buildNameToAttrMapFor<Dielectric>() };
		res.insert(aux.begin(), aux.end());
	}
	return res;
}

void Materials::removeMaterialsNotInList(const NameToAttrMap allowedMaterials)
{

	auto condition = [&allowedMaterials](const Material& mat) {
		return allowedMaterials.count(mat.name) > 0;
		};

	{
		std::vector<PEC> vs;
		std::copy_if(pecs.begin(), pecs.end(), std::back_inserter(vs), condition);
		pecs = vs;
	}
	{
		std::vector<Dielectric> vs;
		std::copy_if(dielectrics.begin(), dielectrics.end(), std::back_inserter(vs), condition);
		dielectrics = vs;
	}
	{
		std::vector<OpenBoundary> vs;
		std::copy_if(openBoundaries.begin(), openBoundaries.end(), std::back_inserter(vs), condition);
		openBoundaries = vs;
	}
}

MaterialId Materials::getMaterialIdFromName(const std::string& name)
{
	std::stringstream ss{ name.substr(name.find("_") + 1) };
	int res;
	ss >> res;
	return res;
}

bool Materials::isDomainMaterial(const std::string& name) const
{
	auto dielectrics = buildNameToAttrMapFor<Dielectric>();
	if (dielectrics.count(name)) {
		return true;
	}
	else {
		return false;
	}

}

bool Materials::hasDielectrics() const
{
	for (const auto& d : dielectrics) {
		if (d.relativePermittivity != VACUUM_RELATIVE_PERMITTIVITY) {
			return true;
		}
	}
	return false;
}

}