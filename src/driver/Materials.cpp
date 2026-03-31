#include "Materials.h"
#include "constants.h"

namespace tulip {


void Materials::removeMaterialsNotInList(const NameToAttrMap allowedMaterials)
{
	// TODO
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