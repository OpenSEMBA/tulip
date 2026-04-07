#include "Domain.h"

namespace tulip {
	
DomainTree::DomainTree(const Domain::IdToDomain& domains)
{
	// Check conductor 0 must be in a single domain.
	bool foundOnce{ false };
	for (const auto& [id, dom] : domains) {
		if (dom.conductorIds.count(0) && !foundOnce) {
			foundOnce = true;
			continue;
		}
		if (dom.conductorIds.count(0) && foundOnce) {
			throw std::runtime_error(
				"Conductor 0 must be only in a single domain."
			);
		}
	}
	if (!foundOnce) {
		throw std::runtime_error("Conductor 0 is not present");
	}
	
	std::multimap<ConductorId, Domain::Id> condIdToDomId;
	for (const auto& [domId, dom] : domains) {
		addVertex(domId);
		for (const auto& condId : dom.conductorIds) {
			condIdToDomId.emplace(condId, domId);
		}
	}

	int prevCondId{ -1 };
	int prevDomId{ -1 };
	for (const auto& [cId, dId] : condIdToDomId) {
		if (cId == prevCondId) {
			addEdge(prevDomId, dId);
		}
		else {
			prevCondId = cId;
			prevDomId = dId;
			addVertex(dId);
		}
	}

	// Post-conditions
	assert(findCycles().size() == 0);
}

}