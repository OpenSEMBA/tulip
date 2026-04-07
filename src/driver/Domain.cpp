#include "Domain.h"

namespace tulip {
	
DomainTree::DomainTree(const Domain::IdToDomain& domains)
{
	bool foundOnce{ false };
	for (const auto& [id, dom] : domains) {
		if (dom.conductorIds.count(Domain::UNDEFINED_GROUND) && !foundOnce) {
			foundOnce = true;
			continue;
		}
		if (dom.conductorIds.count(Domain::UNDEFINED_GROUND) && foundOnce) {
			throw std::runtime_error(
				"Only one domain can have undefined ground."
			);
		}
	}
	if (!foundOnce) {
		throw std::runtime_error("Root domain not found.");
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