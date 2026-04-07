#include "Domain.h"

namespace tulip {
	
DomainTree::DomainTree(const Domain::IdToDomain& domains)
{
	idToDomain = domains;

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

std::set<ConductorId> DomainTree::getConductorsInDomain(Domain::Id id) const
{
	auto it = idToDomain.find(id);
	if (it == idToDomain.end()) {
		throw std::runtime_error("Domain id not found.");
	}
	return it->second.conductorIds;
}

std::set<ConductorId> DomainTree::getConductorsInsideConductor(ConductorId conductorId) const
{
	std::set<ConductorId> res{ conductorId };

	std::map<Domain::Id, std::vector<Domain::Id>> adjacency;
	for (const auto& edge : getEdgesAsPairs()) {
		adjacency[edge.first].push_back(edge.second);
	}

	auto collectConductorsInSubTree = [&](Domain::Id rootId) {
		std::vector<Domain::Id> toVisit{ rootId };
		std::set<Domain::Id> visited;
		while (!toVisit.empty()) {
			auto id = toVisit.back();
			toVisit.pop_back();
			if (!visited.insert(id).second) {
				continue;
			}

			auto it = idToDomain.find(id);
			if (it == idToDomain.end()) {
				continue;
			}
			res.insert(it->second.conductorIds.begin(), it->second.conductorIds.end());

			auto aIt = adjacency.find(id);
			if (aIt != adjacency.end()) {
				toVisit.insert(toVisit.end(), aIt->second.begin(), aIt->second.end());
			}
		}
	};

	bool foundGroundDomain{ false };
	for (const auto& [domId, dom] : idToDomain) {
		if (dom.ground == conductorId) {
			foundGroundDomain = true;
			collectConductorsInSubTree(domId);
		}
	}

	if (foundGroundDomain) {
		return res;
	}

	for (const auto& [domId, dom] : idToDomain) {
		if (dom.ground != Domain::UNDEFINED_GROUND) {
			continue;
		}
		if (!dom.conductorIds.count(conductorId)) {
			continue;
		}
		auto aIt = adjacency.find(domId);
		if (aIt == adjacency.end() || aIt->second.empty()) {
			continue;
		}
		collectConductorsInSubTree(domId);
		break;
	}

	return res;
}

}