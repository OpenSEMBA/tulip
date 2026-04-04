#include "Materials.h"
#include "constants.h"

#include <algorithm>
#include <cmath>

namespace tulip {

void Materials::addConductor(Attribute attribute, ConductorId id, bool isGround)
{
	materials_.push_back(std::make_unique<Conductor>(attribute, id, isGround));
}

void Materials::addDielectric(Attribute attribute, double relativePermittivity)
{
	materials_.push_back(std::make_unique<Dielectric>(attribute, relativePermittivity));
}

void Materials::addOpenBoundary(Attribute attribute)
{
	materials_.push_back(std::make_unique<Open>(attribute));
}

bool Materials::hasDielectrics() const
{
	for (const auto* d : getDielectrics()) {
		if (std::abs(d->getRelativePermittivity() - VACUUM_RELATIVE_PERMITTIVITY) > 1e-12) {
			return true;
		}
	}
	return false;
}

std::list<const Conductor*> Materials::getConductors() const
{
	std::list<const Conductor*> res;
	for (const auto& m : materials_) {
		if (const auto* c = dynamic_cast<const Conductor*>(m.get())) {
			res.push_back(c);
		}
	}

	res.sort([](const Conductor* a, const Conductor* b) {
		return a->getConductorId() < b->getConductorId();
	});

	return res;
}

const Conductor* Materials::getConductorWithId(ConductorId cId) const {
	for (auto c : getConductors()) {
		if (c->getConductorId() == cId) {
			return c;
		}
	}
}

std::list<const Dielectric*> Materials::getDielectrics() const
{
	std::list<const Dielectric*> res;
	for (const auto& m : materials_) {
		if (const auto* d = dynamic_cast<const Dielectric*>(m.get())) {
			res.push_back(d);
		}
	}
	return res;
}

std::list<const Open*> Materials::getOpenBoundaries() const
{
	std::list<const Open*> res;
	for (const auto& m : materials_) {
		if (const auto* b = dynamic_cast<const Open*>(m.get())) {
			res.push_back(b);
		}
	}
	return res;
}

std::list<const Material*> Materials::getAll() const
{
	std::list<const Material*> res;
	for (const auto& m : materials_) {
		res.push_back(m.get());
	}
	return res;
}




}