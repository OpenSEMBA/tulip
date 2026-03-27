#pragma once

#include <map>
#include "FES.h"

namespace pulmtln {

using NameToAttrMap = std::map<std::string, int>;

class AttrToValueMap : public std::map<int,double> {
public:
	AttrToValueMap() = default;
	AttrToValueMap(const std::map<int,double>& attVals) :
		std::map<int,double>{attVals}
	{}
	
	mfem::Array<int> getAttributesAsArray() const
	{
		int dbcSize{ (int) size() };
		mfem::Array<int> bi(dbcSize);
		auto it{ begin() };
		for (int i = 0; i < dbcSize; i++) {
			bi[i] = it->first;
			++it;
		}
		return bi;
	}

	mfem::Vector getValuesAsArray() const
	{
		int dbcSize{ (int) size() };
		mfem::Vector bv(dbcSize);
		auto it{ begin() };
		for (int i = 0; i < dbcSize; i++) {
			bv[i] = it->second;
			++it;
		}
		return bv;
	}
};

}