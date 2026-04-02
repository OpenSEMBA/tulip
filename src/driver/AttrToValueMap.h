#pragma once

#include <map>
#include <string>
#include <vector>

namespace tulip {

using NameToAttrMap = std::map<std::string, int>;

class AttrToValueMap : public std::map<int,double> {
public:
	AttrToValueMap() = default;
	AttrToValueMap(const std::map<int,double>& attVals) :
		std::map<int,double>{attVals}
	{}

	AttrToValueMap(std::initializer_list<std::pair<const int, double>> init) :
		std::map<int, double>{init}
	{}

	AttrToValueMap& operator=(std::initializer_list<std::pair<const int, double>> init)
	{
		this->std::map<int, double>::operator=(init);
		return *this;
	}
	
	std::vector<int> getAttributesAsArray() const
	{
		int dbcSize{ (int) size() };
		std::vector<int> bi(dbcSize);
		auto it{ begin() };
		for (int i = 0; i < dbcSize; i++) {
			bi[i] = it->first;
			++it;
		}
		return bi;
	}

	std::vector<double> getValuesAsArray() const
	{
		int dbcSize{ (int) size() };
		std::vector<double> bv(dbcSize);
		auto it{ begin() };
		for (int i = 0; i < dbcSize; i++) {
			bv[i] = it->second;
			++it;
		}
		return bv;
	}
};

}