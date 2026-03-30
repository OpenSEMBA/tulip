#pragma once

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>
#include <fstream>

using json = nlohmann::json;

static std::string testDataFolder(){ return "./testData/"; }
static std::string casesFolder()   { return testDataFolder(); }
static std::string outFolder()     { return "Results/tests/"; }


static std::string smbCase(const std::string& name)
{
	return casesFolder() + name + "/" + name + ".smb.json";
}

static std::string inputCase(const std::string& name)
{
	return casesFolder() + name + "/" + name + ".pulmtln.in.json";
}

static std::string getCaseName()
{
	return ::testing::UnitTest::GetInstance()->current_test_info()->name();
}

static std::string getTestCaseName() 
{
	std::string suiteName{
		::testing::UnitTest::GetInstance()->current_test_info()->test_suite_name()
	};	
	return suiteName + "." + getCaseName();
}

static json readJSON(const std::string& fn)
{
	std::ifstream stream(fn);
	json j;
	stream >> j;
	return j;
}

static double relError(double expectedVal, double val)
{
	if (expectedVal == 0.0) {
		throw std::runtime_error(
			"Unable to compute relative error of a 0.0 expected value");
	}
	return std::abs(val - expectedVal) / std::abs(expectedVal);
}