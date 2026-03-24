#pragma once

#include <string>

static std::string testDataPath() { return "./testData/"; }

static std::string stepFileFromCaseName(const std::string& caseName) {
    return testDataPath() + caseName + "/" + caseName + ".step";
}
