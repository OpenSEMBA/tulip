#pragma once

#include <string>

static std::string testDataPath() { return "./testData/"; }

static std::string inputFileFromCaseName(const std::string& caseName) {
    return testDataPath() + caseName + "/" + caseName + ".tulip.input.json";
}
