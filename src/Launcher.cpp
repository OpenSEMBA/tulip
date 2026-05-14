#include "Launcher.h"

#include "Adapter.h"
#include "AdaptedInputParser.h"
#include "Driver.h"

#include <filesystem>
#include <gmsh.h>
#include <iostream>
#include <stdexcept>

namespace tulip {

namespace {

bool hasSuffix(const std::string& value, const std::string& suffix)
{
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool isInputJson(const std::string& filename)
{
    return hasSuffix(filename, ".tulip.input.json");
}

bool isAdaptedJson(const std::string& filename)
{
    return hasSuffix(filename, ".tulip.adapted.json");
}

std::string ensureTrailingSlash(const std::string& folder)
{
    if (folder.empty() || folder.back() == '/') {
        return folder;
    }
    return folder + "/";
}

std::string extractCaseName(const std::string& inputFile)
{
    const auto filename = std::filesystem::path(inputFile).filename().string();
    if (hasSuffix(filename, ".tulip.input.json")) {
        return filename.substr(0, filename.size() - std::string(".tulip.input.json").size());
    }
    if (hasSuffix(filename, ".tulip.adapted.json")) {
        return filename.substr(0, filename.size() - std::string(".tulip.adapted.json").size());
    }
    return std::filesystem::path(filename).stem().string();
}

} // namespace

Launcher::Launcher(const std::string& inputFile, const std::string& exportFolder)
    : inputFile_(inputFile),
      exportFolder_(exportFolder)
{
    if (exportFolder_.empty()) {
        exportFolder_ = "./" + std::filesystem::path(inputFile).parent_path().string() + "/";
    }
}

void Launcher::run()
{
    std::cout << "Loading input file: " << inputFile_ << std::endl;
    const std::string outputPrefix = extractCaseName(inputFile_) + ".";
    const std::string driverExportFolder = ensureTrailingSlash(exportFolder_) + outputPrefix;

    if (isAdaptedJson(inputFile_)) {
        auto driver = Driver::loadFromAdaptedFile(inputFile_);
        driver.setExportFolder(driverExportFolder);
        std::cout << "Running Tulip analysis..." << std::endl;
        driver.run();
    }
    else if (isInputJson(inputFile_)) {
        const bool initializedHere = !gmsh::isInitialized();
        if (initializedHere) {
            gmsh::initialize();
        }

        try {
            Adapter adapter(inputFile_);
            AdaptedInputParser parser(inputFile_, adapter.getAdaptedInputJSON());
            Driver driver(parser.readModel(), parser.readDriverOptions());
            driver.setExportFolder(driverExportFolder);
            std::cout << "Running Tulip analysis..." << std::endl;
            driver.run();
        }
        catch (...) {
            if (initializedHere) {
                gmsh::finalize();
            }
            throw;
        }

        if (initializedHere) {
            gmsh::finalize();
        }
    }
    else {
        throw std::runtime_error(
            "Unsupported input file extension: expected .tulip.input.json or .tulip.adapted.json");
    }
    
    std::cout << "-- tulip finished successfully --" << std::endl;
}

} // namespace tulip
