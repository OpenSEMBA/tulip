#include "Launcher.h"

#include "Adapter.h"
#include "AdaptedInputParser.h"
#include "Driver.h"

#include <chrono>
#include <filesystem>
#include <gmsh.h>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace tulip {

namespace {

using Clock = std::chrono::steady_clock;

double elapsedSeconds(Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

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

void printTimingSummary(
    bool hasMeshTiming,
    double meshSeconds,
    const DriverTimings& timings)
{
    std::cout << "Timing summary:" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    if (hasMeshTiming) {
        std::cout << "  mesh:       " << meshSeconds << " s" << std::endl;
    }
    else {
        std::cout << "  mesh:       N/A" << std::endl;
    }
    std::cout << "  solve:      " << timings.solveSeconds << " s" << std::endl;
    if (timings.multipolarComputed) {
        std::cout << "  multipolar: " << timings.multipolarSeconds << " s" << std::endl;
    }
    else {
        std::cout << "  multipolar: N/A" << std::endl;
    }
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
    const std::string caseNamePrefix = extractCaseName(inputFile_) + ".";
    const std::string outputPathPrefix = ensureTrailingSlash(exportFolder_) + caseNamePrefix;
    bool hasMeshTiming = false;
    double meshSeconds = 0.0;

    if (isAdaptedJson(inputFile_)) {
        auto driver = Driver::loadFromAdaptedFile(inputFile_);
        driver.setExportFolder(outputPathPrefix);
        std::cout << "Running Tulip analysis..." << std::endl;
        driver.run();
        printTimingSummary(hasMeshTiming, meshSeconds, driver.getTimings());
    }
    else if (isInputJson(inputFile_)) {
        const bool initializedHere = !gmsh::isInitialized();
        if (initializedHere) {
            gmsh::initialize();
        }

        try {
            const auto meshStart = Clock::now();
            Adapter adapter(inputFile_);
            meshSeconds = elapsedSeconds(meshStart);
            hasMeshTiming = true;
            AdaptedInputParser parser(inputFile_, adapter.getAdaptedInputJSON());
            Driver driver(parser.readModel(), parser.readDriverOptions());
            driver.setExportFolder(outputPathPrefix);
            std::cout << "Running Tulip analysis..." << std::endl;
            driver.run();
            printTimingSummary(hasMeshTiming, meshSeconds, driver.getTimings());
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
