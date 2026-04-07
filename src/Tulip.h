#pragma once

#include <string>

namespace tulip {

/// Minimal launcher for the Tulip transmission line calculator.
class Tulip {
public:
    /// Initialize Tulip with an adapted JSON input file.
    /// @param adaptedJsonFile Path to the adapted JSON file (.tulip.adapted.json)
    /// @param exportFolder Optional export folder for results. If empty, uses input file's folder.
    explicit Tulip(const std::string& adaptedJsonFile, const std::string& exportFolder = "");

    /// Run the transmission line analysis.
    void run();

private:
    std::string adaptedJsonFile_;
    std::string exportFolder_;
};

} // namespace tulip
