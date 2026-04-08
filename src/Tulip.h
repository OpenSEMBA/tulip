#pragma once

#include <string>

namespace tulip {

/// Minimal launcher for the Tulip transmission line calculator.
class Tulip {
public:
    /// Initialize Tulip with an input file.
    /// @param inputFile Path to either .tulip.input.json or .tulip.adapted.json
    /// @param exportFolder Optional export folder for results. If empty, uses input file's folder.
    explicit Tulip(const std::string& inputFile, const std::string& exportFolder = "");

    /// Run the transmission line analysis.
    void run();

private:
    std::string inputFile_;
    std::string exportFolder_;
};

} // namespace tulip
