#include "Tulip.h"
#include "Driver.h"

#include <filesystem>
#include <iostream>

namespace tulip {

Tulip::Tulip(const std::string& adaptedJsonFile, const std::string& exportFolder)
    : adaptedJsonFile_(adaptedJsonFile),
      exportFolder_(exportFolder)
{
    if (exportFolder_.empty()) {
        exportFolder_ = "./" + std::filesystem::path(adaptedJsonFile).parent_path().string() + "/";
    }
}

void Tulip::run()
{
    std::cout << "Loading adapted file: " << adaptedJsonFile_ << std::endl;
    auto driver = Driver::loadFromAdaptedFile(adaptedJsonFile_);
    
    std::cout << "Setting export folder: " << exportFolder_ << std::endl;
    driver.setExportFolder(exportFolder_);
    
    std::cout << "Running Tulip analysis..." << std::endl;
    driver.run();
    
    std::cout << "-- tulip finished successfully --" << std::endl;
}

} // namespace tulip
