#pragma once

#include "DriverOptions.h"
#include "Model.h"

#include <nlohmann/json.hpp>

namespace tulip {

class AdaptedInputParser {
public:
	AdaptedInputParser(const std::string& filename);

	Model readModel() const;
	DriverOptions readDriverOptions() const;

private:
	nlohmann::json json_;
	std::string filename_;
};

}