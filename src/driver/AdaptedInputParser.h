#pragma once

#include "DriverOptions.h"
#include "Model.h"

#include <nlohmann/json.hpp>

namespace tulip {

class Parser {
public:
	Parser(const std::string& filename);

	Model readModel() const;
	DriverOptions readDriverOptions() const;

private:
	nlohmann::json json_;
	std::string filename_;
};

using AdaptedInputParser = Parser;

}