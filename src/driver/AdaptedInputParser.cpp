#include "AdaptedInputParser.h"
#include "constants.h"

#include <filesystem>
#include <fstream>


using json = nlohmann::json;

namespace tulip {

enum class MaterialType {
	Conductor,
	Dielectric,
	OpenBoundary,
	Vacuum
};

const std::map<std::string, MaterialType> LABEL_TO_MATERIAL_TYPE{
	{"conductor", MaterialType::Conductor},
	{"dielectric", MaterialType::Dielectric},
	{"open", MaterialType::OpenBoundary}
};

json readJSON(const std::string& fn)
{
	std::ifstream stream(fn);
	if (!stream.is_open()) {
		throw std::runtime_error("Unable to open file: " + fn);
	}

	return json::parse(stream);
}

Parser::Parser(const std::string& filename) :
	filename_{filename},
	json_(std::move(readJSON(filename)))
{}

Parser::Parser(const std::string& filename, const nlohmann::json& j) :
	filename_{filename},
	json_(j)
{}

Materials readMaterials(const json& j)
{
	Materials res;
	ConductorId nextConductorId = 0;
	auto parseType = [](const std::string& typeLabel) {
		auto typeIt = LABEL_TO_MATERIAL_TYPE.find(typeLabel);
		if (typeIt == LABEL_TO_MATERIAL_TYPE.end()) {
			throw std::runtime_error("Unknown material type: " + typeLabel);
		}
		return typeIt->second;
	};

	auto readAttribute = [](const json& mat) {
		if (mat.contains("attribute")) {
			return mat.at("attribute").get<int>();
		}
		if (mat.contains("tag")) {
			return mat.at("tag").get<int>();
		}
		throw std::runtime_error("Missing material attribute/tag");
	};

	auto parseMaterial = [&](const json& mat) {
		auto type = parseType(mat.at("type").get<std::string>());
		auto attribute = readAttribute(mat);
		switch (type) {
		case MaterialType::Conductor:
		{
			auto conductorId = mat.at("conductorId").get<int>();
			res.addConductor(attribute, conductorId);
			break;
		}
		case MaterialType::OpenBoundary:
			res.addOpenBoundary(attribute);
			break;
		case MaterialType::Vacuum:
			res.addDielectric(attribute, VACUUM_RELATIVE_PERMITTIVITY);
			break;
		case MaterialType::Dielectric:
		{
			double epsR{ VACUUM_RELATIVE_PERMITTIVITY };
			if (mat.contains("relativePermittivity")) {
				epsR = mat.at("relativePermittivity").get<double>();
			}
			res.addDielectric(attribute, epsR);
			break;
		}
		default:
			throw std::runtime_error("Invalid material type");
		}
	};

	if (j.is_array()) {
		for (const auto& mat : j) {
			parseMaterial(mat);
		}
	} else {
		throw std::runtime_error("materials object must be an array.");
	}
	return res;
}

Model Parser::readModel() const
{
	const auto& j = json_.at("model");

	auto gmshFilename =
		(std::filesystem::path{filename_}.parent_path() /
		 j.at("gmshFile").get<std::string>()).generic_string();

	return Model{
		mfem::Mesh::LoadFromFile(gmshFilename),
		readMaterials(j.at("materials"))
	};
}


template <class T>
static void setIfExists(const json& j, T& entry, std::string labelToCheck)
{
	auto const it = j.find(labelToCheck);
	if (it != j.end()) {
		entry = it->get<T>();
	}
}

DriverOptions Parser::readDriverOptions() const
{
	const auto& j = json_.at("driverOptions");
	
	DriverOptions res;
	setIfExists<int>(j,  res.solverOptions.order, "order");
	setIfExists<bool>(j, res.solverOptions.printIterations, "printIterations");
	
	setIfExists<bool>(j, res.exportParaViewSolution, "exportParaViewSolution");
	setIfExists<std::string>(j, res.exportFolder, "exportFolder");

	return res;
}

}