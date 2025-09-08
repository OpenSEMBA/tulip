#include "Parser.h"
#include "constants.h"

#include <filesystem>


using json = nlohmann::json;

namespace pulmtln {

enum class MaterialType {
	PEC,
	Dielectric,
	OpenBoundary,
	Vacuum
};

const std::map<std::string, MaterialType> LABEL_TO_MATERIAL_TYPE{
	{"PEC", MaterialType::PEC},
	{"Dielectric", MaterialType::Dielectric},
	{"OpenBoundary", MaterialType::OpenBoundary},
	{"Vacuum", MaterialType::Vacuum}
};

json readJSON(const std::string& fn)
{
	std::ifstream stream(fn);
	if (!stream.is_open()) {
		std::runtime_error("Unable to open file: " + fn);
	}

	return json::parse(stream);
}

Parser::Parser(const std::string& filename) :
	filename_{filename},
	json_(std::move(readJSON(filename)))
{}

Materials readMaterials(const json& j)
{
	Materials res;
	for (const auto& jMat : j.items()) {
		auto name = jMat.key();
		auto mat = jMat.value();
		auto type{ 
			LABEL_TO_MATERIAL_TYPE.at(
				mat.at("type").get<std::string>()
			)
		};
		auto attribute{ mat.at("tag").get<int>() };
		switch (type) {
		case MaterialType::PEC:
		{
			double area = mat.value("area", 0.0)*1e-6;
			res.pecs.push_back({ name, attribute, area });
			break;
		}
		case MaterialType::OpenBoundary:
			res.openBoundaries.push_back({ name, attribute });
			break;
		case MaterialType::Vacuum:
			res.dielectrics.push_back({ name, attribute, VACUUM_RELATIVE_PERMITTIVITY });
			break;
		case MaterialType::Dielectric:
		{
			double epsR{ mat.at("eps_r").get<double>() };
			res.dielectrics.push_back({ name, attribute, epsR });
			break;
		}
		default:
			throw std::runtime_error("Invalid material type");
		}
	}
	return res;
}

Model Parser::readModel() const
{
	const auto& j = json_.at("model");

	auto directory = 
		"./" + std::filesystem::path{filename_}.parent_path().string() + "/";

	std::string gmshFilename{ 
		directory + j.at("gmshFile").get<std::string>() 
	} ;

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
	const auto& j = json_.at("analysis");
	
	DriverOptions res;
	setIfExists<int>(j,  res.solverOptions.order, "order");
	setIfExists<bool>(j, res.solverOptions.printIterations, "printIterations");
	
	setIfExists<bool>(j, res.exportParaViewSolution, "exportParaviewSolution");
	setIfExists<bool>(j, res.exportVisItSolution, "exportVisItSolution");
	setIfExists<std::string>(j, res.exportFolder, "exportFolder");

	return res;
}

}