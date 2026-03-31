
#include "Driver.h"
#include "Adapter.h"

#include <filesystem>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

// TODO change this to be a minimal launcher which creates an object called Tulip. Include a test with empty_coax.
int main(int argc, char* argv[])
{
	// Parses arguments.
	std::string inputFilename;
	
	po::options_description optionsDescription(
		"-- tulip --\n"
		"Transmission line unit length conductors and in-cell parameters calculator.\n"
		"Visit https://github.com/OpenSEMBA/tulip for more information.\n"
		"Available options"
	);
	optionsDescription.add_options()
		("help,h", "this help message")
		("input,i", po::value(&inputFilename), "input JSON file to use a .step file.")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, optionsDescription), vm);
	po::notify(vm);

	if (vm.count("help") || vm.empty()) {
		std::cout << optionsDescription << std::endl;
		return 0;
	}
	
	std::cout << optionsDescription << std::endl;
	std::cout << "Using input file: " << inputFilename << std::endl;

	// Launcher.
	std::filesystem::path folder{ 
		"./" + std::filesystem::path{inputFilename}.parent_path() + "/"
	};

	auto adapter{ Adapter::loadFromFile(inputFilename) };
	adapter.run();
	
	Driver driver{ 
		adapter.getModel(),
		adapter.getDriverOptions()
	};
	driver.setExportFolder(folder);
	driver.run();

	std::cout << "-- tulip finished successfully --" << std::endl;
}
