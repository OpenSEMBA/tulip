
#include "Driver.h"

#include <filesystem>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
	// Parses arguments.
	std::string inputFilename;
	
	po::options_description desc(
		"-- tulip --\n"
		"Transmission line unit length conductors and in-cell parameters calculator.\n"
		"Visit https://github.com/OpenSEMBA/tulip for more information.\n"
		"Available options"
	);
	desc.add_options()
		("help,h", "this help message")
		("input,i", po::value(&inputFilename), "input filename in JSON format")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help") || vm.empty()) {
		std::cout << desc << std::endl;
		return 0;
	}
	
	std::cout << desc << std::endl;
	std::cout << "Using input file: " << inputFilename << std::endl;


	// Launcher.
	std::string folder{ 
		"./" + std::filesystem::path{inputFilename}.parent_path().string() + "/"
	};

	auto driver{ tulip::Driver::adaptFromFile(inputFilename) };
	driver.setExportFolder(folder);

	driver.run();

	std::cout << "-- tulip finished successfully --" << std::endl;
}
