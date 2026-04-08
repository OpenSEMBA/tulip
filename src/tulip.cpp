
#include "Tulip.h"

#include <boost/program_options.hpp>
#include <iostream>

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
	std::string inputFilename;
	std::string exportFolder;
	
	po::options_description optionsDescription(
		"-- tulip --\n"
		"Transmission line unit length conductors and in-cell parameters calculator.\n"
		"Visit https://github.com/OpenSEMBA/tulip for more information.\n"
		"Available options"
	);
	optionsDescription.add_options()
		("help,h", "this help message")
		("input,i", po::value(&inputFilename), "input JSON file (.tulip.input.json or .tulip.adapted.json)")
		("output,o", po::value(&exportFolder), "optional output folder for results")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, optionsDescription), vm);
	po::notify(vm);

	if (vm.count("help") || vm.empty()) {
		std::cout << optionsDescription << std::endl;
		return 0;
	}

	try {
		tulip::Tulip tulip(inputFilename, exportFolder);
		tulip.run();
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}

	return 0;
}
