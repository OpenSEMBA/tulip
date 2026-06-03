#pragma once

#include <string>
#include <map>

namespace tulip
{
struct AdapterOptions
{
	// Multiplication factor applied to the size of the box created 
	// to determine the inner region used  for unshielded multiwires 
	// bundles with no open boundary defined. 
	double innerRegionBoxScalingFactor = 1.30; 

	// For an open case in which no open boundary is defined a 
	// circle is created enclosing all the defined entities.
	// The diameter is determined using this number.
	double farRegionDiskScalingFactor = 4.0;

	// If empty, it defaults to "CASE_NAME.step", where CASE_NAME is
	// the basename of the input file.
	std::string stepFilename;
	
	// Options which are passed to gmsh. 
	// User specified options modify these options.
	std::map<std::string, double> gmshOptions = {
		{"Mesh.MshFileVersion", 2.2},
		{"Mesh.MeshSizeFromCurvature", 40.0},
		{"Mesh.ElementOrder", 3.0},
		{"Mesh.ScalingFactor", 1e-3},
		{"Mesh.SurfaceFaces", 1.0},
		{"Mesh.MeshSizeMax", 2.5}
	};
};

}