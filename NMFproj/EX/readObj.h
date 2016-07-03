#pragma once

#include <vector>
#include <string>
#include <iostream>
#include "tiny_obj_loader.h"

//shapes[i].mesh.indices[indexoffset+f]
class objReader
{
public:
	int nV;
	int nF;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;

private:


	std::string inputFile;
	std::string err;

public:
	objReader();
	~objReader();
	objReader(objReader &copy);

	bool readObj(std::string inputFile);
	 
};