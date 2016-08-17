#include "readObj.h"

objReader::objReader()
{
}

objReader::~objReader()
{
}

objReader::objReader(objReader & copy)
{
}

bool objReader::readObj(std::string inputFile)
{
	bool ret = tinyobj::LoadObj(shapes, materials, err, inputFile.c_str());

	if (!err.empty())
		std::cerr << err << std::endl;

	if (!ret)
		return ret;
	nF = 0;
	nV = 0;

	int numIndices = 0;
	int numNormals = 0;
	int numTex = 0;
	for (size_t i = 0; i < shapes.size(); i++)
	{
		nF+= shapes[i].mesh.num_vertices.size();
		nV += shapes[i].mesh.positions.size()/3;
		numIndices += shapes[i].mesh.indices.size();
		numNormals += shapes[i].mesh.normals.size()/3;
		numTex += shapes[i].mesh.texcoords.size()/2;
		
		size_t indexOffset = 0;
//		for (size_t n = 0; n < shapes[i].mesh.num_vertices.size(); n++)
//		{
//			int ngon = shapes[i].mesh.num_vertices[n];
//			for (size_t f = 0; f < ngon; f++)
//			{
//				unsigned int v = shapes[i].mesh.indices[indexOffset + f];
//				printf("face [%ld] v[%ld] = (%f, %f, %f)\n", n, f,
//					shapes[i].mesh.positions[3 * v + 0],
//					shapes[i].mesh.positions[3 * v + 1],
//					shapes[i].mesh.positions[3 * v + 2]);
//			}
//			indexOffset += ngon;
//		}
	}

	printf("numVertices = %d\n", nV);
	printf("numFaces = %d\n", nF);

	return ret;
}
