#include <string>
#include <GL/glut.h>
#include <vector>
#ifndef __MODELOBJ_H__
#define __MODELOBJ_H__


typedef unsigned int uint;
struct vector3
{
	float x,y,z;
};
struct TCoord
{
    float t;
    float s;
};

struct Triangle
{
    uint v[3];
    uint vt[3];
    uint vn[3];
};


class ObjLoader
{
public:
	ObjLoader();
	void loadObjModel(char* filename);
	void draw();
	void scale(float scaleX, float scaleY, float scaleZ);
	void translate(float translateX, float translateY, float translateZ);	
	void setMatKs(float _matKs[4]);
	void setMatKd(float _matKd[4]);
	void setMatKa(float _matKa[4]);
	void setMatShininess(float _matShininess);
	float getMaxX();
	float getMinX();
	float getMaxY();
	float getMinY();
	float getMaxZ();
	float getMinZ();

	std::vector<TCoord> tc;
	std::vector<vector3> vert;
	std::vector<vector3> norm;
	std::vector<Triangle> tri;

private:
	bool _loaded;

	GLfloat  matKs[4], matKd[4], matKa[4];
	GLfloat matShininess;
	GLfloat maxX, minX, maxY, minY, maxZ, minZ;
};			 
			 
#endif
