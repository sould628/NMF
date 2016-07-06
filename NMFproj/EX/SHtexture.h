#pragma once

#include <opencv2/opencv.hpp>
#include <GL/glew.h>
#include <GL/glut.h>
#include "cvCtrl.h"

class SHtexture {
private:
	cv::Mat originalNormals[2];
	cv::Mat* SHmaps;
	cv::Mat YlmTex[25];

	int order;


	int oWidth;
	int oHeight;

public:
	SHtexture();
	~SHtexture();
	SHtexture(SHtexture &copy);
	SHtexture(const char* filename, int order=5);

public:
	bool readFile(const char* filename);
	bool createSHtexture(int order);
	void displayMap(int idx, int skip = 0, int destroy = 1);
	void createYlmtex();

	void bindTexture(GLuint *SHtex);
	void bindYlm(GLuint *Ylmtex);
	int getOrder();
};

namespace brdfSH
{
	float BlinnPhong(int l, float exp);
}

namespace SHfunc
{
	float calSHcoef(int l, int m, float val[3]);
}