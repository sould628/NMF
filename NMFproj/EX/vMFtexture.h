#pragma once

#include <iostream>
#include <cstdlib>
#include <opencv2/opencv.hpp>
#include "FreeImage.h"


class vMFtexture {
private:
	float *originalNormals;
	int oWidth, oHeight;
	float ***vMFdata;
	int **vWidth, vHeight;
	int numLobes;
	int mipmapLevel;
	int *vMFmaps; //glGenerate
	cv::Mat cvOriginalNormals;


public:
	vMFtexture();
	vMFtexture(const char* filename, int numLobes = 4, int mipmapLevel = -1);
};


namespace vMFfunc {
	extern FIBITMAP* LoadImage(const char* filename, int &imageWidth, int &imageHeight);
	extern cv::Mat cvLoadImage(const char* filename, int &imageWidth, int &imageHeight);
}