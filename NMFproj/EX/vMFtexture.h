#pragma once

#include<math.h>
#include <iostream>
#include <cstdlib>
#include <opencv2/opencv.hpp>
#include "FreeImage.h"

#define vmfPI 3.141592653589793f

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
	extern float vMF(float normal[3], float mu[3], float kappa);
}