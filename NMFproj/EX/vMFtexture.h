#pragma once

#include <vector>
#include<math.h>
#include <iostream>
#include <cstdlib>
#include <opencv2/opencv.hpp>
#include "FreeImage.h"
#include "cvCtrl.h"

#define vmfPI 3.141592653589793f

class vMFtexture {
private:
	//Vec3f
	cv::Mat originalNormals[2];
	
	int oWidth, oHeight;
	cv::Mat **vMFdata;
	int *vWidth, *vHeight;
	int numLobes;
	int mipmapLevel;
	int *vMFmaps; //glGenerate

	void computeParameters(float *alpha, float **aux, cv::Mat targetRegion, float prevData[4][4]);


public:
	vMFtexture();
	vMFtexture(const char* filename, int numLobes = 4, int mipmapLevel = -1);
	~vMFtexture();
	void showvMFImage(int level, int lobe, int mode=0) const;
	void showOriginalImage(int channel = -1) const;
	void generatevMFmaps();

public:
	int getWidth(int level) const {
		if (level > mipmapLevel - 1) { std::cout << "error\n"; return 0; }
		else return vWidth[level];
	}
	int getHeight(int level) const {
		if (level > mipmapLevel - 1) { std::cout << "error\n"; return 0; }
		else return vHeight[level];
	}
	int getMipmapLevel() const { return mipmapLevel; }
	int getNumLobes() const { return numLobes; }
	float getvMFcompoenent(int level, int lobe, int w, int h, int c) const
	{
		return vMFdata[level][lobe].at<cv::Vec4f>(w, h)[c];
	}
};


namespace vMFfunc {
	extern FIBITMAP* LoadImage(const char* filename, int &imageWidth, int &imageHeight);
	extern cv::Mat cvLoadImage(const char* filename, int &imageWidth, int &imageHeight);
	extern float vMF(float normal[3], float mu[3], float kappa);
}

namespace vectorFunc
{
	void normalize(float input[3]);
	float norm(float input[3]);

}