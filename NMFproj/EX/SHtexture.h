#pragma once

#include <opencv2/opencv.hpp>
#include "cvCtrl.h"

class SHtexture {
private:
	cv::Mat originalNormals[2];
	cv::Mat* SHmaps;

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

};

namespace brdfSH
{
	float BlinnPhong(int l, float val[3]);

}