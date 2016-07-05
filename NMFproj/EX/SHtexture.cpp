#include "SHtexture.h"

SHtexture::SHtexture()
{
	order = 5;
	SHmaps = new cv::Mat[25];
}

SHtexture::~SHtexture()
{
	delete[] SHmaps;
	SHmaps = nullptr;
}

SHtexture::SHtexture(SHtexture & copy)
{
	this->order = copy.order;
	this->oWidth = copy.oWidth;
	this->oHeight = copy.oHeight;
	this->SHmaps = new cv::Mat[order*order];
	for (int i = 0; i < order*order; i++)
	{
		this->SHmaps[i] = copy.SHmaps[i].clone();
	}
}

SHtexture::SHtexture(const char * filename, int order)
{
	this->order = order;
	readFile(filename);
	createSHtexture(order);
}



bool SHtexture::readFile(const char * filename)
{
	this->originalNormals[1] = loadImage(filename, oWidth, oHeight);

	float resizeFactx, resizeFacty;
	resizeFactx = 512.f / (float)oWidth; resizeFacty = 512.f / (float)oHeight;

	oWidth = 512; oHeight = 512;

	cv::resize(originalNormals[1], originalNormals[1], cv::Size(), resizeFactx, resizeFacty);

	std::vector<cv::Mat> change(3); cv::split(this->originalNormals[1], change);
	cv::Mat temp = change[0].clone();
	change[0] = change[2].clone();
	change[2] = temp.clone();

	merge(change, this->originalNormals[1]);


	this->originalNormals[0] = this->originalNormals[1].clone();
	std::vector<cv::Mat> spl(3); cv::split(this->originalNormals[0], spl);
	cv::normalize(spl[0], spl[0], 1.0, -1.0, cv::NORM_MINMAX);
	cv::normalize(spl[1], spl[1], 1.0, -1.0, cv::NORM_MINMAX);
	cv::normalize(spl[2], spl[2], 1.0, 0.0, cv::NORM_MINMAX);

	cv::merge(spl, this->originalNormals[0]);

	return true;
}

bool SHtexture::createSHtexture(int order)
{
	delete[] this->SHmaps;
	int numCoeff = order*order;
	this->SHmaps = new cv::Mat[numCoeff];
	for (int i = 0; i < numCoeff; i++)
	{
		this->SHmaps[i] = cv::Mat::zeros(this->oHeight, this->oWidth, CV_32FC1);
	}
	float coeff;
	for (int h = 0; h < oHeight; h++)
	{
		for (int w = 0; w < oWidth; w++)
		{
			cv::Vec3f curNormal = originalNormals[0].at<cv::Vec3f>(h, w);
			curNormal=cv::normalize(curNormal);
			float x = curNormal[0]; float y = curNormal[1]; float z = curNormal[2];
			float testNorm = cv::norm(curNormal);
			for (int i = 0; i < order; i++)
			{
				switch (i)
				{
				case 0:
				{
					//l=0
					coeff = (1.f / 2.f)*sqrtf(1.f / cvPIf);
					SHmaps[i].at<float>(h, w) = coeff;
					break;
				}
				case 1:
				{
					//l=-1
					coeff = sqrtf(3.f / (4.f*cvPIf));
					coeff *= y;
					SHmaps[i + 0].at<float>(h, w) = coeff;
					//l=0
					coeff = sqrtf(3.f / (4.f*cvPIf));
					coeff *= z;
					SHmaps[i + 1].at<float>(h, w) = coeff;
					//l=1
					coeff = sqrtf(3.f / (4.f*cvPIf));
					coeff *= x;
					SHmaps[i + 2].at<float>(h, w) = coeff;
					break;
				}
				case 2:
				{
					//l=-2
					coeff = (1.f / 2.f)*sqrtf(15.f / cvPIf);
					coeff *= (x*y);
					SHmaps[i*i + 0].at<float>(h, w) = coeff;
					//l=-1
					coeff = (1.f / 2.f)*sqrtf(15.f / cvPIf);
					coeff *= (y*z);
					SHmaps[i*i + 1].at<float>(h, w) = coeff;
					//l=0
					coeff = (1.f / 4.f)*sqrtf(5.f / cvPIf);
					coeff *= ((2.f*z*z)-(x*x)-(y*y));
					SHmaps[i*i + 2].at<float>(h, w) = coeff;
					//l=1
					coeff = (1.f / 2.f)*sqrtf(15.f / cvPIf);
					coeff *= (z*x);
					SHmaps[i*i + 3].at<float>(h, w) = coeff;
					//l=2
					coeff = (1.f / 4.f)*sqrtf(15.f / cvPIf);
					coeff *= ((x*x)-(y*y));
					SHmaps[i*i + 4].at<float>(h, w) = coeff;
					break;
				}
				case 3:
				{
					//l=-3
					coeff = (1.f / 4.f)*sqrtf(35.f / (2.f*cvPIf));
					coeff *= ((3.f * x*x) - (y*y))*y;
					SHmaps[i*i + 0].at<float>(h, w) = coeff;
					//l=-2
					coeff = (1.f / 2.f)*sqrtf(105.f / (cvPIf));
					coeff *= x*y*z;
					SHmaps[i*i + 1].at<float>(h, w) = coeff;
					//l=-1
					coeff = (1.f / 4.f)*sqrtf(21.f / (2.f*cvPIf));
					coeff *= ((4.f * z*z) - (x*x) - (y*y))*y;
					SHmaps[i*i + 2].at<float>(h, w) = coeff;
					//l=0
					coeff = (1.f / 4.f)*sqrtf(7.f / (cvPIf));
					coeff *= ((2.f * z*z) - (3.f*x*x) - (3.f*y*y))*z;
					SHmaps[i*i + 3].at<float>(h, w) = coeff;
					//l=1
					coeff = (1.f / 4.f)*sqrtf(21.f / (2 * cvPIf));
					coeff *= ((4.f * z*z) - (x*x) - (y*y))*x;
					SHmaps[i*i + 4].at<float>(h, w) = coeff;
					//l=2
					coeff = (1.f / 4.f)*sqrtf(105.f / cvPIf);
					coeff *= ((x*x) - (y*y))*z;
					SHmaps[i*i + 5].at<float>(h, w) = coeff;
					//l=3
					coeff = (1.f / 4.f)*sqrt(35.f / (2.f*cvPIf));
					coeff *= ((x*x) - (3.f*y*y))*x;
					SHmaps[i*i + 6].at<float>(h, w) = coeff;
					break;
				}
				case 4:
				{
					//l=-4
					coeff = (3.f / 4.f)*sqrtf(35.f / cvPIf);
					coeff *= ((x*x)-(y*y))*x*y;
					SHmaps[i*i + 0].at<float>(h, w) = coeff;

					//l=-3
					coeff = (3.f / 4.f)*sqrtf(35.f / (2.f*cvPIf));
					coeff *= ((3.f * x*x) - (y*y))*y*z;
					SHmaps[i*i + 1].at<float>(h, w) = coeff;

					//l=-2
					coeff = (3.f / 4.f)*sqrtf(5.f / cvPIf);
					coeff *= ((7.f * z*z) - 1.f)*x*y;
					SHmaps[i*i + 2].at<float>(h, w) = coeff;

					//l=-1
					coeff = (3.f / 4.f)*sqrtf(5.f / (2.f*cvPIf));
					coeff *= ((7.f * z*z) - 3.f)*y*z;
					SHmaps[i*i + 3].at<float>(h, w) = coeff;

					//l=0
					coeff = (3.f / 16.f)*sqrtf(1.f / cvPIf);
					coeff *= ((35.f * z*z*z*z) - (30.f*z*z) + 3.f);
					SHmaps[i*i + 4].at<float>(h, w) = coeff;

					//l=1
					coeff = (3.f / 4.f)*sqrtf(5.f / (2.f*cvPIf));
					coeff *= ((7.f * z*z) - 3.f)*x*z;
					SHmaps[i*i + 5].at<float>(h, w) = coeff;

					//l=2
					coeff = (3.f / 8.f)*sqrtf(5.f / cvPIf);
					coeff *= ((x*x - y*y) * (7.f*z*z - 1.f));
					SHmaps[i*i + 6].at<float>(h, w) = coeff;

					//l=3
					coeff = (3.f / 4.f)*sqrtf(35.f / (2.f*cvPIf));
					coeff *= ((x*x) - (3.f*y*y))*x*z;
					SHmaps[i*i + 7].at<float>(h, w) = coeff;

					//l=4
					coeff = (3.f / 16.f)*sqrtf(35.f / cvPIf);
					coeff *= ((x*x)*(x*x-3.f*y*y) - (y*y)*(3.f*x*x-y*y));
					SHmaps[i*i + 8].at<float>(h, w) = coeff;
					break;
				}
				}//switch
			}//i
		}//w
	}//h
	return true;
}

void SHtexture::displayMap(int idx, int skip, int destroy)
{
	char key = 0;
	cv::namedWindow("SHmap"); cv::namedWindow("SHmap", CV_WINDOW_NORMAL);
	while (key != myESC)
	{
		cv::imshow("SHmap", SHmaps[idx]);
		if (skip == 1)
		{
			cv::waitKey(1);
			key = myESC;
		}
		else {
			key = cv::waitKey();
			switch (key)
			{
			case ',':
			case '<':
				if (idx == 0)
					idx = order*order - 1;
				else idx--;
				break;

			case '.':
			case '>':
				if (idx == (order*order - 1))
					idx = 0;
				else idx++;
				break;
			}
		}
	}
	if (destroy == 1)
		cv::destroyWindow("SHmap");

}

float brdfSH::BlinnPhong(int l, float val[3])
{
	float ret;


	return 0.0f;
}
