#pragma once
#include <vector>
#include <algorithm>
#include <opencv2\opencv.hpp>

#define myESC 27
#define mySpace 
#define myNext '.'
#define myPrev ','




namespace {
	namespace cvCtrl {
		namespace skip { enum Type { yes, no }; };
		namespace destroy { enum Type { yes, no }; };
	}

	void displayImage(const char *windowName, cv::Mat imgToshow, int skip = 0, int flag = 0)
	{

		char key = 0;
		cv::namedWindow(windowName); cv::namedWindow(windowName, CV_WINDOW_NORMAL);
		while (key != myESC)
		{
			cv::imshow(windowName, imgToshow);
			if (skip == 0)
			{
				cv::waitKey(1);
				key = myESC;
			}
			else {
				key = cv::waitKey();
			}
		}
		if (flag == 0)
			cv::destroyWindow(windowName);
	}

	//matType
	//0: CV_32FC1
	//1: CV_32FC3
	cv::Mat elementMul(cv::Mat mat1, cv::Mat mat2, int matType=0)
	{
		int matr, matc, mat2r, mat2c;
		matr = mat1.rows; matc = mat1.cols; mat2r = mat2.rows; mat2c = mat2.cols;
		if ((matr != mat2r) || (matc != mat2c))
		{
			std::cout << "Matrix dimension unmatched" << std::endl;
			return cv::Mat::zeros(1, 1, CV_32FC1);
		}
		else
		{
			if (matType == 0)
			{
				cv::Mat result=cv::Mat::zeros(matr, matc, CV_32FC1);
				for (int i = 0; i < matr; i++)
				{
					for (int j = 0; j < matc; j++)
					{
						result.at<float>(i, j) = mat1.at<float>(i, j)*mat2.at<float>(i, j);
					}
				}
				return result;
			}
			else if (matType == 1)
			{
				cv::Mat result = cv::Mat::zeros(matr, matc, CV_32FC3);
				for (int i = 0; i < matr; i++)
				{
					for (int j = 0; j < matc; j++)
					{
						result.at<cv::Vec3f>(i, j)[0] = mat1.at<cv::Vec3f>(i, j)[0] * mat2.at<cv::Vec3f>(i, j)[0];
						result.at<cv::Vec3f>(i, j)[1] = mat1.at<cv::Vec3f>(i, j)[1] * mat2.at<cv::Vec3f>(i, j)[1];
						result.at<cv::Vec3f>(i, j)[2] = mat1.at<cv::Vec3f>(i, j)[2] * mat2.at<cv::Vec3f>(i, j)[2];

					}
				}
				return result;
			}
		}
		if (matType == 0)
			return cv::Mat::zeros(1, 1, CV_32FC1);
		else
			return cv::Mat::zeros(1, 1, CV_32FC3);
	}

	//mat1/mat2
	//matType
	//0: CV_32FC1
	//1: CV_32FC3
	cv::Mat elementDiv(cv::Mat mat1, cv::Mat mat2, int matType = 0)
	{
		int matr, matc, mat2r, mat2c;
		matr = mat1.rows; matc = mat1.cols; mat2r = mat2.rows; mat2c = mat2.cols;
		if ((matr != mat2r) || (matc != mat2c))
		{
			std::cout << "Matrix dimension unmatched" << std::endl;
			return cv::Mat::zeros(1, 1, CV_32FC1);
		}
		else
		{
			if (matType == 0)
			{
				cv::Mat result = cv::Mat::zeros(matr, matc, CV_32FC1);
				for (int i = 0; i < matr; i++)
				{
					for (int j = 0; j < matc; j++)
					{
						result.at<float>(i, j) = mat1.at<float>(i, j)/mat2.at<float>(i, j);
					}
				}
				return result;
			}
			else if (matType == 1)
			{
				cv::Mat result = cv::Mat::zeros(matr, matc, CV_32FC3);
				for (int i = 0; i < matr; i++)
				{
					for (int j = 0; j < matc; j++)
					{
						result.at<cv::Vec3f>(i, j)[0] = mat1.at<cv::Vec3f>(i, j)[0] / mat2.at<cv::Vec3f>(i, j)[0];
						result.at<cv::Vec3f>(i, j)[1] = mat1.at<cv::Vec3f>(i, j)[1] / mat2.at<cv::Vec3f>(i, j)[1];
						result.at<cv::Vec3f>(i, j)[2] = mat1.at<cv::Vec3f>(i, j)[2] / mat2.at<cv::Vec3f>(i, j)[2];
					}
				}
				return result;
			}
		}
		if (matType == 0)
			return cv::Mat::zeros(1, 1, CV_32FC1);
		else
			return cv::Mat::zeros(1, 1, CV_32FC3);
	}


}