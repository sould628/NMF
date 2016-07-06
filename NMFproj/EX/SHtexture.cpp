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
	createYlmtex();
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

	for (int w = 0; w < oWidth; w++)
	{
		for (int h = 0; h < oHeight; h++)
		{
			cv::Vec3f temp= originalNormals[0].at<cv::Vec3f>(h, w);
			temp = cv::normalize(temp);
			originalNormals[0].at<cv::Vec3f>(h, w) = temp;
		}
	}

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

void SHtexture::createYlmtex()
{
	float center[2] = { (float)oWidth / 2.f, (float)oHeight / 2.f };
	int i = 0;
	float curPos[3] = { 0.f };
	for (int l = 0; l < order; l++)
	{
		for (int m = -l; m <= l; m++)
		{
			YlmTex[i] = cv::Mat::zeros(oHeight, oWidth, CV_32FC1);
			for (int w = 0; w < oWidth; w++)
			{
				for (int h = 0; h < oHeight; h++)
				{
					curPos[0] = ((float)w - center[0]) * 2.f / (float)oWidth;
					curPos[1] = ((float)h - center[1]) * 2.f / (float)oHeight;
					curPos[2] = 0.f;
					float xxpyy = curPos[0] * curPos[0] + curPos[1] * curPos[1];
					if (xxpyy < 1.f)
					{
						curPos[2] = sqrtf(1 - xxpyy);
						YlmTex[i].at<float>(h, w) = SHfunc::calSHcoef(l, m, curPos);
					}
				}//h
			}//w
			i++;
		}//m
	}//l
}

void SHtexture::bindTexture(GLuint &normalizedNMT, GLuint *SHtex)
{
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 0);
	glGenTextures(10, SHtex);
	glGenTextures(1, &normalizedNMT);

	float *shTextureData = new float[4 * oWidth* oHeight];
	float *originalNMT = new float[4 * oWidth*oHeight];
	glActiveTexture(GL_TEXTURE31);
	glBindTexture(GL_TEXTURE_2D, normalizedNMT);
	for (int h = 0; h < oHeight; h++)
	{
		for (int w = 0; w < oWidth; w++)
		{
			originalNMT[4 * (h*oWidth + w) + 0] = this->originalNormals[0].at<cv::Vec3f>(h, w)[0];
			originalNMT[4 * (h*oWidth + w) + 1] = this->originalNormals[0].at<cv::Vec3f>(h, w)[1];
			originalNMT[4 * (h*oWidth + w) + 2] = this->originalNormals[0].at<cv::Vec3f>(h, w)[2];
			originalNMT[4 * (h*oWidth + w) + 3] = 0.f;
		}
	}
	glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, oWidth, oHeight, 0, GL_RGBA, GL_FLOAT, originalNMT);
	GLenum glError = glGetError();


	int idx = 0;
	for (int i = 0; i < order*order; i++)
	{
		if (((i > 0) && (i % 4 == 0))) //Texture Bind
		{
			switch (idx)
			{
			case 0: glActiveTexture(GL_TEXTURE11); break;
			case 1: glActiveTexture(GL_TEXTURE12); break;
			case 2: glActiveTexture(GL_TEXTURE13); break;
			case 3: glActiveTexture(GL_TEXTURE14); break;
			case 4: glActiveTexture(GL_TEXTURE15); break;
			case 5: glActiveTexture(GL_TEXTURE16); break;
			case 6: glActiveTexture(GL_TEXTURE17); break;
			case 7: glActiveTexture(GL_TEXTURE18); break;
			case 8: glActiveTexture(GL_TEXTURE19); break;
			}
			idx++;
			glBindTexture(GL_TEXTURE_2D, SHtex[i]);
			glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, oWidth, oHeight, 0, GL_RGBA, GL_FLOAT, shTextureData);
			GLenum glError = glGetError();
		}
		for (int h = 0; h < oHeight; h++)
		{
			for (int w = 0; w < oWidth; w++)
			{
				shTextureData[4 * (h*oWidth + w) + i % 4] = this->SHmaps[i].at<float>(h, w);

			}//w
		}//h
		if (i == (order*order - 1))
		{
			switch (idx)
			{
			case 0: glActiveTexture(GL_TEXTURE11); break;
			case 1: glActiveTexture(GL_TEXTURE12); break;
			case 2: glActiveTexture(GL_TEXTURE13); break;
			case 3: glActiveTexture(GL_TEXTURE14); break;
			case 4: glActiveTexture(GL_TEXTURE15); break;
			case 5: glActiveTexture(GL_TEXTURE16); break;
			case 6: glActiveTexture(GL_TEXTURE17); break;
			case 7: glActiveTexture(GL_TEXTURE18); break;
			case 8: glActiveTexture(GL_TEXTURE19); break;
			}
			idx++;
			glBindTexture(GL_TEXTURE_2D, SHtex[i]);
			glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, oWidth, oHeight, 0, GL_RGBA, GL_FLOAT, shTextureData);
			GLenum glError = glGetError();
		}
	}
	delete shTextureData;
	glActiveTexture(GL_TEXTURE10);
	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);


}

void SHtexture::bindYlm(GLuint *Ylmtex)
{
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 0);
	glGenTextures(10, Ylmtex);

	float *YlmTextureData = new float[4 * oWidth* oHeight];

	int idx = 0;
	for (int i = 0; i < order*order; i++)
	{
		if (((i > 0) && (i % 4 == 0))) //Texture Bind
		{
			switch (idx)
			{
			case 0: glActiveTexture(GL_TEXTURE21); break;
			case 1: glActiveTexture(GL_TEXTURE22); break;
			case 2: glActiveTexture(GL_TEXTURE23); break;
			case 3: glActiveTexture(GL_TEXTURE24); break;
			case 4: glActiveTexture(GL_TEXTURE25); break;
			case 5: glActiveTexture(GL_TEXTURE26); break;
			case 6: glActiveTexture(GL_TEXTURE27); break;
			case 7: glActiveTexture(GL_TEXTURE28); break;
			case 8: glActiveTexture(GL_TEXTURE29); break;
			}
			idx++;
			glBindTexture(GL_TEXTURE_2D, Ylmtex[i]);
			glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, oWidth, oHeight, 0, GL_RGBA, GL_FLOAT, YlmTextureData);
			GLenum glError = glGetError();
		}
		for (int h = 0; h < oHeight; h++)
		{
			for (int w = 0; w < oWidth; w++)
			{
				YlmTextureData[4 * (h*oWidth + w) + i % 4] = this->YlmTex[i].at<float>(h, w);

			}//w
		}//h
		if (i == (order*order - 1))
		{
			switch (idx)
			{
			case 0: glActiveTexture(GL_TEXTURE21); break;
			case 1: glActiveTexture(GL_TEXTURE22); break;
			case 2: glActiveTexture(GL_TEXTURE23); break;
			case 3: glActiveTexture(GL_TEXTURE24); break;
			case 4: glActiveTexture(GL_TEXTURE25); break;
			case 5: glActiveTexture(GL_TEXTURE26); break;
			case 6: glActiveTexture(GL_TEXTURE27); break;
			case 7: glActiveTexture(GL_TEXTURE28); break;
			case 8: glActiveTexture(GL_TEXTURE29); break;
			}
			idx++;
			glBindTexture(GL_TEXTURE_2D, Ylmtex[i]);
			glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, oWidth, oHeight, 0, GL_RGBA, GL_FLOAT, YlmTextureData);
			GLenum glError = glGetError();
		}
	}
	delete YlmTextureData;
	glActiveTexture(GL_TEXTURE10);
	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);
}

int SHtexture::getOrder()
{
	return this->order;
}

float brdfSH::BlinnPhong(int l, float exp)
{
	float ret=0.f;

	ret = expf(-l*l / (2 * exp));

	return ret;
}

float SHfunc::calSHcoef(int l, int m, float val[3])
{
	float coeff;
	cv::Vec3f temp = cv::Vec3f(val[0], val[1], val[2]);
	temp = cv::normalize(temp);
	float x = temp[0];
	float y = temp[1];
	float z = temp[2];
	switch (l)
	{
	case 0://l=0
	{
		switch (m)
		{
		case 0:
			coeff = (1.f / 2.f)*sqrtf(1.f / cvPIf);
			break;
		}
		break;
	}//l=0
	case 1://l=1
	{
		switch (m)
		{
		case -1:
		{
			coeff = sqrtf(3.f / (4.f*cvPIf));
			coeff *= y;
			break;
		}
		case 0:
		{
			coeff = sqrtf(3.f / (4.f*cvPIf));
			coeff *= z;
			break;
		}
		case 1:
		{
			coeff = sqrtf(3.f / (4.f*cvPIf));
			coeff *= x;
			break;
		}
		}

		break;
	}//l=1
	case 2://l-=2
	{
		switch (m)
		{
		case -2:
		{
			coeff = (1.f / 2.f)*sqrtf(15.f / cvPIf);
			coeff *= (x*y);
			break;
		}
		case -1:
		{
			coeff = (1.f / 2.f)*sqrtf(15.f / cvPIf);
			coeff *= (y*z);
			break;
		}
		case 0:
		{
			coeff = (1.f / 4.f)*sqrtf(5.f / cvPIf);
			coeff *= ((2.f*z*z) - (x*x) - (y*y));
			break;
		}
		case 1:
		{
			coeff = (1.f / 2.f)*sqrtf(15.f / cvPIf);
			coeff *= (z*x);
			break;
		}
		case 2:
		{
			coeff = (1.f / 4.f)*sqrtf(15.f / cvPIf);
			coeff *= ((x*x) - (y*y));
			break;
		}
		}
		break;
	}//l=2
	case 3://l=3
	{
		switch (m)
		{
		case -3:
		{
			coeff = (1.f / 4.f)*sqrtf(35.f / (2.f*cvPIf));
			coeff *= ((3.f * x*x) - (y*y))*y;
			break;
		}
		case -2:
		{
			coeff = (1.f / 2.f)*sqrtf(105.f / (cvPIf));
			coeff *= x*y*z;
			break;
		}
		case -1:
		{
			coeff = (1.f / 4.f)*sqrtf(21.f / (2.f*cvPIf));
			coeff *= ((4.f * z*z) - (x*x) - (y*y))*y;
			break;
		}
		case 0:
		{
			coeff = (1.f / 4.f)*sqrtf(7.f / (cvPIf));
			coeff *= ((2.f * z*z) - (3.f*x*x) - (3.f*y*y))*z;
			break;
		}
		case 1:
		{
			coeff = (1.f / 4.f)*sqrtf(21.f / (2 * cvPIf));
			coeff *= ((4.f * z*z) - (x*x) - (y*y))*x;
			break;
		}
		case 2:
		{
			coeff = (1.f / 4.f)*sqrtf(105.f / cvPIf);
			coeff *= ((x*x) - (y*y))*z;
			break;
		}
		case 3:
		{
			coeff = (1.f / 4.f)*sqrt(35.f / (2.f*cvPIf));
			coeff *= ((x*x) - (3.f*y*y))*x;
			break;
		}
		}
		break;
	}//l=3
	case 4://l=4
	{
		switch (m)
		{
		case -4:
		{
			coeff = (3.f / 4.f)*sqrtf(35.f / cvPIf);
			coeff *= ((x*x) - (y*y))*x*y;
			break;
		}
		case -3:
		{
			coeff = (3.f / 4.f)*sqrtf(35.f / (2.f*cvPIf));
			coeff *= ((3.f * x*x) - (y*y))*y*z;
			break;
		}
		case -2:
		{
			coeff = (3.f / 4.f)*sqrtf(5.f / cvPIf);
			coeff *= ((7.f * z*z) - 1.f)*x*y;
			break;
		}
		case -1:
		{
			coeff = (3.f / 4.f)*sqrtf(5.f / (2.f*cvPIf));
			coeff *= ((7.f * z*z) - 3.f)*y*z;
			break;
		}
		case 0:
		{
			coeff = (3.f / 16.f)*sqrtf(1.f / cvPIf);
			coeff *= ((35.f * z*z*z*z) - (30.f*z*z) + 3.f);
			break;
		}
		case 1:
		{
			coeff = (3.f / 4.f)*sqrtf(5.f / (2.f*cvPIf));
			coeff *= ((7.f * z*z) - 3.f)*x*z;
			break;
		}
		case 2:
		{
			coeff = (3.f / 8.f)*sqrtf(5.f / cvPIf);
			coeff *= ((x*x - y*y) * (7.f*z*z - 1.f));
			break;
		}
		case 3:
		{
			coeff = (3.f / 4.f)*sqrtf(35.f / (2.f*cvPIf));
			coeff *= ((x*x) - (3.f*y*y))*x*z;
			break;
		}
		case 4:
		{
			coeff = (3.f / 16.f)*sqrtf(35.f / cvPIf);
			coeff *= ((x*x)*(x*x - 3.f*y*y) - (y*y)*(3.f*x*x - y*y));
			break;
		}
		}
		break;
	}//l=4
	}//switch l
	return coeff;
}
