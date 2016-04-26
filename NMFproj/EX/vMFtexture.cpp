#include "vMFtexture.h"


//const, dest
vMFtexture::vMFtexture()
{

}
vMFtexture::vMFtexture(const char* filename, int numLobes, int mipmapLevel)
{

	this->numLobes = numLobes;

	this->originalNormals[0] = vMFfunc::cvLoadImage(filename, this->oWidth, this->oHeight);
	this->originalNormals[1] = this->originalNormals[0].clone();
	originalNormals[1] *= 2;
	originalNormals[1] -= 1;
//	std::vector < cv::Mat >spl; cv::split(this->originalNormals[1], spl);
//	zyx[1] *= 2; zyx[2] *= 2;
//	zyx[1] -= 1; zyx[2] -= 1;
//	cv::merge(zyx, originalNormals[1]);


	if (mipmapLevel == -1)
	{
		int wLevel = 0, hLevel = 0;
		int w = this->oWidth; int h = this->oHeight;
		while (w > 1)
		{
			w /= 2;
			wLevel++;
		}
		while (h > 1)
		{
			h /= 2;
			hLevel++;
		}
		this->mipmapLevel = (wLevel < hLevel ? hLevel : wLevel);
	}
	else { this->mipmapLevel = mipmapLevel; }

	vMFmaps = new int[numLobes];

	this->vWidth = new int[this->mipmapLevel]; this->vHeight = new int[this->mipmapLevel];
	int w = this->oWidth; int h = this->oHeight;
	this->vMFdata = new float****[this->mipmapLevel+1];
	for (int m = 0; m <= this->mipmapLevel; m++)
	{
		vMFdata[m] = new float***[numLobes];
		vWidth[m] = w; vHeight[m] = h;
		for (int l = 0; l < numLobes; l++)
		{
			vMFdata[m][l] = new float**[w];
			for (int i = 0; i < w; i++)
			{
				vMFdata[m][l][w] = new float*[h];
				for (int j = 0; j < h; j++)
				{
					vMFdata[m][l][w][h] = new float[4];
				}
			}
		}
		if (w > 1) w /= 2; if (h > 1) h /= 2;
	}
}
vMFtexture::~vMFtexture()
{
	delete[] vMFmaps;
	for (int m = 0; m <= this->mipmapLevel; m++)
	{
		for (int l = 0; l < this->numLobes; l++)
		{
			for (int i = 0; i < this->vWidth[m]; i++)
			{
				for (int j = 0; j < this->vHeight[m]; j++)
				{
					delete[] this->vMFdata[m][l][i][j];
				}
				delete[] this->vMFdata[m][l][i];
			}
			delete[] this->vMFdata[m][l];
		}
		delete[] this->vMFdata[m];
	}
	delete[] this->vMFdata;
	delete[] vWidth; delete[] vHeight;
}


//Parameter Functions


//display Functions
void vMFtexture::showOriginalImage(int channel) const
{
	cv::namedWindow("originalImage"); cv::namedWindow("originalImage", CV_WINDOW_NORMAL);
	char key = 0;
	std::vector<cv::Mat> zyx0;
	std::vector<cv::Mat> zyx1;
	cv::split(this->originalNormals[0], zyx0);
	cv::split(this->originalNormals[1], zyx1);

	int ver = 0;
	while (key != myESC)
	{
		if (channel == -1)
		{
			cv::imshow("originalImage", originalNormals[ver]);
			key = cv::waitKey();
			switch (key)
			{
			case',':case'<':case'.':case'>':
				ver == 0 ? ver = 1 : ver = 0;
				break;
			case 'c':case 'C':
				channel = 0;
				break;
			default:
				break;
			}
		}
		else
		{
			if (ver == 0)
				cv::imshow("originalImage", zyx0[channel]);
			else
				cv::imshow("originalImage", zyx1[channel]);
			key = cv::waitKey();
			switch (key)
			{
			case ',':case'<':
				if (channel == 0) channel = 2;
				else channel--;
				break;
			case '.':case'>':
				if (channel == 2) channel = 0;
				else channel++;
				break;
			case 'c':case 'C':
				channel = -1;
				break;
			default:
				break;
			}
		}
	}


	cv::destroyWindow("originalImage");
}


//vMFfunc
cv::Mat vMFfunc::cvLoadImage(const char* filename, int &imageWidth, int &imageHeight)
{
	cv::Mat load, floatScale;
	load = cv::imread(filename, CV_LOAD_IMAGE_COLOR);
	load.convertTo(floatScale, CV_32FC1, 1.f / 255.f);
	imageWidth = load.cols; imageHeight = load.rows;
	return floatScale;
}
float vMFfunc::vMF(float normal[3], float mu[3], float kappa) {
	double NdotMu = (mu[0] * normal[0]) + (mu[1] * normal[1]) + (mu[2] * normal[2]);
	double Kappa = kappa;
	double result = (Kappa / (4 * vmfPI*sinh(Kappa)))*exp(Kappa*(NdotMu));

	return (float)result;
}
FIBITMAP* vMFfunc::LoadImage(const char* filename, int &imageWidth, int &imageHeight) {
	FREE_IMAGE_FORMAT format = FreeImage_GetFileType(filename);
	RGBQUAD pixel;
	if (format == -1)
	{
		std::cout << "Could not find image: " << filename << " - Aborting." << std::endl;
		system("pause>press anykey to exit");
		exit(-1);
	}
	if (format == FIF_UNKNOWN)
	{
		std::cout << "Couldn't determine file format - attempting to get from file extension..." << std::endl;

		format = FreeImage_GetFIFFromFilename(filename);

		if (!FreeImage_FIFSupportsReading(format))
		{
			std::cout << "Detected image format cannot be read!" << std::endl;
			system("pause>press any key to exit");
			exit(-1);
		}
	}
	FIBITMAP* bitmap = FreeImage_Load(format, filename);
	int bitsPerPixel = FreeImage_GetBPP(bitmap);

	FIBITMAP* bitmap32;

	if (bitsPerPixel == 32)
	{
		std::cout << "Source image has " << bitsPerPixel << " bits per pixel. Skipping conversion." << std::endl;
		bitmap32 = bitmap;
	}
	else
	{
		std::cout << "Source image has " << bitsPerPixel << " bits per pixel. Converting to 32-bit colour." << std::endl;
		bitmap32 = FreeImage_ConvertTo32Bits(bitmap);
	}

	imageWidth = FreeImage_GetWidth(bitmap32);
	imageHeight = FreeImage_GetHeight(bitmap32);

	FreeImage_GetPixelColor(bitmap32, 25, 15, &pixel);
	float r, g, b;
	r = pixel.rgbRed;
	b = pixel.rgbBlue;
	g = pixel.rgbGreen;
	std::cout << "Image: " << filename << " is size: " << imageWidth << "x" << imageHeight << "." << std::endl;

	return bitmap32;

}