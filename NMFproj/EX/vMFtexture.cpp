#include "vMFtexture.h"


//const, dest
vMFtexture::vMFtexture()
{

}
vMFtexture::vMFtexture(const char* filename, int numLobes, int mipmapLevel)
{
	this->mipmapLevel = mipmapLevel;
	this->numLobes = numLobes;

	FIBITMAP* NMdata32bit = vMFfunc::LoadImage(filename, this->oWidth, this->oHeight);
	cvOriginalNormals = vMFfunc::cvLoadImage(filename, this->oWidth, this->oHeight);
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