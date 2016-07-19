

#include "vMFtexture.h"


//const, dest
vMFtexture::vMFtexture()
{

}
vMFtexture::vMFtexture(const char* filename, int numLobes, int mipmapLevel)
{

	this->numLobes = numLobes;

	this->originalNormals[1] = vMFfunc::cvLoadImage(filename, this->oWidth, this->oHeight);

	float resizeFactx, resizeFacty;

	resizeFactx = 256/(float)oWidth; resizeFacty = 256 / (float)oHeight;

	oWidth = 256; oHeight = 256;

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



	if (mipmapLevel == -1)
	{
		int wLevel = 1, hLevel = 1;
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
	this->vMFdata = new cv::Mat*[this->mipmapLevel];
	for (int m = 0; m < this->mipmapLevel; m++)
	{
		vMFdata[m] = new cv::Mat[numLobes];
		vWidth[m] = w; vHeight[m] = h;
		for (int l = 0; l < numLobes; l++)
		{
			vMFdata[m][l] = cv::Mat(h, w, CV_32FC4);
		}
		if (w > 1) w /= 2; if (h > 1) h /= 2;
	}
}
vMFtexture::~vMFtexture()
{
	delete[] vMFmaps;
	for (int m = 0; m < this->mipmapLevel; m++)
	{
		delete[] this->vMFdata[m];
	}
	delete[] this->vMFdata;
	delete[] vWidth; delete[] vHeight;
}


//Parameter Functions
void vMFtexture::generatevMFmaps()
{
	cv::Mat rawData = this->originalNormals[0].clone();
	cv::Mat **vMFdata = this->vMFdata;
	int mipmapLevel = this->mipmapLevel;
	int numLobes = this->numLobes;
	int width = this->oWidth; int height = this->oHeight;

	float **aux= new float*[numLobes];
	float *alpha = new float[numLobes];
	for (int l = 0; l < numLobes; l++)
	{
		aux[l] = new float[3];
	}

	for (int i = 0; i < rawData.rows; i++)
	{
		for (int j = 0; j < rawData.cols; j++)
		{
			cv::Vec3f ntemp = rawData.at<cv::Vec3f>(i, j);
			float nftemp[3] = { ntemp[0], ntemp[1], ntemp[2] };
			vectorFunc::normalize(nftemp);
			rawData.at<cv::Vec3f>(i, j) = cv::Vec3f(nftemp[0], nftemp[1], nftemp[2]);

		}
	}


	//Seed Level
	float seedAlpha = 1.f / numLobes;
	float seedAlpha2 = 1.f / numLobes;
	seedAlpha = 1.f;
	seedAlpha2 = 0.f;
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			cv::Vec4f temp;
			temp[0] = seedAlpha;
			temp[1] = seedAlpha*rawData.at<cv::Vec3f>(j, i)[0];
			temp[2] = seedAlpha*rawData.at<cv::Vec3f>(j, i)[1];
			temp[3] = seedAlpha*rawData.at<cv::Vec3f>(j, i)[2];
			vMFdata[0][0].at<cv::Vec4f>(j, i) = temp;
		}
	}
	for (int l = 1; l < numLobes; l++)
	{
		for (int i = 0; i < width; i++)
		{
			for (int j = 0; j < height; j++)
			{
				cv::Vec4f temp;
				temp[0] = seedAlpha2;
				temp[1] = seedAlpha2*rawData.at<cv::Vec3f>(j, i)[0];
				temp[2] = seedAlpha2*rawData.at<cv::Vec3f>(j, i)[1];
				temp[3] = seedAlpha2*rawData.at<cv::Vec3f>(j, i)[2];
				vMFdata[0][l].at<cv::Vec4f>(j, i) = temp;
			}
		}
	}

	for (int m = 1; m < mipmapLevel; m++)
	{
		width = this->vWidth[m]; height = this->vHeight[m];
		std::cout << "Proceeding mipmap level " << m << " generation (" << width << ", " << height << ")\n";
		int side = pow(2, m);

		for (int w = 0; w < width; w++)
		{
			for (int h = 0; h < height; h++)
			{
				float tKappa = { 0.f };
				float tR[3] = { 0.f };
				float prevData[4][20][4];
				for (int i = 0; i < 4; i++)
				{
					for (int l = 0; l < numLobes; l++)
					{

						cv::Vec4f temp = vMFdata[m - 1][l].at<cv::Vec4f>(h * 2 + (i % 2), w * 2 + (int)(i / 2));
						if (temp[0] != 0)
						{
							temp[1] /= temp[0]; temp[2] /= temp[0]; temp[3] /= temp[0];
							tR[0] = temp[1]; tR[1] = temp[2]; tR[2] = temp[3];
							tKappa = vMFfunc::r2kappa(tR);
							vectorFunc::normalize(tR);
						}
						else
						{
							tKappa = 0.f;
							tR[0] = 0.f; tR[1] = 0.f; tR[2] = 0.f;
						}

						prevData[i][l][0] = tR[0]; prevData[i][l][1] = tR[1]; prevData[i][l][2] = tR[2];
						prevData[i][l][3] = tKappa;
					}
				}
				cv::Rect ROI(w*side, h*side, side, side);
				cv::Mat targetRegion; 
				targetRegion = rawData(ROI).clone();

				//prevData for mu initialization
				this->computeParameters(alpha, aux, targetRegion, prevData);
				
//				if(m>6)
//				vMFfunc::displayvMF(numLobes, alpha, aux, 512, 512, 0, 0);
				//alignment between neighboring pixels of same mipmap level

				for (int l = 0; l < numLobes; l++)
				{

					cv::Vec4f computeData;
					computeData[0] = alpha[l];
					computeData[1] = alpha[l] * aux[l][0];
					computeData[2] = alpha[l] * aux[l][1];
					computeData[3] = alpha[l] * aux[l][2];
					vMFdata[m][l].at<cv::Vec4f>(h, w) = computeData;
				}
			}
		}
	}

	char* outName="out_";
	for (int m = 0; m < mipmapLevel; m++)
	{
		for (int l = 0; l < numLobes; l++)
		{

	
		}
	}

	for (int l = 0; l < numLobes; l++)
	{
		delete[] aux[l];
	}
	delete[] aux; delete[] alpha;
}

void vMFtexture::computeParameters(float *alpha, float **aux, cv::Mat targetRegion, float prevData[4][20][4])
{
	int iteration = 0;

	float mu[20][3], kappa[20];
	int numLobes = this->numLobes;
	int sideX = targetRegion.cols, sideY = targetRegion.rows;
	int area = sideX*sideY;
	float **z;

	z = new float*[numLobes];
	for (int j = 0; j < numLobes; j++)
	{
		z[j] = new float[area];
		for (int i = 0; i < area; i++)
		{
			z[j][i] = 0.f;
		}
	}


	//Initialization method1
	float align[20][3] = { 0.f };
	mu[0][0] = 0.0f; mu[0][1] = 0.0f; mu[0][2] = 1.0f; kappa[0] = 700.f;
	align[0][0] = 0.0f; align[0][1] = 0.0f; align[0][2] = 1.0f;
	for (int i = 0; i < numLobes-1; i++)
	{
		mu[i+1][0] = cos(i*(2 * vmfPI) / (numLobes - 1));
		mu[i+1][1] = sin(i*(2 * vmfPI) / (numLobes - 1));
		mu[i+1][2] = 0.f;
		kappa[i + 1] = 700.f;

		align[i + 1][0] = cos(i*(2 * vmfPI) / (numLobes - 1));
		align[i + 1][1] = sin(i*(2 * vmfPI) / (numLobes - 1));
		align[i + 1][2] = 0.f;
	}

	bool converge = false;
	while (!converge)
	{
		//E
		for (int row = 0; row < sideY; row++)
		{
			for (int col = 0; col < sideX; col++)
			{
				cv::Vec3f Vec3_ni = targetRegion.at<cv::Vec3f>(row, col);
				float float_ni[3] = { Vec3_ni[0], Vec3_ni[1], Vec3_ni[2] };
				double zsum = 0;
				double vMFzij[20];
				for (int j = 0; j < numLobes; j++)
				{
					vMFzij[j] = vMFfunc::vMF(float_ni, mu[j], kappa[j]);
					zsum += vMFzij[j];
				}
				if (zsum < 0.001)
				{
					zsum = 0.001;
				}
				for (int j = 0; j < numLobes; j++)
				{
					z[j][sideX*row + col] = vMFzij[j] / zsum;
					if (isinf<float>(zsum))
					{
						std::cout << "inf zsum detected\n";
						z[j][sideX*row + col] = 0.f;
					}

				}
			}
		}
		//E End


		//M
		float zj[20] = { 0.f };
		cv::Vec3f Vec3_znj[20] = { 0.f };
		for (int j = 0; j < numLobes; j++)
		{
			for (int row = 0; row < sideY; row++)
			{
				for (int col = 0; col < sideX; col++)
				{
					cv::Vec3f Vec3_ni = targetRegion.at<cv::Vec3f>(row, col);
					zj[j] += z[j][sideX*row + col];
					Vec3_znj[j] += z[j][sideX*row + col] * Vec3_ni;
				}
			}
		}

		for (int j = 0; j < numLobes; j++)
		{
			if (zj[j] < 0.01)
				zj[j] = 0.01;
		}
		cv::Vec3f Vec3_aux[20];
		//Alpha
		for (int j = 0; j < numLobes; j++)
		{

			alpha[j] = zj[j] / (float)area;
			Vec3_aux[j] = Vec3_znj[j] / zj[j];

			aux[j][0] = Vec3_aux[j][0];
			aux[j][1] = Vec3_aux[j][1];
			aux[j][2] = Vec3_aux[j][2];
			float lenAux = cv::norm(Vec3_aux[j]);
			kappa[j] = ((3 * lenAux) - (lenAux*lenAux*lenAux)) / (1 - lenAux*lenAux);
			if (kappa[j] < 0.f)
				kappa[j] = 700.f;
			Vec3_aux[j] = cv::normalize(Vec3_aux[j]);
			mu[j][0] = Vec3_aux[j][0];
			mu[j][1] = Vec3_aux[j][1];
			mu[j][2] = Vec3_aux[j][2];
			vectorFunc::normalize(mu[j]);
			float musqr = sqrt(mu[j][0] * mu[j][0] + mu[j][1] * mu[j][1] + mu[j][2] * mu[j][2]);

		}
		if (iteration++ == 100)
			converge = true;
	}
	for (int i = 0; i < numLobes; i++)
	{
		if (isnan<float>(alpha[i]))
			std::cout << alpha[i] << std::endl;
	}
//	std::cout << "alpha: " << alpha[0] << std::endl;
//	std::cout << "zj: " << zj[0] << std::endl;

	//Alignment
	for (int i = 0; i < numLobes; i++)
	{
		float val[20] = { 0.f };
		for (int j = i; j < numLobes; j++)
		{
			val[j] = vectorFunc::dotProd(align[i], aux[j]);
		}
		float highestVal = 0.f;
		int ind;
		for (int k = 0; k < numLobes; k++)
		{
			if (highestVal < val[k])
				highestVal = val[k];
			ind = k;
		}
		float temp[3] = { aux[ind][0], aux[ind][1], aux[ind][2] };
		aux[ind][0] = aux[i][0]; aux[ind][1] = aux[i][1]; aux[ind][2] = aux[i][2];
		aux[i][0] = temp[0]; aux[i][1] = temp[1]; aux[i][2] = temp[2];
	}
	for (int i = 0; i < numLobes; i++)
	{
		delete[] z[i];
	}
	delete[] z;
}

//display Functions
void vMFtexture::showOriginalImage(int channel) const
{
	cv::namedWindow("originalImage"); cv::namedWindow("originalImage", CV_WINDOW_NORMAL);
	char key = 0;
	std::vector<cv::Mat> zyx0(3);
	std::vector<cv::Mat> zyx1(3);
	cv::split(this->originalNormals[0].clone(), zyx0);
	cv::split(this->originalNormals[1].clone(), zyx1);

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
			{
				cv::imshow("originalImage", zyx0[channel]);
			}
			else
			{
				cv::imshow("originalImage", zyx1[channel]);
			}
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
			if ((key == 'v') || (key == 'V'))
			{
				ver == 0 ? ver = 1 : ver = 0;
			}
		}
	}


	cv::destroyWindow("originalImage");

	return;
}
void vMFtexture::showvMFImage(int level, int lobe, int mode) const
{
	cv::namedWindow("vMF"); cv::namedWindow("vMF", CV_WINDOW_NORMAL);
	cv::Mat **data = this->vMFdata;
	char key = 0;

	while (key != myESC)
	{
		cv::imshow("vMF", data[level][lobe]);
		key=cv::waitKey();

		switch (key)
		{
		case '.':case'>':
			if (level == mipmapLevel-1)
				level = 0;
			else
				level++;
			break;
		case ',':case'<':
			if (level == 0)
				level = mipmapLevel-1;
			else
				level--;
			break;
		case 'n':case'N':
			if (lobe == numLobes - 1)
				lobe = 0;
			else
				lobe++;
			break;
		case'b':case'B':
			if (lobe == 0)
				lobe = numLobes - 1;
			else
				lobe--;
			break;
		default:
			break;
		}
	}
	cv::destroyWindow("vMF");
}


//vMFfunc
void vMFfunc::displayvMF(int numLobes, float *alpha, float **aux, int width, int height, int skip, int destroy)
{

	float *a = alpha;
	float **r = aux;

	cv::Mat NMFdisplay = cv::Mat::zeros(width, height, CV_32FC3);

	float center[2] = { (float)height / 2, width / 2 };
	for (int w = 0; w < width; w++)
	{
		for (int h = 0; h < height; h++)
		{
			float val = 0;
			cv::Vec3f value(0.f, 0.f, 0.f);
			for (int i = 0; i < numLobes; i++)
			{
				float mu[3] = { r[i][0], r[i][1], r[i][2] };
				float kappa = 0;
				vectorFunc::normalize(mu);
				kappa = vMFfunc::r2kappa((float*)r[i]);
				float curPos[3] = { 2.f*(float)(w - center[0]) / (float)width,  2.f*(float)(h - center[1]) / (float)height, 0.f };
				if (((curPos[0] * curPos[0]) + (curPos[1] * curPos[1])) > 1.f)
				{
					value = cv::Vec3f(0.f, 30.f, 0.f);
					break;
				}
				else
				{
					curPos[2] = sqrt(1 - ((curPos[0] * curPos[0]) + (curPos[1] * curPos[1])));
					val+= alpha[i] * vMFfunc::vMF(curPos, mu, kappa);
					value = cv::Vec3f(val, val, val);
				}
			}
			NMFdisplay.at<cv::Vec3f>(h, w) = value;
			value = cv::Vec3f(0.f, 0.f, 0.f);
			val = 0.f;
		}
	}
//	normalize(NMFdisplay, NMFdisplay, 1.f, 0.f, cv::NORM_MINMAX);
	displayImage("vMF", NMFdisplay, skip, destroy);

}


void vMFfunc::mukappa2aux(float *aux, float mu[3], float kappa)
{
	double Kappa = kappa;
	float Ak = cosh(Kappa) / sinh(Kappa) - 1 / Kappa;
	aux[0] = mu[0] * Ak; aux[1] = mu[1] * Ak; aux[2] = mu[2] * Ak;
}

float vMFfunc::r2kappa(float r[3])
{
	float result = 0.f;
	float normR = vectorFunc::norm(r);
	result = ((3 * normR) - (normR*normR*normR)) / (1 - (normR*normR));
	if (normR == 1)
		result = 707;
	return result;
}
cv::Mat vMFfunc::cvLoadImage(const char* filename, int &imageWidth, int &imageHeight)
{
	cv::Mat load, floatScale;
	load = cv::imread(filename, CV_LOAD_IMAGE_COLOR);
	load.convertTo(floatScale, CV_32FC1, 1.f / 255.f);
	imageWidth = load.cols; imageHeight = load.rows;
	return floatScale;
}
double vMFfunc::vMF(float normal[3], float mu[3], float kappa) {
	float normalNorm = vectorFunc::norm(normal);
	if ((normalNorm > 1.f) || (normalNorm < 0.99999f))
		vectorFunc::normalize(normal);
	double NdotMu = (mu[0] * normal[0]) + (mu[1] * normal[1]) + (mu[2] * normal[2]);
	double Kappa = kappa;
	if (Kappa == 0)
		Kappa = 0.0001;
	if (Kappa > 700.)
		Kappa = 700.;
	double result = (Kappa / (4 * vmfPI*sinh(Kappa)))*exp(Kappa*(NdotMu));

	if (isinf<double>(result))
	{
		std::cout << "inf vMF value detected\n";

	}
	if (isnan<double>(result))
		std::cout << "wrong vMF value\n";

	return result;
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

//vectorFunc
void vectorFunc::normalize(float input[3])
{
	int counter = 0;
	float norm = vectorFunc::norm(input);
	if (norm == 0)
		return;
	input[0] /= norm;
	input[1] /= norm;
	input[2] /= norm;
	norm = vectorFunc::norm(input);
	while ((norm > 1.00000f)|| (norm < 0.999999f))
	{
		input[0] /= norm;
		input[1] /= norm;
		input[2] /= norm;
		norm = vectorFunc::norm(input);
	}
}
float vectorFunc::norm(float input[3])
{
	return sqrt((input[0] * input[0]) + (input[1] * input[1]) + (input[2] * input[2]));
}

float vectorFunc::dotProd(float a[3], float b[3])
{
	float ret;
	ret = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	return ret;
}

void clusterFunc::doKcluster(Cluster *clusters, int numClusters, std::vector<Sample*> samples)
{
	int numSamples = (int)samples.size();
	Cluster* convChecker = new Cluster[numClusters];

	bool converge = false;
	while (!converge)
	{
		for (int j = 0; j < numSamples; j++)
		{
			int fInd = 0;
			float fDistance = FLT_MAX;
			for (int i = 0; i < numClusters; i++)
			{
				float cent[3] = { clusters[i].centroid[0], clusters[i].centroid[1], clusters[i].centroid[2] };
				float pos[3] = { samples[j]->pos[0],samples[j]->pos[1], samples[j]->pos[2] };

				float distance = sqrt(((cent[0] - pos[0])*(cent[0] - pos[0]))
					+ ((cent[1] - pos[1])*(cent[1] - pos[1]))
					+ ((cent[2] - pos[2])*(cent[2] - pos[2])));
				if (distance < fDistance)
				{
					fInd = i;
					fDistance = distance;
				}
			}
			clusters[fInd].addSample(samples[j]);
		}

		//Check Convergence
		for (int i = 0; i < numClusters; i++)
		{
			int numSampleInClust = clusters[i].samples.size();
			int numChecker = convChecker[i].samples.size();
			if (numSampleInClust != numChecker)
				break;
			for (int j = 0; j < numSampleInClust; j++)
			{
				float newPos[3] = { clusters[i].samples[j]->pos[0],clusters[i].samples[j]->pos[1],clusters[i].samples[j]->pos[2] };
				float prePos[3] = { convChecker[i].samples[j]->pos[0],convChecker[i].samples[j]->pos[1],convChecker[i].samples[j]->pos[2] };
				if ((newPos[0] == prePos[0]) && (newPos[1] == prePos[1]) && (newPos[2] == prePos[2]))
				{
					if ((i == numClusters - 1) && (j== numSampleInClust - 1))
						converge = true;
				}
				else
					break;
			}
		}

		if (!converge)
		{
			for (int i = 0; i < numClusters; i++)
			{
				convChecker[i].clearSamplelist();
				convChecker[i] = clusters[i];
			}

			//Reset Centroids
			for (int i = 0; i < numClusters; i++)
			{
				float newCent[3] = { 0.f };
				int sampSize = clusters[i].samples.size();
				for (int j = 0; j < sampSize; j++)
				{
					newCent[0] += clusters[i].samples[j]->pos[0];
					newCent[1] += clusters[i].samples[j]->pos[1];
					newCent[2] += clusters[i].samples[j]->pos[2];
				}
				newCent[0] /= sampSize; newCent[1] /= sampSize; newCent[2] /= sampSize;
				clusters[i].setCentroid(newCent);
				clusters[i].clearSamplelist();
			}
		}

	}
}
