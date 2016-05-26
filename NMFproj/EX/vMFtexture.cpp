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

	resizeFactx = 256.f/(float)oWidth; resizeFacty = 256.f / (float)oHeight;

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
		std::cout << "Proceeding mipmap level " << m << " generagion (" << width << ", " << height << ")\n";
		int side = pow(2, m);

		for (int w = 0; w < width; w++)
		{
			for (int h = 0; h < height; h++)
			{
				float tKappa = { 0.f };
				float tR[3] = { 0.f };
				float prevData[4][10][4];
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
				targetRegion = originalNormals[0](ROI).clone();

				//prevData for mu initialization
				if (w % 100 == 0)
					std::cout << "(" << w << ", " << h << ")\n";
																																					this->computeParameters(alpha, aux, targetRegion, prevData);

				//alignment between neighboring pixels of same mipmap level

				for (int l = 0; l < numLobes; l++)
				{

					cv::Vec4f computeData;
					computeData[0] = alpha[l];
					computeData[1] = alpha[l]*aux[l][0];
					computeData[2] = alpha[l] * aux[l][1];
					computeData[3] = alpha[l] * aux[l][2];
					vMFdata[m][l].at<cv::Vec4f>(h, w) = computeData;
				}
			}
		}
	}


	for (int l = 0; l < numLobes; l++)
	{
		delete[] aux[l];
	}
	delete[] aux; delete[] alpha;
}

void vMFtexture::computeParameters(float *alpha, float **aux, cv::Mat targetRegion, float prevData[4][10][4])
{
	int iteration = 0;

	float mu[10][3], kappa[10];
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
//
//	//find suitable initial value
//	int numLobes2bUsed[4] = { 0.f };
//	int ind2bUsed[4][10];
//	std::vector<float*> validityStack;
//	for (int i = 0; i < 4; i++)
//	{
//		int count = 0;
//		for (int l = 0; l < this->numLobes; l++)
//		{
//			if (prevData[i][l][3] != 0)
//			{
//				bool valid = true;
//				for (int j = 0; j < validityStack.size(); j++)
//				{
//					if ((validityStack[j][0] == prevData[i][l][0])
//						&& (validityStack[j][1] == prevData[i][l][1])
//						&& (validityStack[j][2] == prevData[i][l][2]))
//						valid = false;
//				}
//				if (valid)
//				{
//					ind2bUsed[i][count] = l;
//					count++;
//					float *newStack = new float[3]{ prevData[i][l][0],prevData[i][l][1] ,prevData[i][l][2] };
//					validityStack.push_back(newStack);
//				}
//			}
//		}
//		numLobes2bUsed[i] = count;
//	}
//	for (int i = 0; i < validityStack.size(); i++)
//		delete[] validityStack[i];
//	validityStack.clear();
//
////	initGraph graph;
////	for (int i = 0; i < 4; i++)
////	{
////		for (int j = 0; j < numLobes2bUsed[i]; j++)
////		{
////			int indUsed = ind2bUsed[i][j];
////			int indkey[2] = { i, indUsed };
////			float pos[3] = { prevData[i][indUsed][0], prevData[i][indUsed][1],prevData[i][indUsed][2] };
////			Node *nodePtr = new Node(indkey, pos);
////			graph.addNode(nodePtr);
////		}
////	}
////	graph.makeComplete();
//	std::vector<Sample*> samples;
//	for (int i = 0; i < 4; i++)
//	{
//		for (int j = 0; j < numLobes2bUsed[i]; j++)
//		{
//			int indKey[2] = { i, ind2bUsed[i][j] };
//			float pos[3] = { prevData[i][indKey[1]][0], prevData[i][indKey[1]][1], prevData[i][indKey[1]][2] };
//			Sample* sample = new Sample(indKey, pos);
//			samples.push_back(sample);
//		}
//	}
//
//	//If numSample=1
//	int numSample = (int)samples.size();
//	int k = numLobes;
//	if (numSample < k)
//	{
//		for (int i = 0; i < numLobes; i++)
//		{
//			if (i < numSample) 
//			{
//				alpha[i] = 1.f;
//				aux[i][0] = samples[i]->pos[0];
//				aux[i][1] = samples[i]->pos[1];
//				aux[i][2] = samples[i]->pos[2];
//			}
//			else
//			{
//				alpha[i] = 0.f;
//				aux[i][0] = samples[numSample - 1]->pos[0];
//				aux[i][1] = samples[numSample - 1]->pos[1];
//				aux[i][2] = samples[numSample - 1]->pos[2];
//			}
//		}
//		return;
//	}
//
//	Cluster *kCluster = new Cluster[7];
//	std::vector<int> indFlag;
//	std::vector<float*> centroids;
//	//Randomly Select k samples
//	for (int i = 0; i < k; i++)
//	{
//		int ind;
//		bool found = true;
//		while (found)
//		{
//			std::cout << "check2";
//			ind = rand() % numSample;
//			for (int ii = 0; ii < indFlag.size(); ii++)
//			{
//				if (ind == indFlag[ii])
//				{
//					found = true;
//					ind = rand() % numSample;
//					ii = 0;
//				}
//				else if (ii == ( indFlag.size() - 1))
//					found = false;
//			}
//			indFlag.push_back(ind);
//			found = false;
//		}
//	
//		float tPos[3] = {
//			samples[ind]->pos[0], samples[ind]->pos[1], samples[ind]->pos[2] };
//		kCluster[i].addSample(samples[ind]);
//		kCluster[i].setCentroid(tPos);
//	}
//	
//	clusterFunc::doKcluster(kCluster, k, samples);


	//Initialization
	mu[0][0] = 0.0f; mu[0][1] = 0.0f; mu[0][2] = 1.0f; kappa[0] = 11.f;
	for (int i = 0; i < numLobes-1; i++)
	{
		mu[i+1][0] = cos(i*(2 * vmfPI) / (numLobes - 1));
		mu[i+1][1] = sin(i*(2 * vmfPI) / (numLobes - 1));
		mu[i+1][2] = 0.f;
		kappa[i + 1] = 1.f;
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
				float zsum = 0;
				for (int j = 0; j < numLobes; j++)
				{
					zsum += vMFfunc::vMF(float_ni, mu[j], kappa[j]);
				}
				if (zsum < 0.0001)
				{
//					std::cout << "too small zsum\n";
				}
				for (int j = 0; j < numLobes; j++)
				{
					z[j][sideX*col + col] = vMFfunc::vMF(float_ni, mu[j], kappa[j]) / zsum;
				}
			}
		}
		//E End

		//M
		float zj[10] = { 0.f };
		cv::Vec3f Vec3_znj[10] = { 0.f };
		for (int j = 0; j < numLobes; j++)
		{
			for (int row = 0; row < sideY; row++)
			{
				for (int col = 0; col < sideX; col++)
				{
					cv::Vec3f Vec3_ni = targetRegion.at<cv::Vec3f>(row, col);
					zj[j] += z[j][sideX*col + col];
					Vec3_znj[j] += z[j][sideX*col + col] * Vec3_ni;
				}
			}
		}

		cv::Vec3f Vec3_aux[10];
		//Alpha
		for (int j = 0; j < numLobes; j++)
		{

			alpha[j] = zj[j] / area;
			Vec3_aux[j] = Vec3_znj[j] / zj[j];
			aux[j][0] = Vec3_aux[j][0];
			aux[j][1] = Vec3_aux[j][1];
			aux[j][2] = Vec3_aux[j][2];
			float lenAux = vectorFunc::norm(aux[j]);
			kappa[j] = ((3 * lenAux) - (lenAux*lenAux*lenAux)) / ((1 - lenAux)*(1 - lenAux));
			mu[j][0] = aux[j][0];
			mu[j][1] = aux[j][1];
			mu[j][2] = aux[j][2];
			vectorFunc::normalize(mu[j]);
		}
		if (iteration++ == 10)
			converge = true;
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
void vMFfunc::displayvMF(int numLobes, float *alpha, float **aux, int width, int height)
{

	float *a = alpha;
	float **r = aux;

	cv::Mat NMFdisplay = cv::Mat::zeros(width, height, CV_32FC1);

	float center[2] = { (float)height / 2, width / 2 };
	for (int w = 0; w < width; w++)
	{
		for (int h = 0; h < height; h++)
		{
			for (int i = 0; i < numLobes; i++)
			{
				float mu[3] = { r[i][0], r[i][1], r[i][2] };
				float kappa = 0;
				vectorFunc::normalize(mu);
				kappa=vMFfunc::r2kappa((float*)r[i]);
				cv::Vec3f value = 0.f;
				float curPos[3] = {2.f*(float)(w - center[0])/(float)width,  2.f*(float)(h - center[1]) / (float)height, 0.f };
				if (((curPos[0] * curPos[0]) + (curPos[1] * curPos[1])) > 1.f)
					continue;
				else
				{
					if ((w == 256) && (h == 256))
						std::cout << "test";
					curPos[2] = sqrt(1 - ((curPos[0] * curPos[0]) + (curPos[1] * curPos[1])));
					NMFdisplay.at<float>(h, w) += alpha[i] * vMFfunc::vMF(curPos, mu, kappa);
				}
			}
		}
	}
	normalize(NMFdisplay, NMFdisplay, 1.f, 0.f, cv::NORM_MINMAX);
	displayImage("vMF", NMFdisplay, cvCtrl::skip::no, cvCtrl::destroy::yes);

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

//vectorFunc
void vectorFunc::normalize(float input[3])
{
	int counter = 0;
	float norm = vectorFunc::norm(input);
	while ((vectorFunc::norm(input) > 1.f)||(counter++>=4))
	{
		input[0] /= norm;
		input[1] /= norm;
		input[2] /= norm;
		norm = vectorFunc::norm(input);
	}
	if (vectorFunc::norm(input)>1.f)
	{
		std::cout << "\nWarning! Normalized vector is larger than 1\n";
	}
}
float vectorFunc::norm(float input[3])
{
	return sqrt((input[0] * input[0]) + (input[1] * input[1]) + (input[2] * input[2]));
}

//graphFunc
extern float graphFunc::calcWeight(Node* n1, Node* n2)
{
	float pos1[3] = { n1->pos[0],n1->pos[1],n1->pos[2] };
	float pos2[3] = { n2->pos[0],n2->pos[1],n2->pos[2] };

	float weight = sqrt(((pos1[0] - pos2[0])*(pos1[0] - pos2[0]))
		+ ((pos1[1] - pos2[1])*(pos1[1] - pos2[1]))
		+ ((pos1[2] - pos2[2])*(pos1[2] - pos2[2]))
		);
	return weight;
}

void initGraph::makeComplete()
{
	this->edgeList.clear();
	this->numEdges = 0;

	int numNodes = this->numNodes;
	for (int i = 0; i < numNodes; i++)
	{
		Node* curNode;
		curNode = (Node*)this->getNode(i);
		curNode->edges.clear();
	}

	std::vector<float> weightList;
	for (int i = 0; i < numNodes; i++)
	{
		Node* curNode;
		curNode = (Node*)this->getNode(i);
		std::vector<float> weightCur;
		for (int j = 0; j < numNodes; j++)
		{
			if (i == j)
				continue;
			Node* neighbor = (Node*)this->getNode(j);

			float weight=graphFunc::calcWeight(curNode, neighbor);

			Edge* edge = new Edge(neighbor, weight);

			curNode->addEdge(edge);
			this->addEdge(edge);
		}
		int numEdges = (int)curNode->edges.size();
		graphFunc::sortEdgeList(curNode->edges, 0, numEdges - 1, numEdges);
	}
	int numEdges = this->edgeList.size();
	graphFunc::sortEdgeList(this->edgeList, 0, numEdges - 1, numEdges);

	std::cout << "making complete graph done\n";
}

void graphFunc::sortEdgeList(std::vector<Edge*> &edgeList,int a, int b, int numEdges)
{

	int numCurrent = b - a + 1;

	if (numCurrent <= 1)
		return;

	int pivot = rand() % numCurrent + a;

	std::vector<Edge*> tempList;

	std::vector<Edge*>::iterator it;
	float pivotWeight = edgeList[pivot]->weight;
	int ii = 0;
	int jj = 0;


	tempList.assign(edgeList.begin(), edgeList.end());

	tempList.erase(tempList.begin()+a);
	tempList.insert(tempList.begin()+a, edgeList[pivot]);
	int p=0;
	for (int i = 0; i < numCurrent; i++)
	{
		if (i != pivot)
		{
			if (pivotWeight < edgeList[i+a]->weight)
			{
				it = tempList.begin()+a+ii;
				tempList.insert(it, edgeList[i+a]);
				it = tempList.begin() + b+1;
				tempList.erase(it);
				ii++; p++;
			}
			else
			{
				it = tempList.begin()+a + p + jj + 1;
				tempList.insert(it, edgeList[i+a]);
				it = tempList.begin() + b+1;
				tempList.erase(it);
				jj++;
			}
		}
	}
	edgeList.swap(tempList);

	int lb, le, rb, re;
	lb = a; le = ii-1; rb = p + 1; re = p + jj;

	graphFunc::sortEdgeList(edgeList, lb, le, numEdges);
	graphFunc::sortEdgeList(edgeList, rb, re, numEdges);

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
