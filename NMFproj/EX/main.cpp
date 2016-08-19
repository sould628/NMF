//kixor.net
//OpenGL Extension registry: http://www.opengl.org/registry
//Realtech VR's OpenGL Extensions Viewer
//#include <fstream>
#include <ctime>
#include <math.h>
#include <algorithm> //min max
#include <float.h> //isinf


#include <GL/glew.h>
#include <GL/glut.h>

#include <stdio.h>
#include <stdlib.h>

#include <map>
#include <string>
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#include "FreeImage.h"

#include "globalVariables.h"
#include "GLSLProgram.h"
#include "Camera.h"
#include "vMFtexture.h"
#include "SHtexture.h"
#include "readObj.h"


#define sysPause system("pause>nul");

objReader objreader;

const int l_index = 4;

int width = 800, height = 800;
int NMwidth, NMheight;


bool renewOBJshader = true;
//Mesh Data
int numVertices, numTexcoord, numNormals, numIndices, numFaces;
float *vertPos, *normals, *texCoord, *tangent, *bitangent;
unsigned int *indices;
float **data;
GLsizeiptr *dataSize;
int *vecSize;
GLuint elementBuffer;


//const char* NMT = "./bricks_normal_map.jpg";

FIBITMAP* NMTdata;
GLubyte* textureData;


GLint button_id;
GLfloat click_pos[2];

float t;

GLuint VAO, VBO[5], TexCoordArray;
GLuint NMbuffer, vMFvertex, vMFtex, vMFnormal, vMFtangent;
GLuint NormalMap, NormalMipMap, normalizedNMT;

GLuint SHmap0, SHmap1, SHmap2, SHmap3, SHmap4, SHmap5, SHmap6;

GLuint vMFmap0, vMFmap1, vMFmap2, vMFmap3, vMFmap4;

GLuint vMFmaps[10];
GLuint SHmaps[10];
GLuint vMFTexMap[10];
GLuint Ylmmaps[10];


GLSLProgram *NMFsh, *NMFvMF, *NMForiginal, *NMFvMFobj;

GLuint windowSH, windowvMF;

time_t sysTime, currentTime;

//vMF
float* alpha;
float** aux;

vMFtexture cVMFtex(NMT, numLobes, alignCtrl);
SHtexture cSHtex(NMT, 4);

static const float exampleData[] =
{
	0.25, -0.25, 0.5, 1.0,
	-0.25, 0.25, 0.5, 1.0,
	0.25, 0.25, 0.5, 1.0
};






void displayCB();
void setMesh();
void generateSHmap(GLubyte* TextureData, int order);
void generatevMFmap(GLubyte* TextureData, int numLobes, float* alpha, float* aux[3], float alignCtrl);

float calculateSHcoeffDelta(int l, int m, float x, float y, float z);
void checkTextureError(GLenum glError);
int vMFparam(float* data[3], float* prev[3], float* tAlpha, float* tAux[3], int numLobes, int mipmapLevel, int maxIteration, float alignCtrl);
float norm(float vector[3]);
void normalize(float* source, float* destination);

inline int mPower(int base, int power)
{
	if (power == 0)
		return 1;
	return mPower(base, power - 1)*base;
}

float norm(float* vector) {
	return sqrt((vector[0] * vector[0]) + (vector[1] * vector[1]) + (vector[2] * vector[2]));
}
void normalize(float* source, float* destination) {
	float value[3];
	value[0] = source[0] / norm(source);
	value[1] = source[1] / norm(source);
	value[2] = source[2] / norm(source);

	destination[0] = value[0];
	destination[1] = value[1];
	destination[2] = value[2];
}



inline float calculateKappa(float* aux)
{
	float kappa;
	kappa = ((3 * norm(aux)) - (norm(aux)*norm(aux)* norm(aux))) / (1 - (norm(aux)*norm(aux)));
	return kappa;
}

//To do:
//mu, kappa initial value should be given
int vMFparam2(float** data, float*** target, int curWidth, int curHeight, int curMipmapWidth, int curMipmapHeight, int dataWidth, int dataHeight, float* tAlpha, float* tAux[3], int numLobes, int mipmapLevel, int maxIteration, float alignCtrl=100.f) { //data=normalized float
	bool bConvergence = false;
	
	int numPixSide = mPower(2, mipmapLevel);
	int numPixel = numPixSide*numPixSide;

	float** mu;
	float* kappa;
	float** z;
	mu = new float*[numLobes];
	kappa = new float[numLobes];
	z = new float*[numLobes];

	float*** iaux;


	iaux = new float**[numLobes];
	int numData = 1;
	for (int i = 0; i < mipmapLevel; i++){
		numData *= 4;
	}
	float **ialpha = new float*[numLobes];

	for (int i = 0; i < numLobes; i++)
	{
		kappa[i] = 10.f;
		z[i] = new float[numData];
		mu[i] = new float[3];
		iaux[i] = new float*[4];
		ialpha[i] = new float[4];
		for (int j = 0; j < 4; j++)
		{
			iaux[i][j] = new float[3];
		}
	}


	
	for (int i = 0; i < numLobes; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			iaux[i][j][0] = target[mipmapLevel - 1][i][(curHeight*curMipmapWidth * 4 + curWidth * 4) * 2 + 1];
			iaux[i][j][1] = target[mipmapLevel - 1][i][(curHeight*curMipmapWidth * 4 + curWidth * 4) * 2 + 2];
			iaux[i][j][2] = target[mipmapLevel - 1][i][(curHeight*curMipmapWidth * 4 + curWidth * 4) * 2 + 3];
			ialpha[i][j] = target[mipmapLevel - 1][i][(curHeight*curMipmapWidth * 4 + curWidth * 4) * 2 + 0];
		}
		normalize(iaux[i][i], mu[i]);
	}



	//Initial Guess Stage (mu, kappa)


	float **targetNormal;
	targetNormal = new float*[numPixel];
	for (int j = 0; j < numPixSide; j++)//UpDown
	{
		for (int i = 0; i < numPixSide; i++)//LeftRight
		{
			targetNormal[j*numPixSide+i] = new float[3];
			targetNormal[j*numPixSide + i][0] = data[(((curHeight*numPixSide) + j)*dataWidth) + (curWidth*numPixSide + i)][0];
			targetNormal[j*numPixSide + i][1] = data[(((curHeight*numPixSide) + j)*dataWidth) + (curWidth*numPixSide + i)][1];
			targetNormal[j*numPixSide + i][2] = data[(((curHeight*numPixSide) + j)*dataWidth) + (curWidth*numPixSide + i)][2];
		}
	}
	//Initial Guess End

	int iteration = 0;
	//EM
	while (!bConvergence)
	{
		//E-STEP
		float* vMFij = new float[numLobes];
		for (int i = 0; i < numPixel; i++)
		{

			float vMFsum = 0.f;
			for (int j = 0; j < numLobes; j++)
			{
				vMFij[j] = vMFfunc::vMF(targetNormal[i], mu[j], kappa[j]);
				vMFsum += vMFij[j];
			}
			for (int j = 0; j < numLobes; j++)
			{
				z[j][i] = vMFij[j] / vMFsum;
			}
		}//E-Step End
		delete[] vMFij;
		float* zsum;
		zsum = new float[numLobes];
		float** expectedNormal;
		expectedNormal = new float*[numLobes];
		for (int jj = 0; jj < numLobes; jj++)
		{
			expectedNormal[jj] = new float[3];
		}
		for (int jj = 0; jj < numLobes; jj++)
		{
			zsum[jj] = 0;
			expectedNormal[jj][0] = 0;
			expectedNormal[jj][1] = 0;
			expectedNormal[jj][2] = 0;
			for (int ii = 0; ii < numData; ii++)
			{
				zsum[jj] += z[jj][ii];

				expectedNormal[jj][0] += z[jj][ii] * data[ii][0];
				expectedNormal[jj][1] += z[jj][ii] * data[ii][1];
				expectedNormal[jj][2] += z[jj][ii] * data[ii][2];
			}
		}
		//M-STEP

		for (int j = 0; j < numLobes; j++)
		{
			int jj = 0;

			tAlpha[j] = zsum[j] / numPixel;

			tAux[j][0] = expectedNormal[j][0] / zsum[j];
			tAux[j][1] = expectedNormal[j][1] / zsum[j];
			tAux[j][2] = expectedNormal[j][2] / zsum[j];

			if ((tAlpha[j] < 0) || (tAlpha[j]>1))
				std::cout << "wrong alpha detected\n";
			if (isinf(tAux[j][0])|| isinf(tAux[j][1])|| isinf(tAux[j][2]))
				std::cout << "inf value detected\n";

//			kappa[j] = (3 * norm(tAux[j]) - ((norm(tAux[j]))*(norm(tAux[j]))*(norm(tAux[j])))) / (1 - (norm(tAux[j]))*(norm(tAux[j])));
			kappa[j] = calculateKappa(tAux[j]);

			normalize(tAux[j], mu[j]);
////////////////////////////////////////////////////
			//Normalize with alignment
////////////////////////////////////////////////////
//			float alignMu[4][3];
//			float alignAlpha[4];
//			for (int qq = 0; qq < 4; qq++)
//			{
//				normalize(iaux[j][qq], alignMu[qq]);
//				alignAlpha[qq] = ialpha[j][qq];
//			}
//			float alignedAux[3] = { 0.f, 0.f, 0.f };
//
//			for (int pp = 0; pp < 3; pp++)
//			{
//				for (int qq = 0;  qq < 4; qq++)
//				{
//
//					alignedAux[pp] += alignAlpha[qq] * alignMu[qq][pp];
//				}
//				alignedAux[pp] *= alignCtrl;
//				alignedAux[pp] += tAux[j][pp];
//			}
//			normalize(tAux[j], mu[j]);

			//Normalize Alignment end

			
		}//M-Step End

		delete[] zsum;
		for (int jj = 0; jj < numLobes; jj++)
		{
			delete[] expectedNormal[jj];
		}
		delete[] expectedNormal;
		//ConvergenceTest
		iteration++;
		static int currentStars = 1;
		if (((float)iteration / (float)MAXITERATION) >(0.1*(float)currentStars))
		{
//			std::cout << "*";
			currentStars++;
		}
		if (iteration == MAXITERATION)
		{
//			std::cout << std::endl;
			bConvergence = true;
		}
	}//EM END

	for (int i = 0; i < numLobes; i++)
	{
		delete[] z[i];
		delete[] mu[i];
	}
	delete[] mu;
	delete[] z;
	return -1;
}


int vMFparam(float* data[3], float* prev[3], float* tAlpha, float* tAux[3], int numLobes, int mipmapLevel, int maxIteration, float alignCtrl) { //data=normalized float
	bool bConvergence = false;

	float** mu;
	float* kappa;
	float** z;
	mu = new float*[numLobes];
	kappa = new float[numLobes];
	z = new float*[numLobes];
	
	int numData = 1;
	for (int i = 0; i < mipmapLevel; i++){
		numData *= 4;
	}
	
	for (int i = 0; i < numLobes; i++)
	{
		z[i] = new float[numData];
		mu[i] = new float[3];
	}



	//Initial Guess Stage (mu, kappa)


	//Initial Guess End

	int iteration = 0;
	//EM
	while (!bConvergence)
	{
		//E-STEP
		float* vMFij = new float[numLobes];
		for (int i = 0; i < numData; i++)
		{

			float vMFsum = 0.f;
			for (int j = 0; j < numLobes; j++)
			{
				vMFij[j] = vMFfunc::vMF(data[i], mu[j], kappa[j]);
				vMFsum += vMFij[j];
			}
			for (int j = 0; j < numLobes; j++)
			{
				z[j][i] = vMFij[j] / vMFsum;
			}
		}//E-Step End
		delete[] vMFij;
		float* zsum;
		zsum = new float[numLobes];
		float** expectedNormal;
		expectedNormal = new float*[numLobes];
		for (int jj = 0; jj < numLobes; jj++)
		{
			expectedNormal[jj] = new float[3];
		}
		for (int jj = 0; jj < numLobes; jj++)
		{
			zsum[jj] = 0;
			expectedNormal[jj][0] = 0;
			expectedNormal[jj][1] = 0;
			expectedNormal[jj][2] = 0;
			for (int ii = 0; ii < numData; ii++)
			{
				zsum[jj] += z[jj][ii];

				expectedNormal[jj][0] += z[jj][ii] * data[ii][0];
				expectedNormal[jj][1] += z[jj][ii] * data[ii][1];
				expectedNormal[jj][2] += z[jj][ii] * data[ii][2];
			}
		}
		//M-STEP

		for (int j = 0; j < numLobes; j++)
		{
			int jj = 0;

			alpha[j] = zsum[j] / (NMwidth*NMheight);

			aux[j][0] = expectedNormal[j][0] / zsum[j];
			aux[j][1] = expectedNormal[j][1] / zsum[j];
			aux[j][2] = expectedNormal[j][2] / zsum[j];

			kappa[j] = (3 * norm(aux[j]) - ((norm(aux[j]))*(norm(aux[j]))*(norm(aux[j])))) / (1 - (norm(aux[j]))*(norm(aux[j])));

			normalize(aux[j], mu[j]);




		}//M-Step End

		delete[] zsum;
		for (int jj = 0; jj < numLobes; jj++)
		{
			delete[] expectedNormal[jj];
		}
		delete[] expectedNormal;
		//ConvergenceTest
		iteration++;
		static int currentStars = 1;
		if (((float)iteration / (float)MAXITERATION) > (0.1*(float)currentStars))
		{
			std::cout << "*";
			currentStars++;
		}
		if (iteration == MAXITERATION)
		{
			std::cout << std::endl;
			bConvergence = true;
		}
	}//EM END

	for (int i = 0; i < numLobes; i++)
	{
		delete[] z[i];
		delete[] mu[i];
	}
	delete[] mu;
	delete[] z;
	return -1;
}

//maxMipmaplevel ??
void generatevMFmap2(GLubyte* TextureData, int numLobes, float alignCtrl){
	int maxMipmapLevel = 0;

	float **normalData;
	normalData = new float*[NMwidth*NMheight];
	for (int i = 0; i < NMwidth*NMheight; i++)
	{
		normalData[i] = new float[3];
	}
	float *prevMu;
	int mipmapWidth = NMwidth;
	int mipmapHeight = NMheight;
	bool bminWidth = false; bool bminHeight = false;
	prevMu = new float[NMwidth*NMheight * numLobes * 4];
	float nn[3] = { 0 };
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 0);
	glGenTextures(1, &vMFmaps[0]); glGenTextures(1, &vMFmaps[1]); glGenTextures(1, &vMFmaps[2]); glGenTextures(1, &vMFmaps[3]); glGenTextures(1, &vMFmaps[4]);
	glGenTextures(1, &vMFmap0);	glGenTextures(1, &vMFmap1);	glGenTextures(1, &vMFmap2);	glGenTextures(1, &vMFmap3); glGenTextures(1, &vMFmap4);
	int tw=NMwidth, th=NMheight;
	while ((tw > 1) || (th > 1))
	{
		if (tw > 1)
			tw /= 2;
		if (th > 1)
			th /= 2;
		maxMipmapLevel++;
	}
	maxMipmapLevel = 6;

	//NormalData
	for (int j = 0; j < NMheight; j++)
	{
		for (int i = 0; i < NMwidth; i++)
		{
			nn[0] = (float)TextureData[j*NMwidth * 4 + i * 4 + 2];
			nn[1] = (float)TextureData[j*NMwidth * 4 + i * 4 + 1];
			nn[2] = (float)TextureData[j*NMwidth * 4 + i * 4 + 0];

			nn[0] = (nn[0] - (255.f / 2.f)) / (255.f / 2.f);
			nn[1] = (nn[1] - (255.f / 2.f)) / (255.f / 2.f);
			nn[2] = (nn[2] - (255.f / 2.f)) / (255.f / 2.f);
						if ((nn[1] * nn[1] + (nn[2] * nn[2]) + (nn[0] * nn[0])) > 1.f)
						{
							normalize(nn, nn);
						}
						if ((nn[1] * nn[1] + (nn[2] * nn[2]) + (nn[0] * nn[0])) > 1.00001f)
						{
							std::cout << "normal value error\n" << std::endl;
							sysPause;
						}
			normalData[j*NMwidth + i][0] = nn[0];
			normalData[j*NMwidth + i][1] = nn[1];
			normalData[j*NMwidth + i][2] = nn[2];
		}//i
	}//j
	//NormalData

	float*** dataBuffer;
	dataBuffer = new float**[maxMipmapLevel+1];
	for (int i = 0; i <= maxMipmapLevel; i++)
	{
		dataBuffer[i] = new float*[numLobes];
		int tWidth, tHeight;
		tWidth = NMwidth; tHeight = NMheight;
		for (int j = 0; j < numLobes; j++)
		{
			dataBuffer[i][j] = new float[tHeight*tWidth * 4];
		}
		if (tHeight>1)tHeight /= 2;
		if (tWidth>1)tWidth /= 2;
	}


	for (int ii = 0; ii <= maxMipmapLevel; ii++)
	{
		switch (ii)
		{
		case 0://MIPMAP Level0
		{
			for (int jj = 0; jj < numLobes; jj++)
			{
				switch (jj)//numLobes
				{
				case 0:
				{
					glActiveTexture(GL_TEXTURE1);
					for (int j = 0; j < NMheight; j++)//j
					{
						for (int i = 0; i < NMwidth; i++)//i
						{
							float aaalpha = 1.f / (float)numLobes;
							dataBuffer[ii][jj][j*NMwidth * 4 + i * 4 + 0] = aaalpha;
							dataBuffer[ii][jj][j*NMwidth * 4 + i * 4 + 1] = aaalpha*normalData[j*NMwidth + i][0];
							dataBuffer[ii][jj][j*NMwidth * 4 + i * 4 + 2] = aaalpha*normalData[j*NMwidth + i][1];
							dataBuffer[ii][jj][j*NMwidth * 4 + i * 4 + 3] = aaalpha*normalData[j*NMwidth + i][2];
						}//i
					}//j
					glBindTexture(GL_TEXTURE_2D, vMFmaps[jj]);
					glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
					//					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
					//					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, maxMipmapLevel+1);
					//					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
					glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NMwidth, NMheight, 0, GL_RGBA, GL_FLOAT, dataBuffer[ii][jj]);
					GLenum glError = glGetError();
					checkTextureError(glError);
					glActiveTexture(GL_TEXTURE10);
					glBindTexture(GL_TEXTURE_2D, 0);
					break;
				}//jj case 0
				default:
				{
					switch (jj)
					{
					case 1:
						glActiveTexture(GL_TEXTURE2);
						break;
					case 2:
						glActiveTexture(GL_TEXTURE3);
						break;
					case 3:
						glActiveTexture(GL_TEXTURE4);
						break;
					case 4:
						glActiveTexture(GL_TEXTURE5);
					}

					for (int j = 0; j < NMheight; j++)//j
					{
						for (int i = 0; i < NMwidth; i++)//i
						{
							float aaalpha = 1.f / (float)numLobes;
							dataBuffer[ii][jj][j*NMwidth * 4 + i * 4 + 0] = aaalpha;
							//							dataBuffer[ii][jj][j*NMwidth * 4 + i * 4 + 1] = normalData[j*NMwidth + i][0];
							//							dataBuffer[ii][jj][j*NMwidth * 4 + i * 4 + 2] = normalData[j*NMwidth + i][1];
							//							dataBuffer[ii][jj][j*NMwidth * 4 + i * 4 + 3] = normalData[j*NMwidth + i][2];
							dataBuffer[ii][jj][j*NMwidth * 4 + i * 4 + 1] = aaalpha*normalData[j*NMwidth + i][0];
							dataBuffer[ii][jj][j*NMwidth * 4 + i * 4 + 2] = aaalpha*normalData[j*NMwidth + i][1];
							dataBuffer[ii][jj][j*NMwidth * 4 + i * 4 + 3] = aaalpha*normalData[j*NMwidth + i][2];
						}//i
					}//j
					glBindTexture(GL_TEXTURE_2D, vMFmaps[jj]);
					glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
					//					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
					//					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, maxMipmapLevel+1);
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
					//					//glTexImage2D(GL_TEXTURE_2D, ii, GL_RGBA, NMwidth, NMheight, 0, GL_RGBA, GL_FLOAT, dataBuffer[ii][jj]);
					glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NMwidth, NMheight, 0, GL_RGBA, GL_FLOAT, dataBuffer[ii][jj]);
					GLenum glError = glGetError();
					checkTextureError(glError);
					glActiveTexture(GL_TEXTURE10);
					glBindTexture(GL_TEXTURE_2D, 0);

					break;
				}//default

				}//switch jj-1
			}//for jj (numLobes)

			break;
		}//ii case 0

		//////////////////////////////////////
		default://ii >0
		{
			if (mipmapHeight > 1)
				mipmapHeight /= 2;
			if (mipmapWidth > 1)
				mipmapWidth /= 2;
			std::cout << "\nmipmap " << ii << " (w,h) = (" << mipmapWidth << ", " << mipmapHeight << ")" << std::endl;


			float* tAlpha, **tAux;
			tAlpha = new float[numLobes];
			tAux = new float*[numLobes];
			for (int zz = 0; zz < numLobes; zz++)
			{
				tAux[zz] = new float[3];
			}
			for (int j = 0; j < mipmapHeight; j++)//j
			{
				for (int i = 0; i < mipmapWidth; i++)//i
				{
					vMFparam2(normalData, dataBuffer, i, j, mipmapWidth, mipmapHeight, NMwidth, NMheight, tAlpha, tAux, numLobes, ii, MAXITERATION);

					for (int jj = 0; jj < numLobes; jj++)
					{
						dataBuffer[ii][jj][j*mipmapWidth * 4 + i * 4 + 0] = tAlpha[jj];
						dataBuffer[ii][jj][j*mipmapWidth * 4 + i * 4 + 1] = tAlpha[jj] * tAux[jj][0];
						dataBuffer[ii][jj][j*mipmapWidth * 4 + i * 4 + 2] = tAlpha[jj] * tAux[jj][1];
						dataBuffer[ii][jj][j*mipmapWidth * 4 + i * 4 + 3] = tAlpha[jj] * tAux[jj][2];
					}//i
				}//j
			}//switch jj-1
			for (int jj = 0; jj < numLobes; jj++)
			{
				switch (jj)//numLobes
				{
				case 0:
					glActiveTexture(GL_TEXTURE1);
					break;
				case 1:
					glActiveTexture(GL_TEXTURE2);
					break;
				case 2:
					glActiveTexture(GL_TEXTURE3);
					break;
				case 3:
					glActiveTexture(GL_TEXTURE4);
					break;
				case 4:
					glActiveTexture(GL_TEXTURE5);
					break;
				}
				glBindTexture(GL_TEXTURE_2D, vMFmaps[jj]);
				//				glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
				//				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
				//				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, maxMipmapLevel);
				//				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
				//				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
				glTexImage2D(GL_TEXTURE_2D, ii, GL_RGBA, mipmapWidth, mipmapHeight, 0, GL_RGBA, GL_FLOAT, dataBuffer[ii][jj]);
				GLenum glError = glGetError();
				checkTextureError(glError);
				glActiveTexture(GL_TEXTURE10);
				glBindTexture(GL_TEXTURE_2D, 0);
			}
		}//for jj (numLobes)
		}

		//switch ii case

	}//for ii
	delete[] dataBuffer;
	delete[] prevMu;

	glDisable(GL_TEXTURE_2D);

}




void generatevMFmap(GLubyte* TextureData, int numLobes, float* alpha, float* aux[3], float alignCtrl){
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 0);
	glGenTextures(10, vMFmaps);
	int nLobes = cVMFtex.getNumLobes();
	int mmlevel = cVMFtex.getMipmapLevel();
	
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, NormalMap);
	glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, cWidth, cHeight, 0, GL_RGBA, GL_FLOAT, vMFtextureData);
	GLenum glError = glGetError();
	checkTextureError(glError);



	for (int l = 0; l < nLobes; l++)
	{
		std::cout << "binding texture of lobe " << l << std::endl;
		for (int m = 0; m < mmlevel; m++)
		{
			int cWidth = cVMFtex.getWidth(m); int cHeight = cVMFtex.getHeight(m);
			float *vMFtextureData = new float[4 * cWidth*cHeight];

			for (int j = 0; j < cHeight; j++)
			{
				for (int i = 0; i < cWidth; i++)
				{
					vMFtextureData[4 * (j*cWidth + i) + 0] = cVMFtex.getvMFcompoenent(m, l, j, i, 0);
					vMFtextureData[4 * (j*cWidth + i) + 1] = cVMFtex.getvMFcompoenent(m, l, j, i, 1);
					vMFtextureData[4 * (j*cWidth + i) + 2] = cVMFtex.getvMFcompoenent(m, l, j, i, 2);
					vMFtextureData[4 * (j*cWidth + i) + 3] = cVMFtex.getvMFcompoenent(m, l, j, i, 3);
				}
			}

			switch (l){
			case 0: glActiveTexture(GL_TEXTURE1); break;
			case 1: glActiveTexture(GL_TEXTURE2); break;
			case 2: glActiveTexture(GL_TEXTURE3); break;
			case 3: glActiveTexture(GL_TEXTURE4); break;
			case 4: glActiveTexture(GL_TEXTURE5); break;
			case 5: glActiveTexture(GL_TEXTURE6); break;
			case 6: glActiveTexture(GL_TEXTURE7); break;
			case 7: glActiveTexture(GL_TEXTURE8); break;
			case 8: glActiveTexture(GL_TEXTURE9); break;
			}

			if (m == 0)
			{
				glBindTexture(GL_TEXTURE_2D, vMFmaps[l]);
				glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			}
			glTexImage2D(GL_TEXTURE_2D, m, GL_RGBA, cWidth, cHeight, 0, GL_RGBA, GL_FLOAT, vMFtextureData);
			GLenum glError = glGetError();
			checkTextureError(glError);



			delete vMFtextureData;
		}
	}
	glActiveTexture(GL_TEXTURE10);
	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);

}


void initSharedMem(){
	NMwidth = 0; NMheight = 0;
	time(&sysTime);
}

float calculateSHcoeffDelta(int l, int m, float x, float y, float z){
	//Real Spherical Harmonics Function
	//https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
	//Y'l,m=spherical harmonics basis function
	//Yl,m=i*sqrt(1/2)(Y'l,m-Y'l,-m) (if m=negative with even absolute)
	//Yl,m=i*sqrt(1/2)(Y'l,m+Y'l,-m) (if m=negative with odd absolute)
	//Yl,m=sqrt(1/2)(Y'l,-m-Y'l,m) (if m=positive with odd absolute)
	//Yl,m=sqrt(1/2)(Y'l,-m+Y'l,m) (if m=positive with even absolute)

	if ((l < 0) || (l>4)){
		std::cout << "Wrong Value of l\n";
		sysPause;
		return -1;
	}

	//x^2+y^2+z^2=1
	float result=0.0;
	switch (l){
	case 0:
	{
		if (m != 0){
			std::cout << "Value of m isn't propriate\n\n";
			sysPause;
			return -1;
		}
		//Y0,0=1/2*(sqrt(1/pi))=0.282095
		return 0.282095;
	}
	case 1:
	{
		//(Y1,-1;Y1,0,Y1,1)=0.4886025119*(y;z;x)
		switch (m)
		{
		case -1: return 0.4886025119*y;
		case 0:	return 0.4886025119*z;
		case 1:	return 0.4886025119*x;
		default:
			std::cout << "Value of m isn't propriate\n\n";
			sysPause;
			return -1;
		}
	}
	case 2:
	{
		//(Y2,-2;Y2,-1,Y2,1)=1.09254843059*(xy;yz;zx)
		//(Y2,0)=1.09254843059*(3*z^2-1)
		//(Y2,2)=0.54627421529*(x^2-y^2)/
		switch (m)
		{
		case -2: return  1.09254843059*x*y;
		case -1: return 1.09254843059*y*z;
		case 0: return 0.31539156525*((z*z) - 1);
		case 1: return 1.09254843059*z*x;
		case 2: return 0.54627421529*((x*x) - (y*y));
		default:
			std::cout << "Value of m isn't propriate\n\n";
			sysPause;
			return -1;
		}
	}
	case 3:
	{
		//(Y3,-3;Y3,3)=0.59004358992*(3*x^2*y-y^3;x^3-3*x*y^2)
		//(Y3,-2)=2.89061144264*(x*y*z)
		//(Y3,-1;Y3,1)=0.45704579946*(y*(4*z^2-x^2-y^2);x*(4*z^2-x^2-y^2))
		//(Y3,0)=0.37317633259*(z*(2*z^2-3*x^2-3*y^2))
		//(Y3,2)=1.44530572132*(x^2*z-y^2*z)
		
		switch (m)
		{
		case -3: return 0.59004358992*((3*x*x*y)-(y*y*y));
		case -2: return 2.89061144264*(x*y*z);
		case -1: return 0.45704579946*(y*((4*z*z)-(x*x)-(y*y)));
		case 0: return 0.37317633259*(z*((2*z*z) - (3*x*x) - (3*y*y)));
		case 1: return 0.45704579946*(x*((4*z*z) - (x*x) - (y*y)));
		case 2: return 1.44530572132*((x*x*z) - (y*y*z));
		case 3: return 0.59004358992*((x*x*x) - (3*x*y*y));
		default:
			std::cout << "Value of m isn't propriate\n\n";
			sysPause;
			return -1;
		}
	}
	case 4:
	{
		//(Y4,-4)=2.5033429418*xy(x^2-y*2)
		//(Y4,-3;Y4,3)=1.77013076978*((3*x^2-y^2)*yz);(x^2-3*y^2)*xz))
		//(Y4,-2)=0.94617469575*(xy*(7*z^2-1))
		//(Y4,-1,Y4,1)=0.66904654355*((yz*(7*z^2-3));(xz*(7*z^2-3))
		//(Y4,0)=0.10578554691*(35*z^4-30^z^2+3)
		//(Y4,2)=0.47308734787*(x^2-y^2)*(7*z^2-1)
		//(Y4,4)=0.62583573544*(x^2*(x^2-3*y^2)-y^2*(3*x^2-y^2))

		switch (m)
		{
		case -4: return 2.5033429418*(x*y*((x*x) - (y*y)));
		case -3: return 1.77013076978*(((3*x*x) - (y*y))*y*z);
		case -2: return 0.94617469575*(x*y*((7*z*z) - 1));
		case -1: return 0.66904654355*(y*z*((7*z*z) - 3));
		case 0: return 0.10578554691*((35 * z*z*z*z) - (30 * z*z) + 3);
		case 1: return 0.66904654355*(x*z*((7*z*z) - 3));
		case 2: return 0.47308734787*((x*x) - (y*y))*((7*z*z) - 1);
		case 3: return 1.77013076978*(((x*x) - (3*y*y))*x*z);
		case 4: return 0.62583573544*(((x*x)*((x*x) - (3*y*y))) - ((y*y)*((3*x*x) - (y*y))));
		default:
			std::cout << "Value of m isn't propriate\n\n";
			sysPause;
			return -1;
		}
	}
	}

	return result;
}

//void generateSHmap(GLubyte* TextureData, int order){
//	float nn[3] = { 0 };
//	int textureCount = 0;
//	GLfloat* coeffTexture;
//	coeffTexture = new GLfloat[NMwidth*NMheight * 4];
//	for (int ll = 0; ll < l_index; ll++)
//	{
//		for (int mm = 0; mm < ll * 2 + 1; mm++)
//		{
//			for (int i = 0; i < NMwidth; i++)
//			{
//				for (int j = 0; j < NMheight; j++)
//				{
//					nn[0] = (float)TextureData[i*NMheight * 4 + j * 4];
//					nn[1] = (float)TextureData[i*NMheight * 4 + j * 4+1];
//					nn[2] = (float)TextureData[i*NMheight * 4 + j * 4+2];
//
//					nn[0] = (nn[0] - 127.f) / 128.f;
//					nn[1] = (nn[1] - 127.f) / 128.f;
//					nn[2] = (nn[2] - 127.f) / 128.f;
//					if ((nn[2] * nn[2] + nn[0] * nn[0])<=1)
//						nn[1] = sqrt(1 - (nn[2] * nn[2] + nn[0] * nn[0]));
//					else
//					{
//						if (nn[2] == 1)
//						{
//							nn[0] = 0.f; nn[1] = 0.f;
//						}
//						else if (nn[0] == 1)
//						{
//							nn[1] = 0.f; nn[2] = 0.f;
//						}
//						else{
//							while ((nn[2] * nn[2] + nn[0] * nn[0]) > 1)
//							{
//								nn[2] = nn[2] * 0.99; nn[0]=nn[0] * 0.99;
//							}
//							nn[1] = sqrt(1 - (nn[2] * nn[2] + nn[0] * nn[0]));
//						}
//					}
//
//					
//					coeffTexture[i*NMheight*4+j*4+textureCount]=calculateSHcoeffDelta(ll, (mm - ll), nn[0], nn[1], nn[2]);
//
//				}
//			}
//			if (textureCount++ == 3)
//			{
//
//				textureCount = 0;
//			}
//		}
//
//	}
//
//	delete[] coeffTexture;
//	coeffTexture = NULL;
//}
//


void drawCube(){

	glColor4f(1, 0, 0, 0.5);
	glBegin(GL_POLYGON);
	glVertex3f(-0.5, -0.5, 0.5); glVertex3f(0.5f, -0.5f, 0.5f); glVertex3f(0.5f, 0.5f, 0.5f); glVertex3f(-0.5f, 0.5f, 0.5f);
	glEnd();
	glColor4f(0, 1, 0, 0.5);
	glBegin(GL_POLYGON);
	glVertex3f(-0.5f, -0.5f, 0.5f); glVertex3f(-0.5f, 0.5f, 0.5f); glVertex3f(-0.5f, 0.5f, -0.5f); glVertex3f(-0.5f, -0.5f, -0.5f);
	glEnd();
	glColor4f(0, 0, 1, 0.5);
	glBegin(GL_POLYGON);
	glVertex3f(-0.5f, -0.5f, -0.5f); glVertex3f(-0.5f, 0.5f, -0.5f); glVertex3f(0.5f, 0.5f, -0.5f); glVertex3f(0.5f, -0.5f, -0.5f);
	glEnd();
	glColor4f(1, 1, 0, 0.5);
	glBegin(GL_POLYGON);
	glVertex3f(0.5f, -0.5f, -0.5f); glVertex3f(0.5f, 0.5f, -0.5f); glVertex3f(0.5f, 0.5f, 0.5f); glVertex3f(0.5f, -0.5f, +0.5f);
	glEnd();
	glColor4f(0, 1, 1, 0.5);
	glBegin(GL_POLYGON);
	glVertex3f(-0.5f, -0.5f, -0.5f); glVertex3f(-0.5f, -0.5f, 0.5f); glVertex3f(0.5f, -0.5f, 0.5f); glVertex3f(0.5f, -0.5f, -0.5f);
	glEnd();
	glColor4f(1, 0, 1, 0.5);
	glBegin(GL_POLYGON);
	glVertex3f(0.5f, 0.5f, -0.5f); glVertex3f(0.5f, 0.5f, 0.5f); glVertex3f(-0.5f, 0.5f, 0.5f); glVertex3f(-0.5f, 0.5f, -0.5f);
	glEnd();
	glColor3f(1, 1, 1);
}

void drawCube2(){
	glEnable(GL_TEXTURE_2D);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, NormalMap);

	GLenum glError = glGetError();
	
	checkTextureError(glError);
	glBegin(GL_POLYGON);
	glTexCoord2f(0.0, 1.0); glVertex3f(-0.5, -0.5, 0.5);
	glTexCoord2f(1.0, 1.0);	 glVertex3f(0.5f, -0.5f, 0.5f);
	glTexCoord2f(1.0, 0.0); glVertex3f(0.5f, 0.5f, 0.5f);
	glTexCoord2f(0.0, 0.0); glVertex3f(-0.5f, 0.5f, 0.5f);
	glEnd();

	glBegin(GL_POLYGON);
	glTexCoord2f(0.0, 1.0); glVertex3f(-0.5f, -0.5f, 0.5f);
	glTexCoord2f(1.0, 1.0);	glVertex3f(-0.5f, 0.5f, 0.5f);
	glTexCoord2f(1.0, 0.0); glVertex3f(-0.5f, 0.5f, -0.5f);
	glTexCoord2f(0.0, 0.0); glVertex3f(-0.5f, -0.5f, -0.5f);
	glEnd();

	glBegin(GL_POLYGON);
	glTexCoord2f(0.0, 1.0); glVertex3f(-0.5f, -0.5f, -0.5f);
	glTexCoord2f(1.0, 1.0);	glVertex3f(-0.5f, 0.5f, -0.5f);
	glTexCoord2f(1.0, 0.0); glVertex3f(0.5f, 0.5f, -0.5f);
	glTexCoord2f(0.0, 0.0); glVertex3f(0.5f, -0.5f, -0.5f);
	glEnd();

	glBegin(GL_POLYGON);
	glTexCoord2f(0.0, 1.0); glVertex3f(0.5f, -0.5f, -0.5f);
	glTexCoord2f(1.0, 1.0);	glVertex3f(0.5f, 0.5f, -0.5f);
	glTexCoord2f(1.0, 0.0); glVertex3f(0.5f, 0.5f, 0.5f);
	glTexCoord2f(0.0, 0.0); glVertex3f(0.5f, -0.5f, +0.5f);
	glEnd();

	glBegin(GL_POLYGON);
	glTexCoord2f(0.0, 1.0); glVertex3f(-0.5f, -0.5f, -0.5f);
	glTexCoord2f(1.0, 1.0);	glVertex3f(-0.5f, -0.5f, 0.5f);
	glTexCoord2f(1.0, 0.0); glVertex3f(0.5f, -0.5f, 0.5f);
	glTexCoord2f(0.0, 0.0); glVertex3f(0.5f, -0.5f, -0.5f);
	glEnd();

	glBegin(GL_POLYGON);
	glTexCoord2f(0.0, 1.0); glVertex3f(0.5f, 0.5f, -0.5f);
	glTexCoord2f(1.0, 1.0);	glVertex3f(0.5f, 0.5f, 0.5f);
	glTexCoord2f(1.0, 0.0); glVertex3f(-0.5f, 0.5f, 0.5f);
	glTexCoord2f(0.0, 0.0); glVertex3f(-0.5f, 0.5f, -0.5f);
	glEnd();
	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);
}



void checkTextureError(GLenum glError){
	if (glError)
	{
		std::cout << "There was an error loading the texture: " << NMT << std::endl;

		switch (glError)
		{
		case GL_INVALID_ENUM:
			std::cout << "Invalid enum." << std::endl;
			break;

		case GL_INVALID_VALUE:
			std::cout << "Invalid value." << std::endl;
			break;

		case GL_INVALID_OPERATION:
			std::cout << "Invalid operation." << std::endl;

		default:
			std::cout << "Unrecognised GLenum." << std::endl;
			break;
		}

		std::cout << "See https://www.opengl.org/sdk/docs/man/html/glTexImage2D.xhtml for further details." << std::endl;
	}

}

#pragma region GLSL
void createProgram() {
	NMFsh = new GLSLProgram("NMF_SH.vert", "NMF_SH.frag");
	NMFvMF = new GLSLProgram("NMF_vMF.vert", "NMF_vMF.frag");
	NMForiginal = new GLSLProgram("NMF_original.vert", "NMF_original.frag");
	NMFvMFobj = new GLSLProgram("NMF_vMF_OBJ.vert", "NMF_vMF_OBJ.frag");
}
void initBuffers() {
	glEnable(GL_TEXTURE_2D);
//	glCreateVertexArrays(1, &VAO);
//	glCreateVertexArrays(1, &TexCoordArray);
	textureData = FreeImage_GetBits(NMTdata);

	glGenTextures(1, &NormalMap);
	glGenTextures(1, &NormalMipMap);

	glDisable(GL_TEXTURE_2D);
	RGBQUAD pixel;
	for (int i = 0; i < NMheight; i++){
		bool found = false;
		for (int j = 0; j < NMwidth; j++){
			FreeImage_GetPixelColor(NMTdata, i, j, &pixel);
			if ((int)pixel.rgbBlue == 130)
			{
				std::cout << "Found rgb=255\nvalue of red = " << (int)pixel.rgbRed << ", value of green= " << (int)pixel.rgbGreen 
					<<", value of blue= " << (int)pixel.rgbBlue << std::endl;
				found = true;
				break;
			}
		}
		if (found)
			break;
	}

	std::cout << "Texture Loaded!\n" << std::endl;


	//vMFparameter
	alpha = new float[numLobes];
	aux = new float*[numLobes];
	for (int i = 0; i < numLobes; i++)
	{
		aux[i] = new float[3];
	}

	///bind texture
	generatevMFmap(textureData, numLobes, alpha, aux, alignCtrl);
//	generatevMFmap2(textureData, numLobes, alignCtrl);
	cSHtex.bindTexture(normalizedNMT, SHmaps);
	cSHtex.bindYlm(Ylmmaps);

	glEnable(GL_TEXTURE_2D);
	glActiveTexture(GL_TEXTURE0);

	glBindTexture(GL_TEXTURE_2D, NormalMap);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NMwidth, NMheight, 0, GL_BGRA, GL_UNSIGNED_BYTE, textureData);
//	gluBuild2DMipmaps(GL_TEXTURE_2D, 4, NMwidth, NMheight, GL_BGRA, GL_UNSIGNED_BYTE, textureData);

	glActiveTexture(GL_TEXTURE10);
	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);

	glEnable(GL_TEXTURE_2D);
	glActiveTexture(GL_TEXTURE10);
	glBindTexture(GL_TEXTURE_2D, NormalMipMap);
	glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NMwidth, NMheight, 0, GL_BGRA, GL_UNSIGNED_BYTE, textureData);

	glActiveTexture(GL_TEXTURE10);
	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);

}


#pragma endregion





void displayCB(){
	glEnable(GL_DEPTH_TEST);
	glClearColor(0.f, 0.f, 0.0f, 1.f);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	
	GLfloat color[] = { (float)sin(t)*0.5f, (float)cos(t)*0.5f, 0.3f, 1.0f };

	cam->glRender();



//	std::cout << "eyePos: " << eyePos[0] << " " << eyePos[1] << " " << eyePos[2] << "\n";

	glGetFloatv(GL_PROJECTION_MATRIX, projectionMatrix);
	glGetFloatv(GL_MODELVIEW_MATRIX, modelviewMatrix);
	

	normalMatrix[0] = modelviewMatrix[0]; normalMatrix[3] = modelviewMatrix[4]; normalMatrix[6] = modelviewMatrix[8];
	normalMatrix[1] = modelviewMatrix[1]; normalMatrix[4] = modelviewMatrix[5]; normalMatrix[7] = modelviewMatrix[9];
	normalMatrix[2] = modelviewMatrix[2]; normalMatrix[5] = modelviewMatrix[6]; normalMatrix[8] = modelviewMatrix[10];

	switch (renderMode){
	case 0:
	{
		//GLSL Spherical Harmonics
		NMFsh->enable();
		NMFsh->SetUniformMatrix4fv("mv_matrix", modelviewMatrix, false);
		NMFsh->SetUniformMatrix4fv("proj_matrix", projectionMatrix, false);
		NMFsh->SetUniformMatrix3fv("normal_matrix", normalMatrix, false);
		NMFsh->setUniform1i("order", cSHtex.getOrder());
		NMFsh->setUniform1f("BPexp", BPexp);
		NMFsh->setUniform1i("renderScene", renderScene);
		NMFsh->setUniform1i("MipMapped", MipMapped);

		//blinphong Coeff
		NMFsh->setUniform1f("b0", brdfSH::BlinnPhong(0, BPexp));
		NMFsh->setUniform1f("b1", brdfSH::BlinnPhong(1, BPexp));
		NMFsh->setUniform1f("b2", brdfSH::BlinnPhong(2, BPexp));
		NMFsh->setUniform1f("b3", brdfSH::BlinnPhong(3, BPexp));
		NMFsh->setUniform1f("b4", brdfSH::BlinnPhong(4, BPexp));
		NMFsh->setUniform1f("b5", brdfSH::BlinnPhong(5, BPexp));
		NMFsh->setUniform1f("b6", brdfSH::BlinnPhong(6, BPexp));

		//		NMFvMF->setUniform3f("eyePos", eyePos[0], eyePos[1], eyePos[2]);
//		glVertexArrayVertexBuffer(VAO, 0, vMFvertex, 0, (GLsizei)(sizeof(float) * 4));
//		glVertexArrayVertexBuffer(VAO, 1, vMFtex, 0, (GLsizei)(sizeof(float) * 4));
//		glVertexArrayVertexBuffer(VAO, 2, vMFnormal, 0, (GLsizei)(sizeof(float) * 4));
//		glVertexArrayVertexBuffer(VAO, 3, vMFtangent, 0, (GLsizei)(sizeof(float) * 4));
//		glVertexArrayAttribFormat(VAO, 0, 4, GL_FLOAT, GL_FALSE, 0);
//		glVertexArrayAttribFormat(VAO, 1, 4, GL_FLOAT, GL_FALSE, 0);
//		glVertexArrayAttribFormat(VAO, 2, 4, GL_FLOAT, GL_TRUE, 0);
//		glVertexArrayAttribFormat(VAO, 3, 4, GL_FLOAT, GL_TRUE, 0);
//		glEnableVertexArrayAttrib(VAO, 0);
//		glEnableVertexArrayAttrib(VAO, 1);
//		glEnableVertexArrayAttrib(VAO, 2);
//		glEnableVertexArrayAttrib(VAO, 3);

		//		glVertexArrayVertexBuffer(TexCoordArray, 0, vMFtex, 0, (GLsizei)(sizeof(float) * 4));
		//		glVertexArrayAttribFormat(TexCoordArray, 0, 4, GL_FLOAT, GL_FALSE, 0);
		//		glEnableVertexArrayAttrib(TexCoordArray, 0);

		glDrawArrays(GL_QUADS, 0, 4);

		NMFsh->disable();
		break;
	}
	case 1://drawObj
	{
		if (renewOBJshader)
		{
			glDetachShader(NMFvMFobj->getProgram(), NMFvMFobj->getFragShader());
			glDetachShader(NMFvMFobj->getProgram(), NMFvMFobj->getVertShader());
			glDeleteShader(NMFvMFobj->getFragShader());
			glDeleteShader(NMFvMFobj->getVertShader());
			delete NMFvMFobj;
			NMFvMFobj = new GLSLProgram();
			NMFvMFobj->GLSLCreateShader("NMF_vMF_OBJ.vert", "NMF_vMF_OBJ.frag");
			glBindAttribLocation(NMFvMFobj->getProgram(), 0, "in_vertices");
			glBindAttribLocation(NMFvMFobj->getProgram(), 1, "in_normalVector");
			glBindAttribLocation(NMFvMFobj->getProgram(), 2, "in_texCoord");
			glBindAttribLocation(NMFvMFobj->getProgram(), 3, "in_tangent");
			glBindAttribLocation(NMFvMFobj->getProgram(), 4, "in_bitangent");
			NMFvMFobj->GLSLLinkShader();
			std::cout << "renewed NMF_vMF_OBJ\n";
			renewOBJshader = false;
		}

		NMFvMFobj->enable();
		NMFvMFobj->SetUniformMatrix4fv("mv_matrix", modelviewMatrix, false);
		NMFvMFobj->SetUniformMatrix4fv("proj_matrix", projectionMatrix, false);
		NMFvMFobj->SetUniformMatrix3fv("normal_matrix", normalMatrix, false);
		NMFvMFobj->setUniform1i("numLobes", numLobes);
		NMFvMFobj->setUniform1f("BPexp", BPexp);
		NMFvMFobj->setUniform1i("renderScene", renderScene);
		NMFvMFobj->setUniform1i("MipMapped", MipMapped);
		NMFvMFobj->setUniform1i("brdfSelect", brdfSelect);
		NMFvMFobj->setUniform1f("MicroSigma", MicroSigma);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementBuffer);
		glDrawElements(GL_TRIANGLES, numIndices, GL_UNSIGNED_INT, (void*)0);



		NMFvMFobj->disable();
//		drawCube2();
//		int indexOffset = 0;
//		for (int f = 0; f < numFaces; f++)
//		{
//			int ngon = objreader.shapes[0].mesh.num_vertices[f];
//			glBegin(GL_TRIANGLES);
//			for (int ind = 0; ind < ngon; ind++)
//			{
//				int v = indices[indexOffset + ind];
//				glVertex3f(objreader.shapes[0].mesh.positions[4 * v + 0],
//					objreader.shapes[0].mesh.positions[4 * v + 1],
//					objreader.shapes[0].mesh.positions[4 * v + 2]);
//			}
//			glEnd();
//			indexOffset += ngon;
//		}


		break;
	}
	case 2://vMF
	{


		break;
	}
	case 3:
	{

		//GLSLvMF
		NMFvMF->enable();
		NMFvMF->SetUniformMatrix4fv("mv_matrix", modelviewMatrix, false);
		NMFvMF->SetUniformMatrix4fv("proj_matrix", projectionMatrix, false);
		NMFvMF->SetUniformMatrix3fv("normal_matrix", normalMatrix, false);
		NMFvMF->setUniform1i("numLobes", numLobes);
		NMFvMF->setUniform1f("BPexp", BPexp);
		NMFvMF->setUniform1i("renderScene", renderScene);
		NMFvMF->setUniform1i("MipMapped", MipMapped);
		NMFvMF->setUniform1i("brdfSelect", brdfSelect);
		NMFvMF->setUniform1f("MicroSigma", MicroSigma);

//		glBufferData(GL_ARRAY_BUFFER, 3 * numVertices*sizeof(GLfloat), objreader.shapes[0].mesh.positions);

		//		NMFvMF->setUniform3f("eyePos", eyePos[0], eyePos[1], eyePos[2]);
//		glVertexArrayVertexBuffer(VAO, 0, vMFvertex, 0, (GLsizei)(sizeof(float) * 4));
//		glVertexArrayVertexBuffer(VAO, 1, vMFtex, 0, (GLsizei)(sizeof(float) * 4));
//		glVertexArrayVertexBuffer(VAO, 2, vMFnormal, 0, (GLsizei)(sizeof(float) * 4));
//		glVertexArrayVertexBuffer(VAO, 3, vMFtangent, 0, (GLsizei)(sizeof(float) * 4));
//		
//		glVertexArrayAttribFormat(VAO, 0, 4, GL_FLOAT, GL_FALSE, 0);
//		glVertexArrayAttribFormat(VAO, 1, 4, GL_FLOAT, GL_FALSE, 0);
//		glVertexArrayAttribFormat(VAO, 2, 4, GL_FLOAT, GL_TRUE, 0);
//		glVertexArrayAttribFormat(VAO, 3, 4, GL_FLOAT, GL_TRUE, 0);
//		glEnableVertexArrayAttrib(VAO, 0);
//		glEnableVertexArrayAttrib(VAO, 1);
//		glEnableVertexArrayAttrib(VAO, 2);
//		glEnableVertexArrayAttrib(VAO, 3);

//		glVertexArrayVertexBuffer(TexCoordArray, 0, vMFtex, 0, (GLsizei)(sizeof(float) * 4));
//		glVertexArrayAttribFormat(TexCoordArray, 0, 4, GL_FLOAT, GL_FALSE, 0);
//		glEnableVertexArrayAttrib(TexCoordArray, 0);

		glDrawArrays(GL_QUADS, 0, 4);

		NMFvMF->disable();
		glFlush();
		break;
	}
	case 5:
	{
		NMForiginal->enable();
		NMForiginal->SetUniformMatrix4fv("mv_matrix", modelviewMatrix, false);
		NMForiginal->SetUniformMatrix4fv("proj_matrix", projectionMatrix, false);
		NMForiginal->SetUniformMatrix3fv("normal_matrix", normalMatrix, false);
		NMForiginal->setUniform1i("numLobes", numLobes);
		NMForiginal->setUniform1f("BPexp", BPexp);
		NMForiginal->setUniform1i("renderScene", renderScene);
		NMForiginal->setUniform1i("MipMapped", MipMapped);
//		NMForiginal->setUniform3f("eyePos", eyePos[0], eyePos[1], eyePos[2]);
//		glVertexArrayVertexBuffer(VAO, 0, vMFvertex, 0, (GLsizei)(sizeof(float) * 4));
//		glVertexArrayVertexBuffer(VAO, 1, vMFtex, 0, (GLsizei)(sizeof(float) * 4));
//		glVertexArrayVertexBuffer(VAO, 2, vMFnormal, 0, (GLsizei)(sizeof(float) * 4));
//		glVertexArrayVertexBuffer(VAO, 3, vMFtangent, 0, (GLsizei)(sizeof(float) * 4));
//		glVertexArrayAttribFormat(VAO, 0, 4, GL_FLOAT, GL_FALSE, 0);
//		glVertexArrayAttribFormat(VAO, 1, 4, GL_FLOAT, GL_FALSE, 0);
//		glVertexArrayAttribFormat(VAO, 2, 4, GL_FLOAT, GL_TRUE, 0);
//		glVertexArrayAttribFormat(VAO, 3, 4, GL_FLOAT, GL_TRUE, 0);
//		glEnableVertexArrayAttrib(VAO, 0);
//		glEnableVertexArrayAttrib(VAO, 1);
//		glEnableVertexArrayAttrib(VAO, 2);
//		glEnableVertexArrayAttrib(VAO, 3);

		//		glVertexArrayVertexBuffer(TexCoordArray, 0, vMFtex, 0, (GLsizei)(sizeof(float) * 4));
		//		glVertexArrayAttribFormat(TexCoordArray, 0, 4, GL_FLOAT, GL_FALSE, 0);
		//		glEnableVertexArrayAttrib(TexCoordArray, 0);

		glDrawArrays(GL_QUADS, 0, 4);

		NMForiginal->disable();

		break;
	}
	case 4://Original NMF
	{
		float *NMFpixels = new float[width*height*3];
		BYTE *NMFbyte = (BYTE*)malloc(width*height * 3);
		BYTE nn[3];
		int pixelIndex[2] = { 0 };
		for (int i = 0; i < width*height * 3; i++){
			NMFpixels[i] = 0.f;
		}
		for (int j = 0; j < NMheight; j++)
		{
			for (int i = 0; i < NMwidth; i++)
			{
				nn[0] = (float)textureData[j*NMwidth * 4 + i * 4 + 2];
				nn[1] = (float)textureData[j*NMwidth * 4 + i * 4 + 1];
				nn[2] = (float)textureData[j*NMwidth * 4 + i * 4 + 0];

				pixelIndex[0] = (int)((float)nn[0] * (float)(width-1) / 255.f);
				pixelIndex[1] = (int)((float)nn[1] * (float)(width-1) / 255.f);;

				NMFpixels[pixelIndex[1] * width * 3 + pixelIndex[0] * 3 + 0] += 1.1f / (float)(NMheight*NMwidth);
				NMFpixels[pixelIndex[1] * width * 3 + pixelIndex[0] * 3 + 1] += 1.1f / (float)(NMheight*NMwidth);
				NMFpixels[pixelIndex[1] * width * 3 + pixelIndex[0] * 3 + 2] += 1.1f / (float)(NMheight*NMwidth);


			}
		}
		float max = 0.f;
		for (int i = 0; i < width*height * 3; i++){
			NMFbyte[i] = 1500000 * NMFpixels[i];
			if (NMFpixels[i]>max)
				max = NMFpixels[i];
		}
		std::cout << "max value: " << max << std::endl;
		glDrawPixels(width, height, GL_RGB, GL_UNSIGNED_BYTE, NMFbyte);
		delete[] NMFpixels;
		delete[] NMFbyte;
		break;
	}
	}



	glFlush();
}
void idleCB(){
	time(&currentTime);
	t = (float)difftime(currentTime, sysTime);

//	glutPostRedisplay();
}

void keyboardCB(unsigned char key, int x, int y){
	if (renderMode == 0)
	{
		switch (key)
		{
		case 'r':case'R':
			delete NMFsh;
			NMFsh = new GLSLProgram("NMF_SH.vert", "NMF_SH.frag");
			std::cout << "Renewed NMF_SH GLSL\n";
			break;
		case 'c':case'C':
			printf("New BPexp: ");
			scanf("%f", &BPexp);
			break;
		}
	}
	if ((renderMode == 3)||(renderMode==1))
	{
		switch (key){
		case 'b':case'B':
			brdfSelect == 0 ? (brdfSelect = 1) : (brdfSelect = 0);
			break;
		case 'c':case'C':
			printf("New BPexp: ");
			scanf("%f", &BPexp);
			break;
		case 'n':case'N':
			printf("New Micro: ");
			scanf("%f", &MicroSigma);
			break;
		case 's':case'S':
			renderScene == 0 ? (renderScene = 1) : (renderScene = 0);
			break;
		case 'm':case'M':
			(MipMapped == 0) ? (MipMapped = 1) : (MipMapped = 0);
//			if (MipMapped == 0) { std::cout << "NotMipMapped\n"; } else { std::cout << "MipMapped\n"; }
			(MipMapped == 0) ? (std::cout << "NotMipMapped\n") : (std::cout << "MipMapped\n");
			break;
		case 'r':case'R':
			if (renderMode == 3)
			{
				delete NMFvMF;
				NMFvMF = new GLSLProgram("NMF_vMF.vert", "NMF_vMF.frag");
				std::cout << "Renewed NMF_vMF Program\n";
			}
			else
			{
				renewOBJshader = true;
			}
			break;
		case 'q':case'Q':
		{
			GLfloat m[16];
			glGetFloatv(GL_MODELVIEW_MATRIX, m);
			std::ofstream camFile("camPos");
			float pos[3], up[3], lookat[3];
			cam->getPosition(pos, lookat, up);
			camFile << pos[0] << " " << pos[1] << " " << pos[2]<<"\n";
			camFile << lookat[0] << " " << lookat[1] << " " << lookat[2] << "\n";
			camFile << up[0] << " " << up[1] << " " << up[2] << "\n";

			camFile.close();
			std::cout << "Saved Cam Pos\n";
			std::cout << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
			break;
		}

		case 'w':
		{
			GLfloat m[16];
			glGetFloatv(GL_MODELVIEW_MATRIX, m);

			float pos[3], lookat[3], up[3];
			std::ifstream camFile("camPos");
			camFile >> pos[0] >> pos[1]  >> pos[2];
			camFile >> lookat[0] >> lookat[1] >> lookat[2];
			camFile >> up[0] >> up[1] >> up[2];

			cam->setPosition(pos, lookat, up);

			std::cout << "loaded cam\n";

			break;
		}

		}

	}
//	glutSetWindow(windowSH);
	glutPostRedisplay();
}
void mouseCB(int button, int state, int x, int y) {
	if ((renderMode == 2) || (renderMode == 4)) {}
	else
	{
		if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
			click_pos[0] = (x - width / 2.0f) / (width / 2.0f);
			click_pos[1] = (y - height / 2.0f) / (height / 2.0f);
		}
		else if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) {
			click_pos[0] = (x - width / 2.0f) / (width / 2.0f);
			click_pos[1] = (y - height / 2.0f) / (height / 2.0f);
			if (glutGetModifiers() == GLUT_ACTIVE_CTRL)
			{
				cam->trans = true;
			}
			else
			{
				cam->trans = false;
			}
		}
		button_id = button;
	}
	glutSetWindow(windowSH);
	glutPostRedisplay();
}
void motionCB(int x, int y){
	GLfloat present[2];
	if ((renderMode == 2) || (renderMode == 4)){}
	else
	{
		switch (button_id) {
		case GLUT_LEFT_BUTTON:
			present[0] = (GLfloat)(x - width / 2.0f) / (GLfloat)(width / 2.0f);
			present[1] = (GLfloat)(y - height / 2.0f) / (GLfloat)(height / 2.0f);

			cam->trackball(click_pos, present);

			click_pos[0] = present[0];
			click_pos[1] = present[1];
			break;
		case GLUT_MIDDLE_BUTTON:
			present[0] = (GLfloat)(x - width / 2.0f) / (GLfloat)(width / 2.0f);
			present[1] = (GLfloat)(y - height / 2.0f) / (GLfloat)(height / 2.0f);
			if (!cam->trans)
			{
				if (present[1] - click_pos[1] < 0)
					cam->dollyin();
				else if (present[1] - click_pos[1] > 0)
					cam->dollyout();

				click_pos[0] = present[0];
				click_pos[1] = present[1];
			}
			else
			{
				glGetFloatv(GL_PROJECTION_MATRIX, projectionMatrix);
				cam->translate(click_pos, present, projectionMatrix);
				click_pos[0] = present[0];
				click_pos[1] = present[1];
			}
			break;
		}
	}
	glutSetWindow(windowSH);
	glutPostRedisplay();
}

void mainMenu(int op) {
	switch (op) {
	case 0:
	{
		delete[] alpha;
		for (int i = 0; i < numLobes; i++)
		{
			delete[] aux[i];

		}
		delete[] aux;
		exit(0);
	}
	}
	glutPostRedisplay();
}
void renderMenu(int op) {
	switch (op) {
	case 0:
		renderMode = 0;
		break;
	case 3:
		renderMode = 3;
		break;
	case 1:
		renderMode = 1;
		break;
	case 2:
		renderMode = 2;
		break;
	case 4:
		renderMode = 4;
		break;
	}
	glutPostRedisplay();
}
void initGLUT(int argc, char** argv){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	windowSH=glutCreateWindow("SH Test");
	glutSetWindow(windowSH);
	glutDisplayFunc(displayCB);
	glutIdleFunc(idleCB);
	glutMouseFunc(mouseCB);
	glutKeyboardFunc(keyboardCB);
	glutMotionFunc(motionCB);

	int RenderMenu = glutCreateMenu(renderMenu);
	glutAddMenuEntry("GLSL", 0);
	glutAddMenuEntry("Render Scene(vMF)", 3);
	glutAddMenuEntry("drawOBJ", 1);
	glutAddMenuEntry("vMF", 2);
	glutAddMenuEntry("Original NDF", 4);
	int MainMenu = glutCreateMenu(mainMenu);
	glutAddSubMenu("RenderMode", RenderMenu);
	glutAddMenuEntry("Clean&Exit", 0);

	glutAttachMenu(GLUT_RIGHT_BUTTON);


}


void initGL(int argc, char** argv){

	if (glewInit() != GLEW_OK){
		printf("GLEW init failed.\n");
		sysPause;
		exit(EXIT_FAILURE);
	}
	if (!GLEW_VERSION_4_5){
		printf("openGL 4.5 not supported\n");
		sysPause;
		exit(EXIT_FAILURE);
	}
	else
		printf("opengl 4.5 supported\n");

	createProgram();
	setMesh();

	cVMFtex.generatevMFmaps();

	initBuffers();

	cam = new Camera();

}

void computeTangentBasis(float* vertices, float* uvs, float* normals, float* tangent, float* bitangent, unsigned int* indices, int numVertices, int numFaces)
{
	std::vector< std::vector<float> > tIndexing;
	std::vector< std::vector<float> > bIndexing;
	tIndexing.resize(numVertices);
	bIndexing.resize(numVertices);
	int idx[3];
	for (int f = 0; f < numFaces; f++)
	{
		idx[0] = indices[3 * f + 0];
		idx[1] = indices[3 * f + 1];
		idx[2] = indices[3 * f + 2];

		float v0[3] = { vertices[idx[0] * 3 + 0], vertices[idx[0] * 3 + 1], vertices[idx[0] * 3 + 2] };
		float v1[3] = { vertices[idx[1] * 3 + 0], vertices[idx[1] * 3 + 1], vertices[idx[1] * 3 + 2] }; 
		float v2[3] = { vertices[idx[2] * 3 + 0], vertices[idx[2] * 3 + 1], vertices[idx[2] * 3 + 2] };
													   
		float uv0[2] = { uvs[idx[0] * 2 + 0], uvs[idx[0] * 2 + 1] };
		float uv1[2] = { uvs[idx[1] * 2 + 0], uvs[idx[1] * 2 + 1] };
		float uv2[2] = { uvs[idx[2] * 2 + 0], uvs[idx[2] * 2 + 1] };

		float dv1[3] = { v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2] };
		float dv2[3] = { v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2] };

		float deltaUV1[2] = { uv1[0] - uv0[0], uv1[1] - uv0[1] };
		float deltaUV2[2] = { uv2[0] - uv0[0], uv2[1] - uv0[1] };

		float r = 1.0f / (deltaUV1[0] * deltaUV2[1] - deltaUV1[1] * deltaUV2[0]);
		
		float t[3] = { 0 }, b[3] = { 0 };
		t[0] = (dv1[0] * deltaUV2[1] - dv2[0] * deltaUV2[1])*r;
		t[1] = (dv1[1] * deltaUV2[1] - dv2[1] * deltaUV2[1])*r;
		t[2] = (dv1[2] * deltaUV2[1] - dv2[2] * deltaUV2[1])*r;

		b[0] = (dv2[0] * deltaUV2[0] - dv1[0] * deltaUV2[0])*r;
		b[1] = (dv2[1] * deltaUV2[0] - dv1[1] * deltaUV2[0])*r;
		b[2] = (dv2[2] * deltaUV2[0] - dv1[2] * deltaUV2[0])*r;
		
		tIndexing[idx[0]].push_back(t[0]); tIndexing[idx[0]].push_back(t[1]); tIndexing[idx[0]].push_back(t[2]);
	//	tIndexing[idx[1]].push_back(t[0]); tIndexing[idx[1]].push_back(t[1]); tIndexing[idx[1]].push_back(t[2]);
	//	tIndexing[idx[2]].push_back(t[0]); tIndexing[idx[2]].push_back(t[1]); tIndexing[idx[2]].push_back(t[2]);

		bIndexing[idx[0]].push_back(b[0]); bIndexing[idx[0]].push_back(b[1]); bIndexing[idx[0]].push_back(b[2]);
	//	bIndexing[idx[1]].push_back(b[0]); bIndexing[idx[1]].push_back(b[1]); bIndexing[idx[1]].push_back(b[2]);
	//	bIndexing[idx[2]].push_back(b[0]); bIndexing[idx[2]].push_back(b[1]); bIndexing[idx[2]].push_back(b[2]);
	}
	for (int i = 0; i < numVertices; i++)
	{
		float t[3] = { 0 };
		float b[3] = { 0 };
		int sizei = tIndexing[i].size() / 3;
		for (int j = 0; j < sizei; j++)
		{
			t[0] += tIndexing[i][3 * j + 0]; t[1] += tIndexing[i][3 * j + 1]; t[2] += tIndexing[i][3 * j + 2];
			b[0] += tIndexing[i][3 * j + 0]; b[1] += tIndexing[i][3 * j + 1]; t[2] += tIndexing[i][3 * j + 2];
		}
		tangent[i * 3 + 0] = t[0]; tangent[i * 3 + 1] = t[1]; tangent[i * 3 + 2] = t[2];
		bitangent[i * 3 + 0] = b[0]; bitangent[i * 3 + 1] = b[1]; bitangent[i * 3 + 2] = b[2];
	}

}

void setMesh()
{
	if (objreader.shapes.size() > 1)
	{
		std::cout << "More than 1 shape exists: " << objreader.shapes.size() << std::endl;
		sysPause;
	}
	numFaces = objreader.shapes[0].mesh.num_vertices.size();
	numVertices = objreader.shapes[0].mesh.positions.size()/3;
	numNormals = objreader.shapes[0].mesh.normals.size()/3;
	numTexcoord = objreader.shapes[0].mesh.texcoords.size()/2;
	numIndices = objreader.shapes[0].mesh.indices.size();
	vertPos = new float[numVertices * 4];
	normals = new float[numNormals * 3];
	texCoord = new float[numTexcoord * 2];
	indices = new unsigned int[numIndices];
	tangent = new float[numVertices * 3];
	bitangent = new float[numVertices * 3];


	for (int i = 0; i < numVertices; i++)
	{
		vertPos[i * 4 + 0] = objreader.shapes[0].mesh.positions[3 * i + 0];
		vertPos[i * 4 + 1] = objreader.shapes[0].mesh.positions[3 * i + 1];
		vertPos[i * 4 + 2] = objreader.shapes[0].mesh.positions[3 * i + 2];
		vertPos[i * 4 + 3] = 1.f;
	}
	for (int i = 0; i < numNormals; i++)
	{
		normals[i * 3 + 0] = objreader.shapes[0].mesh.normals[3 * i + 0];
		normals[i * 3 + 1] = objreader.shapes[0].mesh.normals[3 * i + 1];
		normals[i * 3 + 2] = objreader.shapes[0].mesh.normals[3 * i + 2];
	}
	for (int i = 0; i < numTexcoord; i++)
	{
		texCoord[i * 2 + 0] = objreader.shapes[0].mesh.texcoords[2 * i + 0];
		texCoord[i * 2 + 1] = objreader.shapes[0].mesh.texcoords[2 * i + 1];
	}
	for (int i = 0; i < numIndices; i++)
	{
		indices[i] = objreader.shapes[0].mesh.indices[i];
	}

	computeTangentBasis(vertPos, texCoord, normals, tangent, bitangent, indices, numVertices, numFaces);


	//compute tangents
	

	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);
	glGenBuffers(5, VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
	glBufferData(GL_ARRAY_BUFFER, 4 * numVertices * sizeof(GLfloat), vertPos, GL_STATIC_DRAW);
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
	glBufferData(GL_ARRAY_BUFFER, 3 * numNormals * sizeof(GLfloat), normals, GL_STATIC_DRAW);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, VBO[2]);
	glBufferData(GL_ARRAY_BUFFER, 2 * numTexcoord * sizeof(GLfloat), texCoord, GL_STATIC_DRAW);
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(2);
	glBindBuffer(GL_ARRAY_BUFFER, VBO[3]);
	glBufferData(GL_ARRAY_BUFFER, 3 * numVertices * sizeof(GLfloat), tangent, GL_STATIC_DRAW);
	glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(3);
	glBindBuffer(GL_ARRAY_BUFFER, VBO[4]);
	glBufferData(GL_ARRAY_BUFFER, 3 * numVertices * sizeof(GLfloat), bitangent, GL_STATIC_DRAW);
	glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(4);


	//indices
	glGenBuffers(1, &elementBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, numIndices * sizeof(unsigned int), indices, GL_STATIC_DRAW);

}

int main(int argc, char** argv){
	srand(time(NULL));

	initSharedMem();

//	cSHtex.displayMap(4);

	objreader.readObj(inputObj);

	//float testa = 1.f;
	//float *testaptr = &testa;
	//float testAux[3] = { 0.f, 0.f, 1.f };
	//float *testA1 = testAux;
	//float **testAuxptr = &testA1;

	//vMFfunc::displayvMF(1, testaptr, testAuxptr, 512, 512);


//	cVMFtex.showOriginalImage(0);


//	cVMFtex.showvMFImage(0, 0, 0);
	

//	NMTdata = vMFfunc::LoadImage(NMT, NMwidth, NMheight);



	initGLUT(argc, argv);
	initGL(argc, argv);

	printf("Running Well\n");

	glutMainLoop();
	sysPause;
	return 0;

}