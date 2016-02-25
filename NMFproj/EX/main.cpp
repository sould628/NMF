//kixor.net
//OpenGL Extension registry: http://www.opengl.org/registry
//Realtech VR's OpenGL Extensions Viewer

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

#include "FreeImage.h"
#include "objLoader.h"

#include "globalVariables.h"


#include "GLSLProgram.h"
#include "Camera.h"

#define sysPause system("pause>nul");

#define MAXITERATION 10






const int l_index = 4;

int width = 800, height = 800;
int NMwidth, NMheight;


objLoader *objData;
const char* NMT = "./2.jpg";
//const char* NMT = "./bricks_normal_map.jpg";

FIBITMAP* NMTdata;
GLubyte* textureData;

Camera *cam;
GLint button_id;
GLfloat click_pos[2];

float t;

GLuint VAO, TexCoordArray;
GLuint NMbuffer, vMFvertex, vMFtex, vMFnormal, vMFtangent;
GLuint NormalMap, NormalMipMap;

GLuint SHmap0, SHmap1, SHmap2, SHmap3, SHmap4, SHmap5, SHmap6;

GLuint vMFmap0, vMFmap1, vMFmap2, vMFmap3, vMFmap4;



GLSLProgram *NMF, *NMFvMF, *NMForiginal;


time_t sysTime, currentTime;

//vMF
float* alpha;
float** mu;
float* kappa;
float** aux;
const int numLobes = 5;


static const float exampleData[] =
{
	0.25, -0.25, 0.5, 1.0,
	-0.25, 0.25, 0.5, 1.0,
	0.25, 0.25, 0.5, 1.0
};






void displayCB();
void drawTexture();
void generateSHmap(GLubyte* TextureData);
void generatevMFmap(int numLobes, float* alpha, float* aux[3]);

float calculateSHcoeffDelta(int l, int m, float x, float y, float z);
void checkTextureError(GLenum glError);
int vMFparam(GLubyte* TextureData, float* alpha, float** mu, float* kappa, float* aux[3], int numLobes);
float vMF(float normal[3], float mu[3], float kappa);
float norm(float vector[3]);
void normalize(float* source, float* destination);



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

float vMF(float normal[3], float mu[3], float kappa) {
	float NdotMu = (mu[0] * normal[0]) + (mu[1] * normal[1]) + (mu[2] * normal[2]);
	
	float result = (kappa / (4 * PI*sinh(kappa)))*exp(kappa*(NdotMu));

	//if (result >= 1)
	//{
	//	std::cout << "Error vMF value detected\n";
	//	sysPause;
	//}


	return result;
	
}

int vMFparam(GLubyte* TextureData, float* alpha, float** mu, float* kappa, float* aux[3],  int numLobes) {
	bool bConvergence = false;
	float* normalData;
	normalData = new float[NMwidth*NMheight * 3];
	float nn[3] = { 0 };

	float** z;
	z = new float*[numLobes];
	for (int i = 0; i < numLobes; i++)
	{
		z[i] = new float[NMwidth*NMheight];
	}

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

//			else
//			{
//				if (nn[2] == 1)
//				{
//					nn[0] = 0.f; nn[1] = 0.f;
//				}
//				else if (nn[0] == 1)
//				{
//					nn[1] = 0.f; nn[2] = 0.f;
//				}
//				else {
//					while ((nn[2] * nn[2] + nn[0] * nn[0]) > 1)
//					{
//						nn[2] = nn[2] * 0.99; nn[0] * 0.99;
//					}
//					nn[1] = sqrt(1 - (nn[2] * nn[2] + nn[0] * nn[0]));
//				}
//			}
			normalData[j*NMwidth * 3 + i * 3] = nn[0];
			normalData[j*NMwidth * 3 + i * 3 + 1] = nn[1];
			normalData[j*NMwidth * 3 + i * 3 + 2] = nn[2];
		}//i
	}//j
	//NormalData

	int iteration = 0;
	//EM
	while (!bConvergence)
	{
		//E-STEP
		float* vMFij = new float[numLobes];
		for (int i = 0; i < NMwidth*NMheight; i++)
		{

			float vMFsum = 0.f;
			for (int j = 0; j < numLobes; j++)
			{
				vMFij[j] = vMF(&normalData[(i / NMheight)*NMheight * 3 + (i%NMheight) * 3], mu[j], kappa[j]);
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
			for (int ii = 0; ii < NMwidth*NMheight; ii++)
			{
				zsum[jj] += z[jj][ii];

				expectedNormal[jj][0] += z[jj][ii] * normalData[(ii / NMheight)*NMheight * 3 + ii%NMheight * 3];
				expectedNormal[jj][1] += z[jj][ii] * normalData[(ii / NMheight)*NMheight * 3 + ii%NMheight * 3 + 1];
				expectedNormal[jj][2] += z[jj][ii] * normalData[(ii / NMheight)*NMheight * 3 + ii%NMheight * 3 + 2];
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
	delete[] normalData;
	for (int i = 0; i < numLobes; i++)
	{
		delete[] z[i];
	}
	delete[] z;
	return -1;
}

void generatevMFmap(int numLobes, float* alpha, float* aux[3]){
	

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 0);
	glGenTextures(1, &vMFmap0);	glGenTextures(1, &vMFmap1);	glGenTextures(1, &vMFmap2);	glGenTextures(1, &vMFmap3); glGenTextures(1, &vMFmap4);
	float* dataBuffer;
	dataBuffer = new float[NMwidth*NMheight * 4];
	switch (numLobes)
	{
	case 6:
	{

	}
	case 5:			//주석처리 했을 때 안했을 때 결과가 달라짐?
	{
		for (int j = 0; j < NMheight; j++)
		{
			for (int i = 0; i < NMwidth; i++)
			{
				dataBuffer[j*NMwidth * 4 + i * 4 + 0] = alpha[4];
				dataBuffer[j*NMwidth * 4 + i * 4 + 1] = alpha[4]*aux[4][0];
				dataBuffer[j*NMwidth * 4 + i * 4 + 2] = alpha[4]*aux[4][1];
				dataBuffer[j*NMwidth * 4 + i * 4 + 3] = alpha[4]*aux[4][2];
			}
		}
		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_2D, vMFmap4);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NMwidth, NMheight, 0, GL_RGBA, GL_FLOAT, dataBuffer);
//		gluBuild2DMipmaps(GL_TEXTURE_2D, 5, NMwidth, NMheight, GL_RGBA, GL_FLOAT, dataBuffer);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		GLenum glError = glGetError();
		checkTextureError(glError);
		glActiveTexture(GL_TEXTURE10);
		glBindTexture(GL_TEXTURE_2D, 0);
	}
	case 4:
	{
		for (int j = 0; j < NMheight; j++)
		{
			for (int i = 0; i < NMwidth; i++)
			{
				dataBuffer[j*NMwidth * 4 + i * 4 + 0] = alpha[3];
				dataBuffer[j*NMwidth * 4 + i * 4 + 1] = alpha[3]*aux[3][0];
				dataBuffer[j*NMwidth * 4 + i * 4 + 2] = alpha[3]*aux[3][1];
				dataBuffer[j*NMwidth * 4 + i * 4 + 3] = alpha[3]*aux[3][2];
			}
		}
		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D, vMFmap3);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NMwidth, NMheight, 0, GL_RGBA, GL_FLOAT, dataBuffer);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		GLenum glError = glGetError();
		checkTextureError(glError);
		glActiveTexture(GL_TEXTURE10);
		glBindTexture(GL_TEXTURE_2D, 0);
	}
	case 3:
	{
		for (int j = 0; j < NMheight; j++)
		{
			for (int i = 0; i < NMwidth; i++)
			{
				dataBuffer[j*NMwidth * 4 + i * 4 + 0] = alpha[2];
				dataBuffer[j*NMwidth * 4 + i * 4 + 1] = alpha[2]*aux[2][0];
				dataBuffer[j*NMwidth * 4 + i * 4 + 2] = alpha[2]*aux[2][1];
				dataBuffer[j*NMwidth * 4 + i * 4 + 3] = alpha[2]*aux[2][2];
			}
		}
		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, vMFmap2);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NMwidth, NMheight, 0, GL_RGBA, GL_FLOAT, dataBuffer);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		GLenum glError = glGetError();
		checkTextureError(glError);
		glActiveTexture(GL_TEXTURE10);
		glBindTexture(GL_TEXTURE_2D, 0);
	}
	case 2:
	{
		for (int j = 0; j < NMheight; j++)
		{
			for (int i = 0; i < NMwidth; i++)
			{
				dataBuffer[j*NMwidth * 4 + i * 4 + 0] = alpha[1];
				dataBuffer[j*NMwidth * 4 + i * 4 + 1] = alpha[1]*aux[1][0];
				dataBuffer[j*NMwidth * 4 + i * 4 + 2] = alpha[1]*aux[1][1];
				dataBuffer[j*NMwidth * 4 + i * 4 + 3] = alpha[1]*aux[1][2];
			}												  
		}
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, vMFmap1);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NMwidth, NMheight, 0, GL_RGBA, GL_FLOAT, dataBuffer);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		GLenum glError = glGetError();
		checkTextureError(glError);
		glActiveTexture(GL_TEXTURE10);
		glBindTexture(GL_TEXTURE_2D, 0);
	}
	case 1:
	{
		for (int j = 0; j < NMheight; j++)
		{
			for (int i = 0; i < NMwidth; i++)
			{
				dataBuffer[j*NMwidth * 4 + i * 4 + 0] = alpha[0];
				dataBuffer[j*NMwidth * 4 + i * 4 + 1] = alpha[0]*aux[0][0];
				dataBuffer[j*NMwidth * 4 + i * 4 + 2] = alpha[0]*aux[0][1];
				dataBuffer[j*NMwidth * 4 + i * 4 + 3] = alpha[0]*aux[0][2];
			}
		}
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, vMFmap0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NMwidth, NMheight, 0, GL_RGBA, GL_FLOAT, dataBuffer);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		GLenum glError = glGetError();
		checkTextureError(glError);
		glActiveTexture(GL_TEXTURE10);
		glBindTexture(GL_TEXTURE_2D, 0);
	}
	default:
		break;
	}
	glDisable(GL_TEXTURE_2D);
	delete[] dataBuffer;
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

void generateSHmap(GLubyte* TextureData){
	float nn[3] = { 0 };
	int textureCount = 0;
	GLfloat* coeffTexture;
	coeffTexture = new GLfloat[NMwidth*NMheight * 4];
	for (int ll = 0; ll < l_index; ll++)
	{
		for (int mm = 0; mm < ll * 2 + 1; mm++)
		{
			for (int i = 0; i < NMwidth; i++)
			{
				for (int j = 0; j < NMheight; j++)
				{
					nn[0] = (float)TextureData[i*NMheight * 4 + j * 4];
					nn[1] = (float)TextureData[i*NMheight * 4 + j * 4+1];
					nn[2] = (float)TextureData[i*NMheight * 4 + j * 4+2];

					nn[0] = (nn[0] - 127.f) / 128.f;
					nn[1] = (nn[1] - 127.f) / 128.f;
					nn[2] = (nn[2] - 127.f) / 128.f;
					if ((nn[2] * nn[2] + nn[0] * nn[0])<=1)
						nn[1] = sqrt(1 - (nn[2] * nn[2] + nn[0] * nn[0]));
					else
					{
						if (nn[2] == 1)
						{
							nn[0] = 0.f; nn[1] = 0.f;
						}
						else if (nn[0] == 1)
						{
							nn[1] = 0.f; nn[2] = 0.f;
						}
						else{
							while ((nn[2] * nn[2] + nn[0] * nn[0]) > 1)
							{
								nn[2] = nn[2] * 0.99; nn[0] * 0.99;
							}
							nn[1] = sqrt(1 - (nn[2] * nn[2] + nn[0] * nn[0]));
						}
					}

					
					coeffTexture[i*NMheight*4+j*4+textureCount]=calculateSHcoeffDelta(ll, (mm - ll), nn[0], nn[1], nn[2]);

				}
			}
			if (textureCount++ == 3)
			{

				textureCount = 0;
			}
		}

	}

	delete[] coeffTexture;
	coeffTexture = NULL;
}



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



#pragma region fbostatusfunction

bool checkFrameBufferObjectStatus(){
	GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	switch (status){
	case GL_FRAMEBUFFER_COMPLETE:
		std::cout << "Framebuffer complete." << std::endl;
		return true;

	case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:
		std::cout << "[ERROR] Framebuffer incomplete: Attachment is NOT complete." << std::endl;
		return false;

	case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT:
		std::cout << "[ERROR] Framebuffer incomplete: No image is attached to FBO." << std::endl;
		return false;
	case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER:
		std::cout << "[ERROR] Framebuffer incomplete: Draw buffer." << std::endl;
		return false;

	case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER:
		std::cout << "[ERROR] Framebuffer incomplete: Read buffer." << std::endl;
		return false;

	case GL_FRAMEBUFFER_UNSUPPORTED:
		std::cout << "[ERROR] Framebuffer incomplete: Unsupported by FBO implementation." << std::endl;
		return false;

	default:
		std::cout << "[ERROR] Framebuffer incomplete: Unknown error." << std::endl;
		return false;
	}
}

#pragma endregion

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
	NMF = new GLSLProgram("NMF_SH.vert", "NMF_SH.frag");
	NMFvMF = new GLSLProgram("NMF_vMF.vert", "NMF_vMF.frag");
	NMForiginal = new GLSLProgram("NMF_original.vert", "NMF_original.frag");
}
void initBuffers() {
	glEnable(GL_TEXTURE_2D);
	glCreateVertexArrays(1, &VAO);
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

	generateSHmap(textureData);
	//vMFparameter
	alpha = new float[numLobes];
	mu = new float*[numLobes];
	kappa = new float[numLobes];
	aux = new float*[numLobes];
	for (int i = 0; i < numLobes; i++) 
	{
		mu[i] = new float[3];
		aux[i] = new float[3];
//		mu[i][0] = 0.5-(float)((int)(i / 2));
//		mu[i][1] = 0.5-(float)((int)(i % 2));
//		mu[i][2] = sqrt(1 - ((mu[i][0] * mu[i][0]) + (mu[i][1] * mu[i][1])));
		aux[i][0] = 0.f;
		aux[i][1] = 0.f;
		aux[i][2] = 0.f;
		alpha[i] = 1.f / (float)numLobes;
		kappa[i] = 30.f;
	}
	if (numLobes < 5)
	{
		for (int i = 0; i < numLobes; i++)
		{
//			mu[i][0] = 0.5 - (float)((int)(i / 2));
//			mu[i][1] = 0.5 - (float)((int)(i % 2));
//			mu[i][2] = sqrt(1 - ((mu[i][0] * mu[i][0]) + (mu[i][1] * mu[i][1])));
			mu[0][0] = 0.f;		mu[0][1] = -0.7f;	mu[0][2] = sqrt(1 - ((mu[0][0] * mu[0][0]) + (mu[0][1] * mu[0][1])));;
			mu[1][0] = 0.7f;	mu[1][1] = 0.f;		mu[1][2] = sqrt(1 - ((mu[1][0] * mu[1][0]) + (mu[1][1] * mu[1][1])));
			mu[2][0] = -0.7f;	mu[2][1] = 0.f;		mu[2][2] = sqrt(1 - ((mu[2][0] * mu[2][0]) + (mu[2][1] * mu[2][1])));;
			mu[3][0] = 0.f;		mu[3][1] = 0.7f;	mu[3][2] = sqrt(1 - ((mu[3][0] * mu[3][0]) + (mu[3][1] * mu[3][1])));;

		}
	}
	else if (numLobes == 5)
	{

			mu[0][0] = 0.f;		mu[0][1] = -0.7f;	mu[0][2] = sqrt(1 - ((mu[0][0] * mu[0][0]) + (mu[0][1] * mu[0][1])));;
			mu[1][0] = 0.7f;	mu[1][1] = 0.f;		mu[1][2] = sqrt(1 - ((mu[1][0] * mu[1][0]) + (mu[1][1] * mu[1][1])));
			mu[2][0] = -0.7f;	mu[2][1] = 0.f;		mu[2][2] = sqrt(1 - ((mu[2][0] * mu[2][0]) + (mu[2][1] * mu[2][1])));;
			mu[3][0] = 0.f;		mu[3][1] = 0.7f;	mu[3][2] = sqrt(1 - ((mu[3][0] * mu[3][0]) + (mu[3][1] * mu[3][1])));;
			mu[4][0] = 0.0f;	mu[4][1] = 0.0f;	mu[4][2] = sqrt(1 - ((mu[4][0] * mu[4][0]) + (mu[4][1] * mu[4][1])));;
	}
	for (int i = 0; i < numLobes; i++)
	{
		normalize(mu[i], mu[i]);
	}

	vMFparam(textureData, alpha, mu, kappa, aux, numLobes);

//	delete[] kappa;
//	delete[] alpha;
//	for (int i = 0; i < numLobes; i++)
//	{
//		delete[] aux[i];
//		delete[] mu[i];
//
//	}
//	delete[] mu;
//	delete[] aux;
	
	generatevMFmap(numLobes, alpha, aux);

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
	glActiveTexture(GL_TEXTURE9);
	glBindTexture(GL_TEXTURE_2D, NormalMipMap);
//	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NMwidth, NMheight, 0, GL_BGRA, GL_UNSIGNED_BYTE, textureData);
	gluBuild2DMipmaps(GL_TEXTURE_2D, 3, NMwidth, NMheight, GL_BGRA, GL_UNSIGNED_BYTE, textureData);
	glActiveTexture(GL_TEXTURE10);
	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);

}


#pragma endregion





void drawObj(objLoader *data){
	if (data->faceCount == 0)
	{
		std::cout << "No Obj Exists" << std::endl;
		return;
	}
	obj_vector **vertList=data->vertexList;
	std::cout << "numFaces: " << data->faceCount << std::endl;
	for (int i = 0; i < data->faceCount; i++)
	{
		glColor3f((GLfloat)((float)i/10.f-i/10), 0.5, 0);
		obj_face *flist;
		flist = data->faceList[i];
		int matIndex=flist->material_index;
		glBegin(GL_POLYGON);
		for (int ii = 0; ii < flist->vertex_count; ii++){
			int* vertIndex = flist->vertex_index;
			glVertex3f(vertList[vertIndex[ii]]->e[0], vertList[vertIndex[ii]]->e[1], vertList[vertIndex[ii]]->e[2]);
//			std::cout << "Current Vertex: "<<vertList[vertIndex[ii]]->e[0] << "," << vertList[vertIndex[ii]]->e[1] << "," << vertList[vertIndex[ii]]->e[2] << std::endl;
		}
		glEnd();


	}
}


void drawTexture(){
	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 1.0);
	glVertex3f(0.0, height, -0.1);
	glTexCoord2f(1.0, 1.0);
	glVertex3f(width, height, -0.1);
	glTexCoord2f(1.0, 0.0);
	glVertex3f(width, 0.0, -0.1);
	glTexCoord2f(0.0, 0.0);
	glVertex3f(0.0, 0.0, -0.1);
	glEnd();
	
}

void displayCB(){
	GLfloat attrib[4];
	glEnable(GL_DEPTH_TEST);
	glClearColor(0.f, 0.f, 0.0f, 1.f);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	
	GLfloat color[] = { (float)sin(t)*0.5f, (float)cos(t)*0.5f, 0.3f, 1.0f };

	cam->glRender();

	float eyePos[3] = { 0 };

	cam->getPosition(eyePos);

//	std::cout << "eyePos: " << eyePos[0] << " " << eyePos[1] << " " << eyePos[2] << "\n";

	glGetFloatv(GL_PROJECTION_MATRIX, projectionMatrix);
	glGetFloatv(GL_MODELVIEW_MATRIX, modelviewMatrix);
	

	normalMatrix[0] = modelviewMatrix[0]; normalMatrix[3] = modelviewMatrix[4]; normalMatrix[6] = modelviewMatrix[8];
	normalMatrix[1] = modelviewMatrix[1]; normalMatrix[4] = modelviewMatrix[5]; normalMatrix[7] = modelviewMatrix[9];
	normalMatrix[2] = modelviewMatrix[2]; normalMatrix[5] = modelviewMatrix[6]; normalMatrix[8] = modelviewMatrix[10];



	switch (renderMode){
	case 0:
	{
		//GLSL Region
		NMF->enable();
		glVertexArrayVertexBuffer(VAO, 0, NMbuffer, 0, (GLsizei)(sizeof(float) * 4));
		glVertexArrayAttribFormat(VAO, 0, 4, GL_FLOAT, GL_FALSE, 0);
		glEnableVertexArrayAttrib(VAO, 0);

		GLfloat attrib[] = { (float)sin(t)*0.5f, (float)cos(t)*0.6f, 0.0f, 0.0f };


		glVertexAttrib4fv(10, attrib);
		glVertexAttrib4fv(11, color);
		glDrawArrays(GL_QUADS, 0, 4);

		NMF->disable();
		//GLSL Region
		break;
	}
	case 1:
	{
		drawCube2();
		break;
	}
	case 2://vMF
	{
		float maxtest=0.f;

		for(int j = 0; j < numLobes; j++)
		{
			
			std::cout << "mu[" << j << "] = (" 
				<< mu[j][0] << ", " 
				<< mu[j][1] << ", " 
				<< mu[j][2] << ")" << std::endl;
			std::cout << "kappa[" << j << "] = "
				<< kappa[j] << std::endl;
			std::cout << "alpha[" << j << "] = "
				<< alpha[j] << std::endl << std::endl;
		}
		BYTE *NMFpixels = (BYTE*)malloc(width*height * 3);
		float sum = 0.f;
		float alphasum = 0.f;
		for (int j = 0; j < numLobes; j++)
		{
			alphasum += alpha[j];
		}
		std::cout << "alpha sum=" << alphasum << std::endl;
		for (int j = 0; j < height; j++)
		{
			for (int i = 0; i < width; i++)
			{
				float currentXY[3] = { 0 };
				currentXY[0] = ((float)i / (float)(width / 2)) - 1.f;
				currentXY[1] = ((float)j / (float)(height / 2)) - 1.f;
				NMFpixels[j*height * 3 + i * 3 + 0] = 255;
				if (1.f >= ((currentXY[0] * currentXY[0]) + (currentXY[1] * currentXY[1])))
				{
					currentXY[2] = sqrt(1 - ((currentXY[0] * currentXY[0]) + (currentXY[1] * currentXY[1])));
					float prob = 0.f;
					for (int k = 0; k < numLobes; k++)
					{
						if (maxtest < vMF(currentXY, mu[k], kappa[k]))
							maxtest = vMF(currentXY, mu[k], kappa[k]);
						prob += alpha[k] * vMF(currentXY, mu[k], kappa[k]);
					}
					sum += prob;

					NMFpixels[j*height * 3 + i * 3 + 0] = std::min(255, (int)(prob * 255));
					NMFpixels[j*height * 3 + i * 3 + 1] = std::min(255, (int)(prob * 255));
					NMFpixels[j*height * 3 + i * 3 + 2] = std::min(255, (int)(prob * 255));
				}
				else
				{
					NMFpixels[j*height * 3 + i * 3 + 0] = 0;
					NMFpixels[j*height * 3 + i * 3 + 1] = 0;
					NMFpixels[j*height * 3 + i * 3 + 2] = 200;
				}
			}
		}
		std::cout << "sinh(30)=" << std::sinh(30) << std::endl;
		std::cout << "MAX of prob="<<maxtest<<std::endl;
		glDrawPixels(width, height, GL_RGB, GL_UNSIGNED_BYTE, NMFpixels);
		free(NMFpixels);

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
//		NMFvMF->setUniform3f("eyePos", eyePos[0], eyePos[1], eyePos[2]);
		glVertexArrayVertexBuffer(VAO, 0, vMFvertex, 0, (GLsizei)(sizeof(float) * 4));
		glVertexArrayVertexBuffer(VAO, 1, vMFtex, 0, (GLsizei)(sizeof(float) * 4));
		glVertexArrayVertexBuffer(VAO, 2, vMFnormal, 0, (GLsizei)(sizeof(float) * 4));
		glVertexArrayVertexBuffer(VAO, 3, vMFtangent, 0, (GLsizei)(sizeof(float) * 4));
		glVertexArrayAttribFormat(VAO, 0, 4, GL_FLOAT, GL_FALSE, 0);
		glVertexArrayAttribFormat(VAO, 1, 4, GL_FLOAT, GL_FALSE, 0);
		glVertexArrayAttribFormat(VAO, 2, 4, GL_FLOAT, GL_TRUE, 0);
		glVertexArrayAttribFormat(VAO, 3, 4, GL_FLOAT, GL_TRUE, 0);
		glEnableVertexArrayAttrib(VAO, 0);
		glEnableVertexArrayAttrib(VAO, 1);
		glEnableVertexArrayAttrib(VAO, 2);
		glEnableVertexArrayAttrib(VAO, 3);

//		glVertexArrayVertexBuffer(TexCoordArray, 0, vMFtex, 0, (GLsizei)(sizeof(float) * 4));
//		glVertexArrayAttribFormat(TexCoordArray, 0, 4, GL_FLOAT, GL_FALSE, 0);
//		glEnableVertexArrayAttrib(TexCoordArray, 0);

		glDrawArrays(GL_QUADS, 0, 4);

		NMFvMF->disable();
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
		glVertexArrayVertexBuffer(VAO, 0, vMFvertex, 0, (GLsizei)(sizeof(float) * 4));
		glVertexArrayVertexBuffer(VAO, 1, vMFtex, 0, (GLsizei)(sizeof(float) * 4));
		glVertexArrayVertexBuffer(VAO, 2, vMFnormal, 0, (GLsizei)(sizeof(float) * 4));
		glVertexArrayVertexBuffer(VAO, 3, vMFtangent, 0, (GLsizei)(sizeof(float) * 4));
		glVertexArrayAttribFormat(VAO, 0, 4, GL_FLOAT, GL_FALSE, 0);
		glVertexArrayAttribFormat(VAO, 1, 4, GL_FLOAT, GL_FALSE, 0);
		glVertexArrayAttribFormat(VAO, 2, 4, GL_FLOAT, GL_TRUE, 0);
		glVertexArrayAttribFormat(VAO, 3, 4, GL_FLOAT, GL_TRUE, 0);
		glEnableVertexArrayAttrib(VAO, 0);
		glEnableVertexArrayAttrib(VAO, 1);
		glEnableVertexArrayAttrib(VAO, 2);
		glEnableVertexArrayAttrib(VAO, 3);

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

				NMFpixels[pixelIndex[1] * width * 3 + pixelIndex[0] * 3 + 0] += 1.f / (float)(NMheight*NMwidth);
				NMFpixels[pixelIndex[1] * width * 3 + pixelIndex[0] * 3 + 1] += 1.f / (float)(NMheight*NMwidth);
				NMFpixels[pixelIndex[1] * width * 3 + pixelIndex[0] * 3 + 2] += 1.f / (float)(NMheight*NMwidth);


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
void mouseCB(int button, int state, int x, int y){
	if ((renderMode == 2) || (renderMode == 4)){}
	else
	{
		if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN){
			click_pos[0] = (x - width / 2.0f) / (width / 2.0f);
			click_pos[1] = (y - height / 2.0f) / (height / 2.0f);
		}
		else if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN){
			click_pos[0] = (x - width / 2.0f) / (width / 2.0f);
			click_pos[1] = (y - height / 2.0f) / (height / 2.0f);
		}
		button_id = button;
		glutPostRedisplay();
	}

}
void keyboardCB(unsigned char key, int x, int y){
	if (renderMode == 3)
	{
		switch (key){
		case 's':
			renderScene == 0 ? (renderScene = 1) : (renderScene = 0);
			break;
		case 'm':
			(MipMapped == 0) ? (MipMapped = 1) : (MipMapped = 0);
//			if (MipMapped == 0) { std::cout << "NotMipMapped\n"; } else { std::cout << "MipMapped\n"; }
			(MipMapped == 0) ? (std::cout << "NotMipMapped\n") : (std::cout << "MipMapped\n");
			break;
		}
	}
	glutPostRedisplay();
}
void motionCB(int x, int y){
	GLfloat present[2];
	if ((renderMode == 2) || (renderMode == 4)){}
	else
	{
		switch (button_id){
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
			if (present[1] - click_pos[1] < 0)
				cam->dollyin();
			else if (present[1] - click_pos[1] > 0)
				cam->dollyout();

			click_pos[0] = present[0];
			click_pos[1] = present[1];
			break;
		}
		glutPostRedisplay();
	}
}

void mainMenu(int op) {
	switch (op) {
	case 0:
	{
		delete[] kappa;
		delete[] alpha;
		for (int i = 0; i < numLobes; i++)
		{
			delete[] aux[i];
			delete[] mu[i];

		}
		delete[] mu;
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
	glutCreateWindow("GLEW Test");
	glutDisplayFunc(displayCB);
	glutIdleFunc(idleCB);
	glutMouseFunc(mouseCB);
	glutKeyboardFunc(keyboardCB);
	glutMotionFunc(motionCB);

	int RenderMenu = glutCreateMenu(renderMenu);
	glutAddMenuEntry("GLSL", 0);
	glutAddMenuEntry("Render Scene(vMF)", 3);
	glutAddMenuEntry("Texture", 1);
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



	initBuffers();
	cam = new Camera();

}



FIBITMAP* LoadImage(const char* filename, int &imageWidth, int &imageHeight){
	FREE_IMAGE_FORMAT format = FreeImage_GetFileType(filename);
	RGBQUAD pixel;
	if (format == -1)
	{
		std::cout << "Could not find image: " << filename << " - Aborting." << std::endl;
		sysPause
		exit(-1);
	}
	if (format == FIF_UNKNOWN)
	{
		std::cout << "Couldn't determine file format - attempting to get from file extension..." << std::endl;

		format = FreeImage_GetFIFFromFilename(filename);

		if (!FreeImage_FIFSupportsReading(format))
		{
			std::cout << "Detected image format cannot be read!" << std::endl;
			sysPause
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


int main(int argc, char** argv){
	initSharedMem();

	objData = new objLoader();
	objData->load("./cube.obj");

	NMTdata=LoadImage(NMT, NMwidth, NMheight);

	printf("Number of vertices: %i\n", objData->vertexCount);
	printf("Number of vertex normals: %i\n", objData->normalCount);
	printf("Number of texture coordinates: %i\n", objData->textureCount);
	printf("\n");

	initGLUT(argc, argv);
	initGL(argc, argv);

	printf("Running Well\n");

	glutMainLoop();
	sysPause;
	return 0;

}