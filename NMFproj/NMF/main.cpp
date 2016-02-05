//kixor.net


#include <GL/glew.h>
#include <GL/glut.h>

#include <stdio.h>
#include <stdlib.h>


#include "FreeImage.h"
#include "objLoader.h"


#include "GLSLProgram.h"



#include "Camera.h"

#define sysPause system("pause>nul");



int width = 800, height = 800;
int NMwidth, NMheight;


objLoader *objData;
const char* NMT = "./bricks_normal_map.jpg";
FIBITMAP* NMTdata;

Camera *cam;
GLint button_id;
GLfloat click_pos[2];




GLSLProgram *NMF;
const int TEXTURE_WIDTH = 4096;
const int TEXTURE_HEIGHT = 4096;



void initSharedMem(){
	NMwidth = 0; NMheight = 0;
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



#pragma region GLSL
void createProgram() {
	NMF = new GLSLProgram("NMF_SH.vert", "NMF_SH.frag");
}
void initFBO() {
	
}


#pragma endregion







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


void displayCB(){
	glEnable(GL_DEPTH_TEST);

	

	glClearColor(0.f, 0.f, 0.f, 0.f);


	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	cam->glRender();


	
	glViewport(0, 0, width, height);

	drawObj(objData);
//	drawCube();

	glFlush();
}
void mouseCB(int button, int state, int x, int y){
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
void keyboardCB(unsigned char key, int x, int y){
	glutPostRedisplay();
}
void motionCB(int x, int y){
	GLfloat present[2];
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
		if (present[1] - click_pos[1] <0)
			cam->dollyin();
		else if (present[1] - click_pos[1] >0)
			cam->dollyout();

		click_pos[0] = present[0];
		click_pos[1] = present[1];
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
	glutMouseFunc(mouseCB);
	glutKeyboardFunc(keyboardCB);
	glutMotionFunc(motionCB);
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
//	createProgram();
	cam = new Camera();
}


void LoadImage(FIBITMAP* data, const char* filename, int &imageWidth, int &imageHeight){
	FREE_IMAGE_FORMAT format = FreeImage_GetFileType(filename);
	if (format == -1)
	{
		std::cout << "Could not find image: " << filename << " - Aborting." << std::endl;
		exit(-1);
	}
	if (format == FIF_UNKNOWN)
	{
		std::cout << "Couldn't determine file format - attempting to get from file extension..." << std::endl;

		format = FreeImage_GetFIFFromFilename(filename);

		if (!FreeImage_FIFSupportsReading(format))
		{
			std::cout << "Detected image format cannot be read!" << std::endl;
			exit(-1);
		}
	}
	FIBITMAP* bitmap = FreeImage_Load(format, filename);
	int bitsPerPixel = FreeImage_GetBPP(bitmap);

	if (bitsPerPixel == 32)
	{
		std::cout << "Source image has " << bitsPerPixel << " bits per pixel. Skipping conversion." << std::endl;
		data = bitmap;
	}
	else
	{
		std::cout << "Source image has " << bitsPerPixel << " bits per pixel. Converting to 32-bit colour." << std::endl;
		data = FreeImage_ConvertTo32Bits(bitmap);
	}

	imageWidth = FreeImage_GetWidth(data);
	imageHeight = FreeImage_GetHeight(data);
	std::cout << "Image: " << filename << " is size: " << imageWidth << "x" << imageHeight << "." << std::endl;





}


int main(int argc, char** argv){
	initSharedMem();

	objData = new objLoader();
	objData->load("./cube.obj");

	LoadImage(NMTdata, NMT, NMwidth, NMheight);

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