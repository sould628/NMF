#include <string>
#include "Camera.h"


#define MAXITERATION 10

std::string inputObj = "./cylinder_cloth_0070.obj";


//extern const char* NMT = "Velvet_N.jpg";
extern const char* NMT = "1.jpg";
extern const int numLobes = 8;
extern const int textLobes = 4;
extern int MipMapLevel = 3;
extern float alignCtrl = 0.1f;


Camera *cam;
extern int renderMode = 3;
extern int renderScene = 0;
extern int MipMapped = 1;

extern float BPexp = 100.0f;
extern float MicroSigma = 20.f;

extern int brdfSelect = 0;


extern float modelviewMatrix[16] = { 0 };
extern float projectionMatrix[16] = { 0 };
extern float normalMatrix[9] = { 0 };
