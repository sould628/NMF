#include "Camera.h"

#define MAXITERATION 10

extern const char* NMT = "./1.jpg";
extern const int numLobes = 8;
extern const int textLobes = 4;
extern int MipMapLevel = 3;
extern float alignCtrl = 2.0;


Camera *cam;
extern int renderMode = 3;
extern int renderScene = 0;
extern int MipMapped = 1;

extern const float BPexp = 100.0f;





extern float modelviewMatrix[16] = { 0 };
extern float projectionMatrix[16] = { 0 };
extern float normalMatrix[9] = { 0 };
