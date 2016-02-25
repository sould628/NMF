#pragma once

#include <iostream>
#include <algorithm>
#include <cmath>
#include <GL\glut.h>


#define PI 3.141592f



class Camera{

private:
	void QuaterMultiply(float a[4], float b[4], float result[4])
	{
		result[0] = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
		result[1] = a[1] * b[0] + a[0] * b[1] - a[3] * b[2] + a[2] * b[3];
		result[2] = a[2] * b[0] + a[0] * b[2] + a[3] * b[1] - a[1] * b[3];
		result[3] = a[3] * b[0] + a[0] * b[3] - a[2] * b[1] + a[1] * b[2];
	}
	float QuaterDot(float a[4], float b[4]){
		float result = 0;
		for (int i = 0; i<4; i++){
			result += a[i] * b[i];
		}
		return result;
	}
	void QuaterInverse(float a[4], float b[4]){
		float len = QuaterNorm(a);
		for (int i = 0; i<4; i++){
			if (i == 0) b[i] = a[i] / len;
			else b[i] = -a[i] / len;
		}
	}
	float QuaterNorm(float a[4]){
		return sqrtf(QuaterDot(a, a));
	}
	void QuaterNormalize(float a[4]){
		float magni = QuaterNorm(a);
		for (int i = 0; i<4; i++){
			a[i] = a[i] / magni;
		}
	}
	void AngleToQuater(float angle[4], float quater[4]){
		if (angle[0] == 0){
			quater[0] = 1;
			for (int i = 0; i<3; i++)
				quater[i + 1] = 0;

		}
		else {
			quater[0] = cosf(angle[0] / 2);
			quater[1] = angle[1] * sinf(angle[0] / 2);
			quater[2] = angle[2] * sinf(angle[0] / 2);
			quater[3] = angle[3] * sinf(angle[0] / 2);
		}
	}

private:
	float select_z(float x, float y)
	{
		float num = x*x + y*y;
		if (num < 1) return sqrt(1 - num);
		else return 0;
	}
	void normalize_vec(float point[3])
	{
		float norm = sqrt(vec_dot(point, point));
		for (int i = 0; i<3; i++)
			point[i] /= norm;
	}
	float vec_dot(float vec1[3], float vec2[3])
	{
		float sum = 0;
		for (int i = 0; i<3; i++)
			sum += vec1[i] * vec2[i];
		return sum;
	}
	void vec_cross(float vec1[3], float vec2[3], float result[3])
	{
		result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
		result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
		result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
	}
	void make_quaternion(float quater[4], float vec1[3], float vec2[3], float x_axis[3], float y_axis[3], float z_axis[3])
	{
		float axis_angle[4];
		float temp_vec[3];
		float new_vec[3];

		normalize_vec(vec1);
		normalize_vec(vec2);
		axis_angle[0] = -acos(std::min((float)1.0, vec_dot(vec1, vec2) / sqrt(vec_dot(vec1, vec1) * vec_dot(vec2, vec2))));

		vec_cross(vec1, vec2, temp_vec);
		normalize_vec(temp_vec);
		for (int i = 0; i<3; i++) new_vec[i] = x_axis[i] * temp_vec[0] + y_axis[i] * temp_vec[1] + z_axis[i] * temp_vec[2];
		normalize_vec(new_vec);
		for (int i = 0; i<3; i++) axis_angle[i + 1] = new_vec[i];

		AngleToQuater(axis_angle, quater);
		QuaterNormalize(quater);
	}
	void rotate(float quaternion[4], float ori[3])
	{
		float inverse[4];
		float result[4];
		float ori_quater[4];
		float temp_quaternion[4];
		ori_quater[0] = 0;
		for (int i = 0; i<3; i++) ori_quater[i + 1] = ori[i];

		QuaterInverse(quaternion, inverse);

		QuaterMultiply(quaternion, ori_quater, temp_quaternion);
		QuaterMultiply(temp_quaternion, inverse, result);
		for (int i = 0; i<3; i++) ori[i] = result[i + 1];
	}


public:
	Camera();
	void glRender();
	void getPosition(float* pos);
	void setPosition(float x, float y, float z);
	void setLookat(float x, float y, float z);
	void setUp(float x, float y, float z);
	void setViewAngle(float viewAngle_);
	void setZNear(float zNear_);
	void setZFar(float zFar_);
	void trackball(GLfloat src[2], GLfloat dst[2]);
	void dollyin();
	void dollyout();
private:
	float ori_vec[3];
	float position[3];
	float lookat[3], up[3];
	float viewAngle;
	float zNear, zFar;

};
void Camera::getPosition(float* pos){
	pos[0] = position[0];
	pos[1] = position[1];
	pos[2] = position[2];

}
Camera::Camera(){
	ori_vec[0] = 0;
	ori_vec[1] = 0;
	ori_vec[2] = 3;
	lookat[0] = lookat[1] = lookat[2] = 0;
	up[0] = up[2] = 0;
	up[1] = 1;
	viewAngle = 60;
	zNear = 0.1;
	zFar = 100;
}
void Camera::setPosition(float x, float y, float z){
	position[0] = x;
	position[1] = y;
	position[2] = z;
}
void Camera::setLookat(float x, float y, float z){
	lookat[0] = x;
	lookat[1] = y;
	lookat[2] = z;
}
void Camera::setUp(float x, float y, float z){
	up[0] = x;
	up[1] = y;
	up[2] = z;
}
void Camera::setViewAngle(float viewAngle_){
	viewAngle = viewAngle_;
}
void Camera::setZNear(float zNear_){
	zNear = zNear_;
}
void Camera::setZFar(float zFar_){
	zFar = zFar_;
}

void Camera::glRender(){
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	
	int w = viewport[2];
    int h = viewport[3];
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	gluPerspective(viewAngle,float(w)/float(h), zNear, zFar);
    glMatrixMode(GL_MODELVIEW);
	for(int i=0; i<3; i++) {
		position[i] = lookat[i] + ori_vec[i];
	}
    glLoadIdentity();
    gluLookAt(position[0], position[1], position[2], lookat[0], lookat[1], lookat[2], up[0], up[1], up[2]);
}

void Camera::trackball(GLfloat src[2], GLfloat dst[2]){
	//virtual trackball function
	GLfloat p1[3], p2[3];
	GLfloat quaternion[4];
	for(int i=0; i<2;++i){
		p1[i] = src[i];
		p2[i] = dst[i];
	}	
	
	p1[2] = select_z(p1[0], p1[1]);
	p2[2] = select_z(p2[0], p2[1]);
	
	normalize_vec(p1);
	normalize_vec(p2);
	float x_axis[3], y_axis[3], z_axis[3];
	for(int i=0; i<3; i++){
		z_axis[i] = -ori_vec[i];
		y_axis[i] = up[i];
	}
	normalize_vec(z_axis);
	normalize_vec(y_axis);

	vec_cross(y_axis, z_axis, x_axis);
	normalize_vec(x_axis);
	vec_cross(z_axis, x_axis, y_axis);
	normalize_vec(y_axis);

	make_quaternion(quaternion, p1, p2,x_axis, y_axis, z_axis);

	rotate(quaternion, ori_vec);
	rotate(quaternion, up);

	normalize_vec(up);
}


void Camera::dollyin(){	
	for(int i=0; i<3; i++){
		ori_vec[i] *= 0.9;
	}

}

void Camera::dollyout(){
	for(int i=0; i<3; i++){
		ori_vec[i] *= 1.1;
	}

}