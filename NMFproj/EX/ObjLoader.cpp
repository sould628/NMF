#pragma warning(disable : 4996)

#include <fstream>
#include <iostream>
#include <sstream>
#include "ObjLoader.h"

using namespace std;
ObjLoader::ObjLoader(){
	maxX=maxY=maxZ=-FLT_MAX;
	minX=minY=minZ=FLT_MAX;
}
void ObjLoader::loadObjModel(char* filename)
{
    fstream file;

	string line;
    file.open(filename,fstream::in);
    if(!file.is_open())
    {
    	cout << "Model File Not Found: " << endl;    	
    	return;
    }
    uint n=0;
    vector3 v;
    TCoord vt;
    vector3 vn;
    Triangle t;
    _loaded=true;
	while(!file.eof())
	{
		getline(file,line);
		if(line.length()>2)
		switch(line.at(0))
		{
			case 'v':
					//vertex or tex coord
					if(line.at(1)==' ')
					{
						//vertex
						n=sscanf(line.c_str(),"v %f %f %f",&v.x,&v.y,&v.z);
						if(n!=3)
						{
							v.x=0;
							v.y=0;
							v.z=0;
							_loaded=false;
						}
						if(v.x > maxX) maxX = v.x;
						if(v.x < minX) minX = v.x;
						if(v.y > maxY) maxY = v.y;
						if(v.y < minY) minY = v.y;
						if(v.z > maxZ) maxZ = v.z;
						if(v.z < minZ) minZ = v.z;
						vert.push_back(v);
					}
					else if(line.at(1)=='n')
					{
						//normal
						n=sscanf(line.c_str(),"vn %f %f %f",&vn.x,&vn.y,&vn.z);
						if(n!=3)
						{
							vn.x=0;
							vn.y=0;
							vn.z=0;
							_loaded=false;
						}
						norm.push_back(vn);
					}
					else if(line.at(1)=='t')
					{
						//texture coord
						n=sscanf(line.c_str(),"vt %f %f",&vt.t,&vt.s);
						if(n!=2)
						{
							vt.t=0;
							vt.s=0;
							_loaded=false;
						}
						tc.push_back(vt);
					}
			break;

			case 'f':
				//face
				n=sscanf(line.c_str(),"f %d/%d/%d "
									"%d/%d/%d "
									"%d/%d/%d",
									&t.v[0],&t.vt[0],&t.vn[0],
									&t.v[1],&t.vt[1],&t.vn[1],
									&t.v[2],&t.vt[2],&t.vn[2]
									);
				if(n!=9)
				{
					t.v[0]=0;
					t.v[1]=0;
					t.v[2]=0;
					
					t.vt[0]=0;
					t.vt[1]=0;
					t.vt[2]=0;
					
					t.vn[0]=0;
					t.vn[1]=0;
					t.vn[2]=0;

					_loaded=false;
				}
				for(int i=0;i<3;i++)
				{
					t.v[i]-=1;
					t.vt[i]-=1;
					t.vn[i]-=1;
				}
				tri.push_back(t);
			break;
			
			default:

			break;
		}
	}
    file.close();
}


void ObjLoader::draw()
{	
	glMaterialfv(GL_FRONT, GL_SPECULAR, matKs);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, matKd);
	glMaterialfv(GL_FRONT, GL_AMBIENT, matKa);
	glMaterialfv(GL_FRONT, GL_SHININESS, &matShininess);

	glBegin(GL_TRIANGLES);
		for(int i=0;i<tri.size();i++){
			Triangle &t=tri.at(i);
			vector3 &v0=vert.at(t.v[0]);
			vector3 &v1=vert.at(t.v[1]);
			vector3 &v2=vert.at(t.v[2]);
	
			vector3 &vn0=norm.at(t.vn[0]);
			vector3 &vn1=norm.at(t.vn[1]);
			vector3 &vn2=norm.at(t.vn[2]);
	
			TCoord &vt0=tc.at(t.vt[0]);
			TCoord &vt1=tc.at(t.vt[1]);
			TCoord &vt2=tc.at(t.vt[2]);
	
	
			glTexCoord2f(vt0.t,vt0.s);
			glNormal3f(vn0.x,vn0.y,vn0.z);
			glVertex3f(v0.x,v0.y,v0.z);
	
	
			glTexCoord2f(vt1.t,vt1.s);
			glNormal3f(vn1.x,vn1.y,vn1.z);
			glVertex3f(v1.x,v1.y,v1.z);
	
			glTexCoord2f(vt2.t,vt2.s);
			glNormal3f(vn2.x,vn2.y,vn2.z);
			glVertex3f(v2.x,v2.y,v2.z);
		}	
	glEnd();
	
}

void ObjLoader::scale(float scaleX, float scaleY, float scaleZ){
	for(int i=0; i<vert.size(); i++){
		vert[i].x *= scaleX;
		vert[i].y *= scaleY;
		vert[i].z *= scaleZ;
	}
	maxX *= scaleX;
	minX *= scaleX;
	maxY *= scaleY;
	minY *= scaleY;
	maxZ *= scaleZ;
	minZ *= scaleZ;
}


void ObjLoader::translate(float translateX, float translateY, float translateZ){
	for(int i=0; i<vert.size(); i++){
		vert[i].x += translateX;
		vert[i].y += translateY;
		vert[i].z += translateZ;
	}
	maxX += translateX;
	minX += translateX;
	maxY += translateY;
	minY += translateY;
	maxZ += translateZ;
	minZ += translateZ;
}

void ObjLoader::setMatKs(float _matKs[4]){
	for(int i=0;i <4; i++){
		matKs[i] = _matKs[i];
	}
}
void ObjLoader::setMatKd(float _matKd[4]){
	for(int i=0;i <4; i++){
		matKd[i] = _matKd[i];
	}
}
void ObjLoader::setMatKa(float _matKa[4]){
	for(int i=0;i <4; i++){
		matKa[i] = _matKa[i];
	}
}
void ObjLoader::setMatShininess(float _matShininess){
	matShininess = _matShininess;
}

float ObjLoader::getMaxX(){
	return maxX;
}
float ObjLoader::getMinX(){
	return minX;
}
float ObjLoader::getMaxY(){
	return maxY;
}
float ObjLoader::getMinY(){
	return minY;
}
float ObjLoader::getMaxZ(){
	return maxZ;
}
float ObjLoader::getMinZ(){
	return minZ;
}
