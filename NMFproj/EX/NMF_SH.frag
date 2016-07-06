#version 450 core

layout (binding=0) uniform sampler2D original;
layout (binding=18) uniform sampler2D SHmap8;
layout (binding=17) uniform sampler2D SHmap7;
layout (binding=16) uniform sampler2D SHmap6;
layout (binding=15) uniform sampler2D SHmap5;
layout (binding=14) uniform sampler2D SHmap4;
layout (binding=13) uniform sampler2D SHmap3;
layout (binding=12) uniform sampler2D SHmap2;
layout (binding=11) uniform sampler2D SHmap1;

layout (binding=28) uniform sampler2D Ylmmap8;
layout (binding=27) uniform sampler2D Ylmmap7;
layout (binding=26) uniform sampler2D Ylmmap6;
layout (binding=25) uniform sampler2D Ylmmap5;
layout (binding=24) uniform sampler2D Ylmmap4;
layout (binding=23) uniform sampler2D Ylmmap3;
layout (binding=22) uniform sampler2D Ylmmap2;
layout (binding=21) uniform sampler2D Ylmmap1;

layout (location = 50) uniform float b0;
layout (location = 51) uniform float b1;
layout (location = 52) uniform float b2;
layout (location = 53) uniform float b3;
layout (location = 54) uniform float b4;
layout (location = 55) uniform float b5;
layout (location = 56) uniform float b6;
layout (location = 57) uniform float b7;
layout (location = 58) uniform float b8;
layout (location = 59) uniform float b9;

layout (binding=9) uniform sampler2D originalMipMap;

//Blinn-Phong Coeff


layout (location = 10) uniform mat4 mv_matrix;


layout (location = 100) uniform int order;
layout (location = 101) uniform float BPexp;

layout (location = 1000) uniform int renderScene;
layout (location = 1001) uniform int MipMapped;

const float PI = 3.141592653589793238462643383;


//input from vertex shader
in VS_OUT{
	vec3 lightPos;
	vec3 eyePos;
	vec3 n;
	vec3 t;
	vec3 b;
	vec4 p;
	vec2 texCoord;
}fs_in;
//output into color buffer frame
out vec4 color;


vec4 lightIntensity=vec4(1.0f, 1.0f, 1.f, 1.0f);

vec4 Kd=vec4(0.1f,0.f,0.f,1.0f);
vec4 Ks=vec4(0.5f,0.5f,0.5f,1.0f);
vec4 Ka=vec4(0.0f,0.0f,0.0f,1.0f);

void main(void)
{
	vec4 effBRDF=vec4(0.0,0.0,0.0,1.0);
	vec3 v=vec3(0.0, 0.0, 0.0);
	color=vec4(0.0, 0.0, 0.0, 0.0);
	
	vec3 lightDir=vec3(0.0, 0.0, 0.0);
	vec3 eyeDir=vec3(0.0, 0.0, 0.0);

	//Calculating Vectors
	//lightDir: w_i
	//toeye: w_o
	//h:half

	v.x=dot(fs_in.eyePos, fs_in.t);
	v.y=dot(fs_in.eyePos, fs_in.b);
	v.z=dot(fs_in.eyePos, fs_in.n);

	eyeDir=normalize(v);
	//diriectional (lightPos =(x,x,x,0.0) in vert shader)
//	v.x=dot(fs_in.lightPos.xyz, fs_in.t);
//	v.y=dot(fs_in.lightPos.xyz, fs_in.b);
//	v.z=dot(fs_in.lightPos.xyz, fs_in.n);

	//point (lightPos = (x,x,x,1.0) in vert shader)
	v.x=dot(fs_in.lightPos.xyz-fs_in.p.xyz, fs_in.t);
	v.y=dot(fs_in.lightPos.xyz-fs_in.p.xyz, fs_in.b);
	v.z=dot(fs_in.lightPos.xyz-fs_in.p.xyz, fs_in.n);

	lightDir=normalize(v);
	vec3 h=normalize(-eyeDir+lightDir);
	
	float theta_h = acos(h.z);
	float phi_h = atan(h.y, h.x);

	vec2 ylmCoord=vec2(h.x/2.+0.5, h.y/2.+0.5);

//	lightDir=normalize((mv_matrix*vec4(fs_in.lightPos,1)).xyz-fs_in.p.xyz);
//	lightDir=normalize((vec4(fs_in.lightPos,1)).xyz-fs_in.p.xyz);
	vec3 toeye=normalize(-fs_in.p.xyz);


	vec4 coeffs[10];
	coeffs[0]=texture2D(SHmap1, fs_in.texCoord.xy);
	coeffs[1]=texture2D(SHmap2, fs_in.texCoord.xy);
	coeffs[2]=texture2D(SHmap3, fs_in.texCoord.xy);
	coeffs[3]=texture2D(SHmap4, fs_in.texCoord.xy);
	coeffs[4]=texture2D(SHmap5, fs_in.texCoord.xy);
	coeffs[5]=texture2D(SHmap6, fs_in.texCoord.xy);
	coeffs[6]=texture2D(SHmap7, fs_in.texCoord.xy);
	coeffs[7]=texture2D(SHmap8, fs_in.texCoord.xy);

	vec4 Ylmcoeffs[10];
	Ylmcoeffs[0]=texture2D(Ylmmap1, ylmCoord);
	Ylmcoeffs[1]=texture2D(Ylmmap2, ylmCoord);
	Ylmcoeffs[2]=texture2D(Ylmmap3, ylmCoord);
	Ylmcoeffs[3]=texture2D(Ylmmap4, ylmCoord);
	Ylmcoeffs[4]=texture2D(Ylmmap5, ylmCoord);
	Ylmcoeffs[5]=texture2D(Ylmmap6, ylmCoord);
	Ylmcoeffs[6]=texture2D(Ylmmap7, ylmCoord);
	Ylmcoeffs[7]=texture2D(Ylmmap8, ylmCoord);

	float brdfCoeff[10];
	brdfCoeff[0]=b0;
	brdfCoeff[1]=b1;
	brdfCoeff[2]=b2;
	brdfCoeff[3]=b3;
	brdfCoeff[4]=b4;
	brdfCoeff[5]=b5;
	brdfCoeff[6]=b6;
	brdfCoeff[7]=b7;
	brdfCoeff[8]=b8;

	int idx=0;
	int count=0;
	float YlmCoeff;
	float ndfCoeff;
	float Bs=0.;
	for(int l=0; l<order; l++)
	{

		for(int m=-l; m<=l; m++)
		{
			if(count==4)
			{
				idx+=1;
				count=0;
			}
			YlmCoeff=Ylmcoeffs[idx][count];
			ndfCoeff=coeffs[idx][count];
			Bs+=brdfCoeff[l]*ndfCoeff*YlmCoeff;
			count++;
		}
	}
//	float LdotN=max(dot(lightDir,mu),0.0);
	effBRDF+=((Ks*Bs+Kd));
	switch(renderScene){
		case 0:
		{
			color=vec4(0., 0., 0., 0.);
			color=(lightIntensity*effBRDF);
//			color=(vec4(texture2D(Ylmmap1, fs_in.texCoord.xy)));
//			color=vec4(ylmCoord.x, ylmCoord.y, 0., 0.);
//			color=vec4(1.0, 1.0, 1.0, 1.0);
			break;
		}
		case 1:
		{
			switch(MipMapped)
			{
				case 0:
				{
					color+=texture2D(SHmap1, fs_in.texCoord.xy);
					break;
				}
				case 1:
				{
					color+=texture2D(originalMipMap, fs_in.texCoord.xy);
					break;
				}
			}
		}
	}
//	color.rgb=texture2D(vMFmap1, fs_in.texCoord.xy).gba;
//	color.a=1;
//	color.rgb=textureLod(vMFmap1, fs_in.texCoord.xy, 1).gba;
//	color.a=1;
}