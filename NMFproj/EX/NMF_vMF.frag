#version 450 core

layout (binding=0) uniform sampler2D original;
layout (binding=8) uniform sampler2D vMFmap8;
layout (binding=7) uniform sampler2D vMFmap7;
layout (binding=6) uniform sampler2D vMFmap6;
layout (binding=5) uniform sampler2D vMFmap5;
layout (binding=4) uniform sampler2D vMFmap4;
layout (binding=3) uniform sampler2D vMFmap3;
layout (binding=2) uniform sampler2D vMFmap2;
layout (binding=1) uniform sampler2D vMFmap1;
layout (binding=9) uniform sampler2D originalMipMap;



layout (location = 10) uniform mat4 mv_matrix;


layout (location = 100) uniform int numLobes;
layout (location = 101) uniform float BPexp;
layout (location = 102) uniform float MicroFacet;

layout (location = 1000) uniform int renderScene;
layout (location = 1001) uniform int MipMapped;
layout (location = 1002) uniform int brdfSelect;

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


vec4 lightIntensity=vec4(0.9f, 0.9f, 0.9f, 1.0f);

vec4 Kd=vec4(0.2f,0.f,0.f,1.0f);
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
	

//	lightDir=normalize((mv_matrix*vec4(fs_in.lightPos,1)).xyz-fs_in.p.xyz);
//	lightDir=normalize((vec4(fs_in.lightPos,1)).xyz-fs_in.p.xyz);
	vec3 toeye=normalize(-fs_in.p.xyz);
//	h=normalize((toeye+lightDir)/2);


	vec4 coeffs[10];
	coeffs[0]=texture2D(vMFmap1, fs_in.texCoord.xy);
	coeffs[1]=texture2D(vMFmap2, fs_in.texCoord.xy);
	coeffs[2]=texture2D(vMFmap3, fs_in.texCoord.xy);
	coeffs[3]=texture2D(vMFmap4, fs_in.texCoord.xy);
	coeffs[4]=texture2D(vMFmap5, fs_in.texCoord.xy);
	coeffs[5]=texture2D(vMFmap6, fs_in.texCoord.xy);
	coeffs[6]=texture2D(vMFmap7, fs_in.texCoord.xy);
	coeffs[7]=texture2D(vMFmap8, fs_in.texCoord.xy);

	switch(brdfSelect)
	{
	//BlinnPhong
	case 0:
	{
		for(int i=0; i<numLobes; i++)
		{
	
			float alpha=0.0;
			vec3 aux=vec3(0.0,0.0,0.0);
			float kappa=0.0;
			vec3 mu=vec3(0.0,0.0,0.0);
			float sPrime=0.0;
			float Bs=0.0;
			float r=0.0;
	
	
			alpha=coeffs[i].x;
			aux=coeffs[i].yzw/max(alpha,0.01);
			r=length(aux);
			kappa=((3*r)-(r*r*r))/max(0.001, (1.0-(r*r)));
	
			mu=normalize(aux);
	
			sPrime=(kappa*BPexp)/(kappa+BPexp);
			float HdotMu=max(dot(h,mu),0.0);
			Bs=((sPrime+1.0)/(2.0*PI))*pow(HdotMu, sPrime);
			float LdotMu=max(dot(lightDir,mu),0.0);
			effBRDF+=(alpha*(Ks*Bs+Kd)*LdotMu);
		}
		break;
	}
	//MicroFacet
	case 1:
	{
		for(int i=0; i<numLobes; i++)
		{
			float alpha=0.0;
			vec3 aux=vec3(0.0,0.0,0.0);
			float kappa=0.0;
			vec3 mu=vec3(0.0,0.0,0.0);
			float r=0.0;

		}
		break;
	}

	}
	switch(renderScene){
		case 0:
		{
			color=(lightIntensity*effBRDF);
//			color=vec4(1.0, 1.0, 1.0, 1.0);
			break;
		}
		case 1:
		{

			color=vec4(texture2D(vMFmap1, fs_in.texCoord.st))
			+vec4(texture2D(vMFmap2, fs_in.texCoord.st))
			+vec4(texture2D(vMFmap3, fs_in.texCoord.st))
			+vec4(texture2D(vMFmap4, fs_in.texCoord.st))
			+vec4(texture2D(vMFmap5, fs_in.texCoord.st))
			+vec4(texture2D(vMFmap6, fs_in.texCoord.st))
			+vec4(texture2D(vMFmap7, fs_in.texCoord.st))
			+vec4(texture2D(vMFmap8, fs_in.texCoord.st))
			;
			vec4 temp=vec4(color.yza, color.x);
			color=temp;
			break;

		}
	}
//	color.rgb=texture2D(vMFmap1, fs_in.texCoord.xy).gba;
//	color.a=1;
//	color.rgb=textureLod(vMFmap1, fs_in.texCoord.xy, 1).gba;
//	color.a=1;
}