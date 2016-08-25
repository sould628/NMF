#version 450 core

layout (binding=0) uniform sampler2D original;
layout (binding=9) uniform sampler2D vMFmap9;
layout (binding=8) uniform sampler2D vMFmap8;
layout (binding=7) uniform sampler2D vMFmap7;
layout (binding=6) uniform sampler2D vMFmap6;
layout (binding=5) uniform sampler2D vMFmap5;
layout (binding=4) uniform sampler2D vMFmap4;
layout (binding=3) uniform sampler2D vMFmap3;
layout (binding=2) uniform sampler2D vMFmap2;
layout (binding=1) uniform sampler2D vMFmap1;
layout (binding=19) uniform sampler2D originalMipMap;



layout (location = 10) uniform mat4 mv_matrix;


layout (location = 100) uniform int numLobes;
//layout (location = 101) uniform float BPexp;
//layout (location = 102) uniform float MicroSigma;

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


vec4 lightIntensity=vec4(1.f, 1.f, 1.f, 1.0f);

vec4 Kd=vec4(0.2f,0.f,0.f,1.0f);
vec4 Ks=vec4(0.5f,0.5f,0.5f,1.0f);
vec4 Ka=vec4(0.0f,0.0f,0.0f,1.0f);

float _Schlick(float refIndex, float incAngle);
float Schlick(float refIndex, float LdotH);

float BPexp=300.f;
float MicroSigma=0.0000001;
float refractiveIdx=1.557;	

//1.557


void main(void)
{

	vec4 effBRDF=vec4(0.0,0.0,0.0,1.0);
	vec3 v=vec3(0.0, 0.0, 0.0);
	color=vec4(0.0, 1.0, 0.0, 0.0);
	
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
//	vec3 toeye=normalize(-fs_in.p.xyz);
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

	vec4 orig_n;
	orig_n=texture2D(originalMipMap, fs_in.texCoord.xy);

	vec4 effBRDF_orig=vec4(0., 0., 0., 0.);
	float LdotN=max(dot(lightDir, orig_n.xyz), 0.0);
	float HdotN=max(dot(h, orig_n.xyz), 0.0);
	float Bspec=(BPexp+1.0)/(2.0*PI)*pow(HdotN, BPexp);
	effBRDF_orig=Kd*LdotN+Ks*Bspec;
	float HdotMu=0.f;
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
			HdotMu=max(dot(h,mu),0.0);
			Bs=((sPrime+1.0)/(2.0*PI))*pow(HdotMu, sPrime);
			float LdotMu=max(dot(lightDir,mu),0.0);
			effBRDF+=(alpha*(Ks*Bs+Kd)*LdotMu);
		}
		break;
	}	
	//MicroFacet
	case 1:
	{
		float fresnel=0.;
		float alpha=0.0;
		vec3 aux=vec3(0.0,0.0,0.0);
		float kappa=0.0;
		vec3 mu=vec3(0.0,0.0,0.0);
		float r=0.0;
		float Ms=0.0;
		float mPrime=0.;
		HdotMu=0.;

		for(int i=0; i<numLobes; i++)
		{	
			alpha=coeffs[i].x;
			aux=coeffs[i].yzw/max(alpha,0.0001);
			r=length(aux);
			kappa=min(700., ((3*r)-(r*r*r))/max(0.0001, (1.0-(r*r))));	
			mu=normalize(aux);

			float sigmaPrime=0.;
			sigmaPrime=sqrt((MicroSigma*MicroSigma)+(1.0/(2.0*kappa)));
			
			vec3 localh=vec3(dot(h, fs_in.t), dot(h, fs_in.b), dot(h, fs_in.n));
			localh=normalize(localh);

			float theta_h=acos(h.z);
			HdotMu=max(dot(mu, h), 0.0);
			theta_h=acos(HdotMu);
			float theta_i=acos(lightDir.z);
			float theta_o=acos(-eyeDir.z);

//			theta_h=acos(localh.z);
			theta_i=acos(lightDir.z);
			theta_o=acos(eyeDir.z);

			float LdotH=max(dot(lightDir, h),0.0);
			fresnel=_Schlick(refractiveIdx, lightDir.z);
//			fresnel=Schlick(refractiveIdx, LdotH);

			Ms=(1./(PI*sigmaPrime*sigmaPrime))*exp(-(theta_h/sigmaPrime)*(theta_h/sigmaPrime));
			float LdotMu=max(dot(lightDir,mu),0.0);
			effBRDF+=(alpha*(Kd+(Ks*Ms*fresnel/(4*lightDir.z*(-eyeDir.z))))*LdotMu);
//			effBRDF+=(alpha*(Kd+(Ks*fresnel/(4*lightDir.z*(-eyeDir.z))))*LdotMu);
//			effBRDF=vec4(theta_h/sigmaPrime, 0., 0., 0.);
		}


		break;
	}

	}
	switch(renderScene){
		case 0:
		{
			color=vec4(0., 0., 0., 0.);
			color=(lightIntensity*effBRDF);
//			color.xyz=fs_in.t;
//			color=vec4(1.0, 1.0, 1.0, 1.0);
//			color=vec4(HdotMu, HdotMu, HdotMu,1.0);
//			color=vec4(HdotMu, 0.f, 0.f, 0.f);
//			color=vec4(-eyeDir, 0.f);
			break;
		}
		case 1:
		{
			color=vec4(0., 0., 0., 0.);
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
//			color=vec4(fract(fs_in.texCoord.x), fract(fs_in.texCoord.y), 0., 1.);
			float mipmapLevel=textureQueryLod(vMFmap1,fs_in.texCoord).x;
			switch(int(mipmapLevel))
			{
			case 0: color=vec4(1., 0., 0., 1.); break;
			case 1: color=vec4(0., 1., 0., 1.); break;
			case 2: color=vec4(0., 0., 1., 1.); break;
			case 3: color=vec4(1., 1., 0., 1.); break;
			case 4: color=vec4(1., 0., 1., 1.); break;
			case 5: color=vec4(0., 1., 1., 1.); break;
			case 6: color=vec4(1., 0.5, 1., 1.); break;
			case 7: color=vec4(1., 1., 1., 1.); break;
			case 8: color=vec4(0.5, 0.5, 1., 1.); break;

			}
			color=orig_n;
			color=effBRDF_orig*lightIntensity;
			break;
		}
	}
//	color.rgb=texture2D(vMFmap1, fs_in.texCoord.xy).gba;
//	color.a=1;
//	color.rgb=textureLod(vMFmap1, fs_in.texCoord.xy, 1).gba;
//	color.a=1;
}

float _Schlick(float n, float cosT)
{
	float ret=0.;
	ret=((n-1)*(n-1)+4*n*(1-cosT)*(1-cosT)*(1-cosT)*(1-cosT)*(1-cosT))/((n+1)*(n+1));
	return ret;
}
float Schlick(float refIndex, float LdotH)
{
//f0=((1-n)/(1+n))^2
//F(f0)=f0+(1-f0)(1-LdotH)^5
	float ret=0.;
	float f0=(((1-refIndex)*(1-refIndex))/((1+refIndex)*(1+refIndex)));
	ret = f0+(1-f0)*pow(1-LdotH, 5);

	return ret;
}

