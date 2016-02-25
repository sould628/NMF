#version 450 core

layout (binding=0) uniform sampler2D original;
layout (binding=5) uniform sampler2D vMFmap5;
layout (binding=4) uniform sampler2D vMFmap4;
layout (binding=3) uniform sampler2D vMFmap3;
layout (binding=2) uniform sampler2D vMFmap2;
layout (binding=1) uniform sampler2D vMFmap1;
layout (binding=9) uniform sampler2D originalMipMap;

layout (location = 10) uniform mat4 mv_matrix;


layout (location = 100) uniform int numLobes;
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
	vec4 texCoord;
}fs_in;
//output into color buffer frame
out vec4 color;


vec4 lightIntensity=vec4(10.0f, 10.f, 10.f, 1.0f);

vec4 Kd=vec4(0.4f,0.4f,0.4f,1.0f);
vec4 Ks=vec4(0.3f,0.3f,0.3f,1.0f);

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
//	vec3 lightAngle=normalize((mv_matrix*lightPos).xyz-fs_in.p.xyz);

	vec3 toeye=normalize(-fs_in.p.xyz);
	vec3 h=normalize((toeye+lightDir)/2);

	v.x=dot(fs_in.eyePos, fs_in.t);
	v.y=dot(fs_in.eyePos, fs_in.b);
	v.z=dot(fs_in.eyePos, fs_in.n);

	eyeDir=normalize(v);

	v.x=dot(fs_in.lightPos.xyz, fs_in.t);
	v.y=dot(fs_in.lightPos.xyz, fs_in.b);
	v.z=dot(fs_in.lightPos.xyz, fs_in.n);

	lightDir=normalize(v);

	h=normalize(-eyeDir+lightDir);
	
	vec4 coeffs[6];
	coeffs[0]=texture2D(vMFmap1, fs_in.texCoord.xy);
	coeffs[1]=texture2D(vMFmap2, fs_in.texCoord.xy);
	coeffs[2]=texture2D(vMFmap3, fs_in.texCoord.xy);
	coeffs[3]=texture2D(vMFmap4, fs_in.texCoord.xy);
	coeffs[4]=texture2D(vMFmap5, fs_in.texCoord.xy);

	for(int i=0; i<numLobes; i++)
	{

		float alpha=0.0;
		vec3 aux=vec3(0.0,0.0,0.0);
		float kappa=0.0;
		vec3 mu=vec3(0.0,0.0,0.0);
		float sPrime=0.0;
		float Bs=0.0;
		vec3 halfVec=vec3(0.0,0.0,0.0);
		float r=0.0;
//		alpha=theta.x;
//		aux=theta.yzw/alpha;
//		kappa=((3*length(aux))+(length(aux)*length(aux)*length(aux)));
//		mu=normalize(aux);
//
//		sPrime=(kappa*BPexp)/(kappa+BPexp);
//		float HdotMu=dot(h,mu);
//		Bs=((sPrime+1)/(2*PI))*pow(HdotMu, sPrime);
//
//		effBRDF+=effBRDF+(alpha*(Ks*Bs+Kd)*(dot(lightAngle,mu)));

		alpha=coeffs[i].x;
		mu=2.0*coeffs[i].yzw-1.0;
		mu/=alpha;
		r=length(mu);
		mu=normalize(mu);
		kappa=(3.0*r-r*r*r)/(1.0-r*r);

		sPrime=kappa*BPexp/(kappa+BPexp);
		float norm=(sPrime+1.0)/(2.0*PI);
		Bs=dot(h,mu);
		Bs=pow(Bs,sPrime);
		Bs*=norm;
		float LdotMu=max(dot(lightDir,mu), 0.0);

		effBRDF.xyz+=effBRDF.xyz+alpha*((Ks.xyz*Bs)+(Kd.xyz*LdotMu));
//		color+=vec4(LdotMu,LdotMu,LdotMu,0);
	}

	switch(renderScene){
		case 0:
		{
			color=(lightIntensity*effBRDF);
			break;
		}
		case 1:
		{
			switch(MipMapped)
			{
				case 0:
				{
					color+=texture2D(original, fs_in.texCoord.xy);
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


//	color=vec4(h,0.0);

//	color+=vec4(h,0.0);

//	color=vec4(lightDir, 0.0);
//	color=vec4(numLobes/, 0, 0, 0);
}