#version 450 core

layout (binding=0) uniform sampler2D original;
layout (binding=5) uniform sampler2D vMFmap5;
layout (binding=4) uniform sampler2D vMFmap4;
layout (binding=3) uniform sampler2D vMFmap3;
layout (binding=2) uniform sampler2D vMFmap2;
layout (binding=1) uniform sampler2D vMFmap1;


layout (location = 100) uniform int numLobes;
layout (location = 101) uniform float BPexp;


//input from vertex shader
in VS_OUT{

	vec4 texCoord;
}fs_in;
//output into color buffer frame
out vec4 color;

float lightPos[4]={1.0f, 1.0f, 1.0f, 1.0f};
float lightIntensity[4]={1.0f, 1.0f, 1.0f, 1.0f};



void main(void)
{
	
	vec3 effBRDF=vec3(0.0,0.0,0.0);

	color=vec4(0.0, 0.0, 0.0, 0.0);
	for(int i=0; i<numLobes; i++)
	{
		vec4 theta=vec4(0.0,0.0,0.0,0.0);
		float alpha=0.0;
		vec3 aux=vec3(0.0,0.0,0.0);
		float kappa=0.0;
		vec3 mu=vec3(0.0,0.0,0.0);
		float sPrime=0.0;
		float Bs=0.0;
		vec3 halfVec=vec3(0.0,0.0,0.0);
		vec3 
		switch(i)
		{
			case 0:
			{
				theta=texture2D(vMFmap1, fs_in.texCoord.xy);
				alpha=theta.x;
				aux=theta.yzw/alpha;
				kappa=((3*length(aux))+(length(aux)*length(aux)*length(aux)));
				mu=normalize(aux);

				sPrime=(kappa*BPexp)/(kappa+BPexp);
//				Bs=


				break;
			}
			case 1:
			{
				
			}
		}
	}

	

	color+=texture2D(original, fs_in.texCoord.xy);

//	color=vec4(fs_in.texCoord.xy, 0.0, 0.0);
//	color=vec4(numLobes/, 0, 0, 0);
}