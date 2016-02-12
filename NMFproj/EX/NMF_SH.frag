#version 450 core

layout (binding=5) uniform sampler2D vMFmap5;
layout (binding=4) uniform sampler2D vMFmap4;
layout (binding=3) uniform sampler2D vMFmap3;
layout (binding=2) uniform sampler2D vMFmap2;
layout (binding=1) uniform sampler2D vMFmap1;

//input from vertex shader
in VS_OUT{
vec4 color;
}fs_in;
//output into color buffer frame
out vec4 color;


void main(void)
{

	color=fs_in.color;
	color=200*texture2D(vMFmap1, vec2(0.0, 0.0));

}