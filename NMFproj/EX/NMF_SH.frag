#version 450 core

layout (binding=0) uniform sampler2D s;

//input from vertex shader
in VS_OUT{
vec4 color;
}fs_in;
//output into color buffer frame
out vec4 color;


void main(void)
{

	color=fs_in.color;
	color=texture2D(s, vec2(0.0, 0.0));

}