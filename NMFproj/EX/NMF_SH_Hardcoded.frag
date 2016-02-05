#version 450 core


//input from vertex shader
in VS_OUT{
vec4 color;
}fs_in;
//output into color buffer frame
out vec4 color;


void main(void)
{
	color=fs_in.color;
}