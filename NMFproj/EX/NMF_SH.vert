#version 450 core

layout (location=10) in vec4 offset;
layout (location=11) in vec4 color;

layout (location=0)

//passed to frag.shader

out VS_OUT{
	vec4 color;
} vs_out;

void main(void)
{
	
	gl_VertexID;

	const vec4 vertices[4] = vec4[](vec4(0.2, -0.2, 0.5, 1.0), vec4(-0.2, -0.2, 0.5, 1.0), vec4(-0.2, 0.2, 0.5, 1.0), vec4(0.2, 0.2, 0.5, 1.0));
	//gl_VertexID: indicated by "glDrawArrays(GL_TRIANGLES, 0, 3);" in c++ with gl_VertexID from 0 to 4
	gl_Position=vertices[gl_VertexID]+offset;

	vs_out.color=color;
}