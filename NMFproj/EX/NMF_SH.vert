#version 450 core

layout (location=10) in vec4 offset;
layout (location=11) in vec4 color;

layout (location=0)

//passed to frag.shader

out VS_OUT{
	vec4 color;
	vec2 texCoord[4];
} vs_out;

void main(void)
{
	int i=0;
	gl_VertexID;

	const vec4 vertices[4] = vec4[](vec4(0.2, -0.2, 0.5, 1.0), vec4(-0.2, -0.2, 0.5, 1.0), vec4(-0.2, 0.2, 0.5, 1.0), vec4(0.2, 0.2, 0.5, 1.0));
	vec2 texCoord[4]=vec2[](vec2(0.0, 0.0),vec2(0.0, 1.0),vec2(1.0, 1.0),vec2(0.0, 0.0));
	for(i=0; i<4; i++)
	{
	 vs_out.texCoord[i]=texCoord[i];
	}
	//gl_VertexID: indicated by "glDrawArrays(GL_TRIANGLES, 0, 3);" in c++ with gl_VertexID from 0 to 4
	gl_Position=vertices[gl_VertexID]+offset;

	vs_out.color=color;
}