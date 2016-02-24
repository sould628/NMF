#version 450 core

layout (location = 10) uniform mat4 mv_matrix;
layout (location = 11) uniform mat4 proj_matrix;

layout (location = 0) in vec4 position;
layout (location = 1) in vec4 texCoord;



//passed to frag.shader

out VS_OUT{

	vec4 texCoord;
} vs_out;

void main(void)
{
	const vec4 vertices[4] = vec4[](vec4(0.5, -0.5, 0.5, 1.0), vec4(-0.5, -0.5, 0.5, 1.0), vec4(-0.5, 0.5, 0.5, 1.0), vec4(0.5, 0.5, 0.5, 1.0));
	const vec4 texCoord[4] = vec4[](vec4(0.0, 0.0, 0.0, 0.0), vec4(0.0, 5.0, 0.0, 0.0), vec4(5.0, 5.0, 0.0, 0.0), vec4(5.0, 0.0, 0.0, 0.0));
	//gl_VertexID: indicated by "glDrawArrays(GL_TRIANGLES, 0, 3);" in c++ with gl_VertexID from 0 to 4
	
	vs_out.texCoord=texCoord[gl_VertexID];
	
	gl_Position=proj_matrix*mv_matrix*vertices[gl_VertexID];
//	gl_Position=vertices[gl_VertexID];

}