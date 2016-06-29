#version 450 core

layout (location = 10) uniform mat4 mv_matrix;
layout (location = 11) uniform mat4 proj_matrix;
layout (location = 12) uniform mat3 normal_matrix;

layout (location = 0) in vec4 position;
layout (location = 1) in vec4 texCoord;
layout (location = 2) in vec4 normalVector;
layout (location = 3) in vec4 tangent;


//passed to frag.shader



out VS_OUT{
	vec3 lightPos;
	vec3 eyePos;
	vec3 n;
	vec3 t;
	vec3 b;
	vec4 p;
	vec2 texCoord;
} vs_out;

void main(void)
{
	
	double texSize=50.0;

	double modelSize=0.5;
	vec4 lightPos=vec4(0.0, 0.0f, 10.0f, 1.0f);
	const vec4 vertices[4] = vec4[](vec4(modelSize, modelSize, 0.0, 1.0), vec4(-modelSize, modelSize, 0.0, 1.0), vec4(-modelSize, -modelSize, 0.0, 1.0), vec4(modelSize, -modelSize, 0.0, 1.0));

	const vec3 normalVector[4] = vec3[](vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0),vec3(0.0, 0.0, 1.0),vec3(0.0, 0.0, 1.0)); 
	const vec2 texCoord[4] = vec2[](vec2(0.0, 0.0), vec2(0.0, texSize), vec2(texSize, texSize), vec2(texSize, 0.0));
	const vec3 tangent[4]= vec3[](vec3(1.0, 0.0, 0.0),vec3(0.0, 1.0, 0.0),vec3(-1.0, 0.0, 0.0),vec3(0.0, -1.0, 0.0));
	//gl_VertexID: indicated by "glDrawArrays(GL_TRIANGLES, 0, 3);" in c++ with gl_VertexID from 0 to 4
	
	vs_out.n=normalize(normal_matrix*normalVector[gl_VertexID]).xyz;
	vs_out.t=normalize(normal_matrix*tangent[gl_VertexID]).xyz;
	vs_out.b=normalize(cross(vs_out.n.xyz,vs_out.t.xyz)).xyz;
	vs_out.texCoord=texCoord[gl_VertexID];
	
	vs_out.eyePos=vec3((mv_matrix*vertices[gl_VertexID]).xyz);
	vs_out.lightPos=vec3(mv_matrix*lightPos);
	vs_out.lightPos=lightPos.xyz;

	vs_out.p=mv_matrix*vertices[gl_VertexID];
	gl_Position=proj_matrix*mv_matrix*vertices[gl_VertexID];
//	gl_Position=vertices[gl_VertexID];

}