#version 450 core



layout (location = 10) uniform mat4 mv_matrix;
layout (location = 11) uniform mat4 proj_matrix;
layout (location = 12) uniform mat3 normal_matrix;

layout (location = 0) in vec4 in_vertices;
layout (location = 1) in vec2 in_texCoord;
layout (location = 2) in vec3 in_normalVector;
layout (location = 3) in vec3 in_tangent;
layout (location = 4) in vec3 in_bitangent;

//passed to frag.shader



out VS_OUT{
	vec3 lightPos;
	vec3 eyePos;
	vec3 n;
	vec3 t;
	vec3 b;
	vec4 p;
	vec2 texCoord;
	vec3 origNormals;
} vs_out;

void main(void)
{
	
	double texSize=100.0;

	double modelSize=0.5;
	vec4 lightPos=vec4(0.0, 100.0f, 0.0f, 1.0f);

	vs_out.n=normalize(normal_matrix*in_normalVector).xyz;
	vs_out.t=normalize(normal_matrix*in_tangent).xyz;
	vs_out.b=normalize(cross(vs_out.n.xyz,vs_out.t.xyz)).xyz;
	vs_out.texCoord=in_texCoord;
	vs_out.origNormals=in_normalVector;

	vs_out.eyePos=vec3((mv_matrix*in_vertices).xyz);
	vs_out.lightPos=vec3(mv_matrix*lightPos);
	vs_out.lightPos=lightPos.xyz;

	vs_out.p=mv_matrix*in_vertices;
	gl_Position=proj_matrix*mv_matrix*in_vertices;
//	gl_Position=vertices[gl_VertexID];

}