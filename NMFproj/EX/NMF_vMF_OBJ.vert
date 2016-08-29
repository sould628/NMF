#version 450 core



layout (location = 10) uniform mat4 mv_matrix;
layout (location = 11) uniform mat4 proj_matrix;
layout (location = 12) uniform mat3 normal_matrix;
layout (location = 16) uniform float texModifier; 

layout (location = 13) uniform float lightPosX;
layout (location = 14) uniform float lightPosY;
layout (location = 15) uniform float lightPosZ;

layout (location = 0) in vec4 in_vertices;
layout (location = 1) in vec3 in_normalVector;
layout (location = 2) in vec2 in_texCoord;
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
	mat3 tbnMatrix;
} vs_out;

void main(void)
{
	

//	vec4 lightPos=vec4(0.0, 30.0f, -100.0f, 1.0f);

	vec4 lightPos = vec4(lightPosX, lightPosY, lightPosZ, 1.f);
	mat3 tbnMatrix;
	tbnMatrix[0]=normalize(in_tangent);
	tbnMatrix[1]=normalize(in_bitangent);
	tbnMatrix[2]=normalize(in_normalVector);

//	tbnMatrix[0]=in_tangent;
//	tbnMatrix[1]=in_bitangent;
//	tbnMatrix[2]=in_normalVector;

	vs_out.tbnMatrix=tbnMatrix;

	vec3 n=normalize(cross(in_bitangent,in_tangent));

	vs_out.n=normalize(mv_matrix*vec4(normalize(in_normalVector), 0.)).xyz;
	vs_out.t=normalize(mv_matrix*vec4(normalize(in_tangent),0.)).xyz;
	vs_out.b=normalize(mv_matrix*vec4(normalize(in_bitangent),0)).xyz;
	vs_out.b=normalize(mv_matrix*vec4(cross(vs_out.n, vs_out.t),0.)).xyz;
	vs_out.t=normalize(mv_matrix*vec4(cross(vs_out.b, vs_out.n), 0.)).xyz;
	vs_out.texCoord=in_texCoord*texModifier;
	vs_out.origNormals=normalize(vec3(in_texCoord, 0.0));

	vs_out.eyePos=vec3((mv_matrix*in_vertices).xyz);
	vs_out.lightPos=vec3(mv_matrix*lightPos);
//	vs_out.lightPos=lightPos.xyz;

	vs_out.p=mv_matrix*in_vertices;
	gl_Position=proj_matrix*mv_matrix*in_vertices;
//	gl_Position=vertices[gl_VertexID];

}