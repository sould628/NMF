#pragma once
#include<iostream>

#include <GL/glew.h>
#include <GL/glut.h>


class GLSLProgram {
private:
	GLuint program, vertShader, tcShader, fragShader;


public:

	char* readFromFile(const char* filename) {
		FILE* fp = fopen(filename, "rt");
		if (fp == NULL) { std::cout << "no file exists" << std::endl; return nullptr; }
		fseek(fp, 0, SEEK_END);
		size_t count = ftell(fp);
		rewind(fp);
		if (count > 0)
		{
			char* text = (char*)malloc(count + 1);
			count = (int)fread(text, 1, count, fp);
			text[count] = '\0';
			fclose(fp);
			return text;
		}
	}

	GLuint createShader(const char* src, GLenum type) {
		if (src == NULL)
			return -1;
		GLuint shader = glCreateShader(type);
		glShaderSource(shader, 1, &src, NULL);
		glCompileShader(shader);
		printShaderInfoLog(shader);
		return shader;
	}

	void printShaderInfoLog(GLuint obj) {
		int len = 0, charsWritten = 0;
		char* infoLog;
		glGetShaderiv(obj, GL_INFO_LOG_LENGTH, &len);
		if (len > 0) {
			infoLog = (char*)malloc(len);
			glGetShaderInfoLog(obj, len, &charsWritten, infoLog);
			printf("%s\n", infoLog);
			free(infoLog);
		}
	}
	void printProgramInfoLog(GLuint obj) {
		int len = 0, charsWritten = 0;
		char* infoLog;
		glGetProgramiv(obj, GL_INFO_LOG_LENGTH, &len);
		if (len > 0) {
			infoLog = (char*)malloc(len);
			glGetProgramInfoLog(obj, len, &charsWritten, infoLog);
			printf("%s\n", infoLog);
			free(infoLog);
		}
	}


public:
	GLSLProgram() {}
	GLSLProgram(char *vertfileName, char *fragfileName) {

		char* vert = readFromFile(vertfileName);
		char* frag = readFromFile(fragfileName);

		vertShader = createShader(vert, GL_VERTEX_SHADER);
		fragShader = createShader(frag, GL_FRAGMENT_SHADER);

		program = glCreateProgram();
		glAttachShader(program, vertShader);
		glAttachShader(program, fragShader);
		glLinkProgram(program);
		printProgramInfoLog(program);
		free(vert);
		free(frag);
	}
	void GLSLCreateShader(char *vertfileName, char*fragfileName)
	{
		char* vert = readFromFile(vertfileName);
		char* frag = readFromFile(fragfileName);

		vertShader = createShader(vert, GL_VERTEX_SHADER);
		fragShader = createShader(frag, GL_FRAGMENT_SHADER);

		program = glCreateProgram();
		glAttachShader(program, vertShader);
		glAttachShader(program, fragShader);
		free(vert);
		free(frag);
	}
	void GLSLLinkShader()
	{
		glLinkProgram(program);
		printProgramInfoLog(program);
	}
public:
	GLuint getProgram() { return program; }
	GLuint getVertShader() { return vertShader; }
	GLuint getFragShader() { return fragShader; }
	void enable() { glUseProgram(program); }
	void disable() { glUseProgram(0); }
	void setUniform1i(const GLchar *name, GLint x) {
		GLint loc = glGetUniformLocation(program, name);
		glUniform1i(loc, x);
	}
	void setUniform1f(const GLchar *name, GLfloat x) {
		GLint loc = glGetUniformLocation(program, name);
		glUniform1f(loc, x);
	}
	void setUniform2f(const GLchar *name, GLfloat x, GLfloat y) {
		GLint loc = glGetUniformLocation(program, name);
		glUniform2f(loc, x, y);
	}
	void setUniform3f(const char *name, float x, float y, float z) {
		GLint loc = glGetUniformLocation(program, name);
		glUniform3f(loc, x, y, z);
	}
	void setUniform4f(const char *name, float x, float y, float z, float w) {
		GLint loc = glGetUniformLocation(program, name);
		glUniform4f(loc, x, y, z, w);
	}
	void setUniform4fv(const char *name, float *v) {
		GLint loc = glGetUniformLocation(program, name);
		glUniform4fv(loc, 4, v);
	}
	void setUniformfv(const GLchar *name, GLfloat *v, int elementSize, int count = 1) {
		GLint loc = glGetUniformLocation(program, name);
		glUniform1fv(loc, elementSize, v);
	}
	void SetUniformMatrix4fv(const GLchar *name, GLfloat *m, bool transpose) {
		GLint loc = glGetUniformLocation(program, name);
		glUniformMatrix4fv(loc, 1, transpose, m);
	}

	void SetUniformMatrix3fv(const GLchar *name, GLfloat *m, bool transpose) {
		GLint loc = glGetUniformLocation(program, name);
		glUniformMatrix3fv(loc, 1, transpose, m);
	}

	void bindTexture(const char *name, GLuint tex, GLenum target, GLint unit) {
		glActiveTexture(GL_TEXTURE0_ARB + unit);
		glBindTexture(target, tex);
		GLint loc = glGetUniformLocation(program, name);
		glUniform1i(loc, unit);
	}

};