#ifndef SHADERS_H
#define SHADERS_H

#include "GL/glew.h"
#include <string>

GLuint CreateShaderProgram(const std::string & vsSource,
	const std::string & fsSource);

GLuint CreateShaderProgram(const std::string & vsSource,
	const std::string & gsSource,
	const std::string & fsSource);

bool checkCompileStatus(GLint shaderID);
bool checkLinkStatus(GLint programID);

std::string loadShaderStringFromFile(const std::string & filePath);













#endif