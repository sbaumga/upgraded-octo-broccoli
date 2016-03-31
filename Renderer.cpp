#include "Renderer.h"

mat4 winRatio(1.f);

unsigned int WIN_WIDTH;
unsigned int WIN_HEIGHT;

void updateWinRatio()
{
	float minDimension = std::min((float)WIN_WIDTH, (float)WIN_HEIGHT);
	winRatio[0][0] = minDimension / ((float)WIN_WIDTH);
	winRatio[1][1] = minDimension / ((float)WIN_HEIGHT);
}

void glErrorCheck(const char* location)
{
	GLenum code = glGetError();

	switch (code)
	{
	case GL_INVALID_ENUM:
		cout << "GL_INVALID_ENUM - " << location << endl;
		break;
	case GL_INVALID_VALUE:
		cout << "GL_INVALID_VALUE - " << location << endl;
		break;
	case GL_INVALID_OPERATION:
		cout << "GL_INVALID_OPERATION - " << location << endl;
		break;
	case GL_INVALID_FRAMEBUFFER_OPERATION:
		cout << "GL_INVALID_FRAMEBUFFER_OPERATION - " << location << endl;
		break;
	case GL_OUT_OF_MEMORY:
		cout << "GL_OUT_OF_MEMORY - " << location << endl;
		break;
	}
}

Renderer::Renderer(GLFWwindow* window) :window(window), light(normalize(vec3(1.f, 1.f, 1.f))), modelVertices(NULL), modelNormals(NULL), modelIndices(NULL), drawingVertices(NULL), drawingColors(NULL)
{
	generateIDs();
	setupVAOs();
	loadShaders();

	GLint vp[4];
	glGetIntegerv(GL_VIEWPORT, vp);
	WIN_WIDTH = vp[2];
	WIN_HEIGHT = vp[3];
}

void Renderer::generateIDs()
{
	glErrorCheck("Begin generate IDs");

	// load IDs given from OpenGL
	glGenVertexArrays(VAO::COUNT, vao);
	glGenBuffers(VBO::COUNT, vbo);

	glErrorCheck("End generate IDs");
}

void Renderer::loadShaders()
{
	string vertSource = loadShaderStringFromFile("./src/drawing.vert");
	string fragSource = loadShaderStringFromFile("./src/drawing.frag");
	shader[Shader::DRAWING] = CreateShaderProgram(vertSource, fragSource);

	vertSource = loadShaderStringFromFile("./src/model.vert");
	fragSource = loadShaderStringFromFile("./src/model.frag");
	shader[Shader::MODEL] = CreateShaderProgram(vertSource, fragSource);
}

void Renderer::setupVAOs()
{
	glErrorCheck("Pre-setup VAO");

	//Drawing vao
	glBindVertexArray(vao[VAO::DRAWING]);
	//Vertex array
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, vbo[VBO::DRAWING_VERTS]);
	glVertexAttribPointer(
		0,				//Attribute
		2,				//Size
		GL_FLOAT,		//Type
		GL_FALSE,		//Normalized
		sizeof(vec2),	//Stride
		(void*)0		//Offset
		);
	//Normal array
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(
		1,				//Attribute
		3,				//Size
		GL_FLOAT,		//Type
		GL_FALSE,		//Normalized
		sizeof(vec3),	//Stride
		(void*)0		//Offset
		);

	//Model vao
	glBindVertexArray(vbo[VBO::MODEL_VERTS]);
	//Vertex array
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, vbo[VBO::MODEL_VERTS]);
	glVertexAttribPointer(
		0,				//Attribute
		3,				//Size
		GL_FLOAT,		//Type
		GL_FALSE,		//Normalized
		sizeof(vec3),	//Stride
		(void*)0		//Offset
		);
	//Normal array
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(
		1,				//Attribute
		3,				//Size
		GL_FLOAT,		//Type
		GL_FALSE,		//Normalized
		sizeof(vec3),	//Stride
		(void*)0		//Offset
		);

	glErrorCheck("Setup VAO");
}

void Renderer::loadModelUniforms()
{
	GLuint uniformLocation = glGetUniformLocation(shader[Shader::MODEL], "transform");
	glUniformMatrix4fv(uniformLocation, 1, false, &transform[0][0]);

	uniformLocation = glGetUniformLocation(shader[Shader::MODEL], "objectTransform");
	glUniformMatrix4fv(uniformLocation, 1, false, &modelview[0][0]);

	uniformLocation = glGetUniformLocation(shader[Shader::MODEL], "light");
	glUniform3f(uniformLocation, light.x, light.y, light.z);

	uniformLocation = glGetUniformLocation(shader[Shader::MODEL], "viewPos");
	glUniform3f(uniformLocation, cam->getPos().x, cam->getPos().y, cam->getPos().z);

}

void Renderer::loadDrawingUniforms()
{
	GLuint uniformLocation = glGetUniformLocation(shader[Shader::MODEL], "transform");
}

void Renderer::loadModelBuffer(vector<vec3>* vertices, vector<vec3>* normals, vector<unsigned int>* indices)
{
	modelVertices = vertices;
	modelNormals = normals;
	modelIndices = indices;

	glBindVertexArray(vao[VAO::MODEL]);
	glBindBuffer(GL_ARRAY_BUFFER, vbo[VBO::MODEL_VERTS]);

	glBufferData(GL_ARRAY_BUFFER,
		sizeof(vec3)*vertices->size(),
		vertices->data(),
		GL_STATIC_DRAW
		);

	glBindBuffer(GL_ARRAY_BUFFER, vbo[VBO::MODEL_NORMALS]);
	glBufferData(GL_ARRAY_BUFFER,
		sizeof(vec3)*normals->size(),
		normals->data(),
		GL_STATIC_DRAW
		);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[VBO::MODEL_INDICES]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		sizeof(unsigned int)*indices->size(),
		indices->data(),
		GL_STATIC_DRAW
		);

	glBindVertexArray(0);
}

void Renderer::loadDrawBuffer(vector<vec2>* vertices, vector<vec3>* color)
{
	drawingVertices = vertices;
	drawingColors = color;

	glBindVertexArray(vao[VAO::MODEL]);
	glBindBuffer(GL_ARRAY_BUFFER, vbo[VBO::MODEL_VERTS]);

	glBufferData(GL_ARRAY_BUFFER,
		sizeof(vec2)*vertices->size(),
		vertices->data(),
		GL_STATIC_DRAW
		);

	glBindBuffer(GL_ARRAY_BUFFER, vbo[VBO::MODEL_NORMALS]);
	glBufferData(GL_ARRAY_BUFFER,
		sizeof(vec3)*color->size(),
		color->data(),
		GL_STATIC_DRAW
		);

	glBindVertexArray(0);
}



void Renderer::loadCamera(Camera* _cam)
{
	cam = _cam;
}

void Renderer::render2DView()
{
	glErrorCheck("Begin render 2D");

	glUseProgram(shader[Shader::DRAWING]);

	loadDrawingUniforms();

	glBindVertexArray(vao[VAO::DRAWING]);

	glDrawArrays(GL_LINES, 0, drawingVertices->size());

	glErrorCheck("End render 2D");
}

void Renderer::render3DView()
{
	glErrorCheck("Begin render 3D");

	glUseProgram(shader[Shader::MODEL]);

	loadModelUniforms();

	glBindVertexArray(vao[VAO::MODEL]);

	glDrawElements(GL_TRIANGLES, modelIndices->size(), GL_UNSIGNED_INT, modelIndices);

	glErrorCheck("End render 3D");
}

void resizeEvent(GLFWwindow* window, int width, int height)
{

	WIN_WIDTH = width;
	WIN_HEIGHT = height;

	glViewport(0, 0, width, height);

	updateWinRatio();
}