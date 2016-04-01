#ifndef RENDERER_H
#define RENDERER_H

#include "GL/glew.h"
#include "GLFW/glfw3.h"
#include "glm/glm.hpp"
#include <vector>

#include "noise.h"
#include "camera.h"
#include "shaders.h"

using namespace glm;

struct VAO
{
	enum{ DRAWING, MODEL, COUNT};
};

struct VBO
{
	enum { DRAWING_VERTS, DRAWING_COLOR, MODEL_VERTS, MODEL_NORMALS, MODEL_INDICES, COUNT};
};

struct Shader
{
	enum { DRAWING, MODEL, COUNT };
};

class Renderer
{
private:
	//Members
	GLFWwindow* window;
	Camera* cam;
	mat4 transform;
	mat4 modelview;

	vector<vec3>* modelVertices;
	vector<vec3>* modelNormals;
	vector<unsigned int>* modelIndices;

	vector<vec2>* drawingVertices;
	vector<vec3>* drawingColors;

	GLuint vao[VAO::COUNT];
	GLuint vbo[VBO::COUNT];
	GLuint shader[Shader::COUNT];

	//Functions
	void generateIDs();	
	void deleteIDs();

	void setupVAOs();
	void loadShaders();

	void loadModelUniforms();
	void loadDrawingUniforms();

	void updateTransform();

public:
	//Members
	vec3 light;
	mat4 projection;

	//Functions
	Renderer(GLFWwindow* window);

	void loadModelBuffer(vector<vec3>* vertices, vector<vec3>* normals, vector<unsigned int>* indices);

	void loadDrawBuffer(vector<vec2>* vertices, vector<vec3>* color);

	void loadCamera(Camera* cam);

	void render2DView();
	void render3DView();

};



void resizeEvent(GLFWwindow* window, int width, int height);

void generatePlane(unsigned int widthSegments, unsigned int depthSegments, float width, float depth, vector<vec3>* points, vector<vec3>* normals, vector<unsigned int>* indices);


#endif