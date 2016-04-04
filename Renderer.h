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

	mat4 transform2D;
	
	mat4 panning2D;
	mat4 scaling2D;

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
	void loadDrawingUniforms(vec3 color);

public:
	//Members
	vec3 light;
	mat4 projection;

	//Functions
	Renderer() {}		//Temporary, not sufficient
	Renderer(GLFWwindow* window);

	void loadModelBuffer(vector<vec3>* vertices, vector<vec3>* normals, vector<unsigned int>* indices);

	void loadDrawBuffer(vector<vec2>* vertices);

	void loadCamera(Camera* cam);

	void updateTransform();

	void render2DView(vec3 color);
	void render3DView();

	mat4 getWinRatio();

};

mat4 perspectiveMatrix(float n, float f, float fov);

void resizeEvent(GLFWwindow* window, int width, int height);

void generatePlane(unsigned int widthSegments, unsigned int depthSegments, float width, float depth, vector<vec3>* points, vector<vec3>* normals, vector<vec2>* uvs, vector<unsigned int>* indices);


#endif