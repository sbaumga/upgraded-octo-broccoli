#include <iostream>
#include <GL/glew.h>  
#include <GLFW/glfw3.h>  
#include <glm/glm.hpp>
#include <stdio.h>  
#include <stdlib.h>  
#include <vector>
#include <cmath>

#define _USE_MATH_DEFINES
#include <math.h>
using namespace std;
using namespace glm;

enum
{
	GROUND,
	WATER,
	TREES,
};

GLFWwindow *window;
int w, h;
double mouseX, mouseY;

vector < int > lineTypes;
vector < vector < vec2 > > lines;
vector < vec2 > points;
bool drawing;
int drawType = GROUND;


void renderDrawing() {
	glBegin(GL_LINES);
	

	if (lines.size() >= 1) {
		for (int i = 0; i < lines.size(); i++) {
			vector < vec2 > line = lines[i];
			int type = lineTypes[i];
			if (type == GROUND) {
				// brown
				glColor3f(0.4, 0.2, 0.0);
			}
			else if (type == WATER) {
				glColor3f(0.0, 0.0, 1.0);
			}
			else if (type == TREES) {
				glColor3f(0.0, 1.0, 0.0);
			}
			for (int j = 0; j < line.size() - 1; j++) {
				glVertex2d(line[j].x, line[j].y);
				glVertex2d(line[j + 1].x, line[j + 1].y);
			}
		}
	}

	if (points.size() > 1) {
		if (drawType == GROUND) {
			// brown
			glColor3f(0.4, 0.2, 0.0);
		}
		else if (drawType == WATER) {
			glColor3f(0.0, 0.0, 1.0);
		}
		else if (drawType == TREES) {
			glColor3f(0.0, 1.0, 0.0);
		}
		for (int i = 0; i < points.size() - 1; i++) {
			glVertex2d(points[i].x, points[i].y);
			glVertex2d(points[i + 1].x, points[i + 1].y);
		}
	}

	glEnd();
}

/*
 * Renders and calculates what to render
 */
void render() {
	glEnable(GL_DEPTH_TEST);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//Functions for changing transformation matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//Functions for changing projection matrix
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1, 1, -1, 1, -1, 1);
	//gluPerspective (fov, aspect ratio, near plane, far plane)
	//glFrustum

	renderDrawing();
}

GLFWwindow* initializeWindow()
{

	//Initialize GLFW  
	if (!glfwInit())
	{
		exit(EXIT_FAILURE);
	}

	//Set the GLFW window creation hints - these are optional  
	//glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); //Request a specific OpenGL version  
	//glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3); //Request a specific OpenGL version  
	//glfwWindowHint(GLFW_SAMPLES, 4); //Request 4x antialiasing  
	//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  

	//Declare a window object  
	GLFWwindow* window;

	//Create a window and create its OpenGL context  
	window = glfwCreateWindow(750, 750, "B-Spline", NULL, NULL);

	//If the window couldn't be created  
	if (!window)
	{
		fprintf(stderr, "Failed to open GLFW window.\n");
		glfwTerminate();
		exit(EXIT_FAILURE);
	}

	//This function makes the context of the specified window current on the calling thread.   
	glfwMakeContextCurrent(window);

	//Initialize GLEW  
	GLenum err = glewInit();

	//If GLEW hasn't initialized  
	if (err != GLEW_OK)
	{
		fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
		return NULL;
	}

	return window;

}

void keyboard(GLFWwindow *sender, int key, int scancode, int action, int mods) {
	if (key == GLFW_KEY_G && action == GLFW_PRESS) {
		drawType = GROUND;
		if (drawing) {
			lineTypes.back() = drawType;
		}
	}
	else if (key == GLFW_KEY_W && action == GLFW_PRESS) {
		drawType = WATER;
		if (drawing) {
			lineTypes.back() = drawType;
		}
	}
	else if (key == GLFW_KEY_T && action == GLFW_PRESS) {
		drawType = TREES;
		if (drawing) {
			lineTypes.back() = drawType;
		}
	}
}

void mouseClick(GLFWwindow *sender, int button, int action, int mods) {
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		drawing = true;
		lineTypes.push_back(drawType);
	}
	else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
		drawing = false;
		lines.push_back(points);
		points.clear();
	}
}

void mousePos(GLFWwindow *sender, double x, double y) {	
	mouseX = (2 * x / w) - 1;
	mouseY = (-2 * y / h) + 1;
	
	if (drawing) {
		if (points.size() > 0) {
			cout << x << ", " << points.back().x << "\n";
			double distanceX = abs(mouseX - points.back().x);
			double distanceY = abs(mouseY - points.back().y);
			cout << distanceX << "\n";
			if (distanceX >= 0.01 || distanceY >= 0.01) {
				points.push_back(vec2(mouseX, mouseY));
			}
		}
		else {
			points.push_back(vec2(mouseX, mouseY));
		}
	}
}

void scroll(GLFWwindow *sender, double x, double y) {
	
}

int main() {
	window = initializeWindow();
	if (window == NULL) {
		return -1;
	}

	drawing = false;

	glfwSetKeyCallback(window, keyboard);
	glfwSetMouseButtonCallback(window, mouseClick);
	glfwSetCursorPosCallback(window, mousePos);
	glfwSetScrollCallback(window, scroll);
	while (!glfwWindowShouldClose(window)) {
		glfwGetFramebufferSize(window, &w, &h);
		glViewport(0, 0, w, h);

		render();

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}
