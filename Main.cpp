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
vector < vector < vec2 > > controlPoints;
vector < vec2 > points;
vector < vector < double > > knots;
vector < vector < vec2 > > plottedPoints;
int order = 4;
bool drawing;
int drawType = GROUND;


/*
* Finds delta for the given u using knots[lineNum]
*/
int findDelta(double u, int lineNum) {
	for (int i = 0; i < knots[lineNum].size() - order; i++) {
		if (u >= knots[lineNum][i] && u < knots[lineNum][i + 1]) {
			return i;
		}
	}
	return knots[lineNum].size() - order;
}


/*
* Calculates the x and y coords of a B-spline curve at point u and pushes them to vectors
* Uses deBoors algorithm
* Interpolated from pseudo-code shown in class
*/
void calculateBSpline(double u, int lineNum) {
	int delta = findDelta(u, lineNum);
	vector<double> x, y;

	for (int i = 0; i <= order - 1; i++) {
		x.push_back(controlPoints[lineNum][delta - i].x);
		y.push_back(controlPoints[lineNum][delta - i].y);
	}

	for (int r = order; r >= 2; r--) {
		int i = delta;
		for (int s = 0; s <= r - 2; s++) {
			double denom = (knots[lineNum][i + r - 1] - knots[lineNum][i]);
			double omega;
			if (denom == 0) {
				omega = 0;
			}
			else {
				omega = (u - knots[lineNum][i]) / denom;
			}
			x[s] = omega*x[s] + (1 - omega) * x[s + 1];
			y[s] = omega*y[s] + (1 - omega) * y[s + 1];
			i--;
		}
	}

	plottedPoints[lineNum].push_back(vec2(x[0], y[0]));
}

/*
* Recalculates a standard knot sequence whenever a control point is added or removed
*/
void recalculateKnots(int lineNum) {
	knots[lineNum].clear();

	if (order <= controlPoints[lineNum].size()) {
		int numKnots = controlPoints[lineNum].size() - (order - 2);
		if (numKnots != 0) {
			// Setting up standard knot sequence
			for (int i = 0; i < order - 1; i++) {
				knots[lineNum].push_back(0);

			}

			for (int i = 0; i < numKnots; i++) {
				knots[lineNum].push_back((double)i / (numKnots - 1));
			}

			// Setting up standard knot sequence
			for (int i = 0; i < order - 1; i++) {
				knots[lineNum].push_back(1);
			}
		}
	}
}

void renderDrawing() {
	glBegin(GL_LINES);
	

	if (plottedPoints.size() >= 1) {
		for (int i = 0; i < plottedPoints.size(); i++) {
			vector < vec2 > line = plottedPoints[i];
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
		controlPoints.push_back(points);

		vector <double> temp;
		knots.push_back(temp);
		recalculateKnots(controlPoints.size() - 1);
		vector <vec2> temp2;
		plottedPoints.push_back(temp2);
		for (double i = 0; i < 1; i += 0.01) {
			calculateBSpline(i, controlPoints.size() - 1);
		}
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
