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
double uStep = 1000;
int order = 4;
bool drawing;
int drawType = GROUND;
double selectDistance = 0.05;
float rotateX, rotateZ = 0;

bool renderTerrain = false;


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

	glEnd();

	// Draw circle when drawing near end point of line
	// If cursor is in circle, last point of line will be the first point
	// to create a fully connected loop
	if (points.size() > 0) {
		if (abs(points.front().x - mouseX) < selectDistance && abs(points.front().y - mouseY) < selectDistance) {
			glBegin(GL_POLYGON);
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
			for (double j = 0; j < 2 * M_PI; j += M_PI / 25) {
				double pointX = sin(j) * 0.05 + points.front().x;
				double pointY = cos(j) * 0.05 + points.front().y;
				glVertex2d(pointX, pointY);
			}

			glEnd();
		}
	}


	glBegin(GL_LINES);
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

void renderGround(vector<vec2> points) {
}

void renderWater(vector<vec2> points){
	// Index of the middle point of the curve
	// double so if the number of points is even
	// mid index will be x.5
	double midIndex = (double)points.size() / 2.0;

	// The line the curve gets revolved around
	vector<vec3> axisOfRevolution;
	// Normals for each circle
	// They are normalized
	vector<vec3> normals;
	// A normalized vector from the centre of the circle to
	// any point on the circumference
	// Uses vector from axis of revolution to point on plotted line
	vector<vec3> circDirection;
	// The radius for each circle
	vector<double> radii;

	/*
	 * Loop through half the points
	 * If there  are an odd number of points, the middle one is never reached
	 */
	for (int i = 0; i < midIndex; i++) {
		vec3 axis;
		// Saves the midpoint between the ith points coming from both ends of the curve
		axis.x = (points[i].x + points[points.size() - 1 - i].x) / 2;
		axis.y = (points[i].y + points[points.size() - 1 - i].y) / 2;
		// z = y, might need to change this later
		axis.z = axis.y;
		axisOfRevolution.push_back(axis);
		vec3 temp;
		temp.x = points[i].x;
		temp.y = points[i].y;
		temp.z = temp.y;
		// Radius is the distance from the axis to the associated point on the curve 
		radii.push_back(length(axis - temp));
		// CircDirection is normalized radius
		circDirection.push_back(normalize(axis - temp));
	}

	// Normal calculations probably not right
	// Normal is just distance between first axis point to next
	for (int i = 0; i < axisOfRevolution.size() - 1; i++) {
		normals.push_back(normalize(axisOfRevolution[i + 1] - axisOfRevolution[i]));
	}
	// Last normal is just copy of second last as there is no next axis point to use
	normals.push_back(normals.back());

	vector < vector < vec3 > > waterPoints;
	for (int i = 0; i < axisOfRevolution.size(); i++) {
		vector<vec3> setOfpoints;
		// Starts at PI because only want bottom of lake, not top as well
		// Creates a circle from the ith point to the (numPoints - i)th point
		for (double j = M_PI; j < 2 * M_PI; j += M_PI / 25) {
			setOfpoints.push_back((float)(radii[i] * cos(j)) * circDirection[i] + (float)(radii[i] * sin(j))*cross(normals[i], circDirection[i]) + axisOfRevolution[i]);
		}
		waterPoints.push_back(setOfpoints);
	}

	for (int i = 0; i < waterPoints.size() - 1; i++) {
		glBegin(GL_QUAD_STRIP);
		glColor3f(0.0f, 0.0f, 1.0f);
		for (int j = 0; j < waterPoints[i].size(); j++) {
			glVertex3d(waterPoints[i][j].x, waterPoints[i][j].y, waterPoints[i][j].z);
			glVertex3d(waterPoints[i + 1][j].x, waterPoints[i + 1][j].y, waterPoints[i + 1][j].z);
		}
		glEnd();
	}
}

void renderTrees(vector<vec2> points) {

}

void render3D() {
	for (int i = 0; i < plottedPoints.size(); i++) {
		if (lineTypes[i] == GROUND) {
			renderGround(plottedPoints[i]);
		}
		else if (lineTypes[i] == WATER) {
			renderWater(plottedPoints[i]);
		}
		else {
			renderTrees(plottedPoints[i]);
		}
	}
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
	if (renderGround) {
		glRotatef(rotateX, 1.0f, 0.0f, 0.0f);
		glRotatef(rotateZ, 0.0f, 0.0f, 1.0f);
	}

	//Functions for changing projection matrix
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1, 1, -1, 1, -1, 1);
	//gluPerspective (fov, aspect ratio, near plane, far plane)
	//glFrustum


	if (!renderTerrain) {
		renderDrawing();
	}
	else {
		render3D();
	}
}

void finishLine() {
	if (points.size() > 0) {
		if (abs(points.front().x - mouseX) < selectDistance && abs(points.front().y - mouseY) < selectDistance) {
			if (points.back().x != points.front().x && points.back().y != points.front().y) {
				points.push_back(vec2(points.front().x, points.front().y));
			}
		}
	}

	controlPoints.push_back(points);

	vector <double> temp;
	knots.push_back(temp);
	recalculateKnots(controlPoints.size() - 1);
	vector <vec2> temp2;
	plottedPoints.push_back(temp2);
	for (double u = 0; u <= 1; u += 1.0 / uStep) {
		calculateBSpline(u, controlPoints.size() - 1);
	}
	points.clear();
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
	if (!renderTerrain) {
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
		else if (key == GLFW_KEY_ENTER && action == GLFW_PRESS) {
			renderTerrain = true;
			if (drawing) {
				finishLine();
			}
		}
	}
	else {
		if (key == GLFW_KEY_ENTER && action == GLFW_PRESS) {
			renderTerrain = false;
			rotateX = 0.f;
			rotateZ = 0.f;
		}
		else if (key == GLFW_KEY_UP && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			if (renderGround) {
				rotateX += 1.0f;
			}
		}
		else if (key == GLFW_KEY_DOWN && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			if (renderGround) {
				rotateX -= 1.0f;
			}
		}
		else if (key == GLFW_KEY_LEFT && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			if (renderGround) {
				rotateZ -= 1.0f;
			}
		}
		else if (key == GLFW_KEY_RIGHT && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			if (renderGround) {
				rotateZ += 1.0f;
			}
		}
	}
}

void mouseClick(GLFWwindow *sender, int button, int action, int mods) {
	if (!renderTerrain) {
		if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
			drawing = true;
			lineTypes.push_back(drawType);
		}
		else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
			drawing = false;

			finishLine();
		}
	}
}

void mousePos(GLFWwindow *sender, double x, double y) {	
	mouseX = (2 * x / w) - 1;
	mouseY = (-2 * y / h) + 1;
	
	if (!renderTerrain) {
		if (drawing) {
			if (points.size() > 0) {
				double distanceX = abs(mouseX - points.back().x);
				double distanceY = abs(mouseY - points.back().y);
				if (distanceX >= 0.01 || distanceY >= 0.01) {
					points.push_back(vec2(mouseX, mouseY));
				}
			}
			else {
				points.push_back(vec2(mouseX, mouseY));
			}
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
