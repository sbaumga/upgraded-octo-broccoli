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


/*
 * Written by Shayne Baumgartner, 10098339
 * A program to create a B-spline or NURBS curve by placing control points
 *
 * Controls:
 * Left click on an empty part on screen to place a control point
 * Left click and hold on a control point to move it
 * Left click and hold on a control point in NURBS mode and use the scroll wheel to change the weight
 *     The selection radius of the points is dependent on their weight, so be careful as making points weigh less will make them harder to select again
 * Right click to delete a point
 *
 * Left/Right arrow keys to increase/decrease order
 * Up/Down arrow keys to increase/decrease number of lines drawn in the curve
 * Space bar to toggle showing the deBoor's algorithm
 *     A/D to choose the point on the line for the deBoor's algorithm display
 *
 * N to switch between NURBs and B-Spline curves
 *     All set weights cleared when switched
 *
 * R to revolve a B-spline curve around the Y axis
 * W/S to rotate the surface about the X axis
 * Q/E to rotate the surface about the Z axis
 *     The line can be modified to dynamically modify the surface, but only when the rotation is 0
 * Backspace to reset the rotation
 */

GLFWwindow *window;
int w, h;
double mouseX, mouseY;

int selected = -1;
double selectDistance = 0.05;
vector<double> pointsX, pointsY, knots, plottedPointsX, plottedPointsY, nurbsWeights;
int order = 2;
bool needToRender = true;
double numPoints = 100;
bool showDeBoor = false;
int selectedPoint = 0;
bool nurbs = false;
bool revolveMode = false;
float rotateX, rotateZ = 0;
vector<vector<double>> revolveYX;
vector < vector<double>> revolveYZ;

/*
 * Finds delta for the given u using knots[]
 */
int findDelta(double u) {
	for (int i = 0; i < knots.size() - order; i++) {
		if (u >= knots[i] && u < knots[i+1]) {
			return i;
		}
	}
	return knots.size() - order;
}

/*
 * Recursively calculates the B-spline basis function for the given u, i (control point), and k (order)
 * Only used in NURBS calculcation
 */
double bSplineBasis(double u, int i, int k) {
	if (k == 1) {
		if (knots[i] <= u && u < knots[i + 1]) {
			return 1;
		}
		else {
			return 0;
		}
	}
	else {
		double mod1;
		double denom1 = knots[i + k - 1] - knots[i];
		if (denom1 == 0) {
			mod1 = 0;
		}
		else {
			mod1 = (u - knots[i]) / denom1;
		}
		double mod2;
		double denom2 = knots[i + k] - knots[i + 1];
		if (denom2 == 0) {
			mod2 = 0;
		}
		else {
			mod2 = (knots[i + k] - u) / denom2;
		}
		
		return mod1 * bSplineBasis(u, i, k - 1) + mod2 * bSplineBasis(u, i + 1, k - 1);
	}
}

/*
 * Calculates the x and y coords of a NURBs curve at point u and pushes them to vectors
 */
void calculateNurbs(double u) {
	double topSumX = 0;
	double topSumY = 0;
	double botSum = 0;
	for (int i = 0; i < pointsX.size(); i++) {
		double basis = bSplineBasis(u, i, order);

		botSum += basis * nurbsWeights[i];
		topSumX += basis * nurbsWeights[i] * pointsX[i];
		topSumY += basis * nurbsWeights[i] * pointsY[i];
	}

	plottedPointsX.push_back(topSumX / botSum);
	plottedPointsY.push_back(topSumY / botSum);
}

/*
* Calculates the x and y coords of a B-spline curve at point u and pushes them to vectors
* Uses deBoors algorithm
* Interpolated from pseudo-code shown in class
*/
void calculateBSpline(double u) {
	int delta = findDelta(u);
	vector<double> x, y;

	for (int i = 0; i <= order - 1; i++) {
		x.push_back(pointsX[delta - i]);
		y.push_back(pointsY[delta - i]);
	}

	for (int r = order; r >= 2; r--) {
		int i = delta;
		for (int s = 0; s <= r - 2; s++) {
			double denom = (knots[i + r - 1] - knots[i]);
			double omega;
			if (denom == 0) {
				omega = 0;
			}
			else {
				omega = (u - knots[i]) / denom;
			}
			x[s] = omega*x[s] + (1 - omega) * x[s + 1];
			y[s] = omega*y[s] + (1 - omega) * y[s + 1];
			i--;
		}
	}

	plottedPointsX.push_back(x[0]);
	plottedPointsY.push_back(y[0]);
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
	if (revolveMode) {
		glRotatef(rotateX, 1.0f, 0.0f, 0.0f);
		glRotatef(rotateZ, 0.0f, 0.0f, 1.0f);
	}

	//Functions for changing projection matrix
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1, 1, -1, 1, -1, 1);
	//gluPerspective (fov, aspect ratio, near plane, far plane)
	//glFrustum


	// Rendering the control points
	for (int i = 0; i < pointsX.size(); i++) {
		glBegin(GL_POLYGON); 
		glColor3f(1.0f, 1.0f, 1.0f);
		// Size of the points is modified by the NURBs weight of that point
		for (double j = 0; j < 2 * M_PI; j += M_PI / 25) {
			double pointX = sin(j) * (0.01 * nurbsWeights[i]) + pointsX[i];
			double pointY = cos(j) * (0.01 * nurbsWeights[i]) + pointsY[i];
			glVertex2f(pointX, pointY);
		}
		
		glEnd();
	}

	// Calculating the curve
	if (order <= pointsX.size()) {
		plottedPointsX.clear();
		plottedPointsY.clear();
		for (double u = 0; u < 1; u += (1 / numPoints)) {
			if (nurbs) {
				calculateNurbs(u);
			}
			else {
				calculateBSpline(u);
			}
		}
	}

	// Showing deBoor's algorithm
	if (showDeBoor && !revolveMode) {
		// Rendering the selected point as a circle
		glBegin(GL_POLYGON);
		glColor3f(0.0f, 1.0f, 0.0f);
		for (double j = 0; j < 2 * M_PI; j += M_PI / 25) {
			double pointX = sin(j) * 0.01 + plottedPointsX[selectedPoint];
			double pointY = cos(j) * 0.01 + plottedPointsY[selectedPoint];
			glVertex2f(pointX, pointY);
		}

		glEnd();


		glBegin(GL_LINES);
		glColor3f(0.0f, 0.0f, 1.0f);

		// Finding the u value of the selected point
		double u = (double)(selectedPoint / numPoints);
		int delta = findDelta(u);
		vector<double> x, y;


		// Uses same method a calculateBSpline, but renders some intermediary steps
		for (int i = 0; i <= order - 1; i++) {
			x.push_back(pointsX[delta - i]);
			y.push_back(pointsY[delta - i]);
		}

		for (int r = order; r >= 2; r--) {
			int i = delta;
			for (int s = 0; s <= r - 2; s++) {
				glVertex2f(x[s], y[s]);
				glVertex2f(x[s + 1], y[s + 1]);
			}
			for (int s = 0; s <= r - 2; s++) {
				double denom = (knots[i + r - 1] - knots[i]);
				double omega;
				if (denom == 0) {
					omega = 0;
				}
				else {
					omega = (u - knots[i]) / denom;
				}
				x[s] = omega*x[s] + (1 - omega) * x[s + 1];
				y[s] = omega*y[s] + (1 - omega) * y[s + 1];
				i--;
			}
		}

		glVertex2f(x[0], y[0]);
		glVertex2f(x[1], y[1]);

		glEnd();
	}

	// Rendering the curve
	glBegin(GL_LINES);
	glColor3f(1.0f, 0.0f, 0.0f);

	for (int i = 1; i < plottedPointsX.size(); i++) {
		glVertex2f(plottedPointsX[i - 1], plottedPointsY[i - 1]);
		glVertex2f(plottedPointsX[i], plottedPointsY[i]);
	}

	glEnd();

	// Revolving the curve
	if (revolveMode && !nurbs) {
		revolveYX.clear();
		revolveYZ.clear();

		// Calculating the X and Z points by creating a circle with a radius from 
		// the y-axis to the x point of a point of the curve
		for (int i = 0; i < plottedPointsY.size(); i++) {
			vector<double> shapePointsX, shapePointsZ;
			for (double j = 0; j < 2 * M_PI; j += M_PI / 25) {
				double pointX = sin(j) * abs(plottedPointsX[i]);
				shapePointsX.push_back(pointX);
				double pointZ = cos(j) * abs(plottedPointsX[i]);
				shapePointsZ.push_back(pointZ);
			}
			revolveYX.push_back(shapePointsX);
			revolveYZ.push_back(shapePointsZ);
		}

		// Rendering the surface using the calculated X and Z points
		// with the Y value being taken from the corresponding point on the curve
		for (int i = 0; i < revolveYX.size() - 2; i++) {
			glBegin(GL_QUAD_STRIP);
			glColor3f(0.0f, 1.0f, 0.0f);
			for (int j = 0; j < revolveYX[i].size(); j++) {
				glVertex3f(revolveYX[i][j], plottedPointsY[i], revolveYZ[i][j]);
				glVertex3f(revolveYX[i + 1][j], plottedPointsY[i + 1], revolveYZ[i + 1][j]);
			}
			// Ensures that a solid ring is created
			glVertex3f(revolveYX[i][0], plottedPointsY[i], revolveYZ[i][0]);
			glVertex3f(revolveYX[i + 1][0], plottedPointsY[i + 1], revolveYZ[i + 1][0]);
			glEnd();
		}
	}
}

/*
 * Recalculates a standard knot sequence whenever a control point is added or removed
 */
void recalculateKnots() {
	knots.clear();

	if (order <= pointsX.size()) {
		int numKnots = pointsX.size() - (order - 2);
		if (numKnots != 0) {
			// Setting up standard knot sequence
			for (int i = 0; i < order - 1; i++) {
				knots.push_back(0);

			}

			for (int i = 0; i < numKnots; i++) {
				knots.push_back((double)i / (numKnots - 1));
			}

			// Setting up standard knot sequence
			for (int i = 0; i < order - 1; i++) {
				knots.push_back(1);
			}
		}
	}
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
	if (key == GLFW_KEY_LEFT && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		if (order > 2) {
			order--;
			recalculateKnots();
			needToRender = true;
		}
		cout << "Order is now " << order << "\n";
	}
	else if (key == GLFW_KEY_RIGHT && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		if (order + 1 <= pointsX.size()) {
			order++;
			recalculateKnots();
			needToRender = true;
		}
		cout << "Order is now " << order << "\n";
	}
	else if (key == GLFW_KEY_DOWN && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		if (numPoints > 2) {
			numPoints--;
			needToRender = true;
		}
		cout << "Number of lines drawn is now " << numPoints << "\n";
	}
	else if (key == GLFW_KEY_UP && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		numPoints++;
		cout << "Number of lines drawn is now " << numPoints << "\n";
		needToRender = true;
	}
	else if (key == GLFW_KEY_SPACE && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		if (nurbs) {
			showDeBoor = false;
			cout << "DeBoor's algorithm cannot be used for NURBs curves\n";
		}
		else {
			if (showDeBoor) {
				showDeBoor = false;
				cout << "Disabling deBoor geometry display\n";
			}
			else {
				showDeBoor = true;
				cout << "Enabling deBoor geometry display\n";
			}
			needToRender = true;
		}
	}
	else if (key == GLFW_KEY_A && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		if (showDeBoor) {
			if (selectedPoint > 0) {
				selectedPoint -= 1;
				needToRender = true;
			}
			cout << "Selected point is " << selectedPoint << "\n";
		}
	}
	else if (key == GLFW_KEY_D && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		if (showDeBoor && !revolveMode && !nurbs) {
			if (selectedPoint < numPoints - 1) {
				selectedPoint += 1;
				needToRender = true;
			}
			cout << "Selected point is " << selectedPoint << "\n";
		}
	}
	else if (key == GLFW_KEY_N && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		if (revolveMode) {
			cout << "NURBs not allowed with revolving\n";
			nurbs = false;
		}
		else {
			needToRender = true;
			if (nurbs) {
				nurbs = false;
				for (int i = 0; i < nurbsWeights.size(); i++) {
					nurbsWeights[i] = 1;
				}
				cout << "Displaying B-spline curve\n";
			}
			else {
				nurbs = true;
				cout << "Displaying NURBS curve\n";
				if (showDeBoor) {
					showDeBoor = false;
					cout << "DeBoor's algorithm cannot be used for NURBs curves\n";
				}
			}
		}
	}
	else if (key == GLFW_KEY_R && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		if (nurbs) {
			cout << "Revolving not allowed on NURBs curves\n";
		}
		else {
			if (revolveMode) {
				cout << "Revolving canceled\n";
				revolveMode = false;
				needToRender = true;
				rotateX = 0.0f;
				rotateZ = 0.0f;
			}
			else {
				cout << "Revolving shown\n";
				revolveMode = true;
				needToRender = true;
			}
		}
	} 
	else if (key == GLFW_KEY_W && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		if (revolveMode) {
			rotateX += 1.0f;
			needToRender = true;
			cout << "X rotation: " << rotateX << "\n";
		}
		else {
			cout << "Rotating only allowed in revolve mode\n";
		}
	}
	else if (key == GLFW_KEY_S && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		if (revolveMode) {
			rotateX -= 1.0f;
			needToRender = true;
			cout << "X rotation: " << rotateX << "\n";
		}
		else {
			cout << "Rotating only allowed in revolve mode\n";
		}
	}
	else if (key == GLFW_KEY_Q && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		if (revolveMode) {
			rotateZ -= 1.0f;
			needToRender = true;
			cout << "Z rotation: " << rotateZ << "\n";
		}
		else {
			cout << "Rotating only allowed in revolve mode\n";
		}
	}
	else if (key == GLFW_KEY_E && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		if (revolveMode) {
			rotateZ += 1.0f;
			needToRender = true;
			cout << "Z rotation: " << rotateZ << "\n";
		}
		else {
			cout << "Rotating only allowed in revolve mode\n";
		}
	}
	else if (key == GLFW_KEY_BACKSPACE && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		if (revolveMode) {
			rotateX = 0;
			rotateZ = 0;
			needToRender = true;
			cout << "Reset rotation\n";
		}
	}
}

void mouseClick(GLFWwindow *sender, int button, int action, int mods) {
	if (rotateX == 0 && rotateZ == 0) {
		if (action == GLFW_PRESS) {
			for (int i = 0; i < pointsX.size(); i++) {
				if (abs(pointsX[i] - mouseX) < selectDistance * nurbsWeights[i] && abs(pointsY[i] - mouseY) < selectDistance * nurbsWeights[i]) {
					selected = i;
				}
			}

			if (button == GLFW_MOUSE_BUTTON_LEFT && selected == -1) {
				pointsX.push_back(mouseX);
				pointsY.push_back(mouseY);
				nurbsWeights.push_back(1.0);
				recalculateKnots();
			}
			if (button == GLFW_MOUSE_BUTTON_RIGHT && selected > -1) {
				pointsX.erase(pointsX.begin() + selected);
				pointsY.erase(pointsY.begin() + selected);
				nurbsWeights.erase(nurbsWeights.begin() + selected);
				recalculateKnots();
				selected = -1;
			}

		}
		if (action == GLFW_RELEASE) {
			selected = -1;
		}

		needToRender = true;
	}
	else {
		cout << "Editing curve can only bo done when no rotation is applied to prevent" <<
			" unexpected behaviour.\n X rotation: " << rotateX << " Z rotation: " << rotateZ << "\n";
	}
}

void mousePos(GLFWwindow *sender, double x, double y) {
	mouseX = (2 * x / w) - 1;
	mouseY = (-2 * y / h) + 1;

	if (selected > -1) {
		pointsX[selected] = mouseX;
		pointsY[selected] = mouseY;
		needToRender = true;
	}
}

void scroll(GLFWwindow *sender, double x, double y) {
	if (nurbs && selected > -1) {
		if (nurbsWeights[selected] + 0.1 * y > 0) {
			nurbsWeights[selected] = nurbsWeights[selected] + 0.1 * y;
			needToRender = true;
		}
		else {
			cout << "Weights cannot be reduced below 0\n";
		}
	}
}

int main() {
	window = initializeWindow();
	if (window == NULL) {
		return -1;
	}

	glfwSetKeyCallback(window, keyboard);
	glfwSetMouseButtonCallback(window, mouseClick);
	glfwSetCursorPosCallback(window, mousePos);
	glfwSetScrollCallback(window, scroll);
	while (!glfwWindowShouldClose(window)) {
		if (needToRender) {
			glfwGetFramebufferSize(window, &w, &h);
			glViewport(0, 0, w, h);

			render();

			glfwSwapBuffers(window);
			needToRender = false;
		}
		glfwPollEvents();
	}

	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}
