#include <iostream>
#include <GL/glew.h>  
#include <GLFW/glfw3.h>  
#include <glm/glm.hpp>
#include <stdio.h>  
#include <stdlib.h>  
#include <vector>
#include <cmath>
#include <functional>

#include "Renderer.h"
#include "camera.h"
#include "mapLoader.h"

#define _USE_MATH_DEFINES
#include <math.h>
using namespace std;
using namespace glm;

const int PLANE_DIM_X = 400;
const int PLANE_DIM_Z = 400;

enum
{
	GROUND,
	WATER,
	TREES,
	MOUNTAIN,
};

enum
{
	NEW_START_AT_START,
	NEW_START_AT_END,
	NEW_END_AT_START,
	NEW_END_AT_END,
};

const unsigned int NUMBER_OF_SLICES = 50;


unsigned int OCTAVES = 5;
float PERSISTENCE = 0.4f;
unsigned int RANGE = 40;

unsigned int LAND_OCTAVES = 7;
float LAND_PERSISTENCE = 0.6;
unsigned int LAND_RANGE = 8;

PerlinNoise mountainNoise(OCTAVES, PERSISTENCE, RANGE);
PerlinNoise landNoise(LAND_OCTAVES, LAND_PERSISTENCE, LAND_RANGE);

GLFWwindow *window;
int w, h;
double mouseX, mouseY;

vector < int > lineTypes;
vector < int > lineTypes2;
vector < vector < vec2 > > controlPoints;
vector < vec2 > points;
vector < vector < double > > knots;
vector < vector < vec2 > > plottedPoints;

vector<vec2> mountainCenters;
vector<float> mountainRadii;

double uStep = 1000;
int order = 4;
bool drawing;
int drawType = GROUND;
double selectDistance = 0.06;

bool erasing = false;
float erasingRadius = 0.05f;

float rotateX, rotateZ = 0;

bool renderTerrain = false;
bool debugging = false;

//3D manipulation
Camera lookat;
vec2 lastPos;		//Last mouse position
bool rotatingCamera = false;

//2D mouse manipulation
bool panningCamera = false;

//Terrain
vector<vec3> terrainPoints;
vector<vec3> terrainNormals;
vector<vec2> terrainUVs;
vector<unsigned int> terrainIndices;

float canvasWidth = 10.f;
float canvasHeight = 10.f;

//Rendering
Renderer r;
bool wireframe = false;

float xTrans, yTrans, zTrans = 0.0f;
float scale = std::max(1.f/(canvasWidth*1.1f), 1.f/(canvasHeight*1.1f));

// Define Infinite (Using INT_MAX caused overflow problems)
#define INF 10000

// Given three colinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool onSegment(vec2 p, vec2 q, vec2 r)
{
	if (q.x <= std::max(p.x, r.x) && q.x >= std::min(p.x, r.x) &&
		q.y <= std::max(p.y, r.y) && q.y >= std::min(p.y, r.y))
		return true;
	return false;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(vec2 p, vec2 q, vec2 r)
{
	int val = (q.y - p.y) * (r.x - q.x) -
		(q.x - p.x) * (r.y - q.y);

	if (val == 0) return 0;  // colinear
	return (val > 0) ? 1 : 2; // clock or counterclock wise
}


// The function that returns true if line segment 'p1q1'
// and 'p2q2' intersect
bool doIntersect(vec2 p1, vec2 q1, vec2 p2, vec2 q2) {
	// Find the four orientations needed for general and
	// special cases
	int o1 = orientation(p1, q1, p2);
	int o2 = orientation(p1, q1, q2);
	int o3 = orientation(p2, q2, p1);
	int o4 = orientation(p2, q2, q1);

	// General case
	if (o1 != o2 && o3 != o4)
		return true;

	// Special Cases
	// p1, q1 and p2 are colinear and p2 lies on segment p1q1
	if (o1 == 0 && onSegment(p1, p2, q1)) return true;

	// p1, q1 and p2 are colinear and q2 lies on segment p1q1
	if (o2 == 0 && onSegment(p1, q2, q1)) return true;

	// p2, q2 and p1 are colinear and p1 lies on segment p2q2
	if (o3 == 0 && onSegment(p2, p1, q2)) return true;

	// p2, q2 and q1 are colinear and q1 lies on segment p2q2
	if (o4 == 0 && onSegment(p2, q1, q2)) return true;

	return false; // Doesn't fall in any of the above cases
}

bool isInside(vector<vec2> line, vec2 point) {
	if (line.size() < 3) {
		return false;
	}

	vec2 extreme = vec2(INF, point.y);
	
	// Count intersections of the extreme line with sides of line
	int count = 0;
	int i = 0;
	do {
		int next = (i + 1) % line.size();

		// Check if the line segment from 'point' to 'extreme' intersects
		// with the line segment from 'line[i]' to 'line[next]'
		if (doIntersect(line[i], line[next], point, extreme)){
			// If the point is colinear with line segment 'i-next',
			// then check if it lies on segment. If it lies, return true
			// otherwise false
			if (orientation(line[i], point, line[next]) == 0) {
				return onSegment(line[i], point, line[next]);
			}

			count++;
		}
		i = next;
	} while (i != 0);

	// Return true if count is odd, false otherwise
	return count & 1;
}

void getBoundingCirclesFromMountainLines(vector<vec2>* centers, vector<float>* radii)
{
	centers->clear();
	radii->clear();

	for (unsigned int i = 0; i < controlPoints.size(); i++)
	{
		if (lineTypes[i] == MOUNTAIN)
		{
			vec2 leftMost(1000.f, 0.f);
			vec2 rightMost(-1000.f, 0.f);
			vec2 topMost(0.f, -1000.f);

			for (unsigned int j = 0; j < controlPoints[i].size(); j++)
			{
				//Find top point, bottom left point, and bottom right point
				if (controlPoints[i][j].x < leftMost.x)
					leftMost = controlPoints[i][j];
				if (controlPoints[i][j].x > rightMost.x)
					rightMost = controlPoints[i][j];
				if (controlPoints[i][j].y > topMost.y)
					topMost = controlPoints[i][j];
			}

			vec2 center = (leftMost + rightMost + topMost) / 3.f;
			float radius = length(leftMost - center);

			centers->push_back(center);
			radii->push_back(radius);

		}
	}
}

/*
 * Returns a number which represents the how many bodies of water the square is surrounded by
 * 0 = not surrounded at all
 * 1 = one corner of the square is surrounded by one body of water
 * 2 = two corners are surrounded by one body of water
 * OR two corners are each surrounded by their own body of water
 * OR one corner is surrounded by two bodies of water
 * 3 and up are similar
 *
 * All involved methods are adapted from http://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/
 */
int isSurroundedByWater(vec2 topLeft, float width, float height){
	int waterPointCounter = 0;
	float xMin = topLeft.x;
	float xMax = topLeft.x + width;
	float yMin = topLeft.y;
	float yMax = topLeft.y + height;
	for (int i = 0; i < plottedPoints.size(); i++) {
		// If line is water and is a closed curve
		if (lineTypes[i] == WATER && plottedPoints[i][0] == plottedPoints[i][plottedPoints[i].size() - 1]) {
			if (isInside(plottedPoints[i], vec2(xMin, yMin))) {
				waterPointCounter++;
			}
			if (isInside(plottedPoints[i], vec2(xMax, yMin))) {
				waterPointCounter++;
			}
			if (isInside(plottedPoints[i], vec2(xMin, yMax))) {
				waterPointCounter++;
			}
			if (isInside(plottedPoints[i], vec2(xMax, yMax))) {
				waterPointCounter++;
			}
		}
	}
	return waterPointCounter;
}

/*
 * Returns the number of water points inside the given square
 */
int containsWaterBorders(vec2 topLeft, float width, float height) {
	int waterPointCounter = 0;
	float xMin = topLeft.x;
	float xMax = topLeft.x + width;
	float yMin = topLeft.y;
	float yMax = topLeft.y + height;
	for (int i = 0; i < plottedPoints.size(); i++) {
		if (lineTypes[i] == WATER) {
			for (int j = 0; j < plottedPoints[i].size(); j++) {
				float x = plottedPoints[i][j].x;
				float y = plottedPoints[i][j].y;
				if ((x >= xMin && x <= xMax) && (y >= yMin && y <= yMax)) {
					waterPointCounter++;
				}
			}
		}
	}
	return waterPointCounter;
}

/*
* Returns the number of ground points inside the given square
*/
int containsGroundBorders(vec2 topLeft, float width, float height) {
	int groundPointCounter = 0;
	float xMin = topLeft.x;
	float xMax = topLeft.x + width;
	float yMin = topLeft.y;
	float yMax = topLeft.y + height;
	for (int i = 0; i < plottedPoints.size(); i++) {
		if (lineTypes[i] == GROUND) {
			for (int j = 0; j < plottedPoints[i].size(); j++) {
				float x = plottedPoints[i][j].x;
				float y = plottedPoints[i][j].y;
				if ((x >= xMin && x <= xMax) && (y >= yMin && y <= yMax)) {
					groundPointCounter++;
				}
			}
		}
	}
	return groundPointCounter;
}



/*
 * Returns a vector containing the two points that make up each line that crosses either border of the horizontal slice
 * defined by the two y values
 * Only returns Edges for lines of the given line type
 *
 * topY should be the higher y value (i.e. closest to the top of the screen)
 * and bottomY should be the lower y value
 */
vector<Edge> findEdges(int lineType, float topY, float bottomY) {
	vector<Edge> crossingEdges;
	for (int lineNum = 0; lineNum < controlPoints.size(); lineNum++) {
		if (lineTypes[lineNum] == lineType) {
			bool inside;
			if (controlPoints[lineNum].front().y <= topY && controlPoints[lineNum].front().y > bottomY) {
				inside = true;
			}
			else {
				inside = false;
			}

			// Starts at 1 because we already know where the first point lie from above
			// Also because an edge needs two points to be defined
			for (int i = 1; i < controlPoints[lineNum].size(); i++) {
				if (inside) {
					// If inside, check if point lies outside
					if (controlPoints[lineNum][i].y > topY || controlPoints[lineNum][i].y <= bottomY) {
						inside = false;
						Edge temp;
						temp.p1 = &controlPoints[lineNum][i - 1];
						temp.p2 = &controlPoints[lineNum][i];

						crossingEdges.push_back(temp);
					}
				}
				else {
					// If outside, check if point lies inside
					if (controlPoints[lineNum][i].y <= topY && controlPoints[lineNum][i].y > bottomY) {
						inside = true;
						Edge temp;
						temp.p1 = &controlPoints[lineNum][i - 1];
						temp.p2 = &controlPoints[lineNum][i];

						crossingEdges.push_back(temp);
					}
				}
			}

			//Additional case to add all internal, non-crossing edges
			for (int i = 1; i < controlPoints[lineNum].size(); i++)
			{
				if ((controlPoints[lineNum][i - 1].y < topY) &&
					(controlPoints[lineNum][i - 1].y > bottomY) &&
					(controlPoints[lineNum][i].y < topY) &&
					(controlPoints[lineNum][i].y > bottomY))
				{
					Edge temp;
					temp.p1 = &controlPoints[lineNum][i - 1];
					temp.p2 = &controlPoints[lineNum][i];
					crossingEdges.push_back(temp);
				}
			}
		}
	}

	return crossingEdges;
}

void partitionEdgesIntoSlices(int lineType, float minY, float maxY, unsigned int numSlices, vector<vector<Edge>>* partition)
{
	partition->clear();

	float increment = (maxY - minY) / (float)(numSlices - 1);
	float low = minY;

	for (unsigned int i = 0; i < numSlices; i++)
	{
		partition->push_back(findEdges(lineType, low + increment, low));
		low += increment;
	}
}

/*
* Finds delta for the given u using knots[lineNum]
*/
int findDelta(double u, int lineNum) {
	for (int i = 0; i < (int)knots[lineNum].size() - order; i++) {
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
	glClearColor(0.5f, 0.5f, 0.5f, 1.f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glBegin(GL_QUADS);
		//Shadow
		glColor3f(0.f, 0.f, 0.f);
		glVertex3f(-canvasWidth*0.55, -canvasHeight*0.55, -1.f);
		glVertex3f(canvasWidth*0.45, -canvasHeight*0.55, -1.f);
		glVertex3f(canvasWidth*0.45, canvasHeight*0.45, -1.f);
		glVertex3f(-canvasWidth*0.55, canvasHeight*0.45, -1.f);

		//Canvas
		glColor3f(1.f, 1.f, 1.f);
		glVertex3f(-canvasWidth*0.5, -canvasHeight*0.5, -0.99f);
		glVertex3f(canvasWidth*0.5, -canvasHeight*0.5, -0.99f);
		glVertex3f(canvasWidth*0.5, canvasHeight*0.5, -0.99f);
		glVertex3f(-canvasWidth*0.5, canvasHeight*0.5, -0.99f);
	glEnd();
	
	
	if (debugging)
	{
		// Debugging edge detection
		float sliceStart = 0.f;
		float sliceSize = 0.1f;

		vector<Edge> edgesInSlice = findEdges(GROUND, sliceStart + sliceSize, sliceStart);
		vector<vector<Edge>> partition;

		partitionEdgesIntoSlices(GROUND, canvasHeight*-0.5, canvasHeight*0.5, NUMBER_OF_SLICES, &partition);

		printf("Segments in slice = %d\n", partition.size());

		glBegin(GL_LINES);

		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex2f(-10.0f, sliceStart);
		glVertex2f(10.0f, sliceStart);

		glVertex2f(-10.0f, sliceStart + sliceSize);
		glVertex2f(10.0f, sliceStart + sliceSize);


		/*	glLineWidth(5.f);
			glColor3f(1.f, 0.f, 0.f);
			for (unsigned int i = 0; i < partition.size(); i++)
			{
			if (i % 2 == 0)
			{
			glColor3f(
			1.f - (float)i / (float)(partition.size() - 1),
			0.f,
			(float)i / (float)(partition.size() - 1));

			for (unsigned int j = 0; j < partition[i].size(); j++)
			{
			glVertex2f(partition[i][j].p1->x + 0.1f, partition[i][j].p1->y);
			glVertex2f(partition[i][j].p2->x + 0.1f, partition[i][j].p2->y);
			}
			}

			}

			glEnd();

			glLineWidth(1.f);*/

		float increment = canvasHeight / 29.f;

		glBegin(GL_POINTS);
		for (unsigned int i = 0; i < 30; i++)
		{
			for (unsigned int j = 0; j < 30; j++)
			{
				vec2 point(
					increment*(float)i - canvasWidth*0.5f,
					increment*(float)j - canvasHeight*0.5f);

				if (inPartitionedPolygon(point, &partition, -canvasHeight*0.5f, canvasHeight*0.5f))
					glVertex2f(point.x, point.y);
			}
		}
		glEnd();

		//Debugging
		getBoundingCirclesFromMountainLines(&mountainCenters, &mountainRadii);

		glPointSize(2.f);
		glBegin(GL_POINTS);
		glColor3f(1.f, 0.f, 0.f);
		for (unsigned int i = 0; i < mountainCenters.size(); i++)
		{
			glVertex2f(mountainCenters[i].x, mountainCenters[i].y);
		}
		glEnd();

		glBegin(GL_LINES);
		glColor3f(0.f, 1.f, 1.f);
		for (unsigned int i = 0; i < mountainCenters.size(); i++)
		{
			glVertex2f(mountainCenters[i].x, mountainCenters[i].y);
			glVertex2f(mountainCenters[i].x, mountainCenters[i].y + mountainRadii[i]);
		}
		glEnd();

		glPointSize(2.f);
	}
	
	

	glBegin(GL_LINES);
	
	// Draw finished lines
	if (plottedPoints.size() >= 1) {
		for (int i = 0; i < plottedPoints.size(); i++) {
			vector < vec2 > line = plottedPoints[i];
			int type = lineTypes2[i];
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
			else if (type == MOUNTAIN) {
				glColor3f(0.5f, 0.5f, 0.5f);
			}

			for (int j = 0; j < (int)(line.size() - 1); j++) {
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
			else if (drawType == MOUNTAIN) {
				glColor3f(0.5f, 0.5f, 0.5f);
			}
			for (double j = 0; j < 2 * M_PI; j += M_PI / 25) {
				double pointX = sin(j) * 0.05 + points.front().x;
				double pointY = cos(j) * 0.05 + points.front().y;
				glVertex2d(pointX, pointY);
			}

			glEnd();
		}
	}

	// Draw a circle if near an end point for a line of the same draw type
	for (int i = 0; i < controlPoints.size(); i++) {
		if (drawType == lineTypes[i]) {
			if (abs(controlPoints[i].front().x - mouseX) < selectDistance && abs(controlPoints[i].front().y - mouseY) < selectDistance) {
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
				else if (drawType == MOUNTAIN) {
					glColor3f(0.5f, 0.5f, 0.5f);
				}
				for (double j = 0; j < 2 * M_PI; j += M_PI / 25) {
					double pointX = sin(j) * 0.05 + controlPoints[i].front().x;
					double pointY = cos(j) * 0.05 + controlPoints[i].front().y;
					glVertex2d(pointX, pointY);
				}

				glEnd();
			}
			else if (abs(controlPoints[i].back().x - mouseX) < selectDistance && abs(controlPoints[i].back().y - mouseY) < selectDistance) {
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
				else if (drawType == MOUNTAIN) {
					glColor3f(0.5f, 0.5f, 0.5f);
				}
				for (double j = 0; j < 2 * M_PI; j += M_PI / 25) {
					double pointX = sin(j) * 0.05 + controlPoints[i].back().x;
					double pointY = cos(j) * 0.05 + controlPoints[i].back().y;
					glVertex2d(pointX, pointY);
				}

				glEnd();
			}
		}
	}

	// Draw current line
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
		else if (drawType == MOUNTAIN) {
			glColor3f(0.5f, 0.5f, 0.5f);
		}
		for (int i = 0; i < points.size() - 1; i++) {
			glVertex2d(points[i].x, points[i].y);
			glVertex2d(points[i + 1].x, points[i + 1].y);
		}
	}

	glEnd();


	if (erasing) {
		// Draw erasing circle
		glBegin(GL_LINE_STRIP);

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
		else if (drawType == MOUNTAIN) {
			glColor3f(0.5f, 0.5f, 0.5f);
		}

		for (double j = 0; j < 2 * M_PI; j += M_PI / 25) {
			double pointX = sin(j) * (erasingRadius / scale) + mouseX;
			double pointY = cos(j) * (erasingRadius / scale) + mouseY;
			glVertex2d(pointX, pointY);
		}

		double pointX = sin(0) * (erasingRadius / scale) + mouseX;
		double pointY = cos(0) * (erasingRadius / scale) + mouseY;
		glVertex2d(pointX, pointY);

		glEnd();
	}

	

	//Range
	float stepSize = canvasHeight / (float)(RANGE - 1);
	float offset = -canvasHeight*0.5;
	glColor3f(0.8f, 0.8f, 0.8f);
	glBegin(GL_LINES);
	//Rows
	for (unsigned int i = 0; i < RANGE; i++)
	{
		glVertex3f(offset, offset + stepSize*(float)i, -0.9f);
		glVertex3f(-offset, offset + stepSize*(float)i, -0.9f);
	}
	//Columns
	for (unsigned int i = 0; i < RANGE; i++)
	{
		glVertex3f(offset + stepSize*(float)i, offset, -0.9f);
		glVertex3f(offset + stepSize*(float)i, -offset, -0.9f);
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
	/*for (int i = 0; i < plottedPoints.size(); i++) {
		if (lineTypes[i] == GROUND) {
			renderGround(plottedPoints[i]);
		}
		else if (lineTypes[i] == WATER) {
			renderWater(plottedPoints[i]);
		}
		else {
			renderTrees(plottedPoints[i]);
		}
	}*/
	glClearColor(1.f, 1.f, 1.f, 1.f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	r.loadModelBuffer(&terrainPoints, &terrainNormals, &terrainIndices);
	r.updateTransform();
	r.render3DView();

	/*vector<vec3> verts, normals;
	vector<vec2> uvs;
	vector<unsigned int> indices;
	generatePlane(2, 2, canvasWidth, canvasHeight, &verts, &normals, &uvs, &indices);
	r.loadModelBuffer(&verts, &normals, &indices);
	r.render3DView();*/
}

void erase() {
	// These lines will be added to the plottedPoints vector after all have been considered
	vector < vector <vec2> > linesAfterEreasing;
	// These are the indices of of the lines that the erasing effects
	// The lines corresponding to these indices will be erased after all have been considered
	vector<int> lineIndices;
	for (int lineNum = 0; lineNum < controlPoints.size(); lineNum++) {
		
		// Only erases lines of the correct type
		if (lineTypes[lineNum] == drawType) {
			bool inLineIndices = false;
			bool startPointInside = false;
			// Tracks if the intersection will be moving out of or into the circle
			bool inside = false;
			
			// These are the indices of the points just before the line intersects with the circle
			vector<int> pointsBeforeIntersections;

			for (int j = 0; j < controlPoints[lineNum].size(); j++) {
				// Uses parametric equation of circle (centreX^2 + centreY^2 = r^2)
				// if x^2 + y^2 <= 1 the point is in or on the circle
				float x = controlPoints[lineNum][j].x;
				float y = controlPoints[lineNum][j].y;
				if ((x - mouseX) * (x - mouseX) + (y - mouseY) * (y - mouseY) <= (erasingRadius / scale) * (erasingRadius / scale)) {
					// If last point was outside, probably intersection
					if (!inside) {
						if (!inLineIndices) {
							// Push back line num, as this line will be affected by the erasing
							inLineIndices = true;
							lineIndices.push_back(lineNum);
						}

						// Now inside the circle
						inside = true;
						if (j == 0) {
							// If this is the first point, then no intersection has taken place
							// But this will mean that the points going from odd to even intersection numbers
							// will be kept, instead of erased like normal
							startPointInside = true;
						}
						else {
							// Pushing point j - 1 as we want to keep track of the points outside of the circle
							// and j is outside, but j - 1 was not
							pointsBeforeIntersections.push_back(j - 1);
						}
					}
				}
				// Otherwise point outside of circle
				else {
					// If last point was inside, intersection
					if (inside) {
						// Now outside
						inside = false;

						// Pushing back j because we want points outside of the circle and this time j is outside
						pointsBeforeIntersections.push_back(j);
					}
				}
			}
			// Whole line has been checked now

			// Might need to add case where no intersections, but whole line is in circle
			if (pointsBeforeIntersections.size() > 0) {
				// Determines wether to starting storing points from the first point of plottedPoints
				// or to start from the first intersection point
				int startingIntersection;
				if (!startPointInside) {
					// Normal behaivour
					// Start from firstPoint in plottedPoints
					// Store points from even intersections to odd
					startingIntersection = -1;
				}
				else {
					// Start from first intersection
					// Store points from odd intersections to even
					startingIntersection = 0;
				}



				for (int i = startingIntersection; i < ((int)pointsBeforeIntersections.size()); i = i + 2) {
					vector<vec2> line;

					// Indicates the starting index of plottedPoints
					int startPoint;
					if (i == -1) {
						startPoint = 0;
					}
					else {
						startPoint = pointsBeforeIntersections[i];
					}

					// Indicates the ending index of plottedPoints
					int endPoint;
					if (i == pointsBeforeIntersections.size() - 1) {
						// If the line does not end in the erasing zone
						endPoint = controlPoints[lineNum].size() - 1;
					}
					else {
						endPoint = pointsBeforeIntersections[i + 1];
					}

					// Loops through and saves the correct points
					for (int j = startPoint; j <= endPoint; j++) {
						line.push_back(controlPoints[lineNum][j]);
					}

					// All points for these intersections now added to line
					linesAfterEreasing.push_back(line);
				}
			}
		}
		// End of draw type check
	}

	vector < vector <vec2> > linesAfterEreasing2;
	// These are the indices of of the lines that the erasing effects
	// The lines corresponding to these indices will be erased after all have been considered
	vector<int> lineIndices2;
	for (int lineNum = 0; lineNum < plottedPoints.size(); lineNum++) {

		// Only erases lines of the correct type
		if (lineTypes2[lineNum] == drawType) {
			bool inLineIndices = false;
			bool startPointInside = false;
			// Tracks if the intersection will be moving out of or into the circle
			bool inside = false;

			// These are the indices of the points just before the line intersects with the circle
			vector<int> pointsBeforeIntersections;

			for (int j = 0; j < plottedPoints[lineNum].size(); j++) {
				// Uses parametric equation of circle (centreX^2 + centreY^2 = r^2)
				// if x^2 + y^2 <= 1 the point is in or on the circle
				float x = plottedPoints[lineNum][j].x;
				float y = plottedPoints[lineNum][j].y;
				if ((x - mouseX) * (x - mouseX) + (y - mouseY) * (y - mouseY) <= (erasingRadius / scale) * (erasingRadius / scale)) {
					// If last point was outside, probably intersection
					if (!inside) {
						if (!inLineIndices) {
							// Push back line num, as this line will be affected by the erasing
							inLineIndices = true;
							lineIndices2.push_back(lineNum);
						}

						// Now inside the circle
						inside = true;
						if (j == 0) {
							// If this is the first point, then no intersection has taken place
							// But this will mean that the points going from odd to even intersection numbers
							// will be kept, instead of erased like normal
							startPointInside = true;
						}
						else {
							// Pushing point j - 1 as we want to keep track of the points outside of the circle
							// and j is outside, but j - 1 was not
							pointsBeforeIntersections.push_back(j - 1);
						}
					}
				}
				// Otherwise point outside of circle
				else {
					// If last point was inside, intersection
					if (inside) {
						// Now outside
						inside = false;

						// Pushing back j because we want points outside of the circle and this time j is outside
						pointsBeforeIntersections.push_back(j);
					}
				}
			}
			// Whole line has been checked now

			// Might need to add case where no intersections, but whole line is in circle
			if (pointsBeforeIntersections.size() > 0) {
				// Determines wether to starting storing points from the first point of plottedPoints
				// or to start from the first intersection point
				int startingIntersection;
				if (!startPointInside) {
					// Normal behaivour
					// Start from firstPoint in plottedPoints
					// Store points from even intersections to odd
					startingIntersection = -1;
				}
				else {
					// Start from first intersection
					// Store points from odd intersections to even
					startingIntersection = 0;
				}



				for (int i = startingIntersection; i < ((int)pointsBeforeIntersections.size()); i = i + 2) {
					vector<vec2> line;

					// Indicates the starting index of plottedPoints
					int startPoint;
					if (i == -1) {
						startPoint = 0;
					}
					else {
						startPoint = pointsBeforeIntersections[i];
					}

					// Indicates the ending index of plottedPoints
					int endPoint;
					if (i == pointsBeforeIntersections.size() - 1) {
						// If the line does not end in the erasing zone
						endPoint = plottedPoints[lineNum].size() - 1;
					}
					else {
						endPoint = pointsBeforeIntersections[i + 1];
					}

					// Loops through and saves the correct points
					for (int j = startPoint; j <= endPoint; j++) {
						line.push_back(plottedPoints[lineNum][j]);
					}

					// All points for these intersections now added to line
					linesAfterEreasing2.push_back(line);
				}
			}
		}
		// End of draw type check
	}

	// All lines have now been checked

	// Removing now incorrect lines
	// Removing from back to front so as not to affect placement of other lines to be dealt with
	for (int i = lineIndices.size() - 1; i >= 0; i--) {
		controlPoints.erase(controlPoints.begin() + lineIndices[i]);
		//knots.erase(knots.begin() + lineIndices[i]);
		lineTypes.erase(lineTypes.begin() + lineIndices[i]);
	}

	for (int i = lineIndices2.size() - 1; i >= 0; i--) {
		plottedPoints.erase(plottedPoints.begin() + lineIndices2[i]);
		lineTypes2.erase(lineTypes2.begin() + lineIndices2[i]);
	}

	// Adding new lines
	for (int i = 0; i < linesAfterEreasing.size(); i++) {
		if (linesAfterEreasing[i].size() >= order) {
			controlPoints.push_back(linesAfterEreasing[i]);
			lineTypes.push_back(drawType);

			/*
			vector <double> temp;
			knots.push_back(temp);
			recalculateKnots(controlPoints.size() - 1);
			vector <vec2> temp2;
			plottedPoints.push_back(temp2);
			for (double u = 0; u <= 1; u += 1.0 / uStep) {
				calculateBSpline(u, controlPoints.size() - 1);
			}
			calculateBSpline(0.99999999999, controlPoints.size() - 1);
			*/
		}
	}

	for (int i = 0; i < linesAfterEreasing2.size(); i++) {
		if (linesAfterEreasing2[i].size() > 1) {
			plottedPoints.push_back(linesAfterEreasing2[i]);
			lineTypes2.push_back(drawType);
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

	glTranslatef(xTrans, yTrans, zTrans);
	if (renderGround) {
		glRotatef(rotateX, 1.0f, 0.0f, 0.0f);
		glRotatef(rotateZ, 0.0f, 0.0f, 1.0f);
	}
	glScalef(scale, scale, 1.f);

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
	if (points.size() >= order) {
		if (points.size() > 0) {
			if (abs(points.front().x - mouseX) < selectDistance && abs(points.front().y - mouseY) < selectDistance) {
				if (points.back().x != points.front().x && points.back().y != points.front().y) {
					points.push_back(vec2(points.front().x, points.front().y));
				}
			}
			else {
				for (int i = 0; i < controlPoints.size(); i++) {
					if (lineTypes[i] == drawType) {
						if (abs(controlPoints[i].front().x - mouseX) < selectDistance && abs(controlPoints[i].front().y - mouseY) < selectDistance) {
							// If within selecting distance of start of line of same type
							// Make end of drawn line = start of line i
							points.push_back(controlPoints[i].front());
						}
						else if (abs(controlPoints[i].back().x - mouseX) < selectDistance && abs(controlPoints[i].back().y - mouseY) < selectDistance) {
							// If within selecting distance of start of line of same type
							// Make end of drawn line = end of line i
							points.push_back(controlPoints[i].back());
						}
					}
				}
			}
		}

		bool mergeLines = false;
		bool mergeAtStart = false;
		bool mergeAtEnd = false;
		int mergeTypeStart = -1;
		int mergeTypeEnd = -1;
		int mergeIndexStart = -1;
		int mergeIndexEnd = -1;
		
		for (int i = 0; i < controlPoints.size(); i++) {
			if (lineTypes[i] == drawType) {
				if (points.front() == controlPoints[i].front()) {
					// If the first point of the drawing line matches exactly with the
					// first point of another line of the same type
					// means that the drawing line starts at the start of the matching line
					mergeLines = true;
					mergeAtStart = true;
					mergeTypeStart = NEW_START_AT_START;
					mergeIndexStart = i;
				}
				else if (points.front() == controlPoints[i].back()) {
					// New line starts at the end of the ith line
					mergeLines = true;
					mergeAtStart = true;
					mergeTypeStart = NEW_START_AT_END;
					mergeIndexStart = i;
				}
				else if (points.back() == controlPoints[i].front()) {
					// New line ends at start of the ith line
					mergeLines = true;
					mergeAtEnd = true;
					mergeTypeEnd = NEW_END_AT_START;
					mergeIndexEnd = i;
				}
				else if (points.back() == controlPoints[i].back()) {
					// New line ends at end of the ith line
					mergeLines = true;
					mergeAtEnd = true;
					mergeTypeEnd = NEW_END_AT_END;
					mergeIndexEnd = i;
				}
			}
		}

		if (!mergeLines) {
			// Normal behaviours
			controlPoints.push_back(points);
		}
		else {
			vector <vec2> mergedLine;
			if (mergeTypeStart == NEW_START_AT_START && mergeTypeEnd == -1) {
				// <--d--o-----> goes to --d-->----->
				// example of line
				// d is the drawn line, o is the joining point for the two lines

				// Push back drawn line in reverse order, EXCEPT for the last point
				// Last point of this line, is the same as the last point for the second line
				// Do not want to dupelicate it and give it more weight than all other points
				for (int i = (int)points.size() - 1; i > 0; i--) {
					mergedLine.push_back(points[i]);
				}

				// Push back points of other line in same order
				for (int i = 0; i < controlPoints[mergeIndexStart].size(); i++) {
					mergedLine.push_back(controlPoints[mergeIndexStart][i]);
				}

			}
			else if (mergeTypeStart == NEW_START_AT_END && mergeTypeEnd == -1) {
				// ----->o--d--> goes to ----->--d-->

				// Again, do not want to duplicate shared point
				for (int i = 0; i < controlPoints[mergeIndexStart].size() - 1; i++) {
					mergedLine.push_back(controlPoints[mergeIndexStart][i]);
				}

				for (int i = 0; i < points.size(); i++) {
					mergedLine.push_back(points[i]);
				}
			}
			else if (mergeTypeEnd == NEW_END_AT_START && mergeTypeStart == -1) {
				// --d-->o-----> goes to --d-->----->

				for (int i = 0; i < points.size() - 1; i++) {
					mergedLine.push_back(points[i]);
				}

				for (int i = 0; i < controlPoints[mergeIndexEnd].size(); i++) {
					mergedLine.push_back(controlPoints[mergeIndexEnd][i]);
				}
			}
			else if (mergeTypeEnd == NEW_END_AT_END && mergeTypeStart == -1) {
				// ----->o<--d-- goes to ----->--d-->

				for (int i = 0; i < controlPoints[mergeIndexEnd].size() - 1; i++) {
					mergedLine.push_back(controlPoints[mergeIndexEnd][i]);
				}

				for (int i = (int)points.size() - 1; i > 0; i--) {
					mergedLine.push_back(points[i]);
				}
			}
			// Below are cases where the drawn line starts and ends at end points of other lines
			else if (mergeTypeStart == NEW_START_AT_START && mergeTypeEnd == NEW_END_AT_START) {
				// <-----o<--d--o-----> goes to ----->--d-->----->

				for (int i = controlPoints[mergeIndexEnd].size() - 1; i > 0; i--) {
					mergedLine.push_back(controlPoints[mergeIndexEnd][i]);
				}

				for (int i = (int)points.size() - 1; i > 0; i--) {
					mergedLine.push_back(points[i]);
				}

				if (mergeIndexStart != mergeIndexEnd) {
					// If not merging a broken line into a loop
					// Push back points of other line in same order
					for (int i = 0; i < controlPoints[mergeIndexStart].size(); i++) {
						mergedLine.push_back(controlPoints[mergeIndexStart][i]);
					}
				}
				else {
					// Creating a loop

					// Pushing back the last point of points
					// As it was not pushed before
					mergedLine.push_back(points.back());

					// This will prevent duplicating the already completed line
				}
			}
			else if (mergeTypeStart == NEW_START_AT_START && mergeTypeEnd == NEW_END_AT_END) {
				// ----->o<--d--o-----> goes to ----->--d-->----->

				for (int i = 0; i < controlPoints[mergeIndexEnd].size() - 1; i++) {
					mergedLine.push_back(controlPoints[mergeIndexEnd][i]);
				}

				for (int i = (int)points.size() - 1; i > 0; i--) {
					mergedLine.push_back(points[i]);
				}

				if (mergeIndexStart != mergeIndexEnd) {
					// If not merging a broken line into a loop
					// Push back points of other line in same order
					for (int i = 0; i < controlPoints[mergeIndexStart].size(); i++) {
						mergedLine.push_back(controlPoints[mergeIndexStart][i]);
					}
				}
				else {
					// Creating a loop

					// Pushing back the last point of points
					// As it was not pushed before
					mergedLine.push_back(points.back());
				}
			}
			else if (mergeTypeStart == NEW_START_AT_END && mergeTypeEnd == NEW_END_AT_START) {
				// ----->o--d-->o-----> goes to ----->--d-->----->

				for (int i = 0; i < controlPoints[mergeIndexStart].size() - 1; i++) {
					mergedLine.push_back(controlPoints[mergeIndexStart][i]);
				}

				for (int i = 0; i < points.size() - 1; i++) {
					mergedLine.push_back(points[i]);
				}

				if (mergeIndexStart != mergeIndexEnd) {
					// If not merging a broken line into a loop
					// Push back points of other line in same order
					for (int i = 0; i < controlPoints[mergeIndexEnd].size(); i++) {
						mergedLine.push_back(controlPoints[mergeIndexEnd][i]);
					}
				}
				else {
					// Creating a loop

					// Pushing back the last point of points
					// As it was not pushed before
					mergedLine.push_back(points.back());
				}
			}
			else if (mergeTypeStart == NEW_START_AT_END && mergeTypeEnd == NEW_END_AT_END) {
				// ----->o--d-->o<-----

				for (int i = 0; i < controlPoints[mergeIndexStart].size() - 1; i++) {
					mergedLine.push_back(controlPoints[mergeIndexStart][i]);
				}

				for (int i = 0; i < points.size() - 1; i++) {
					mergedLine.push_back(points[i]);
				}

				if (mergeIndexStart != mergeIndexEnd) {
					// If not merging a broken line into a loop
					// Push back points of other line in same order
					for (int i = controlPoints[mergeIndexEnd].size() - 1; i >= 0; i--) {
						mergedLine.push_back(controlPoints[mergeIndexEnd][i]);
					}
				}
				else {
					// Creating a loop

					// Pushing back the last point of points
					// As it was not pushed before
					mergedLine.push_back(points.back());
				}
			}
			
			if (mergeIndexStart != -1) {
				// Must erease modified line everywhere it is referenced
				plottedPoints.erase(plottedPoints.begin() + mergeIndexStart);
				controlPoints.erase(controlPoints.begin() + mergeIndexStart);
				lineTypes.erase(lineTypes.begin() + mergeIndexStart);
				lineTypes2.erase(lineTypes2.begin() + mergeIndexStart);
				knots.erase(knots.begin() + mergeIndexStart);

				if (mergeIndexEnd != -1 && mergeIndexStart != mergeIndexEnd) {
					if (mergeIndexEnd > mergeIndexStart) {
						// Decrement to accommodate the changes made by the first removal
						mergeIndexEnd--;
					}

					plottedPoints.erase(plottedPoints.begin() + mergeIndexEnd);
					controlPoints.erase(controlPoints.begin() + mergeIndexEnd);
					lineTypes.erase(lineTypes.begin() + mergeIndexEnd);
					lineTypes2.erase(lineTypes2.begin() + mergeIndexEnd);
					knots.erase(knots.begin() + mergeIndexEnd);
				}
			}
			else if (mergeIndexEnd != -1) {
				plottedPoints.erase(plottedPoints.begin() + mergeIndexEnd);
				controlPoints.erase(controlPoints.begin() + mergeIndexEnd);
				lineTypes.erase(lineTypes.begin() + mergeIndexEnd);
				lineTypes2.erase(lineTypes2.begin() + mergeIndexEnd);
				knots.erase(knots.begin() + mergeIndexEnd);
			}
			
			// Correct line type was already pushed in correct location when line 
			// was started
			// So just need to push the merged line
			controlPoints.push_back(mergedLine);


		}

		vector <double> temp;
		knots.push_back(temp);
		recalculateKnots(controlPoints.size() - 1);
		vector <vec2> temp2;
		lineTypes2.push_back(lineTypes.back());
		plottedPoints.push_back(temp2);
		for (double u = 0; u <= 1; u += 1.0 / uStep) {
			calculateBSpline(u, controlPoints.size() - 1);
		}
		calculateBSpline(0.99999999999, controlPoints.size() - 1);
		points.clear();
	}
}


void createRandomizedPlane()
{
	cout << "Create new plane";

	getBoundingCirclesFromMountainLines(&mountainCenters, &mountainRadii);

	vector<vector<Edge>> partition;

	partitionEdgesIntoSlices(GROUND, canvasHeight*-0.5, canvasHeight*0.5, NUMBER_OF_SLICES, &partition);

	//noise = PerlinNoise(OCTAVES, PERSISTENCE, RANGE, canvasWidth, canvasHeight, mountainCenters, mountainRadii);
	mountainNoise.generateMountainNoise(OCTAVES, PERSISTENCE, RANGE, canvasWidth, canvasHeight, mountainCenters, mountainRadii);

	landNoise.generateLandmassNoise(LAND_OCTAVES, LAND_PERSISTENCE, LAND_RANGE, canvasWidth, canvasHeight, -0.5f*canvasHeight, &partition);

	//Initialize Plane
	generatePlane(PLANE_DIM_X, PLANE_DIM_Z, canvasHeight, canvasWidth, &terrainPoints, &terrainNormals, &terrainUVs, &terrainIndices);

	float scale = 1.f*std::min(canvasHeight, canvasHeight) / (float)RANGE;
	float landScale = 0.125f*std::min(canvasHeight, canvasHeight) / (float)LAND_RANGE;

	//Add perlin noise to plane
	for (unsigned int i = 0; i < terrainPoints.size(); i++)
	{
		vec2 uv = terrainUVs[i];
		vec3 point = terrainPoints[i];

		float noiseContribution = landNoise.get(-uv.x + 1.f, uv.y)*landScale + mountainNoise.get(-uv.x + 1.f, uv.y)*scale;

		terrainPoints[i] = vec3(
			point.x, 
			std::max(noiseContribution, 0.f),
			point.z);
	}

	//Calculate normals
	for (int i = 0; i < PLANE_DIM_Z; i++)
	{
		for (int j = 0; j < PLANE_DIM_X; j++)
		{
			int i0 = std::max(i - 1, 0);
			int i1 = std::min(i + 1, PLANE_DIM_Z-1);

			int j0 = std::max(j - 1, 0);
			int j1 = std::min(j + 1, PLANE_DIM_X-1);

			vec3 xTangent = terrainPoints[i*PLANE_DIM_X + j1] 
				- terrainPoints[i*PLANE_DIM_X + j0];
			vec3 zTangent = terrainPoints[i1*PLANE_DIM_X + j]
				- terrainPoints[i0*PLANE_DIM_X + j1];
			terrainNormals[i*PLANE_DIM_X + j] = cross(xTangent, zTangent);
		}
	}

	r.loadModelBuffer(&terrainPoints, &terrainNormals, &terrainIndices);
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
			/*
			Removed as will not be doing tree stuff
			drawType = TREES;
			if (drawing) {
			lineTypes.back() = drawType;
			}
			*/
		}
		else if (key == GLFW_KEY_M && action == GLFW_PRESS) {
			drawType = MOUNTAIN;
			if (drawing) {
				lineTypes.back() = drawType;
			}
		}

		else if (key == GLFW_KEY_E && action == GLFW_PRESS) {
			if (!drawing) {
				erasing = !erasing;
			}
		}
		else if (key == GLFW_KEY_D && action == GLFW_PRESS){
			debugging = !debugging;
		}

		else if (key == GLFW_KEY_ENTER && action == GLFW_PRESS) {
			renderTerrain = true;
			//createRandomizedPlane();
			if (drawing) {
				finishLine();
			}
		}
		else if (key == GLFW_KEY_S && action == GLFW_PRESS) {
			if (saveMap("maps/saved.map", &controlPoints, &lineTypes))
				printf("Map saved\n");

		}
		else if (key == GLFW_KEY_L && action == GLFW_PRESS){
			if (loadMap("maps/load.map", &controlPoints, &lineTypes))
			{
				printf("Map loaded\n");

				knots.clear();
				lineTypes2.clear();
				plottedPoints.clear();

				for (unsigned int i = 0; i < controlPoints.size(); i++)
				{
					vector <double> temp;
					knots.push_back(temp);
					recalculateKnots(i);
					vector <vec2> temp2;
					lineTypes2.push_back(lineTypes[i]);
					plottedPoints.push_back(temp2);
					for (double u = 0; u <= 1; u += 1.0 / uStep) {
						calculateBSpline(u, i);
					}
					calculateBSpline(0.99999999999, i);
				}
			

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
		else if ((key == GLFW_KEY_SPACE) && (action == GLFW_PRESS)) {
			createRandomizedPlane();
		}
		else if ((key == GLFW_KEY_1) && (action == GLFW_PRESS)) {
			if (wireframe)
			{
				wireframe = false;
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				glDisable(GL_CULL_FACE);
			}
			else
			{
				wireframe = true;
				glPolygonMode(GL_FRONT, GL_LINE);
				glEnable(GL_CULL_FACE);
				glCullFace(GL_BACK);
			}
		}
	}
	if (key == GLFW_KEY_I && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		yTrans += 0.01f;
	}
	else if (key == GLFW_KEY_K && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		yTrans -= 0.01f;
	}
	else if (key == GLFW_KEY_L && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		xTrans += 0.01f;
	}
	else if (key == GLFW_KEY_J && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		xTrans -= 0.01f;
	}
	else if (key == GLFW_KEY_O && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		zTrans += 0.01f;
	}
	else if (key == GLFW_KEY_U && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		zTrans -= 0.01f;
	}
}


void mouseClick(GLFWwindow *sender, int button, int action, int mods) {
	if (!renderTerrain) {
		if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
			drawing = true;

			if (!erasing) {
				for (int i = 0; i < controlPoints.size(); i++) {
					if (lineTypes[i] == drawType) {
						if (abs(controlPoints[i].front().x - mouseX) < selectDistance && abs(controlPoints[i].front().y - mouseY) < selectDistance) {
							// If within selecting distance of start of line of same type
							// Make start of drawn line = start of line i
							points.push_back(controlPoints[i].front());
						}
						else if (abs(controlPoints[i].back().x - mouseX) < selectDistance && abs(controlPoints[i].back().y - mouseY) < selectDistance) {
							// If within selecting distance of end of line of same type
							// Make start of drawn line = end of line i
							points.push_back(controlPoints[i].back());
						}
					}
				}
				lineTypes.push_back(drawType);
			}
			else {
				erase();
			}
		}
		else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
			drawing = false;

			if (!erasing) {
				finishLine();
			}
			else {
				knots.clear();
				plottedPoints.clear();
				lineTypes2.clear();

				for (int i = 0; i < controlPoints.size(); i++) {
					vector <double> temp;
					knots.push_back(temp);
					recalculateKnots(i);
					vector <vec2> temp2;
					lineTypes2.push_back(lineTypes[i]);
					plottedPoints.push_back(temp2);
					for (double u = 0; u <= 1; u += 1.0 / uStep) {
						calculateBSpline(u, i);
					}
					calculateBSpline(0.99999999999, i);
					points.clear();
				}
			}

			cout << "Num of lines: " << plottedPoints.size() << "\n";
		}
		else if (button == GLFW_MOUSE_BUTTON_MIDDLE && ((action == GLFW_PRESS) || (action == GLFW_RELEASE)))
		{
			panningCamera = !panningCamera;
		}
	}
	else
	{
		if (button == GLFW_MOUSE_BUTTON_LEFT){
			if (action == GLFW_PRESS){
				rotatingCamera = true;
			}
			else if (action == GLFW_RELEASE){
				rotatingCamera = false;
			}
				

		}
	}
}


void mousePos(GLFWwindow *sender, double x, double y) {	
	mouseX = (2 * x / w) - 1 - xTrans;
	mouseY = (-2 * y / h) + 1 - yTrans;

	mouseX = mouseX / scale;
	mouseY = mouseY / scale;

	//cout << "mouseX: " << mouseX << " mouseY: " << mouseY << "\n";

	vec2 diff = vec2((2 * x / w) - 1, (-2 * y / h) + 1) - lastPos;
	lastPos = vec2((2 * x / w) - 1, (-2 * y / h) + 1);
	
	if (!renderTerrain) {
		if (drawing) {
			if (!erasing) {
				if (points.size() > 1) {
					double distanceX = abs(points[points.size() - 2].x - points.back().x);
					double distanceY = abs(points[points.size() - 2].y - points.back().y);
					if (distanceX >= selectDistance || distanceY >= selectDistance) {
						points.push_back(vec2(mouseX, mouseY));
					}
					else {
						points[points.size() - 1] = vec2(mouseX, mouseY);
					}
				}
				else {
					points.push_back(vec2(mouseX, mouseY));
				}
			}
			else {
				erase();
			}
		}
		else if (panningCamera)
		{
			yTrans += diff.y;
			xTrans += diff.x;
		}
	}
	else
	{
		if (rotatingCamera)
		{
			lookat.rotateViewAround(-diff.x, -diff.y);
		}
	}
}

void scroll(GLFWwindow *sender, double x, double y) {
	
	if (!renderTerrain)
	{
		scale = std::max(scale + 0.1f * y, 0.1);
		cout << scale << "\n";
	}
	else
	{
		if (y > 0)
			lookat.zoom(1/1.1f);
		else
			lookat.zoom(1.1f);
	}
}



void init()
{
	//Setup renderer
	lookat = Camera(vec3(0.f, 0.f, -1.f), vec3(0.f, 1.f, 0.f), vec3(0.f, 0.f, 2.f), MODELVIEWER_CAMERA);

	r.projection = perspectiveMatrix(0.1f, 30.f, 80.f);
	r.loadCamera(&lookat);

	glClearColor(0.5f, 0.5f, 0.5f, 1.f);

	createRandomizedPlane();


}


int main() {
	window = initializeWindow();
	if (window == NULL) {
		return -1;
	}

	r = Renderer(window);

	drawing = false;

	glfwSetKeyCallback(window, keyboard);
	glfwSetMouseButtonCallback(window, mouseClick);
	glfwSetCursorPosCallback(window, mousePos);
	glfwSetScrollCallback(window, scroll);

	init();

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
