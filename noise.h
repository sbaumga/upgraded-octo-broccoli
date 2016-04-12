#ifndef NOISE_H
#define NOISE_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <glm/glm.hpp>
 
#include <vector>
#include <ctime>
#include <math.h>

#define M_PI 3.14159265359

using namespace std;
using namespace glm;

const int range = 3;

struct Edge {
	vec2* p1;
	vec2* p2;
};

class Octave
{
private:
	int dimension;
	float amplitude;

public:
	int freq;
	float persistence;
	vector<float> noise;

	Octave();

	Octave(int _freq, float _persistence);

	Octave(int _freq, float _persistence, int _range);

	Octave(int _freq, float _persistence, int _range, float xWidth, float yWidth, const vector<vec2>& centers, const vector<float>& radii);

	float getValueAt(int x, int y);

	float getValueAt(float x, float y);
};

class PerlinNoise
{
public:
	vector<Octave> octaves;

	PerlinNoise();

	PerlinNoise(int _maxFreq, float _persistence);

	PerlinNoise(int _maxFreq, float _persistence, int _range);

	PerlinNoise(int _maxFreq, float _persistence, int _range, float xWidth, float yWidth, const vector<vec2>& centers, const vector<float>& radii);		//Constructor for circle detection

	void generateMountainNoise(int _maxFreq, float _persistence, int _range, float xWidth, float yWidth, const vector<vec2>& centers, const vector<float>& radii);

	float get(float x, float y);

	glm::vec3 getNormal(float x, float y);
};

bool withinRadius(vec2 center, float radius, const vector<vec2>& centers, const vector<float>& radii);
float areaCovered(vec2 center, float radius, const vector<vec2>& centers, const vector<float>& radii);

#endif
