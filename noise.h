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

const int range = 5;

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

	float getValueAt(int x, int y);

	float getValueAt(float x, float y);
};

class PerlinNoise
{
public:
	vector<Octave> octaves;

	PerlinNoise();

	PerlinNoise(int _maxFreq, float _persistence);

	float get(float x, float y);

	glm::vec3 getNormal(float x, float y);
};


#endif
