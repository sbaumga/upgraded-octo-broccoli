#ifndef NOISE_CPP
#define NOISE_CPP

#include "noise.h"


float interpolation(float a, float b, float u)
{
    float angle = u*M_PI;
    u = (1 - cos(angle))*0.5;

    return (1-u)*a + u*b;
}

Octave::Octave():
    freq(0), persistence(0) {}


Octave::Octave(int _freq, float _persistence):
    freq(_freq), persistence(_persistence), 
    dimension(pow(2, _freq)*range), amplitude(pow(_persistence, _freq))
{

    for(int i=0; i<dimension; i++)
    {
        for(int j=0; j<dimension; j++)
        {
            float randomValue = ((float)rand()/RAND_MAX)*2.f- 0.5f;
            noise.push_back(randomValue*amplitude);
        }
    }
}

Octave::Octave(int _freq, float _persistence, int _range) :
freq(_freq), persistence(_persistence),
dimension(pow(2, _freq)*_range), amplitude(pow(_persistence, _freq))
{

	for (int i = 0; i<dimension; i++)
	{
		for (int j = 0; j<dimension; j++)
		{
			if (abs(i - dimension / 2) + abs(j - dimension / 2) <= dimension / 4)
			{
				float randomValue = ((float)rand() / RAND_MAX)*2.f - 0.5f;
				noise.push_back(randomValue*amplitude);
			}
			else
				noise.push_back(0.f);
		}
	}
}

//Perlin noise constructor for mountain
Octave::Octave(int _freq, float _persistence, int _range, float xWidth, float yWidth, const vector<vec2>& centers, const vector<float>& radii):
freq(_freq), persistence(_persistence),
dimension(pow(2, _freq)*_range), amplitude(pow(_persistence, _freq))
{
	float radius = std::max(xWidth, yWidth) / ((float)(dimension - 1));
	float width = xWidth/ ((float)(dimension - 1));
	float height = yWidth/ ((float)(dimension - 1));
	float xOffset = -xWidth*0.5f;
	float yOffset = -yWidth*0.5f;

	for (int i = 0; i<dimension; i++)
	{
		for (int j = 0; j<dimension; j++)
		{
			vec2 center(
				xOffset + width*(float)j,
				yOffset + height*(float)i);

			/*if (withinRadius(center, radius, centers, radii))
			{
				float randomValue = ((float)rand() / RAND_MAX)*2.f - 0.5f;
				noise.push_back(randomValue*amplitude);
			}
			else
				noise.push_back(0.f);*/

			float coverage = areaCovered(center, radius, centers, radii);

			float randomValue = ((float)rand() / RAND_MAX)*2.f - 0.5f;
			noise.push_back(randomValue*amplitude*coverage);
		}
	}
}

//Constructor for landmass
Octave::Octave(int _freq, float _persistence, int _range, float xWidth, float yWidth, float minY, vector<vector<Edge>>* partitions):freq(_freq), persistence(_persistence),
dimension(pow(2, _freq)*_range), amplitude(pow(_persistence, _freq))
{
	float radius = std::max(xWidth, yWidth) / ((float)(dimension - 1));
	float width = xWidth / ((float)(dimension - 1));
	float height = yWidth / ((float)(dimension - 1));
	float xOffset = -xWidth*0.5f;
	float yOffset = -yWidth*0.5f;

	int comparisons = std::max(5-freq, 1);

	

	for (int i = 0; i<dimension; i++)
	{
		for (int j = 0; j<dimension; j++)
		{
			vec2 center(
				xOffset + width*(float)j,
				yOffset + height*(float)i);

			float coverage; 
/*			if (inPartitionedPolygon(center, partitions, minY, minY + yWidth))
				coverage = 1.f;
			else
				coverage = -1.f;*/

			coverage = polygonCoverage(center, partitions, minY, minY + yWidth, width, height, comparisons);

			float randomValue = ((float)rand() / RAND_MAX);
			if (coverage >= 0.99f)
				noise.push_back((randomValue*amplitude*0.5f + 0.5f)*coverage);
			else
				noise.push_back(-randomValue*amplitude);
		}
	}
}



bool withinRadius(vec2 center, float radius, const vector<vec2>& centers, const vector<float>& radii)
{

	for (unsigned int i = 0; i < centers.size(); i++)
	{
		vec2 diff = center - centers[i];
		if (dot(diff, diff) < (radius + radii[i])*(radius + radii[i]))
			return true;
	}

	return false;
}

float clamp(float x, float low, float high)
{
	return std::max(std::min(x, high), low);
}

//Will return an approximation of the covered area
//	Spheres are modelled as intersecting squares lying on the same axis. 
//	The diagonal of the square is equal to the diameter
float areaCovered(vec2 center, float radius, const vector<vec2>& centers, const vector<float>& radii)
{
	float totalArea = 2.f*radius*radius;		//Approximation given by square with 2*radius as the diagonal edge
	float intersectedArea = 0.f;

	for (unsigned int i = 0; i < centers.size(); i++)
	{
		vec2 diff = center - centers[i];
		float d = length(diff);
		float diagonal = clamp(radii[i] + radius - d, 0.f, 2.f*radius);

		intersectedArea += diagonal*diagonal*0.5f;
	}

	return std::min(intersectedArea / totalArea, 1.f);
}

float Octave::getValueAt(int x, int y)
{
    x = std::max(std::min(x, dimension-1), 0);
    y = std::max(std::min(y, dimension-1), 0);

    //if(noise[x*dimension + y] < 0.0)
      //  printf("Error = %f\n", noise[x*dimension+y]);

    return noise[x*dimension + y];
}

float Octave::getValueAt(float x, float y)
{
    float f_dimension = (float)dimension;

    int x_i = (int)(f_dimension*x);
    int y_i = (int)(f_dimension*y);

    float x_r = x*dimension -(float)x_i;
    float y_r = y*dimension - (float)y_i;

    float x0y0 = getValueAt(x_i, y_i);
    float x1y1 = getValueAt(x_i+1, y_i+1);
    float x0y1 = getValueAt(x_i, y_i+1);
    float x1y0 = getValueAt(x_i+1, y_i);

    return ( (1-x_r)*interpolation(x0y0, x0y1, y_r) + x_r*interpolation(x1y0, x1y1, y_r)
            + (1-y_r)*interpolation(x0y0, x1y0, x_r) + y_r*interpolation(x0y1, x1y1, x_r) )*0.5; 
}

PerlinNoise::PerlinNoise()
{
	srand(time(0));
}


PerlinNoise::PerlinNoise(int _maxFreq, float _persistence)
{
	srand(time(0));

    for(int i=0; i<_maxFreq; i++)
    {
        octaves.push_back(Octave(i, _persistence));
    }
}

PerlinNoise::PerlinNoise(int _maxFreq, float _persistence, int _range)
{
	srand(time(0));

	for (int i = 0; i<_maxFreq; i++)
	{
		octaves.push_back(Octave(i, _persistence, _range));
	}
}

PerlinNoise::PerlinNoise(int _maxFreq, float _persistence, int _range, float xWidth, float yWidth, const vector<vec2>& centers, const vector<float>& radii)
{
	srand(time(0));

	for (int i = 0; i<_maxFreq; i++)
	{
		octaves.push_back(Octave(i, _persistence, _range, xWidth, yWidth, centers, radii));
	}
}

void PerlinNoise::generateMountainNoise(int _maxFreq, float _persistence, int _range, float xWidth, float yWidth, const vector<vec2>& centers, const vector<float>& radii)
{
	srand(time(0));
	octaves.clear();
	octaves.clear();

	for (int i = 0; i<_maxFreq; i++)
	{
		octaves.push_back(Octave(i, _persistence, _range, xWidth, yWidth, centers, radii));
	}
}

void PerlinNoise::generateLandmassNoise(int _maxFreq, float _persistence, int _range, float xWidth, float yWidth, float minY, vector<vector<Edge>>* partitions)
{
	srand(time(0));
	octaves.clear();
	octaves.clear();

	for (int i = 0; i < _maxFreq; i++)
	{
		octaves.push_back(Octave(i, _persistence, _range, xWidth, yWidth, minY, partitions));
	}
}

float PerlinNoise::get(float x, float y)
{
    float value = 0.0;

    for(int i=0; i<octaves.size(); i++)
    {
        value += octaves[i].getValueAt(x, y);
    }

    return value;
}

glm::vec3 PerlinNoise::getNormal(float x, float y)
{
	float inc = 0.0000001;

	glm::vec3 base(0.f, get(x, y), 0.f);
	glm::vec3 u(inc, get(x + inc, y), 0.f);
	glm::vec3 v(0.f, get(x, y + inc), inc);

	return normalize(cross((v - base), (u - base)));
}

float polygonCoverage(vec2 point, vector<vector<Edge>>* partitions, float minY, float maxY, float width, float height, unsigned int comparisons)
{
	float coverage = 0.f;
	float comparisonsMade = 0.f;

	for (unsigned int i = 0; i < comparisons; i++)
	{
		for (unsigned int j = 0; j < comparisons; j++)
		{
			vec2 offset(
				width*((float)j) / (float)(comparisons - 1) - width*0.5f,
				height*((float)i) / (float)(comparisons - 1) - height*0.5f);

			if (inPartitionedPolygon(point + offset, partitions, minY, maxY))
				coverage++;
			comparisonsMade++;
		}
	}

	coverage = coverage / comparisonsMade++;

	return coverage;
}

//Detecting if within polygon
bool inPolygon(vec2 point, vector<Edge>* edges)
{
	unsigned int intersections = 0;

	for (unsigned int i = 0; i < edges->size(); i++)
	{
		vec2 p1 = *edges->at(i).p1;
		vec2 dir = *edges->at(i).p2 - p1;

		float t = (point.y - p1.y) / (dir.y);
		vec2 intersectedPoint = p1 + dir*t;

		if ((intersectedPoint.x > point.x) && (t > 0.f) && (t < 1.f))
			intersections++;
	}

	return (intersections & 1);
}

bool inPartitionedPolygon(vec2 point, vector<vector<Edge>>* partition, float minY, float maxY)
{
	float increment = (maxY - minY) / (float)(partition->size() - 1);

	unsigned int index = clamp((point.y - minY) / increment, 0, partition->size()-1);

	return inPolygon(point, &partition->at(index));
}


#endif
