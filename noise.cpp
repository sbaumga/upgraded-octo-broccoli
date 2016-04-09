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

Octave::Octave(int _freq, float _persistence, int _range, float xWidth, float yWidth, const vector<vec2>& centers, const vector<float>& radii):
freq(_freq), persistence(_persistence),
dimension(pow(2, _freq)*_range), amplitude(pow(_persistence, _freq))
{
	float radius = ((float)(dimension - 1)) / std::max(xWidth, yWidth);
	float width = ((float)(dimension - 1)) / xWidth;
	float height = ((float)(dimension - 1)) / yWidth;
	float xOffset = -xWidth*0.5f;
	float yOffset = -yWidth*0.5f;

	for (int i = 0; i<dimension; i++)
	{
		for (int j = 0; j<dimension; j++)
		{
			vec2 center(
				xOffset + width*(float)j,
				yOffset + height*(float)i);

			if (withinRadius(center, radius, centers, radii))
			{
				float randomValue = ((float)rand() / RAND_MAX);
				noise.push_back(randomValue*amplitude);
			}
			else
				noise.push_back(0.f);
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


#endif
