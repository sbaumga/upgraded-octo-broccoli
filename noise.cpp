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
            float randomValue = ((float)rand()/RAND_MAX);
            noise.push_back(randomValue*amplitude);
        }
    }
}

float Octave::getValueAt(int x, int y)
{
    x = max(min(x, dimension-1), 0);
    y = max(min(y, dimension-1), 0);

    if(noise[x*dimension + y] < 0.0)
        printf("Error = %f\n", noise[x*dimension+y]);

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
    
}


PerlinNoise::PerlinNoise(int _maxFreq, float _persistence)
{
    for(int i=0; i<_maxFreq; i++)
    {
        octaves.push_back(Octave(i, _persistence));
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


#endif
