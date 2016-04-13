#ifndef MAPLOADER_H
#define MAPLOADER_H

#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include "glm/glm.hpp"

using namespace std;
using namespace glm;

bool loadMap(string fileName, vector<vector<vec2>>* lineSegments, vector<int>* lineTypes);

bool saveMap(string fileName, vector < vector<vec2>>* lineSegments, vector<int>* lineTypes);















#endif