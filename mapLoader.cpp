#include "mapLoader.h"


bool saveMap(string fileName, vector < vector<vec2>>* lineSegments, vector<int>* lineTypes)
{
	std::ofstream f(fileName);
	if (!f.is_open())
	{
		printf("Could not save file\n");
		return false;
	}

	for (unsigned int i = 0; i < lineSegments->size(); i++)
	{
		f << "Line " << lineTypes->at(i) << endl;

		for (unsigned int j = 0; j < lineSegments->at(i).size(); j++)
		{
			f << "p " << lineSegments->at(i)[j].x << " "
				<< lineSegments->at(i)[j].y << endl;
		}
	}

	f.close();

	return true;
}

bool loadMap(string fileName, vector<vector<vec2>>* lineSegments, vector<int>* lineTypes)
{

	std::ifstream f(fileName);
	if (!f.is_open())
	{
		printf("File could not be opened\n");
		return false;
	}

	lineSegments->clear();
	lineTypes->clear();

	std::string line;
	vec2 point;
	int lineType = 0;

	while (getline(f, line))
	{
		if (sscanf(line.c_str(), "Line %d", &lineType) == 1)
		{
			printf("Linetype = %d\n", lineType);

			//Add new line
			lineTypes->push_back(lineType);
			vector<vec2> newLineSeg;
			lineSegments->push_back(newLineSeg);
		}
		else if (sscanf(line.c_str(), "p %f %f", &point.x, &point.y) == 2)
		{
			lineSegments->back().push_back(point);
		}
	}

	f.close();

	return true;
}