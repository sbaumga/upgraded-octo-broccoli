#version 330 core
layout(location = 0) in vec3 vertexPosition;
layout(location = 1) in vec3 vertexColor;

uniform mat4 transform;

out vec3 color;

void main()
{	
	color - vertexColor;
 	gl_Position =  transform*vec4(vertexPosition, 1);
}
