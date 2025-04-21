#version 430 core

in vec3 a_position;

uniform mat4 ModelViewMatrix;
uniform mat4 ProjectionMatrix;

out vec3 worldPos;

void main(void)
{
    mat4 MVP    = ProjectionMatrix * ModelViewMatrix;
    gl_Position = MVP * (vec4(a_position, 1.0));
    worldPos    = a_position;
} 
