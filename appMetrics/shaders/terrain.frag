#version 430 core

uniform vec2 u_worldMin;
uniform vec2 u_worldSize;
uniform sampler2D u_texture;
uniform float u_cursorRadius;
uniform vec2  u_cursorWorld;
uniform vec4  u_cursorColor;

in vec3 worldPos;

out vec4 fragment;

void main()
{
    vec2 uv = (worldPos.xy - u_worldMin)/u_worldSize;
    vec4 c = texture(u_texture, uv);
    if (u_cursorRadius > 0) {
        float d = length(u_cursorWorld - worldPos.xy);
        float t = smoothstep(u_cursorRadius, 0.5*u_cursorRadius, d);
        vec4 ccursor = (1 - u_cursorColor.a)*c + u_cursorColor.a*u_cursorColor;
        c = mix(c, ccursor, t);
    }
    fragment = vec4(c.rgb, 1.0);
}
