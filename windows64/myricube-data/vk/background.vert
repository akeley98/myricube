#version 460
#extension GL_GOOGLE_include_directive : enable

// Dumb vertex shader for filling the whole screen, when drawn with 4
// vertices as triangle strip.
layout(location=0) out vec2 normalized_screen_xy;
const float MAX_DEPTH = 1.0;

void main() {
    switch (gl_VertexIndex) {
        case 0: gl_Position = vec4(+1, +1, MAX_DEPTH, 1.0); break;
        case 1: gl_Position = vec4(+1, -1, MAX_DEPTH, 1.0); break;
        case 2: gl_Position = vec4(-1, +1, MAX_DEPTH, 1.0); break;
        default:gl_Position = vec4(-1, -1, MAX_DEPTH, 1.0); break;
    }
    normalized_screen_xy = gl_Position.xy;
}
