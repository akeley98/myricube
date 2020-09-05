































// Dumb vertex shader for filling the whole screen, when drawn with 4
// vertices as triangle strip.
out vec2 normalized_screen_xy;
const float MAX_DEPTH = 0.9999998807907104;

void main() {
    if (gl_VertexID == 0) gl_Position = vec4(-1, -1, MAX_DEPTH, 1.0);
    if (gl_VertexID == 1) gl_Position = vec4(+1, -1, MAX_DEPTH, 1.0);
    if (gl_VertexID == 2) gl_Position = vec4(-1, +1, MAX_DEPTH, 1.0);
    if (gl_VertexID == 3) gl_Position = vec4(+1, +1, MAX_DEPTH, 1.0);
    normalized_screen_xy = gl_Position.xy;
}
