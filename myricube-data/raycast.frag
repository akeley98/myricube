


















in vec3 aabb_residue_coord;
in float border_fade;
flat in ivec3 aabb_low;
flat in ivec3 aabb_high;
flat in vec3 floor_ceil_fudge;
uniform vec3 eye_relative_group_origin;
uniform sampler3D chunk_blocks;
uniform bool chunk_debug;
uniform int far_plane_squared;
out vec4 color;
void main() {
    vec3 init_coord = aabb_residue_coord + floor_ceil_fudge;
    if (chunk_debug) {
        int x_floor = int(floor(init_coord.x));
        int y_floor = int(floor(init_coord.y));
        int z_floor = int(floor(init_coord.z));
        color = vec4(x_floor & 1, y_floor & 1, z_floor & 1, 1);
        return;
    }
    const float d = BORDER_WIDTH;
    float x0 = eye_relative_group_origin.x;
    float y0 = eye_relative_group_origin.y;
    float z0 = eye_relative_group_origin.z;
    vec3 slope = aabb_residue_coord - eye_relative_group_origin;
    float xm = slope.x;
    float ym = slope.y;
    float zm = slope.z;
    float rcp = 1.0/GROUP_SIZE;
    float best_t = 1.0 / 0.0;
    vec4 best_color = vec4(0,0,0,0);
    vec3 best_coord = vec3(0,0,0);
    int iter = 0;
    int x_init = int(xm > 0 ? ceil(init_coord.x)
                            : floor(init_coord.x));
    int x_end = xm > 0 ? aabb_high.x : aabb_low.x;
    int x_step = xm > 0 ? 1 : -1;
    float x_fudge = xm > 0 ? .25 : -.25;
    for (int x = x_init; x != x_end; x += x_step) {
        if (iter++ >= 255) { color = vec4(1,0,1,1); return; }
        float t = (x - x0) / xm;
        float y = y0 + ym * t;
        float z = z0 + zm * t;
        if (y < aabb_low.y || y > aabb_high.y) break;
        if (z < aabb_low.z || z > aabb_high.z) break;
        vec3 texcoord = vec3(x + x_fudge, y, z) * rcp;
        vec4 lookup_color = textureLod(chunk_blocks, texcoord, 0);
        if (lookup_color.a > 0 && t > 0) {
            if (best_t > t) {
                best_t = t;
                best_color = lookup_color;
                best_coord = vec3(x,y,z);
                if (y - floor(y + d) < d || z - floor(z + d) < d) {
                    best_color.rgb *= border_fade;
                }
            }
            break;
        }
    }
    int y_init = int(ym > 0 ? ceil(init_coord.y)
                            : floor(init_coord.y));
    int y_end = ym > 0 ? aabb_high.y : aabb_low.y;
    int y_step = ym > 0 ? 1 : -1;
    float y_fudge = ym > 0 ? .25 : -.25;
    for (int y = y_init; y != y_end; y += y_step) {
        if (iter++ >= 255) { color = vec4(1,0,1,1); return; }
        float t = (y - y0) / ym;
        float x = x0 + xm * t;
        float z = z0 + zm * t;
        if (x < aabb_low.x || x > aabb_high.x) break;
        if (z < aabb_low.z || z > aabb_high.z) break;
        vec3 texcoord = vec3(x, y + y_fudge, z) * rcp;
        vec4 lookup_color = textureLod(chunk_blocks, texcoord, 0);
        if (lookup_color.a > 0 && t > 0) {
            if (best_t > t) {
                best_t = t;
                best_color = lookup_color;
                best_coord = vec3(x,y,z);
                if (x - floor(x + d) < d || z - floor(z + d) < d) {
                    best_color.rgb *= border_fade;
                }
            }
            break;
        }
    }
    int z_init = int(zm > 0 ? ceil(init_coord.z)
                            : floor(init_coord.z));
    int z_end = zm > 0 ? aabb_high.z : aabb_low.z;
    int z_step = zm > 0 ? 1 : -1;
    float z_fudge = zm > 0 ? .25 : -.25;
    for (int z = z_init; z != z_end; z += z_step) {
        if (iter++ >= 255) { color = vec4(1,0,1,1); return; }
        float t = (z - z0) / zm;
        float x = x0 + xm * t;
        float y = y0 + ym * t;
        if (x < aabb_low.x || x > aabb_high.x) break;
        if (y < aabb_low.y || y > aabb_high.y) break;
        vec3 texcoord = vec3(x, y, z + z_fudge) * rcp;
        vec4 lookup_color = textureLod(chunk_blocks, texcoord, 0);
        if (lookup_color.a > 0 && t > 0) {
            if (best_t > t) {
                best_t = t;
                best_color = lookup_color;
                best_coord = vec3(x,y,z);
                if (x - floor(x + d) < d || y - floor(y + d) < d) {
                    best_color.rgb *= border_fade;
                }
            }
            break;
        }
    }
    if (best_color.a == 0) discard;
    vec3 disp = best_coord - eye_relative_group_origin;
    float dist_squared = dot(disp, disp);
    float raw_fog_fade = FOG_SCALAR * (1 - dist_squared/far_plane_squared);
    float fog_fade = clamp(raw_fog_fade, 0, 1);
    color = vec4(best_color.rgb * fog_fade, 1);
}
