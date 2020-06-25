


















in vec3 residue_coord;
in float border_fade;
flat in ivec3 aabb_low;
flat in ivec3 aabb_high;
// TODO: Add bias (based on normal vector) to deal with floor/ceil
// rounding errors. Also hide this mess in a glsl file somewhere.
uniform vec3 eye_relative_group_origin;
uniform sampler3D chunk_blocks;
uniform bool chunk_debug;
uniform int far_plane_squared;
out vec4 color;
void main() {
    if (chunk_debug) {
        int x_floor = int(floor(residue_coord.x));
        int y_floor = int(floor(residue_coord.y));
        int z_floor = int(floor(residue_coord.z));
        color = vec4(x_floor & 1, y_floor & 1, z_floor & 1, 1);
        return;
    }
    const float d = BORDER_WIDTH;
    float x0 = eye_relative_group_origin.x;
    float y0 = eye_relative_group_origin.y;
    float z0 = eye_relative_group_origin.z;
    vec3 slope = vec3(residue_coord) - eye_relative_group_origin;
    float xm = slope.x;
    float ym = slope.y;
    float zm = slope.z;
    float rcp = 1.0/GROUP_SIZE;
    float best_t = 1.0 / 0.0;
    vec4 best_color = vec4(0,0,0,0);
    vec3 best_coord = vec3(0,0,0);
    int iter = 0;
    int x_init = int(xm > 0 ? ceil(residue_coord.x) 
                            : floor(residue_coord.x));
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
        vec4 lookup_color = texture(chunk_blocks, texcoord);
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
    int y_init = int(ym > 0 ? ceil(residue_coord.y) 
                            : floor(residue_coord.y));
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
        vec4 lookup_color = texture(chunk_blocks, texcoord);
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
    int z_init = int(zm > 0 ? ceil(residue_coord.z) 
                            : floor(residue_coord.z));
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
        vec4 lookup_color = texture(chunk_blocks, texcoord);
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
    float fog_fade = clamp(2 * (1 - dist_squared/far_plane_squared), 0, 1);
    color = vec4(best_color.rgb * fog_fade, 1);
}

