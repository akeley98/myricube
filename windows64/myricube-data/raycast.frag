





























in vec3 aabb_residue_coord;
flat in ivec3 aabb_low;
flat in ivec3 aabb_high;
flat in vec3 aabb_half_normal;
uniform vec3 eye_relative_group_origin;
uniform bool chunk_debug;
uniform bool use_old_raycast;

out vec4 out_color;

vec4 read_group_voxel(ivec3 residue);

vec3 fog_color_from_world_direction(vec3 world_direction);

vec4 fog_border_color(
    vec3 base_color,
    float dist_squared,
    vec2 uv,
    vec3 fog_color);

void old_raycast(
    vec3 frag_residue_coord,
    out vec4 voxel_color,
    out vec3 disp,
    out vec2 uv);

void raycast(
    vec3 frag_residue_coord,
    out vec4 voxel_color,
    out vec3 disp,
    out vec2 uv);

void main()
{
    vec3 frag_residue_coord = aabb_residue_coord;

    vec4 voxel_color;
    vec3 disp;
    vec2 uv;

    if (chunk_debug) {
        ivec3 residue_floor =
            ivec3(floor(frag_residue_coord + aabb_half_normal));
        voxel_color = vec4(residue_floor & ivec3(1, 1, 1), 1);
        disp = frag_residue_coord - eye_relative_group_origin;
        uv = vec2(0.5, 0.5);
    }
    else if (use_old_raycast) {
        old_raycast(frag_residue_coord, voxel_color, disp, uv);
    }
    else {
        raycast(frag_residue_coord, voxel_color, disp, uv);
    }

    if (voxel_color.a == 0) discard;

    float dist_squared = dot(disp, disp);
    vec3 fog_color = fog_color_from_world_direction(disp);
    out_color = fog_border_color(
        voxel_color.rgb, dist_squared, uv, fog_color);
}

// Given the residue coordinate of the current fragment (position
// relative to currently-drawn chunk group's origin), raycast through
// the current chunk (positioned by aabb_low/aabb_high; also residue
// coordinates), and return
//
// voxel_color: color of first visible voxel hit. a=0 if no hit.
//
// disp: vector from camera to the hit voxel. Undefined if voxel_color.a=0
//
// uv: [0,1] x [0,1] coordinate showing, if we hit the face of a visible
// voxel, where on that face the raycast intersection was.
void raycast(
    vec3 frag_residue_coord,
    out vec4 voxel_color,
    out vec3 disp,
    out vec2 uv)
{
    // Compute the parametric equation for the ray being cast: v0 + m*t.
    // At the actual coordinate of the fragment when t=0.
    vec3 v0 = frag_residue_coord;
    vec3 m = v0 - eye_relative_group_origin;
    if (m.x == 0) m.x = 1e-7;
    if (m.y == 0) m.y = 1e-7;
    if (m.z == 0) m.z = 1e-7;
    vec3 m_rcp = vec3(1.0) / m;

    // Which way the ray is travelling on the x/y/z axis. +1 positive,
    // -1 negative.
    ivec3 isign = ivec3(sign(m));

    // Like isign but +1 for positive, 0 for negative.
    ivec3 i1_or_0 = (isign + abs(isign)) >> 1;

    // Which voxel (in residue coordinates) the ray is "currently" in.
    // This is moved along the ray until a hit is found or we exit the AABB.
    // The aabb_half_normal (vectors of 1/2 magnitude pointing outwards
    // of the AABB box) is needed to safely sink the coordinate into the
    // interior of a voxel before taking the dangerous floor.
    ivec3 current_residue =
        ivec3(floor(frag_residue_coord - aabb_half_normal));
    float t = 0;

    // AABB gives half-open range of possible visible voxels; I want a
    // closed range.
    ivec3 aabb_high_inclusive = aabb_high - ivec3(1, 1, 1);

    // We also need to keep track of which kind of face (constant
    // x/y/z) we passed in order to figure out the uv eventually.
    // Recycle X/Y/Z SHIFT as useful constants.
    // TODO Get this directly from the vertex shader?
    int last_face_passed =
        aabb_half_normal.x != 0 ? X_SHIFT :
        aabb_half_normal.y != 0 ? Y_SHIFT : Z_SHIFT;

    // iter is needed to defend against accidental infinite loop
    // (crahses entire GPU).
    for (int iter = 0; iter < CHUNK_SIZE * 3; ++iter) {
        // Read voxel here and return if visible voxel hit.
        voxel_color = read_group_voxel(current_residue);
        if (voxel_color.a > 0) {
            vec3 position = v0 + m * t;
            vec3 uv_candidates = position - floor(position);
            disp = position - eye_relative_group_origin;

            switch (last_face_passed) {
                case X_SHIFT: uv = uv_candidates.yz; break;
                case Y_SHIFT: uv = uv_candidates.xz; break;
                case Z_SHIFT: uv = uv_candidates.xy; break;
                default: voxel_color = vec4(1, .5, 0, 0); break;
            }

            return;
        }

        // Otherwise, figure out which "back face" of the voxel this
        // ray exits from. This tells us which voxel the ray enters next.
        // v = v0 + m*t
        // Given target x coord, t = (x - x0) * reciprocal(m). Same for y/z

        // t when ray passes plane defined by back x/y/z face.
        vec3 ts = (current_residue + i1_or_0 - v0) * m_rcp;

        // Exited x back face. (Passed its plane first).
        if (ts.x < ts.y && ts.x < ts.z) {
            current_residue.x += isign.x;
            t = ts.x;
            last_face_passed = X_SHIFT;
        }
        // Exited y back face.
        else if (ts.y < ts.z) {
            current_residue.y += isign.y;
            t = ts.y;
            last_face_passed = Y_SHIFT;
        }
        // Exited z back face.
        else {
            current_residue.z += isign.z;
            t = ts.z;
            last_face_passed = Z_SHIFT;
        }

        // AABB exit check; report no hit.
        ivec3 clamped = clamp(current_residue, aabb_low, aabb_high_inclusive);
        if (clamped != current_residue) {
            voxel_color = vec4(0, 0, 0, 0);
            return;
        }
    }

    // At this point, we reached to loop limit; return debug color.
    voxel_color = vec4(1, 0, 1, 1);
    disp = vec3(1, 0, 0);
    uv = vec2(0, 0);
}

void old_raycast(
    vec3 frag_residue_coord,
    out vec4 best_color,
    out vec3 disp,
    out vec2 uv)
{
    vec3 init_coord = frag_residue_coord + aabb_half_normal;
    float x0 = eye_relative_group_origin.x;
    float y0 = eye_relative_group_origin.y;
    float z0 = eye_relative_group_origin.z;
    vec3 slope = aabb_residue_coord - eye_relative_group_origin;
    float xm = slope.x;
    float ym = slope.y;
    float zm = slope.z;
    float rcp = 1.0/GROUP_SIZE;
    float best_t = 1.0 / 0.0;
    best_color = vec4(0,0,0,0);
    vec3 best_coord = vec3(0,0,0);
    int iter = 0;

    int x_init = int(xm > 0 ? ceil(init_coord.x)
                            : floor(init_coord.x));
    int x_end = xm > 0 ? aabb_high.x : aabb_low.x;
    int x_step = xm > 0 ? 1 : -1;
    float x_bias = xm > 0 ? .25 : -.25;
    for (int x = x_init; x != x_end; x += x_step) {
        if (iter++ >= 255) { best_color = vec4(1,0,1,1); return; }
        float t = (x - x0) / xm;
        float y = y0 + ym * t;
        float z = z0 + zm * t;
        if (y < aabb_low.y || y >= aabb_high.y) break;
        if (z < aabb_low.z || z >= aabb_high.z) break;
        ivec3 texcoord = ivec3(floor(vec3(x + x_bias, y, z)));
        vec4 lookup_color = read_group_voxel(texcoord);
        if (lookup_color.a > 0 && t > 0) {
            if (best_t > t) {
                best_t = t;
                best_color = lookup_color;
                best_coord = vec3(x,y,z);
                uv.x = y;
                uv.y = z;
            }
            break;
        }
    }
    int y_init = int(ym > 0 ? ceil(init_coord.y)
                            : floor(init_coord.y));
    int y_end = ym > 0 ? aabb_high.y : aabb_low.y;
    int y_step = ym > 0 ? 1 : -1;
    float y_bias = ym > 0 ? .25 : -.25;
    for (int y = y_init; y != y_end; y += y_step) {
        if (iter++ >= 255) { best_color = vec4(1,0,1,1); return; }
        float t = (y - y0) / ym;
        float x = x0 + xm * t;
        float z = z0 + zm * t;
        if (x < aabb_low.x || x >= aabb_high.x) break;
        if (z < aabb_low.z || z >= aabb_high.z) break;
        ivec3 texcoord = ivec3(floor(vec3(x, y + y_bias, z)));
        vec4 lookup_color = read_group_voxel(texcoord);
        if (lookup_color.a > 0 && t > 0) {
            if (best_t > t) {
                best_t = t;
                best_color = lookup_color;
                best_coord = vec3(x,y,z);
                uv.x = x;
                uv.y = z;
            }
            break;
        }
    }
    int z_init = int(zm > 0 ? ceil(init_coord.z)
                            : floor(init_coord.z));
    int z_end = zm > 0 ? aabb_high.z : aabb_low.z;
    int z_step = zm > 0 ? 1 : -1;
    float z_bias = zm > 0 ? .25 : -.25;
    for (int z = z_init; z != z_end; z += z_step) {
        if (iter++ >= 255) { best_color = vec4(1,0,1,1); return; }
        float t = (z - z0) / zm;
        float x = x0 + xm * t;
        float y = y0 + ym * t;
        if (x < aabb_low.x || x >= aabb_high.x) break;
        if (y < aabb_low.y || y >= aabb_high.y) break;
        ivec3 texcoord = ivec3(floor(vec3(x, y, z + z_bias)));
        vec4 lookup_color = read_group_voxel(texcoord);
        if (lookup_color.a > 0 && t > 0) {
            if (best_t > t) {
                best_t = t;
                best_color = lookup_color;
                best_coord = vec3(x,y,z);
                uv.x = x;
                uv.y = y;
            }
            break;
        }
    }
    disp = best_coord - eye_relative_group_origin;
}
