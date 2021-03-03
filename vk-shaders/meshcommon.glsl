#ifndef MYRICUBE_MESHCOMMON_GLSL_
#define MYRICUBE_MESHCOMMON_GLSL_

// Bit assignments for unpacking packed voxels.
#define POS_X_FACE_BIT (1 << 24)
#define NEG_X_FACE_BIT (1 << 25)
#define POS_Y_FACE_BIT (1 << 26)
#define NEG_Y_FACE_BIT (1 << 27)
#define POS_Z_FACE_BIT (1 << 28)
#define NEG_Z_FACE_BIT (1 << 29)

#define X_SHIFT 0
#define Y_SHIFT 8
#define Z_SHIFT 16

#define RED_SHIFT 0
#define GREEN_SHIFT 8
#define BLUE_SHIFT 16

#ifndef __cplusplus
// Unit cube model, used in shaders only.

// Position coordinate of unit box:
vec3 unit_box_verts[36] = vec3[] (
    vec3(0, 1, 1),    vec3(0, 1, 0),    vec3(0, 0, 1),
    vec3(0, 1, 0),    vec3(0, 0, 0),    vec3(0, 0, 1), // -x face

    vec3(1, 0, 0),    vec3(1, 1, 0),    vec3(1, 0, 1),
    vec3(1, 1, 0),    vec3(1, 1, 1),    vec3(1, 0, 1), // +x face

    vec3(0, 0, 0),    vec3(1, 0, 1),    vec3(0, 0, 1),
    vec3(0, 0, 0),    vec3(1, 0, 0),    vec3(1, 0, 1), // -y face

    vec3(1, 1, 1),    vec3(1, 1, 0),    vec3(0, 1, 0),
    vec3(0, 1, 1),    vec3(1, 1, 1),    vec3(0, 1, 0), // +y face

    vec3(0, 0, 0),    vec3(0, 1, 0),    vec3(1, 0, 0),
    vec3(1, 1, 0),    vec3(1, 0, 0),    vec3(0, 1, 0), // -z face

    vec3(0, 1, 1),    vec3(1, 0, 1),    vec3(1, 1, 1),
    vec3(1, 0, 1),    vec3(0, 1, 1),    vec3(0, 0, 1));// +z face

// Bit corresponding to the face being drawn of the unit box; for
// example, if we're drawing the +x face, this is POS_X_FACE_BIT.
int unit_box_face_bits[36] = int[] (
    NEG_X_FACE_BIT,    NEG_X_FACE_BIT,    NEG_X_FACE_BIT,
    NEG_X_FACE_BIT,    NEG_X_FACE_BIT,    NEG_X_FACE_BIT,

    POS_X_FACE_BIT,    POS_X_FACE_BIT,    POS_X_FACE_BIT,
    POS_X_FACE_BIT,    POS_X_FACE_BIT,    POS_X_FACE_BIT,

    NEG_Y_FACE_BIT,    NEG_Y_FACE_BIT,    NEG_Y_FACE_BIT,
    NEG_Y_FACE_BIT,    NEG_Y_FACE_BIT,    NEG_Y_FACE_BIT,

    POS_Y_FACE_BIT,    POS_Y_FACE_BIT,    POS_Y_FACE_BIT,
    POS_Y_FACE_BIT,    POS_Y_FACE_BIT,    POS_Y_FACE_BIT,

    NEG_Z_FACE_BIT,    NEG_Z_FACE_BIT,    NEG_Z_FACE_BIT,
    NEG_Z_FACE_BIT,    NEG_Z_FACE_BIT,    NEG_Z_FACE_BIT,

    POS_Z_FACE_BIT,    POS_Z_FACE_BIT,    POS_Z_FACE_BIT,
    POS_Z_FACE_BIT,    POS_Z_FACE_BIT,    POS_Z_FACE_BIT);

// Texture coordinate (for now, just used for the color border).
vec2 uv_array[36] = vec2[] (
    vec2(1, 1),    vec2(1, 0),    vec2(0, 1),
    vec2(1, 0),    vec2(0, 0),    vec2(0, 1), // -x face

    vec2(0, 0),    vec2(1, 0),    vec2(0, 1),
    vec2(1, 0),    vec2(1, 1),    vec2(0, 1), // +x face

    vec2(0, 0),    vec2(1, 1),    vec2(0, 1),
    vec2(0, 0),    vec2(1, 0),    vec2(1, 1), // -y face

    vec2(1, 1),    vec2(1, 0),    vec2(0, 0),
    vec2(0, 1),    vec2(1, 1),    vec2(0, 0), // +y face

    vec2(0, 0),    vec2(0, 1),    vec2(1, 0),
    vec2(1, 1),    vec2(1, 0),    vec2(0, 1), // -z face

    vec2(0, 1),    vec2(1, 0),    vec2(1, 1),
    vec2(1, 0),    vec2(0, 1),    vec2(0, 0));// +z face
#endif

#endif
