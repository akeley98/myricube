#ifndef MYRICUBE_SRGB_GLSL_
#define MYRICUBE_SRGB_GLSL_

// https://github.com/nvpro-samples/vk_compute_mipmaps/blob/main/shaders/srgb.h

uint srgb8_from_linear_bias(float arg, float bias)
{
    float srgb = arg <= 0.0031308F ? (323.0F / 25.0F) * arg :
                                     1.055F * pow(arg, 1.0F / 2.4F) - 0.055F;
    return uint(clamp(srgb * 255.0F + bias, 0.F, 255.F));
}

// Convert float linear red/green/blue value to 8-bit sRGB component.
uint srgb8_from_linear(float arg)
{
    return srgb8_from_linear_bias(arg, 0.5);
}

float linear_from_srgb8(uint arg)
{
    arg = min(arg, 255u);
    float u = float(arg) * (1.0F / 255.0F);
    return u <= 0.04045F ? u * (25.0F / 323.0F) :
                           pow((200.0F * u + 11.0F) * (1.0F / 211.0F), 2.4F);

}
#endif
