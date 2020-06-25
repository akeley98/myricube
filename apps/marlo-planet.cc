#include <iostream>
#include <string>
#include <fstream>
#include "app.hh"
#include "FastNoise.cpp"

using namespace myricube;

void marlo(int Radius, VoxelWorld& world)
{
    constexpr double box = 14.0;
    auto relu = [] (double n) -> double
    {
        return std::max(0.0, n);
    };
    FastNoise noise;
    noise.SetNoiseType(FastNoise::Perlin);
    auto noiser = [&noise] (glm::dvec3 v) -> double
    {
        float result = noise.GetNoise(v.x * 200, v.y * 200, v.z * 200);
        return result;
    };
    auto l2norm = [] (glm::dvec3 v) -> double
    {
        return sqrt(glm::dot(v,v));
    };

    struct color_control_point
    {
        glm::vec3 color;
        float radius;
    };

    constexpr int color_count = 8;
    static const color_control_point ccp[color_count] =
    {
        { glm::vec3(1.0, 1.0, 1.0), 0.0f },
        { glm::vec3(1.0, 1.0, 0.0), 7.1f },
        { glm::vec3(1.0, 0.0, 0.0), 7.5f },
        { glm::vec3(0.0, 0.8, 0.2), 8.5f },
        { glm::vec3(0.0, 0.0, 1.0), 11.5f },
        { glm::vec3(0.0, 1.0, 1.0), 12.5f },
        { glm::vec3(1.0, 0.0, 1.0), 12.9f },
        { glm::vec3(1.0, 1.0, 0.0), 100.0f },
    };
    auto voxel_from_norm = [] (float n) -> Voxel
    {
        for (int i = 0; i < color_count-1; ++i) {
            if (ccp[i].radius <= n && n < ccp[i+1].radius) {
                float segment_length = ccp[i+1].radius-ccp[i].radius;
                float interpolant = (n - ccp[i].radius) / segment_length;
                glm::vec3 float_color = glm::mix(
                    ccp[i].color, ccp[i+1].color, interpolant);
                uint8_t red(float_color.x * 255);
                uint8_t green(float_color.g * 255);
                uint8_t blue(float_color.b * 255);
                return Voxel(red, green, blue);
            }
        }
        return Voxel(255, 255, 255);
    };
    
    for (int i = -Radius; i <= +Radius; ++i) {
        for (int j = -Radius; j <= +Radius; ++j) {
            for (int k = -Radius; k <= +Radius; ++k) {
                auto pos = glm::dvec3(i, j, k) * double(box / Radius);
                auto norm = l2norm(pos);
                auto shell1 = relu(-std::abs(norm - 8) + 0.5);
                auto shell2 = relu(-std::abs(norm - 12) + 0.5);
                auto ret = shell1 + shell2;
                ret += relu((noiser(pos*0.5) + 1)*.8
                     - std::abs(2-std::abs(norm-10)))*(std::abs(norm-10)<2.0);
                auto noise1 = noiser(-1020.0 + pos*0.2);
                ret += (3.0-std::abs(norm-10.2)) <= 0 ? 0 : relu(-0.06 + noiser(pos*0.194));
                ret *= relu(l2norm((pos - glm::dvec3(12.2,0,0)) * glm::dvec3(2,1,1))
                            -3.2*(noise1+1) ) * relu(l2norm((pos-glm::dvec3(-8,0,0)) * glm::dvec3(2,1,1)) - 4*(noise1+1));
                if (ret > 0) {
                    world.set(glm::ivec3(i, j, k), voxel_from_norm(norm));
                }
            }
        }
    }
    string line;
    int j = 0;
    while (getline (file, line) {
        for(std::string::size_type i = 0; i < str.size(); ++i) {
            if (str[i]!=' '){
	        world.set(glm:ivec3(i,j,0), voxel_from_norm(100));  
	    }
        }
     	j +=1;	
    }
}

namespace myricube {

VoxelWorld* p_world;
std::mt19937 rng;

void app_init(VoxelWorld& world, Window& window)
{
    marlo(300, world);
}

void app_update(VoxelWorld& world)
{

}

} // end namespace
