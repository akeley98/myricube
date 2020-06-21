// Copied bad code.
#include "myricube.hh"

#include <random>
#include <stdint.h>
#include <time.h>

#include "app.hh"

static constexpr int Size = 365;
static constexpr int Border = 8;
static constexpr double Load = 0.4;

using namespace myricube;

template <int N, int Border>
class congestion_model
{
    constexpr static uint16_t
        red_category = 1,
        green_category = 2,
        blue_category = 3,
        blank = 0;

    struct tile
    {
        int red_next_ptrdiff = 0;
        int green_next_ptrdiff = 0;
        int blue_next_ptrdiff = 0;
        int x = 0, y = 0, z = 0;
        uint16_t old_category = 0, current_category = 0;
        uint8_t red = 0, green = 0, blue = 0;
    };
    std::array<tile, N*N*6> tile_array;
    uint16_t tick_color = blue_category;

  public:
    congestion_model(double probability, VoxelWorld& world)
    {
        std::mt19937 rng(1337);
        tile* pos_x_face = &tile_array[N*N*0];
        tile* pos_y_face = &tile_array[N*N*1];
        tile* pos_z_face = &tile_array[N*N*2];
        tile* neg_x_face = &tile_array[N*N*3];
        tile* neg_y_face = &tile_array[N*N*4];
        tile* neg_z_face = &tile_array[N*N*5];

        tile* red_next = nullptr;
        tile* green_next = nullptr;
        tile* blue_next = nullptr;

        uint32_t thresh = uint32_t(2147483648.0 * probability);
        auto maybe_fill_tile = [&] (tile& t, uint16_t c0, uint16_t c1)
        {
            assert(t.current_category == blank);
            uint32_t r = rng();
            uint16_t set_category = blank;
            if (r < thresh) {
                t.current_category = c0;
                set_category = c0;
            }
            if (r > ~thresh) {
                t.current_category = c1;
                set_category = c1;
            }
            uint8_t red = 0;
            uint8_t green = 0;
            uint8_t blue = 0;
            if (set_category != blank) {
                auto bump0 = rng() % 72;
                auto bump1 = rng() % 40;
                if (bump0 <= 15 && bump1 <= 8) {
                    switch (set_category) {
                        default: assert(0);
                        break; case red_category:
                            red = 255;
                            green = 136;
                            blue = 160;
                        break; case green_category:
                            red = 255;
                            green = 255;
                            blue = 56;
                        break; case blue_category:
                            red = 0;
                            green = 120;
                            blue = 255;
                    }
                }
                else {
                    switch (set_category) {
                        default: assert(0);
                        break; case red_category:
                            red = 128 + bump0;
                            green = 32;
                            blue = 56 + bump1;
                        break; case green_category:
                            red = 64 + bump1;
                            green = 224 - bump0;
                            blue = 0;
                        break; case blue_category:
                            red = 152 + bump0;
                            green = red;
                            blue = 255 - bump1;
                    }
                }
                world.set(glm::ivec3(t.x, t.y, t.z), Voxel(red, green, blue));
                t.red = red;
                t.green = green;
                t.blue = blue;
            }
        };

        for (int y = 0; y < N; ++y) {
            for (int z = 0; z < N; ++z) {
                tile& t = neg_x_face[y*N + z];

                green_next = z == 0 ?
                            &neg_z_face[0*N + y] : &neg_x_face[y*N + z-1];
                blue_next = y == 0 ?
                            &neg_y_face[0*N + z] : &neg_x_face[(y-1)*N + z];
                t.green_next_ptrdiff = int(green_next - &t);
                t.blue_next_ptrdiff  = int(blue_next  - &t);
                t.x = -1;
                t.y = y;
                t.z = z;
                if (Border <= y and y < N - Border
                and Border <= z and z < N - Border) {
                    maybe_fill_tile(t, green_category, blue_category);
                }
            }
        }

        for (int y = 0; y < N; ++y) {
            for (int z = 0; z < N; ++z) {
                tile& t = pos_x_face[y*N + z];

                green_next = z == N-1 ?
                            &pos_z_face[(N-1)*N + y] : &pos_x_face[y*N + z+1];
                blue_next = y == N-1 ?
                            &pos_y_face[(N-1)*N + z] : &pos_x_face[(y+1)*N + z];
                t.green_next_ptrdiff = int(green_next - &t);
                t.blue_next_ptrdiff  = int(blue_next  - &t);
                t.x = N;
                t.y = y;
                t.z = z;
                if (Border <= y and y < N - Border
                and Border <= z and z < N - Border) {
                    maybe_fill_tile(t, green_category, blue_category);
                }
            }
        }

        for (int x = 0; x < N; ++x) {
            for (int z = 0; z < N; ++z) {
                tile& t = neg_y_face[x*N + z];

                red_next = z == 0 ?
                          &neg_z_face[x*N + 0] : &neg_y_face[x*N + z-1];
                blue_next = x == N-1 ?
                          &pos_x_face[0*N + z] : &neg_y_face[(x+1)*N + z];

                t.red_next_ptrdiff  = int(red_next - &t);
                t.blue_next_ptrdiff = int(blue_next  - &t);
                t.x = x;
                t.y = -1;
                t.z = z;
                if (Border <= x and x < N - Border
                and Border <= z and z < N - Border) {
                    maybe_fill_tile(t, red_category, blue_category);
                }
            }
        }

        for (int x = 0; x < N; ++x) {
            for (int z = 0; z < N; ++z) {
                tile& t = pos_y_face[x*N + z];

                red_next = z == N-1 ?
                          &pos_z_face[x*N + N-1] : &pos_y_face[x*N + z+1];
                blue_next = x == 0 ?
                          &neg_x_face[(N-1)*N + z] : &pos_y_face[(x-1)*N + z];

                t.red_next_ptrdiff  = int(red_next - &t);
                t.blue_next_ptrdiff = int(blue_next  - &t);
                t.x = x;
                t.y = N;
                t.z = z;
                if (Border <= x and x < N - Border
                and Border <= z and z < N - Border) {
                    maybe_fill_tile(t, red_category, blue_category);
                }
            }
        }

        for (int x = 0; x < N; ++x) {
            for (int y = 0; y < N; ++y) {
                tile& t = neg_z_face[x*N + y];

                red_next = y == N-1 ?
                          &pos_y_face[x*N + 0] : &neg_z_face[x*N + y+1];
                green_next = x == N-1 ?
                          &pos_x_face[y*N + 0] : &neg_z_face[(x+1)*N + y];

                t.red_next_ptrdiff  = int(red_next - &t);
                t.green_next_ptrdiff = int(green_next - &t);
                t.x = x;
                t.y = y;
                t.z = -1;
                if (Border <= x and x < N - Border
                and Border <= y and y < N - Border) {
                    maybe_fill_tile(t, red_category, green_category);
                }
            }
        }

        for (int x = 0; x < N; ++x) {
            for (int y = 0; y < N; ++y) {
                tile& t = pos_z_face[x*N + y];

                red_next = y == 0 ?
                          &neg_y_face[x*N + N-1] : &pos_z_face[x*N + y-1];
                green_next = x == 0 ?
                          &neg_x_face[y*N + N-1] : &pos_z_face[(x-1)*N + y];

                t.red_next_ptrdiff  = int(red_next - &t);
                t.green_next_ptrdiff = int(green_next - &t);
                t.x = x;
                t.y = y;
                t.z = N;
                if (Border <= x and x < N - Border
                and Border <= y and y < N - Border) {
                    maybe_fill_tile(t, red_category, green_category);
                }
            }
        }
    }

    template <uint16_t ColorCategory>
    void update_tile(tile& t, VoxelWorld& world) {
        auto relptr = ColorCategory == red_category ? t.red_next_ptrdiff :
                      ColorCategory == green_category ? t.green_next_ptrdiff :
                      ColorCategory == blue_category ? t.blue_next_ptrdiff : 0;
        tile& next = (&t)[relptr];
        assert(&*tile_array.begin() <= &next && &next < &*tile_array.end());
        if (next.old_category == blank && t.old_category == ColorCategory) {
            next.current_category = ColorCategory;
            t.current_category = blank;
            auto red = next.red = t.red;
            auto green = next.green = t.green;
            auto blue = next.blue = t.blue;
            world.set(glm::ivec3(t.x, t.y, t.z), Voxel());
            Voxel v = Voxel(red, green, blue);
            world.set(glm::ivec3(next.x, next.y, next.z), v);
        }
    }

    void update(VoxelWorld& world)
    {
        for (tile& t : tile_array) {
            t.old_category = t.current_category;
        }
        switch (tick_color) {
          default: assert(0);
          break; case red_category:
            tick_color = green_category;
            for (tile& t : tile_array) {
                update_tile<red_category>(t, world);
            }
          break; case green_category:
            tick_color = blue_category;
            for (tile& t : tile_array) {
                update_tile<green_category>(t, world);
            }
          break; case blue_category:
            tick_color = red_category;
            for (tile& t : tile_array) {
                update_tile<blue_category>(t, world);
            }
        }
    }
};

static congestion_model<Size, Border>* ptr_congestion_model = nullptr;

namespace myricube {

void app_init(VoxelWorld& world, Window& window)
{
    ptr_congestion_model = new congestion_model<Size, Border>(Load, world);

    KeyTarget skip_100;
    skip_100.down = [&world] (KeyArg a)
    {
        if (a.repeat) return false;
        for (int i = 0; i < 100; ++i) {
            ptr_congestion_model->update(world);
        }
        return true;
    };
    window.add_key_target("do_it", skip_100);
}

static double previous_update = -1.0/0.0;

void app_update(VoxelWorld& world)
{
    struct timespec ts;
    clock_gettime(CLOCK_BOOTTIME, &ts);
    double ns = ts.tv_sec + 1e-9 * ts.tv_nsec;
    double dt = ns - previous_update;
    constexpr double interval = 0.075;
    if (dt > interval * 1.5) {
        previous_update = ns;
    }
    else if (dt >= interval) {
        previous_update += interval;
    }
    else {
        return;
    }

    ptr_congestion_model->update(world);
}

}
