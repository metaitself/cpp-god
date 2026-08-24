#include <iostream>
#include <cmath>
#include <string>

static constexpr double M = 1.0;
static constexpr double RS = 2.0 * M;
static constexpr double ISCO = 6.0 * M;
static constexpr double DISK_MAX = 12.0 * M;
static constexpr double DT = 0.05;

struct vec3 {
    double x, y, z;
    vec3 operator+(vec3 o) const { return {x + o.x, y + o.y, z + o.z}; }
    vec3 operator-(vec3 o) const { return {x - o.x, y - o.y, z - o.z}; }
    vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
};

static inline double dot(vec3 a, vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline double len(vec3 v) { return std::sqrt(dot(v, v)); }
static inline vec3 norm(vec3 v) { double l = len(v); return l > 0 ? v * (1.0 / l) : v; }
static inline vec3 cross(vec3 a, vec3 b) {
    return { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x };
}

static inline vec3 accel(vec3 p, vec3 v) {
    double r = len(p);
    if (r <= RS) return {0, 0, 0};
    vec3 L = cross(p, v);
    double h2 = dot(L, L);
    return p * (-1.5 * RS * h2 / (r * r * r * r * r));
}

static void rk4(vec3 &p, vec3 &v, double dt) {
    vec3 p1 = v;
    vec3 v1 = accel(p, v);

    vec3 p2 = v + v1 * (dt * 0.5);
    vec3 v2 = accel(p + p1 * (dt * 0.5), p2);

    vec3 p3 = v + v2 * (dt * 0.5);
    vec3 v3 = accel(p + p2 * (dt * 0.5), p3);

    vec3 p4 = v + v3 * dt;
    vec3 v4 = accel(p + p3 * dt, p4);

    p = p + (p1 + p2 * 2.0 + p3 * 2.0 + p4) * (dt / 6.0);
    v = v + (v1 + v2 * 2.0 + v3 * 2.0 + v4) * (dt / 6.0);
}

static double trace(vec3 ro, vec3 rd) {
    vec3 p = ro;
    vec3 v = rd;

    for (int i = 0; i < 1500; ++i) {
        double r = len(p);
        if (r <= RS * 1.02) return 0.0;

        vec3 old_p = p;
        rk4(p, v, DT);

        if (old_p.z * p.z <= 0.0) {
            double rd_disk = std::sqrt(p.x * p.x + p.y * p.y);
            if (rd_disk >= ISCO && rd_disk <= DISK_MAX) {
                return std::pow(ISCO / rd_disk, 1.5);
            }
        }

        if (r > 35.0 && dot(p, v) > 0.0) {
            return 0.05;
        }
    }
    return 0.0;
}

int main() {
    const int W = 80;
    const int H = 40;
    const double fov = 1.05;

    vec3 cam = {0.0, -18.0, 4.2};
    vec3 target = {0.0, 0.0, 0.0};
    vec3 fwd = norm(target - cam);
    vec3 right = norm(cross(fwd, {0.0, 0.0, 1.0}));
    vec3 up = cross(right, fwd);

    const char* ramp = " .:-=+*#%@";

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            double u = ((double)x / W - 0.5) * fov * ((double)W / H);
            double v = (0.5 - (double)y / H) * fov;

            vec3 dir = norm(fwd + right * u + up * v);
            double val = trace(cam, dir);

            if (val <= 0.001) {
                std::cout << ' ';
            } else {
                int idx = (int)(val * 9.0);
                if (idx > 9) idx = 9;
                std::cout << ramp[idx];
            }
        }
        std::cout << '\n';
    }

    return 0;
}
