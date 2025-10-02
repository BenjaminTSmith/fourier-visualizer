#include <raylib.h>
#include <raymath.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <nanosvg.h>

#define TRAIL_POINT_COUNT 1000
#define DRAW_TIME 20
#define VECTOR_SCALE 1
#define COEFFICIENTS 250

#define TRANS_GRAY ((Color){ 0xD3, 0xD3, 0xDE, 90 })

void DrawVector(Vector2 vec, Vector2 origin, float width) {
    float length = Vector2Length(vec);
    float tip_length = length * 0.1;
    vec = Vector2Scale(vec, 0.9);
    Vector2 tip = Vector2Scale(Vector2Normalize(vec), tip_length);
    Vector2 v1 = Vector2Add(Vector2Add(vec, tip), origin);
    Vector2 v2 = Vector2Add(Vector2Add(vec, Vector2Rotate(tip, -2.472)), origin);
    Vector2 v3 = Vector2Add(Vector2Add(vec, Vector2Rotate(tip, 2.472)), origin);
    DrawLineEx(origin, Vector2Add(origin, vec), width, WHITE);
    DrawTriangle(v1, v2, v3, WHITE);
}

float lerp(float a, float b, float t) {
    return (1 - t) * a + t * b;
}

// only defined between 0 and 2
float *path;
int npts;
int width;
int height;
float scale_factor;
float x_offset;
float y_offset;
complex f(float x) {
    float n = (npts - 1) * x / 3;
    float t = n - (int)n;
    int i = (int)n * 6;
    
    float x0, y0, cpx1, cpy1, cpx2, cpy2, x1, y1;
    x0 = path[i];
    y0 = path[i+1];
    cpx1 = path[i+2];
    cpy1 = path[i+3];
    cpx2 = path[i+4];
    cpy2 = path[i+5];
    x1 = path[i+6];
    y1 = path[i+7];

    float xq0 = lerp(x0, cpx1, t);
    float yq0 = lerp(y0, cpy1, t);

    float xq1 = lerp(cpx1, cpx2, t);
    float yq1 = lerp(cpy1, cpy2, t);

    float xq2 = lerp(cpx2, x1, t);
    float yq2 = lerp(cpy2, y1, t);

    float xr0 = lerp(xq0, xq1, t);
    float yr0 = lerp(yq0, yq1, t);

    float xr1 = lerp(xq1, xq2, t);
    float yr1 = lerp(yq1, yq2, t);

    float real = lerp(xr0, xr1, t) * scale_factor - 800 + x_offset;
    float imag = lerp(yr0, yr1, t) * scale_factor - 450 + y_offset;

    return real + I*imag;
}

int main(int argc, char **argv) {
    InitWindow(1600, 900, "Fourier Series Visualization");

    NSVGimage *svg;
    if (argc > 1) {
        svg = nsvgParseFromFile(argv[1], "px", 1);
    } else {
        svg = nsvgParseFromFile("fourier.svg", "px", 1);
    }

    width = svg->width;
    height = svg->height;
    if (width/1600.f > height/900.f) {
        scale_factor = 1600.f / width * 0.9f;
    } else {
        scale_factor = 900.f / height * 0.9f;
    }
    x_offset = (1600 - scale_factor * width) / 2;
    y_offset = (900 - scale_factor * height) / 2;

    NSVGpath *paths = svg->shapes->paths;
    npts = 0;
    while (paths != NULL) {
        npts += paths->npts;
        paths = paths->next;
    }
    int j = 0;
    path = (float *)malloc(sizeof(float) * npts * 2);

    paths = svg->shapes->paths;
    while (paths != NULL) {
        for (int i = 0; i < paths->npts * 2; i++) {
            path[j] = paths->pts[i]; 
            j++;
        }
        paths = paths->next;
    }

    Vector2 *points = (Vector2 *)malloc(sizeof(Vector2) * npts);
    for (int i = 0; i < npts * 2 - 1; i+=2) {
        points[i/2].x = path[i];
        points[i/2].y = path[i+1];
    }

    complex coefficients[COEFFICIENTS * 2 + 1];
    int a = 0;
    int b = 1;
    int n = 10000;
    for (int i = -COEFFICIENTS; i <= COEFFICIENTS; i++) {
        // Calculate real part
        complex real = (f(a)*cos(PI*i*a) + f(b)*cos(PI*i*b)) / 2;
        for (int k = 1; k <= n - 1; k++) {
            real += f(a + k * (b - a) / (float)n) * cos(2*PI*i*(a + k * (b - a) / (float)n));
        }
        real *= (b - a) / (float)n;

        // Calculate complex part
        complex comp = (f(a)*sin(PI*i*a) + f(b)*sin(PI*i*b)) / 2;
        for (int k = 1; k <= n - 1; k++) {
            comp += f(a + k * (b - a) / (float)n) * sin(2*PI*i*(a + k * (b - a) / (float)n));
        }
        comp *= (b - a) / (float)n;

        coefficients[i + COEFFICIENTS] = real  - I*comp;
    }

    for (int i = 0; i < COEFFICIENTS * 2 + 1; i++) {
        printf("k%d: %f %fi\n", i - COEFFICIENTS, creal(coefficients[i]), cimag(coefficients[i]));
    }

    Vector2 origin = { 800, 450 };

    // Calculate initial point
    Vector2 initial_origin = origin;
    complex phasor = coefficients[COEFFICIENTS];
    float x = creal(phasor) * VECTOR_SCALE;
    float y = cimag(phasor) * VECTOR_SCALE;
    Vector2 vec = { x, y };
    initial_origin = Vector2Add(vec, initial_origin);

    for (int k = 1; k <= COEFFICIENTS; k++) {
        int i_positive = k + COEFFICIENTS;
        int i_negative = -k + COEFFICIENTS;
        phasor = coefficients[i_positive];
        x = creal(phasor) * VECTOR_SCALE;
        y = cimag(phasor) * VECTOR_SCALE;
        Vector2 vec1 = { x, y };
        initial_origin = Vector2Add(vec1, initial_origin);
        phasor = coefficients[i_negative];
        x = creal(phasor) * VECTOR_SCALE;
        y = cimag(phasor) * VECTOR_SCALE;
        Vector2 vec2 = { x, y };
        initial_origin = Vector2Add(vec2, initial_origin);
    }

    Vector2 trail_points[TRAIL_POINT_COUNT];
    for (int i = 0; i < TRAIL_POINT_COUNT; i++) {
        trail_points[i] = initial_origin;
    }

    int current_trail_point = 0;

    int printed = 0;

    float t = 0;
    float trail_accumulator = 0;
    while (!WindowShouldClose()) {
        float dt = GetFrameTime();
        t += dt / DRAW_TIME;
        trail_accumulator += dt / DRAW_TIME;

        BeginDrawing();

        ClearBackground(BLACK);

        Vector2 current_origin = origin;
        complex phasor = coefficients[COEFFICIENTS];
        float x = creal(phasor) * VECTOR_SCALE;
        float y = cimag(phasor) * VECTOR_SCALE;
        Vector2 vec = { x, y };
        float vector_width = 4;
        float ring_width = 1.5;
        // Do I want to draw the vector with f = 0?
        // DrawVector(vec, current_origin, vector_width);
        float radius = sqrtf(x*x + y*y);
        // DrawRing(current_origin, radius, radius - ring_width, 0, 360, 100, TRANS_GRAY);
        current_origin = Vector2Add(vec, current_origin);

        for (int k = 1; k <= COEFFICIENTS; k++) {
            int i_positive = k + COEFFICIENTS;
            int i_negative = -k + COEFFICIENTS;
            vector_width *= 0.9f;
            // ring_width *= 0.8f;
            if (vector_width <= 0.5f) vector_width = 0.5f;
            // if (ring_width <= 0.5f) ring_width = 0.5f;

            phasor = coefficients[i_positive] * (cos(2 * PI * k * t) + I * sin(2 * PI * k * t));
            x = creal(phasor) * VECTOR_SCALE;
            y = cimag(phasor) * VECTOR_SCALE;
            Vector2 vec1 = { x, y };
            DrawVector(vec1, current_origin, vector_width);
            radius = sqrtf(x*x + y*y);
            DrawRing(current_origin, radius, radius - ring_width, 0, 360, 100, TRANS_GRAY);
            current_origin = Vector2Add(vec1, current_origin);

            phasor = coefficients[i_negative] * (cos(2 * PI * -k * t) + I * sin(2 * PI * -k * t));
            x = creal(phasor) * VECTOR_SCALE;
            y = cimag(phasor) * VECTOR_SCALE;
            Vector2 vec2 = { x, y };
            DrawVector(vec2, current_origin, vector_width);
            radius = sqrtf(x*x + y*y);
            DrawRing(current_origin, radius, radius - ring_width, 0, 360, 100, TRANS_GRAY);
            current_origin = Vector2Add(vec2, current_origin);
        }

        while (trail_accumulator >= (float)1 / TRAIL_POINT_COUNT) {
            trail_points[current_trail_point] = current_origin;
            current_trail_point++;
            if (current_trail_point >= TRAIL_POINT_COUNT) current_trail_point = 0;
            trail_accumulator -= (float)1 / TRAIL_POINT_COUNT;
        }

        /*for (int i = 0; i < TRAIL_POINT_COUNT; i++) {
            DrawCircleV(trail_points[i], 5, (Color){ 255, 0, 0, 255 });
        }*/

        int i = current_trail_point;
        int last_i = i - 1;
        if (last_i < 0) last_i = TRAIL_POINT_COUNT - 1;
        while (i != last_i) {
            int next_i = i + 1;
            if (next_i >= TRAIL_POINT_COUNT) next_i = 0;

            int distance_from_start;
            if (i > last_i) {
                distance_from_start = i - current_trail_point;
            } else {
                distance_from_start = i + (TRAIL_POINT_COUNT - current_trail_point);
            }
            int alpha = (float)distance_from_start / TRAIL_POINT_COUNT * 255;

            DrawLineEx(trail_points[i], trail_points[next_i], 3, (Color){ 253, 249, 0, alpha });
            i = next_i;
        }

        // DrawSplineBezierCubic(points, npts, 4, WHITE);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
