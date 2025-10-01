#include <raylib.h>
#include <raymath.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <nanosvg.h>

#define TRAIL_POINT_COUNT 1000
#define VECTOR_SCALE 1
#define COEFFICIENTS 500

void DrawVector(Vector2 vec, Vector2 origin, float width) {
    DrawLineEx(origin, Vector2Add(origin, vec), width, (Color){ 207, 159, 255, 255 });
}

float lerp(float a, float b, float t) {
    return (1 - t) * a + t * b;
}

// only defined between 0 and 2
float *path;
int npts;
complex f(float x) {
    // return cos(2*PI*x) + I * sin(2*PI*x);
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

    float real = lerp(xr0, xr1, t) - 800;
    float imag = lerp(yr0, yr1, t) - 450;

    return real + I*imag;
}

int main() {
    InitWindow(1600, 900, "Fourier Series Visualization");

    NSVGimage *fourier_svg = nsvgParseFromFile("fourier.svg", "px", 1);
    path = fourier_svg->shapes[0].paths[0].pts;
    npts = fourier_svg->shapes[0].paths[0].npts;
    Vector2 *points = (Vector2 *)malloc(sizeof(Vector2) * npts);
    for (int i = 0; i < npts * 2 - 1; i+=2) {
        points[i/2].x = path[i];
        points[i/2].y = path[i+1];
    }

    complex coefficients[COEFFICIENTS * 2];
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

    for (int i = 0; i < COEFFICIENTS * 2; i++) {
        printf("k%d: %f %fi\n", i - COEFFICIENTS, creal(coefficients[i]), cimag(coefficients[i]));
    }

    Vector2 origin = { 800, 450 };

    Vector2 trail_points[TRAIL_POINT_COUNT];
    int trail_opacities[TRAIL_POINT_COUNT];
    memset(trail_points, 0, sizeof(trail_points));
    memset(trail_opacities, 0, sizeof(trail_opacities));

    int current_trail_point = 0;

    int printed = 0;

    float t = 0;
    float trail_accumulator = 0;
    float fade_accumulator = 0;
    while (!WindowShouldClose()) {
        float dt = GetFrameTime();
        t += dt / 20;
        trail_accumulator += dt;
        fade_accumulator += dt;

        BeginDrawing();

        ClearBackground(BLACK);

        for (int i = 0; i < TRAIL_POINT_COUNT; i++) {
            DrawCircleV(trail_points[i], 5, (Color){ 255, 0, 0, trail_opacities[i] });
        }

        Vector2 current_origin = origin;
        complex phasor = coefficients[COEFFICIENTS];
        float x = creal(phasor) * VECTOR_SCALE;
        float y = cimag(phasor) * VECTOR_SCALE;
        Vector2 vec = { x, y };
        float width = 10;
        DrawVector(vec, current_origin, width);
        float radius = sqrtf(x*x + y*y);
        DrawRing(current_origin, radius, radius - 5, 0, 360, 100, BLUE);
        current_origin = Vector2Add(vec, current_origin);

        for (int k = 1; k <= COEFFICIENTS; k++) {
            int i_positive = k + COEFFICIENTS;
            int i_negative = -k + COEFFICIENTS;
            width *= 0.9f;
            
            phasor = coefficients[i_positive] * (cos(PI * k * t) + I * sin(PI * k * t));
            x = creal(phasor) * VECTOR_SCALE;
            y = cimag(phasor) * VECTOR_SCALE;
            Vector2 vec1 = { x, y };
            DrawVector(vec1, current_origin, width);
            radius = sqrtf(x*x + y*y);
            DrawRing(current_origin, radius, radius - 5, 0, 360, 100, BLUE);
            current_origin = Vector2Add(vec1, current_origin);

            phasor = coefficients[i_negative] * (cos(PI * -k * t) + I * sin(PI * -k * t));
            x = creal(phasor) * VECTOR_SCALE;
            y = cimag(phasor) * VECTOR_SCALE;
            Vector2 vec2 = { x, y };
            DrawVector(vec2, current_origin, width);
            radius = sqrtf(x*x + y*y);
            DrawRing(current_origin, radius, radius - 5, 0, 360, 100, BLUE);
            current_origin = Vector2Add(vec2, current_origin);
        }

        if (trail_accumulator >= 0.01f) {
            trail_points[current_trail_point] = current_origin;
            trail_opacities[current_trail_point] = 255;
            current_trail_point++;
            if (current_trail_point >= TRAIL_POINT_COUNT) current_trail_point = 0;
            trail_accumulator = 0;
        }

        if (fade_accumulator >= 0.00784313725) {
            for (int i = 0; i < TRAIL_POINT_COUNT; i++) {
                trail_opacities[i] -= 1;
                if (trail_opacities[i] < 0) trail_opacities[i] = 0;
            }
            fade_accumulator = 0;
        }

        // DrawSplineBezierCubic(points, npts, 3, WHITE);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
