#include <raylib.h>
#include <raymath.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <nanosvg.h>

#define TRAIL_POINT_COUNT 1000
#define VECTOR_SCALE 500
#define COEFFICIENTS 10

void DrawVector(Vector2 vec, Vector2 origin) {
    DrawLineEx(origin, Vector2Add(origin, vec), 10, (Color){ 207, 159, 255, 255 });
}

// only defined between 0 and 2
float f(float x) {
    if (x <= 1) return 1;
    return 0;
}

int main() {
    InitWindow(1600, 900, "Fourier Series Visualization");

    NSVGimage *fourier_svg = nsvgParseFromFile("fourier.svg", "px", 1);
    float *path = fourier_svg->shapes[0].paths[0].pts;
    int npts = fourier_svg->shapes[0].paths[0].npts;
    Vector2 *points = (Vector2 *)malloc(sizeof(Vector2) * npts);
    for (int i = 0; i < npts * 2 - 1; i+=2) {
        points[i/2].x = path[i];
        points[i/2].y = path[i+1];
    }

    printf("%f %f\n", points[npts * 3-3].x, points[npts * 3-3].y);
    printf("%f %f\n", points[npts * 3-2].x, points[npts * 3-2].y);
    printf("%f %f\n", points[npts * 3-1].x, points[npts * 3-1].y);

    complex coefficients[COEFFICIENTS * 2];
    int a = 0;
    int b = 2;
    int n = 1000;
    for (int i = -COEFFICIENTS; i <= COEFFICIENTS; i++) {
        // Calculate real part
        float real = (f(a)*cos(PI*i*a) + f(b)*cos(PI*i*b)) / 2;
        for (int k = 1; k <= n - 1; k++) {
            real += f(a + k * (b - a) / (float)n) * cos(PI*i*(a + k * (b - a) / (float)n));
        }
        real *= (b - a) / (float)n;
        real /= 2;

        // Calculate complex part
        float comp = (f(a)*sin(PI*i*a) + f(b)*sin(PI*i*b)) / 2;
        for (int k = 1; k <= n - 1; k++) {
            comp += f(a + k * (b - a) / (float)n)*sin(PI*i*(a + k * (b - a) / (float)n));
        }
        comp *= (b - a) / (float)n;
        comp /= 2;

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
        t += dt / 5;
        trail_accumulator += dt;
        fade_accumulator += dt;

        BeginDrawing();

        ClearBackground(BLACK);

        // DrawVector(end_point, origin);

        for (int i = 0; i < TRAIL_POINT_COUNT; i++) {
            DrawCircleV(trail_points[i], 5, (Color){ 255, 0, 0, trail_opacities[i] });
        }

        Vector2 current_origin = origin;
        complex phasor = coefficients[COEFFICIENTS];
        float x = creal(phasor) * VECTOR_SCALE;
        float y = cimag(phasor) * VECTOR_SCALE;
        Vector2 vec = { x, y };
        DrawVector(vec, current_origin);
        float radius = sqrtf(x*x + y*y);
        DrawRing(current_origin, radius, radius - 5, 0, 360, 100, BLUE);
        current_origin = Vector2Add(vec, current_origin);

        for (int k = 1; k <= COEFFICIENTS; k++) {
            int i_positive = k + COEFFICIENTS;
            int i_negative = -k + COEFFICIENTS;
            
            phasor = coefficients[i_positive] * (cos(PI * k * t) + I * sin(PI * k * t));
            x = creal(phasor) * VECTOR_SCALE;
            y = cimag(phasor) * VECTOR_SCALE;
            Vector2 vec1 = { x, y };
            DrawVector(vec1, current_origin);
            radius = sqrtf(x*x + y*y);
            DrawRing(current_origin, radius, radius - 5, 0, 360, 100, BLUE);
            current_origin = Vector2Add(vec1, current_origin);

            phasor = coefficients[i_negative] * (cos(PI * -k * t) + I * sin(PI * -k * t));
            x = creal(phasor) * VECTOR_SCALE;
            y = cimag(phasor) * VECTOR_SCALE;
            Vector2 vec2 = { x, y };
            DrawVector(vec2, current_origin);
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

        DrawSplineBezierCubic(points, npts, 3, WHITE);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
