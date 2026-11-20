/*
 * blackhole.c: renders a Kerr black hole with a thin accretion disk
 * to a PNG. The physics lives in physics.h; live.c is the interactive viewer.
 *
 * Build: make blackhole
 * Run:   ./blackhole            (renders blackhole.png)
 *        ./blackhole --check    (compares the integrator against known results)
 *        ./blackhole --help
 */
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "physics.h"
#include "png.h"

#define STILL_STEP 0.02

typedef struct {
    int width, height, ss, threads;
    double exposure;
    const char *out;
    View view;
} Settings;

static Settings S = {
    .width = 1920, .height = 1080, .ss = 3, .threads = 0, .exposure = 1.0,
    .out = "blackhole.png",
    .view = {.r_cam = 80.0, .incl_deg = 84.0, .azim_deg = 0.0, .fov_deg = 20.0,
             .t_peak = 9000.0, .r_out = 30.0, .spin = 0.9},
};

/* ---------- checks against known results ---------- */

/* Ray in the equatorial plane at radius r0, heading inward, with l = L/E. */
static Ray equatorial_ray(double r0, double l) {
    Ray ray = {.l = l, .eta = 0, .e_obs = 1};
    double u = 1 / r0;
    ray.s = (Geo){u, sqrt(fmax(u_potential(u, l, 0), 0)), 0, 0, 0};
    return ray;
}

/* Bisects along screen direction psi for the shadow edge; returns the edge ray. */
static Ray edge_ray(const Camera *c, double psi) {
    double lo = 0, hi = 1;
    for (int k = 0; k < 50; k++) {
        double mid = 0.5 * (lo + hi);
        vec3 n = normalize(add(c->fwd, add(mul(c->right, mid * cos(psi)), mul(c->up, mid * sin(psi)))));
        Ray ray = camera_ray_dir(c, n);
        if (trace(&ray, 0, R_ESCAPE, STILL_STEP, STILL_STEP).kind == HIT_HORIZON) lo = mid; else hi = mid;
    }
    vec3 n = normalize(add(c->fwd, add(mul(c->right, hi * cos(psi)), mul(c->up, hi * sin(psi)))));
    return camera_ray_dir(c, n);
}

static int check(const char *name, double got, double want, double tol) {
    double err = fabs(got - want) / fmax(fabs(want), 1e-300);
    printf("%-34s %14.9f  expected %14.9f  rel err %.1e  %s\n", name, got, want, err, err < tol ? "ok" : "FAIL");
    return err < tol;
}

static int run_checks(void) {
    int ok = 1;
    bb_init();
    stars_init();

    /* Schwarzschild limit */
    spacetime_init(0);
    double lo = 4.5, hi = 6.0;
    for (int i = 0; i < 48; i++) {
        double b = 0.5 * (lo + hi);
        Ray ray = equatorial_ray(1e4, b);
        if (trace(&ray, 0, 2e4, STILL_STEP, STILL_STEP).kind == HIT_HORIZON) lo = b; else hi = b;
    }
    ok &= check("a=0 critical impact parameter", 0.5 * (lo + hi), sqrt(27.0), 1e-6);

    double b = 100, R0 = 1e7;
    Ray ray = equatorial_ray(R0, b);
    Hit h = trace(&ray, 0, R0, STILL_STEP, STILL_STEP);
    double defl = fabs(h.end.phi) - (M_PI - asin(b / R0) - asin(b * h.end.u));
    double series = 4 / b + 15 * M_PI / (4 * b * b) + 128 / (3 * b * b * b) + 3465 * M_PI / (64 * pow(b, 4)) +
                    3584 / (5 * pow(b, 5));
    ok &= check("a=0 deflection at b=100", defl, series, 1e-6);

    disk_init(1, 60);
    ok &= check("a=0 disk flux peak radius", nt_rpeak, 9.5525, 2e-3);
    ok &= check("a=0 radiative efficiency", disk_efficiency(), 1 - sqrt(8.0 / 9.0), 1e-6);

    /* Kerr: an equatorial observer sees the shadow's left and right edges at the
       l of the co- and counter-rotating circular photon orbits (Bardeen 1973). */
    for (int k = 0; k < 2; k++) {
        double a = k == 0 ? 0.5 : 0.9;
        View v = {.r_cam = 50, .incl_deg = 90, .fov_deg = 30, .t_peak = 4500, .r_out = 30, .spin = a};
        physics_sync(&v);
        Camera c = camera_make(v, 100, 100, STILL_STEP, STILL_STEP);
        double l1 = edge_ray(&c, 0).l, l2 = edge_ray(&c, M_PI).l;
        double co = fmax(l1, l2), counter = fmin(l1, l2);
        char name[64];
        snprintf(name, sizeof name, "a=%.1f shadow edge, co-rotating", a);
        ok &= check(name, co, photon_orbit_l(K.r_ph_min), 1e-5);
        snprintf(name, sizeof name, "a=%.1f shadow edge, counter-rotating", a);
        ok &= check(name, counter, photon_orbit_l(K.r_ph_max), 1e-5);

        /* the ISCO formula sits at the minimum of the circular-orbit energy */
        double best = 0, emin = 1e9;
        for (double r = K.r_ph_min + 0.05; r < 12; r += 1e-5) {
            double e = orbit_E(r);
            if (e < emin) { emin = e; best = r; }
        }
        snprintf(name, sizeof name, "a=%.1f ISCO vs min of E(r)", a);
        ok &= check(name, best, K.r_isco, 1e-4);

        /* a ray skimming the photon shell keeps (du/dlambda)^2 = U(u) */
        Ray edge = edge_ray(&c, M_PI / 3);
        Hit e = trace(&edge, 0, R_ESCAPE, STILL_STEP, STILL_STEP);
        double drift = fabs(e.end.vu * e.end.vu - u_potential(e.end.u, edge.l, edge.eta));
        snprintf(name, sizeof name, "a=%.1f constraint drift (edge ray)", a);
        printf("%-34s %14.3e  limit 1e-8  %s\n", name, drift, drift < 1e-8 ? "ok" : "FAIL");
        ok &= drift < 1e-8;
    }

    /* Lensing of a distant star: seen from distance D, a star at angle beta from
       the hole has two images that solve beta = theta - alpha(b), b = D sin(theta),
       with the bending angle alpha = 4M/b + (15 pi / 4) M^2/b^2 + (128/3) M^3/b^3. To first order
       this is the Einstein point lens, theta_E^2 = 4M/D. */
    {
        View v = {.r_cam = 1e5, .incl_deg = 84, .fov_deg = 1, .t_peak = 4500, .r_out = 30, .spin = 0};
        physics_sync(&v);
        Camera cam = camera_make(v, 100, 100, STILL_STEP, STILL_STEP);
        double theta_e = sqrt(4 / v.r_cam), beta = 2 * theta_e;
        for (int img = 0; img < 2; img++) {
            /* t: signed screen offset along cam.right. Both brackets see the star
               short of beta at lo and past it at hi. */
            double lo = img == 0 ? tan(beta) : -tan(theta_e), hi = img == 0 ? tan(3 * beta) : -tan(0.3 * theta_e);
            for (int k = 0; k < 60; k++) {
                double mid = 0.5 * (lo + hi);
                Ray ray = camera_ray_dir(&cam, normalize(add(cam.fwd, mul(cam.right, mid))));
                Hit h = trace(&ray, 0, 1e9, STILL_STEP, STILL_STEP);
                double seen = atan2(dot(h.dir, cam.right), dot(h.dir, cam.fwd));
                if (seen < beta) lo = mid; else hi = mid;
            }
            double theta = atan(0.5 * (lo + hi));
            double wlo = img == 0 ? beta : -theta_e, whi = img == 0 ? 3 * beta : -0.3 * theta_e;
            for (int k = 0; k < 100; k++) {
                double w = 0.5 * (wlo + whi), b = v.r_cam * fabs(sin(w));
                double bend = 4 / b + 15 * M_PI / (4 * b * b) + 128 / (3 * b * b * b);
                if (w - (w > 0 ? bend : -bend) < beta) wlo = w; else whi = w;
            }
            ok &= check(img == 0 ? "lensed star, outer image" : "lensed star, inner image", theta,
                        0.5 * (wlo + whi), 1e-5);
        }
    }

    /* Sky orientation: Orion as seen from Earth, not its mirror image. Betelgeuse
       is north-east of Rigel, so with north up it must appear up and to the left. */
    {
        View v = S.view;
        physics_sync(&v);
        Camera cam = camera_make(v, 100, 100, STILL_STEP, STILL_STEP);
        vec3 bet = v3(0, 0, 0), rig = v3(0, 0, 0);
        for (int i = 0; i < BSC5_COUNT; i++) {
            if (BSC5[i].hr == 2061) bet = radec_to_scene(BSC5[i].ra, BSC5[i].dec);
            if (BSC5[i].hr == 1713) rig = radec_to_scene(BSC5[i].ra, BSC5[i].dec);
        }
        double bx = dot(bet, cam.right) / dot(bet, cam.fwd), by = dot(bet, cam.up) / dot(bet, cam.fwd);
        double rx = dot(rig, cam.right) / dot(rig, cam.fwd), ry = dot(rig, cam.up) / dot(rig, cam.fwd);
        int sky_ok = bx < rx && by > ry;
        printf("%-34s Betelgeuse (%+.3f, %+.3f) vs Rigel (%+.3f, %+.3f)  %s\n", "Orion orientation", bx, by, rx, ry,
               sky_ok ? "ok" : "FAIL");
        ok &= sky_ok;
        ok &= check("sky aim: Orion behind the hole", dot(radec_to_scene(SKY_RA0, SKY_DEC0), cam.fwd), 1.0, 1e-12);
    }

    vec3 c = xyz_to_srgb(bb_xyz(6504));
    double m = fmax(c.x, fmax(c.y, c.z));
    int white = c.x / m > 0.9 && c.y / m > 0.9 && c.z / m > 0.9;
    printf("%-34s %.3f %.3f %.3f  %s\n", "6504 K blackbody in sRGB", c.x / m, c.y / m, c.z / m, white ? "ok" : "FAIL");
    ok &= white;

    printf(ok ? "all checks passed\n" : "CHECKS FAILED\n");
    return ok ? 0 : 1;
}

/* ---------- main ---------- */

static void usage(void) {
    printf("usage: blackhole [options]\n"
           "  --size WxH      image size (default 1920x1080)\n"
           "  --ss N          N x N samples per pixel (default 3)\n"
           "  --dist R        camera distance in M (default 80)\n"
           "  --incl DEG      angle between spin axis and line of sight (default 84)\n"
           "  --azim DEG      camera rotation around the spin axis (default 0)\n"
           "  --fov DEG       vertical field of view (default 20)\n"
           "  --temp K        peak disk temperature (default 9000)\n"
           "  --rout R        disk outer radius in M (default 30)\n"
           "  --spin A        black hole spin a/M, -0.998 to 0.998 (default 0.9)\n"
           "  --exposure X    exposure multiplier (default 1)\n"
           "  --threads N     worker threads (default: all cores)\n"
           "  --out FILE      output PNG (default blackhole.png)\n"
           "  --check         verify the integrator and color model, then exit\n");
}

int main(int argc, char **argv) {
    View *v = &S.view;
    for (int i = 1; i < argc; i += 2) {
        const char *a = argv[i], *val = i + 1 < argc ? argv[i + 1] : NULL;
        if (!strcmp(a, "--check")) return run_checks();
        if (!strcmp(a, "--help") || !strcmp(a, "-h")) { usage(); return 0; }
        if (!val) { usage(); return 1; }
        if (!strcmp(a, "--size")) sscanf(val, "%dx%d", &S.width, &S.height);
        else if (!strcmp(a, "--ss")) S.ss = atoi(val);
        else if (!strcmp(a, "--dist")) v->r_cam = atof(val);
        else if (!strcmp(a, "--incl")) v->incl_deg = atof(val);
        else if (!strcmp(a, "--azim")) v->azim_deg = atof(val);
        else if (!strcmp(a, "--fov")) v->fov_deg = atof(val);
        else if (!strcmp(a, "--temp")) v->t_peak = atof(val);
        else if (!strcmp(a, "--rout")) v->r_out = atof(val);
        else if (!strcmp(a, "--spin")) v->spin = atof(val);
        else if (!strcmp(a, "--exposure")) S.exposure = atof(val);
        else if (!strcmp(a, "--threads")) S.threads = atoi(val);
        else if (!strcmp(a, "--out")) S.out = val;
        else { usage(); return 1; }
    }
    if (fabs(v->spin) > 0.998 || v->r_cam <= 3.5 || S.width <= 0 || S.height <= 0 || S.ss <= 0) {
        fprintf(stderr, "need |spin| <= 0.998, camera distance > 3.5 M and positive sizes\n");
        return 1;
    }

    bb_init();
    stars_init();
    physics_sync(v);
    Camera cam = camera_make(*v, S.width, S.height, STILL_STEP, STILL_STEP);

    int threads = S.threads > 0 ? S.threads : (int)sysconf(_SC_NPROCESSORS_ONLN);
    if (threads < 1) threads = 1;
    size_t pixels = (size_t)S.width * S.height;
    float *image = malloc(sizeof(float) * 3 * pixels);

    fprintf(stderr, "rendering %dx%d, %d spp, %d threads\n", S.width, S.height, S.ss * S.ss, threads);
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    render_rows(&cam, image, NULL, S.width, S.height, S.ss, 0, S.height, threads);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    fprintf(stderr, "done in %.1f s\n", (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9);

    double exposure = S.exposure * EXPOSURE;
    uint8_t *out = malloc(3 * pixels);
    for (size_t i = 0; i < pixels; i++) {
        double c[3];
        tonemap(image + 3 * i, exposure, c);
        for (int k = 0; k < 3; k++) out[3 * i + k] = (uint8_t)lround(255 * c[k]);
    }
    if (!write_png(S.out, out, S.width, S.height)) {
        fprintf(stderr, "could not write %s\n", S.out);
        return 1;
    }
    fprintf(stderr, "wrote %s\n", S.out);
    return 0;
}
