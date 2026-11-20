/*
 * physics.h: Kerr ray tracing, Novikov-Thorne disk, blackbody color.
 * Shared by the still renderer (blackhole.c) and the live viewer (live.c).
 *
 * Units: G = c = M = 1. Spin a in (-1, 1); a = 0 is Schwarzschild.
 *
 * Each pixel traces a light ray backward from a camera through Kerr spacetime.
 * Rays are integrated in Boyer-Lindquist coordinates with Mino time, where
 * Carter's separation gives two decoupled equations with polynomial
 * right-hand sides. In u = 1/r and mu = cos(theta):
 *     u'' = U'(u) / 2      U(u) = u^4 R(1/u) = (1 + c u^2)^2 - K u^2 (1 - 2u + a^2 u^2)
 *     mu'' = M'(mu) / 2    M(mu) = eta - (eta + l^2 - a^2) mu^2 - a^2 mu^4
 * with c = a^2 - a l, K = eta + (l - a)^2, l = L/E and eta = Q/E^2 (Carter's
 * constant). The second-order form passes turning points without sign
 * bookkeeping, and u keeps the far field well conditioned: at a = 0 the u
 * equation is the classic u'' + u = 3u^2 (scaled).
 *
 * The disk is optically thick and radiates as a local blackbody with the
 * Novikov-Thorne (Page-Thorne) flux profile. Because I_nu / nu^3 is invariant,
 * a blackbody at temperature T seen with frequency ratio g looks exactly like a
 * blackbody at g*T. g includes gravitational redshift, Doppler shift from the
 * orbiting gas, and the camera's own frame.
 */
#ifndef PHYSICS_H
#define PHYSICS_H

#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define R_ESCAPE 1.0e4

/* ---------- vectors ---------- */

typedef struct { double x, y, z; } vec3;

static inline vec3 v3(double x, double y, double z) { return (vec3){x, y, z}; }
static inline vec3 add(vec3 a, vec3 b) { return v3(a.x + b.x, a.y + b.y, a.z + b.z); }
static inline vec3 sub(vec3 a, vec3 b) { return v3(a.x - b.x, a.y - b.y, a.z - b.z); }
static inline vec3 mul(vec3 a, double s) { return v3(a.x * s, a.y * s, a.z * s); }
static inline double dot(vec3 a, vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline double len(vec3 a) { return sqrt(dot(a, a)); }
static inline vec3 normalize(vec3 a) { return mul(a, 1.0 / len(a)); }
static inline vec3 cross(vec3 a, vec3 b) {
    return v3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

/* ---------- spacetime ---------- */

typedef struct {
    double a;                  /* spin; negative means the disk orbits against the spin */
    double r_plus;             /* outer event horizon */
    double r_isco;             /* innermost stable circular orbit of the disk */
    double r_ph_min, r_ph_max; /* equatorial circular photon orbits, co- and counter-rotating */
    double r_capture;          /* inside this radius no ray from outside can turn around */
} Spacetime;

static Spacetime K;

static void spacetime_init(double a) {
    double s = fabs(a);
    double z1 = 1 + cbrt(1 - s * s) * (cbrt(1 + s) + cbrt(1 - s));
    double z2 = sqrt(3 * s * s + z1 * z1);
    K.a = a;
    K.r_plus = 1 + sqrt(1 - a * a);
    K.r_isco = 3 + z2 + (a >= 0 ? -1 : 1) * sqrt((3 - z1) * (3 + z1 + 2 * z2)); /* Bardeen+ 1972 */
    K.r_ph_min = 2 * (1 + cos(2.0 / 3.0 * acos(-s)));
    K.r_ph_max = 2 * (1 + cos(2.0 / 3.0 * acos(s)));
    /* Turning points of rays from outside lie beyond the photon region, which
       starts at r_ph_min, so a ray halfway between it and the horizon is lost. */
    K.r_capture = 0.5 * (K.r_plus + K.r_ph_min);
}

/* Circular equatorial orbits in the disk's sense of rotation (Bardeen, Press &
   Teukolsky 1972): angular velocity, u^t, energy and angular momentum. */
static double orbit_den(double r) { return pow(r, 0.75) * sqrt(r * sqrt(r) - 3 * sqrt(r) + 2 * K.a); }
static double orbit_omega(double r) { return 1 / (r * sqrt(r) + K.a); }
static double orbit_ut(double r) { return (r * sqrt(r) + K.a) / orbit_den(r); }
static double orbit_E(double r) { return (r * sqrt(r) - 2 * sqrt(r) + K.a) / orbit_den(r); }
static double orbit_L(double r) { return (r * r - 2 * K.a * sqrt(r) + K.a * K.a) / orbit_den(r); }

/* l = L/E of the circular photon orbit at r, from the photon-region formulas
   (Bardeen 1973). The shadow's equatorial edges have exactly these l. */
static inline double photon_orbit_l(double r) {
    double a = K.a;
    return -(r * r * r - 3 * r * r + a * a * r + a * a) / (a * (r - 1));
}

/* ---------- geodesics ---------- */

/* u = 1/r, mu = cos(theta), phi and the Mino-time rates of u and mu.
   Rays are traced backward from the camera, so phi runs backward too. */
typedef struct { double u, vu, mu, vmu, phi; } Geo;

typedef struct {
    Geo s;
    double l, eta; /* L/E and Q/E^2 */
    double e_obs;  /* E of the photon when the camera measures energy 1 */
} Ray;

/* R(r): the ray can only be where R >= 0; its zeros are radial turning points. */
static inline double radial_potential(double r, double l, double eta) {
    double a = K.a, P = r * r + a * a - a * l, delta = r * r - 2 * r + a * a;
    return P * P - delta * (eta + (l - a) * (l - a));
}

static inline double u_potential(double u, double l, double eta) {
    double a = K.a, c = a * a - a * l, Kc = eta + (l - a) * (l - a), q = 1 + c * u * u;
    return q * q - Kc * u * u * (1 - 2 * u + a * a * u * u);
}

static inline Geo geo_rates(const Geo *s, double l, double eta) {
    double a = K.a, a2 = a * a, u = s->u, mu = s->mu;
    double c = a2 - a * l, Kc = eta + (l - a) * (l - a);
    double sin2 = fmax(1 - mu * mu, 1e-12);
    /* phi rate a P / Delta - a + l / sin^2(theta), with P/Delta written in u */
    double P_over_delta = (1 + c * u * u) / (1 - 2 * u + a2 * u * u);
    Geo d;
    d.u = s->vu;
    d.vu = 2 * c * u * (1 + c * u * u) - Kc * u * (1 - 3 * u + 2 * a2 * u * u); /* U'(u) / 2 */
    d.mu = s->vmu;
    d.vmu = -(eta + l * l - a2) * mu - 2 * a2 * mu * mu * mu;                     /* M'(mu) / 2 */
    d.phi = -(a * P_over_delta - a + l / sin2);                                   /* backward in time */
    return d;
}

static inline Geo geo_axpy(const Geo *s, const Geo *d, double h) {
    return (Geo){s->u + d->u * h, s->vu + d->vu * h, s->mu + d->mu * h, s->vmu + d->vmu * h, s->phi + d->phi * h};
}

static Geo rk4(const Geo *s, double l, double eta, double h) {
    Geo k1 = geo_rates(s, l, eta);
    Geo t = geo_axpy(s, &k1, h / 2);
    Geo k2 = geo_rates(&t, l, eta);
    t = geo_axpy(s, &k2, h / 2);
    Geo k3 = geo_rates(&t, l, eta);
    t = geo_axpy(s, &k3, h);
    Geo k4 = geo_rates(&t, l, eta);
    Geo o = *s;
    o.u += h / 6 * (k1.u + 2 * k2.u + 2 * k3.u + k4.u);
    o.vu += h / 6 * (k1.vu + 2 * k2.vu + 2 * k3.vu + k4.vu);
    o.mu += h / 6 * (k1.mu + 2 * k2.mu + 2 * k3.mu + k4.mu);
    o.vmu += h / 6 * (k1.vmu + 2 * k2.vmu + 2 * k3.vmu + k4.vmu);
    o.phi += h / 6 * (k1.phi + 2 * k2.phi + 2 * k3.phi + k4.phi);
    return o;
}

/* Scene coordinates: y is the spin axis and the disk lies in y = 0. Boyer-Lindquist
   Cartesian (X, Y, Z) = (x, -z, y), so phi grows in the disk's rotation sense. */
static inline vec3 bl_to_scene(double r, double mu, double phi) {
    double rho = sqrt(r * r + K.a * K.a), s = sqrt(fmax(0, 1 - mu * mu));
    return v3(rho * s * cos(phi), r * mu, -rho * s * sin(phi));
}

/* Direction of motion in scene coordinates, from the Mino-time rates. */
static vec3 geo_direction(const Geo *g, const Geo *d) {
    double r = 1 / g->u, vr = -g->vu * r * r;
    double rho = sqrt(r * r + K.a * K.a), s = sqrt(fmax(1e-18, 1 - g->mu * g->mu));
    double drho = r * vr / rho, ds = -g->mu * g->vmu / s, c = cos(g->phi), sn = sin(g->phi);
    double radial = drho * s + rho * ds;
    double X = radial * c - rho * s * sn * d->phi, Y = radial * sn + rho * s * c * d->phi;
    double Z = vr * g->mu + r * g->vmu;
    return normalize(v3(X, Z, -Y));
}

enum { HIT_HORIZON, HIT_DISK, HIT_SKY };

/* order counts how often the ray crossed the disk plane outside the disk before
   it ended: 0 is the direct image, 1 the secondary image, 2+ the photon ring. */
typedef struct { int kind, order; vec3 pos, dir; Geo end; } Hit;

/* Optional record of a ray's positions, for plotting a single geodesic. */
typedef struct { vec3 *pts; int n, cap; } Path;

/* The disk lies in the equatorial plane between K.r_isco and r_disk_out (0
   disables it). step bounds the change per RK4 step: about step * r in r and
   step in mu and phi. Beyond r = 20 rays are nearly straight, so a live preview
   can let the step grow with r up to step_max; stills pass step_max = step. */
static Hit trace_path(const Ray *ray, double r_disk_out, double r_escape, double step, double step_max,
                      Path *path) {
    Geo s = ray->s;
    double l = ray->l, eta = ray->eta, u_capture = 1 / K.r_capture, u_escape = 1 / r_escape;
    int order = 0;
    for (int i = 0; i < 400000; i++) {
        if (path && path->n < path->cap) path->pts[path->n++] = bl_to_scene(1 / s.u, s.mu, s.phi);
        if (s.u > u_capture) return (Hit){.kind = HIT_HORIZON, .order = order, .end = s};
        Geo d = geo_rates(&s, l, eta);
        if (s.u < u_escape && s.vu < 0)
            return (Hit){.kind = HIT_SKY, .order = order, .pos = bl_to_scene(1 / s.u, s.mu, s.phi),
                         .dir = geo_direction(&s, &d), .end = s};
        double rate = fabs(s.vu) / s.u + fabs(s.vmu) + fabs(d.phi) + 1e-12;
        double h = fmin(step * fmax(1.0, 1 / (20.0 * s.u)), step_max) / rate;
        Geo n = rk4(&s, l, eta, h);
        if (r_disk_out > 0 && (s.mu > 0) != (n.mu > 0)) {
            /* A linear estimate of the crossing skips the bisection when the ray
               clearly misses the disk; otherwise bisect the RK4 step. */
            double t = s.mu / (s.mu - n.mu), r_est = 1 / (s.u + (n.u - s.u) * t);
            double margin = 0.1 * r_est + 0.5;
            if (r_est < K.r_isco - margin || r_est > r_disk_out + margin) {
                s = n;
                order++;
                continue;
            }
            double lo = 0, hi = h;
            Geo c = n;
            for (int k = 0; k < 24; k++) {
                double mid = 0.5 * (lo + hi);
                c = rk4(&s, l, eta, mid);
                if ((c.mu > 0) == (s.mu > 0)) lo = mid; else hi = mid;
            }
            double rc = 1 / c.u;
            if (rc >= K.r_isco && rc <= r_disk_out) {
                if (path && path->n < path->cap) path->pts[path->n++] = bl_to_scene(rc, 0, c.phi);
                return (Hit){.kind = HIT_DISK, .order = order, .pos = bl_to_scene(rc, 0, c.phi), .end = c};
            }
            order++;
        }
        s = n;
    }
    return (Hit){.kind = HIT_HORIZON, .order = order, .end = s};
}

static Hit trace(const Ray *ray, double r_disk_out, double r_escape, double step, double step_max) {
    return trace_path(ray, r_disk_out, r_escape, step, step_max, NULL);
}

/* ---------- blackbody color ---------- */

/* CIE 1931 2-degree color matching functions, multi-lobe fit from
   Wyman, Sloan & Shirley, "Simple Analytic Approximations to the CIE XYZ
   Color Matching Functions" (JCGT 2013). */
static double lobe(double x, double mu, double s1, double s2) {
    double t = (x - mu) / (x < mu ? s1 : s2);
    return exp(-0.5 * t * t);
}
static vec3 cmf(double nm) {
    return v3(1.056 * lobe(nm, 599.8, 37.9, 31.0) + 0.362 * lobe(nm, 442.0, 16.0, 26.7) -
                  0.065 * lobe(nm, 501.1, 20.4, 26.2),
              0.821 * lobe(nm, 568.8, 46.9, 40.5) + 0.286 * lobe(nm, 530.9, 16.3, 31.1),
              1.217 * lobe(nm, 437.0, 11.8, 36.0) + 0.681 * lobe(nm, 459.0, 26.0, 13.8));
}

#define BB_N 4096
#define BB_TMIN 300.0
#define BB_TMAX 300000.0
static vec3 bb_table[BB_N];

/* XYZ of Planck's law at temperature T (arbitrary absolute scale). */
static vec3 planck_xyz(double T) {
    const double c2 = 1.438776877e-2; /* hc/k in m*K */
    vec3 acc = v3(0, 0, 0);
    for (double nm = 360; nm <= 830; nm += 1.0) {
        double lam = nm * 1e-9;
        double b = 1.0 / (pow(lam * 1e6, 5) * expm1(c2 / (lam * T)));
        acc = add(acc, mul(cmf(nm), b));
    }
    return acc;
}

static void bb_init(void) {
    for (int i = 0; i < BB_N; i++) {
        double T = BB_TMIN * pow(BB_TMAX / BB_TMIN, i / (double)(BB_N - 1));
        bb_table[i] = planck_xyz(T);
    }
}

static vec3 bb_xyz(double T) {
    if (T <= BB_TMIN) return v3(0, 0, 0);
    double f = log(T / BB_TMIN) / log(BB_TMAX / BB_TMIN) * (BB_N - 1);
    if (f >= BB_N - 1) return mul(bb_table[BB_N - 1], T / BB_TMAX); /* Rayleigh-Jeans */
    int i = (int)f;
    double t = f - i;
    return add(mul(bb_table[i], 1 - t), mul(bb_table[i + 1], t));
}

static vec3 xyz_to_srgb(vec3 c) {
    return v3(3.2404542 * c.x - 1.5371385 * c.y - 0.4985314 * c.z,
              -0.9692660 * c.x + 1.8760108 * c.y + 0.0415560 * c.z,
              0.0556434 * c.x - 0.2040259 * c.y + 1.0572252 * c.z);
}

/* ---------- Novikov-Thorne disk ---------- */

/* Page & Thorne (1974) flux for circular geodesics:
     F(r) = -Mdot/(4 pi r) * Omega' / (E - Omega L)^2 * Int_isco^r (E - Omega L) L' dr
   evaluated numerically for the current spin. */
#define NT_N 8192
static double nt_rmin, nt_rmax, nt_rpeak, nt_temp[NT_N];

static double nt_deriv(double (*f)(double), double r) {
    double h = 1e-5 * r;
    return (f(r + h) - f(r - h)) / (2 * h);
}

static void disk_init(double t_peak, double r_out) {
    static double flux[NT_N];
    nt_rmin = K.r_isco;
    nt_rmax = fmax(r_out, 60.0);
    double dr = (nt_rmax - nt_rmin) / (NT_N - 1), integral = 0, prev = 0, fmax_ = 0;
    int imax = 0;
    for (int i = 0; i < NT_N; i++) {
        double r = nt_rmin + i * dr;
        double e_ol = orbit_E(r) - orbit_omega(r) * orbit_L(r);
        double integrand = i == 0 ? 0 : e_ol * nt_deriv(orbit_L, r);
        if (i > 0) integral += 0.5 * (prev + integrand) * dr;
        prev = integrand;
        flux[i] = -nt_deriv(orbit_omega, r) / (r * e_ol * e_ol) * integral;
        if (flux[i] > fmax_) { fmax_ = flux[i]; imax = i; }
    }
    nt_rpeak = nt_rmin + imax * dr;
    for (int i = 0; i < NT_N; i++) nt_temp[i] = t_peak * pow(fmax(flux[i], 0) / fmax_, 0.25);
}

static double disk_temp(double r) {
    double f = (r - nt_rmin) / (nt_rmax - nt_rmin) * (NT_N - 1);
    if (f <= 0) return 0;
    if (f >= NT_N - 1) return nt_temp[NT_N - 1];
    int i = (int)f;
    double t = f - i;
    return nt_temp[i] * (1 - t) + nt_temp[i + 1] * t;
}

/* Fraction of rest mass radiated by gas spiralling down to the ISCO. */
static double disk_efficiency(void) { return 1 - orbit_E(K.r_isco); }

/* ---------- star field: the real sky ---------- */

#include "bsc5.h"

/* The sky is the Yale Bright Star Catalog: 9096 stars down to V = 6.5 at their
   J2000 positions, magnitudes and B-V colors. It is turned so that Orion sits
   behind the hole for the default camera (inclination 84°, azimuth 0) with
   celestial north up; orbiting the camera then moves across the real sky. */
#define SKY_RA0 (5.583 * M_PI / 12)  /* 5h35m, between Orion's belt and sword */
#define SKY_DEC0 (-1.5 * M_PI / 180)
#define SKY_INCL0 (84.0 * M_PI / 180)
#define STAR_SCALE 4.0 /* peak brightness of a V = 0 star, relative to the disk peak */

/* Equatorial basis at the aim point (east, north, center) and where it lands in
   the scene. Both are right-handed, so the map is a rotation, not a mirror. */
static struct { vec3 e_eq, n_eq, c_eq, e_sc, n_sc, c_sc; } SKY;

static void sky_init(void) {
    double a = SKY_RA0, d = SKY_DEC0;
    SKY.c_eq = v3(cos(d) * cos(a), cos(d) * sin(a), sin(d));
    SKY.n_eq = v3(-sin(d) * cos(a), -sin(d) * sin(a), cos(d));
    SKY.e_eq = v3(-sin(a), cos(a), 0);
    /* default camera: forward, and screen up; east is screen left */
    SKY.c_sc = v3(0, -cos(SKY_INCL0), sin(SKY_INCL0));
    vec3 right = normalize(cross(SKY.c_sc, v3(0, 1, 0)));
    SKY.n_sc = cross(right, SKY.c_sc);
    SKY.e_sc = mul(right, -1);
}

static vec3 radec_to_scene(double ra, double dec) {
    vec3 q = v3(cos(dec) * cos(ra), cos(dec) * sin(ra), sin(dec));
    return add(add(mul(SKY.e_sc, dot(SKY.e_eq, q)), mul(SKY.n_sc, dot(SKY.n_eq, q))),
               mul(SKY.c_sc, dot(SKY.c_eq, q)));
}

/* Right ascension (0 to 2 pi) and declination of a scene direction. */
static void scene_to_radec(vec3 d, double *ra, double *dec) {
    vec3 q = add(add(mul(SKY.e_eq, dot(SKY.e_sc, d)), mul(SKY.n_eq, dot(SKY.n_sc, d))),
                 mul(SKY.c_eq, dot(SKY.c_sc, d)));
    *ra = fmod(atan2(q.y, q.x) + 2 * M_PI, 2 * M_PI);
    *dec = asin(fmax(-1, fmin(1, q.z)));
}

/* Stars are bucketed in a 3D grid over the unit sphere, so a lookup only
   checks the 27 cells around the ray's final direction. */
#define SKY_N 48
#define SKY_CELLS (2 * SKY_N) * (2 * SKY_N) * (2 * SKY_N)
typedef struct { vec3 dir, xyz; double temp, bright; int hr; } Star;
static Star stars[BSC5_COUNT];
static int *cell_start;

static int cell_index(int ix, int iy, int iz) {
    if (ix < 0 || iy < 0 || iz < 0 || ix >= 2 * SKY_N || iy >= 2 * SKY_N || iz >= 2 * SKY_N) return -1;
    return (ix * 2 * SKY_N + iy) * 2 * SKY_N + iz;
}

static int cell_of(vec3 d) {
    return cell_index((int)floor(d.x * SKY_N) + SKY_N, (int)floor(d.y * SKY_N) + SKY_N,
                      (int)floor(d.z * SKY_N) + SKY_N);
}

/* Effective temperature from B-V (Ballesteros 2012). */
static double bv_temperature(double bv) {
    bv = fmax(-0.4, fmin(2.0, bv));
    return 4600 * (1 / (0.92 * bv + 1.7) + 1 / (0.92 * bv + 0.62));
}

static void stars_init(void) {
    sky_init();
    static Star raw[BSC5_COUNT];
    cell_start = calloc(SKY_CELLS + 1, sizeof(int));
    for (int i = 0; i < BSC5_COUNT; i++) {
        raw[i].dir = radec_to_scene(BSC5[i].ra, BSC5[i].dec);
        raw[i].temp = bv_temperature(BSC5[i].bv);
        raw[i].bright = STAR_SCALE * pow(10, -0.4 * BSC5[i].vmag);
        raw[i].hr = BSC5[i].hr;
        cell_start[cell_of(raw[i].dir) + 1]++;
    }
    for (int c = 0; c < SKY_CELLS; c++) cell_start[c + 1] += cell_start[c];
    int *fill = calloc(SKY_CELLS, sizeof(int));
    for (int i = 0; i < BSC5_COUNT; i++) {
        int c = cell_of(raw[i].dir);
        stars[cell_start[c] + fill[c]++] = raw[i];
    }
    free(fill);
}

/* Star colors as seen by the camera: g_cam blueshifts starlight falling down to it. */
static void stars_tint(double g_cam) {
    static double last = -1;
    if (g_cam == last) return;
    last = g_cam;
    for (int i = 0; i < BSC5_COUNT; i++)
        stars[i].xyz = mul(bb_xyz(stars[i].temp * g_cam), stars[i].bright / bb_xyz(stars[i].temp).y);
}

/* Star image: a sharp core plus a faint halo three times wider, like the
   point-spread function of a real camera, so bright stars read larger than faint
   ones. sigma is the core width in radians; gain keeps a star's total light the
   same when sigma changes with resolution. Lensing acts on these source-plane
   shapes, so magnified images stretch into arcs. */
#define HALO_WIDTH 3.0
#define HALO_LEVEL 0.04
static vec3 sky(vec3 d, double sigma, double gain) {
    double sh = HALO_WIDTH * sigma, cutoff = 9 * sh * sh;
    double inv = 1 / (2 * sigma * sigma), inv_h = 1 / (2 * sh * sh);
    vec3 out = v3(0, 0, 0);
    int ix = (int)floor(d.x * SKY_N) + SKY_N, iy = (int)floor(d.y * SKY_N) + SKY_N,
        iz = (int)floor(d.z * SKY_N) + SKY_N;
    for (int a = -1; a <= 1; a++)
        for (int b = -1; b <= 1; b++)
            for (int c = -1; c <= 1; c++) {
                int cell = cell_index(ix + a, iy + b, iz + c);
                if (cell < 0) continue;
                for (int k = cell_start[cell]; k < cell_start[cell + 1]; k++) {
                    vec3 e = sub(stars[k].dir, d);
                    double d2 = dot(e, e);
                    if (d2 < cutoff)
                        out = add(out, mul(stars[k].xyz, exp(-d2 * inv) + HALO_LEVEL * exp(-d2 * inv_h)));
                }
            }
    return mul(out, gain);
}

/* ---------- camera ---------- */

typedef struct {
    double r_cam;    /* distance from the hole, in M */
    double incl_deg; /* angle between the spin axis and the line of sight */
    double azim_deg; /* rotation around the spin axis */
    double fov_deg;  /* vertical field of view */
    double t_peak;   /* peak disk temperature, K */
    double r_out;    /* disk outer radius, in M */
    double spin;     /* a / M */
} View;

/* Brings the spacetime and disk tables in line with v. Call from one thread. */
static void physics_sync(const View *v) {
    static double spin = NAN, t_peak = NAN, r_out = NAN;
    if (v->spin != spin) spacetime_init(v->spin);
    if (v->spin != spin || v->t_peak != t_peak || v->r_out != r_out) disk_init(v->t_peak, v->r_out);
    spin = v->spin;
    t_peak = v->t_peak;
    r_out = v->r_out;
}

/* The camera is a zero-angular-momentum observer (ZAMO). Its local axes line up
   with the Boyer-Lindquist r, theta and phi directions. */
typedef struct {
    vec3 pos, fwd, right, up;      /* scene */
    vec3 e_r, e_th, e_ph;          /* local axes, Boyer-Lindquist Cartesian */
    double r, mu, phi;             /* position */
    double lapse, omega;           /* sqrt(Sigma Delta / A) and frame-dragging rate */
    double k_ph, k_th, k_r, delta; /* sqrt(A/Sigma) sin(theta), sqrt(Sigma), sqrt(Sigma/Delta) */
    double tan_half, aspect, star_sigma, star_gain, y_norm, r_out, step, step_max;
} Camera;

static Camera camera_make(View v, int width, int height, double step, double step_max) {
    Camera c;
    double incl = v.incl_deg * M_PI / 180, azim = v.azim_deg * M_PI / 180, a = K.a;
    c.pos = v3(v.r_cam * sin(incl) * sin(azim), v.r_cam * cos(incl), -v.r_cam * sin(incl) * cos(azim));
    c.fwd = mul(normalize(c.pos), -1);
    vec3 world_up = fabs(cos(incl)) > 0.9999 ? v3(sin(azim), 0, -cos(azim)) : v3(0, 1, 0);
    /* right x up = -forward, as for a real observer, so the sky is not mirrored */
    c.right = normalize(cross(c.fwd, world_up));
    c.up = cross(c.right, c.fwd);

    /* r_cam is the Boyer-Lindquist r; theta and phi follow from the direction. */
    c.r = v.r_cam;
    c.mu = cos(incl);
    vec3 B = v3(c.pos.x, -c.pos.z, c.pos.y);
    c.phi = atan2(B.y, B.x);
    double st = sin(incl), cp = cos(c.phi), sp = sin(c.phi);
    c.e_r = v3(st * cp, st * sp, c.mu);
    c.e_th = v3(c.mu * cp, c.mu * sp, -st);
    c.e_ph = v3(-sp, cp, 0);

    double r = c.r, sigma = r * r + a * a * c.mu * c.mu, delta = r * r - 2 * r + a * a;
    double A = (r * r + a * a) * (r * r + a * a) - a * a * delta * st * st;
    c.lapse = sqrt(sigma * delta / A);
    c.omega = 2 * a * r / A;
    c.k_ph = sqrt(A / sigma) * st;
    c.k_th = sqrt(sigma);
    c.k_r = sqrt(sigma / delta);
    c.delta = delta;

    stars_tint(1 / c.lapse);
    c.tan_half = tan(v.fov_deg * M_PI / 360);
    c.aspect = (double)width / height;
    c.star_sigma = 0.9 * 2 * c.tan_half / height;
    c.star_gain = (height / 1080.0) * (height / 1080.0);
    c.y_norm = bb_xyz(v.t_peak).y;
    c.r_out = v.r_out;
    c.step = step;
    c.step_max = step_max;
    return c;
}

/* Ray seen along scene direction n from the camera. */
static Ray camera_ray_dir(const Camera *c, vec3 n) {
    vec3 N = v3(n.x, -n.z, n.y);
    /* The photon travels along -n. With local energy 1, its momentum: */
    double d_r = -dot(N, c->e_r), d_th = -dot(N, c->e_th), d_ph = -dot(N, c->e_ph);
    double L = c->k_ph * d_ph;
    double E = c->lapse + c->omega * L;
    double l = L / E, p_th = c->k_th * d_th / E, p_r = c->k_r * d_r / E;
    double mu = c->mu, sin2 = fmax(1 - mu * mu, 1e-12), a = K.a;
    Ray ray;
    ray.l = l;
    ray.eta = p_th * p_th + mu * mu * (l * l / sin2 - a * a);
    ray.e_obs = E;
    /* Mino-time rates are dr = Delta p_r (so du = -u^2 Delta p_r) and
       dmu = -sin(theta) p_theta; tracing backward flips both. */
    double u = 1 / c->r;
    ray.s = (Geo){u, u * u * c->delta * p_r, mu, sqrt(sin2) * p_th, c->phi};
    return ray;
}

/* Direction n through screen point (sx, sy) in [-1, 1]^2, and its ray. */
static Ray camera_ray(const Camera *c, double sx, double sy, vec3 *n_out) {
    vec3 n = normalize(add(c->fwd, add(mul(c->right, sx * c->tan_half * c->aspect),
                                       mul(c->up, sy * c->tan_half))));
    if (n_out) *n_out = n;
    return camera_ray_dir(c, n);
}

/* Frequency ratio g = nu_obs / nu_emit for light from disk radius r. */
static inline double disk_redshift(const Ray *ray, double r) {
    return 1 / (ray->e_obs * orbit_ut(r) * (1 - orbit_omega(r) * ray->l));
}

/* What a pixel saw, for overlays and plots. */
typedef struct {
    float r, phi_deg, g, t_obs; /* disk hits */
    float ra_deg, dec_deg;      /* sky hits: where on the real sky the ray ends up */
    uint8_t kind, order;
} Aux;

/* XYZ seen through screen point (sx, sy). Fills aux when it is not NULL. */
static vec3 shade_ex(const Camera *c, double sx, double sy, Aux *aux) {
    Ray ray = camera_ray(c, sx, sy, NULL);
    Hit h = trace(&ray, c->r_out, R_ESCAPE, c->step, c->step_max);
    if (aux) *aux = (Aux){.kind = (uint8_t)h.kind, .order = (uint8_t)(h.order > 255 ? 255 : h.order)};
    if (h.kind == HIT_DISK) {
        double r = 1 / h.end.u, g = disk_redshift(&ray, r), t = g * disk_temp(r);
        if (aux) {
            aux->r = (float)r;
            aux->phi_deg = (float)(atan2(h.pos.z, h.pos.x) * 180 / M_PI);
            aux->g = (float)g;
            aux->t_obs = (float)t;
        }
        return mul(bb_xyz(t), 1.0 / c->y_norm);
    }
    if (h.kind == HIT_SKY) {
        if (aux) {
            double ra, dec;
            scene_to_radec(h.dir, &ra, &dec);
            aux->ra_deg = (float)(ra * 180 / M_PI);
            aux->dec_deg = (float)(dec * 180 / M_PI);
        }
        return sky(h.dir, c->star_sigma, c->star_gain);
    }
    return v3(0, 0, 0);
}

/* ---------- rendering ---------- */

static double hash01(uint32_t a, uint32_t b, uint32_t c) {
    uint32_t h = a * 0x9E3779B1u ^ b * 0x85EBCA77u ^ c * 0xC2B2AE3Du;
    h ^= h >> 16; h *= 0x7FEB352Du; h ^= h >> 15; h *= 0x846CA68Bu; h ^= h >> 16;
    return h * (1.0 / 4294967296.0);
}

typedef struct {
    const Camera *cam;
    float *out; /* linear sRGB, 3 floats per pixel */
    Aux *aux;   /* optional, filled from the pixel-center sample when ss == 1 */
    int width, height, ss, y_end;
    atomic_int next_row;
} RenderJob;

static void *render_worker(void *arg) {
    RenderJob *job = arg;
    int W = job->width, H = job->height, ss = job->ss;
    for (int y; (y = atomic_fetch_add(&job->next_row, 1)) < job->y_end;) {
        for (int x = 0; x < W; x++) {
            vec3 acc = v3(0, 0, 0);
            for (int j = 0; j < ss; j++)
                for (int i = 0; i < ss; i++) {
                    uint32_t k = (uint32_t)(j * ss + i);
                    double fx = ss == 1 ? 0.5 : (i + hash01(x, y, 2 * k)) / ss;
                    double fy = ss == 1 ? 0.5 : (j + hash01(x, y, 2 * k + 1)) / ss;
                    double sx = 2 * (x + fx) / W - 1, sy = 1 - 2 * (y + fy) / H;
                    Aux *a = job->aux && ss == 1 ? &job->aux[(size_t)y * W + x] : NULL;
                    acc = add(acc, shade_ex(job->cam, sx, sy, a));
                }
            vec3 rgb = xyz_to_srgb(mul(acc, 1.0 / (ss * ss)));
            float *p = job->out + 3 * ((size_t)y * W + x);
            p[0] = (float)rgb.x; p[1] = (float)rgb.y; p[2] = (float)rgb.z;
        }
    }
    return NULL;
}

/* Render rows [y0, y1) of a width x height image with ss x ss samples per pixel.
   aux (optional) gets per-pixel hit data, written only when ss == 1. */
static void render_rows(const Camera *cam, float *out, Aux *aux, int width, int height, int ss,
                        int y0, int y1, int threads) {
    RenderJob job = {.cam = cam, .out = out, .aux = aux, .width = width, .height = height, .ss = ss,
                     .y_end = y1};
    atomic_init(&job.next_row, y0);
    pthread_t th[64];
    if (threads > 64) threads = 64;
    for (int i = 0; i < threads; i++) pthread_create(&th[i], NULL, render_worker, &job);
    for (int i = 0; i < threads; i++) pthread_join(th[i], NULL);
}

/* ---------- shadow edge ---------- */

typedef struct {
    const Camera *cam;
    double *t; /* edge radius per direction, in tan(angle) units from the view center */
    int n;
    atomic_int next;
} EdgeJob;

static bool captured(const Camera *c, double t, double psi) {
    vec3 n = normalize(add(c->fwd, add(mul(c->right, t * cos(psi)), mul(c->up, t * sin(psi)))));
    Ray ray = camera_ray_dir(c, n);
    return trace(&ray, 0, R_ESCAPE, 0.02, 0.5).kind == HIT_HORIZON;
}

static void *edge_worker(void *arg) {
    EdgeJob *job = arg;
    for (int i; (i = atomic_fetch_add(&job->next, 1)) < job->n;) {
        double psi = 2 * M_PI * i / job->n, lo = 0, hi = 0.5;
        if (!captured(job->cam, 0, psi)) { job->t[i] = 0; continue; }
        while (captured(job->cam, hi, psi) && hi < 64) { lo = hi; hi *= 2; }
        for (int k = 0; k < 26; k++) {
            double mid = 0.5 * (lo + hi);
            if (captured(job->cam, mid, psi)) lo = mid; else hi = mid;
        }
        job->t[i] = 0.5 * (lo + hi);
    }
    return NULL;
}

/* The shadow's edge (critical curve) seen by the camera: for n directions
   around the view center, bisect between captured and escaping rays. */
static inline void shadow_edge(const Camera *c, double *t, int n, int threads) {
    EdgeJob job = {.cam = c, .t = t, .n = n};
    atomic_init(&job.next, 0);
    pthread_t th[64];
    if (threads > 64) threads = 64;
    for (int i = 0; i < threads; i++) pthread_create(&th[i], NULL, edge_worker, &job);
    for (int i = 0; i < threads; i++) pthread_join(th[i], NULL);
}

/* ---------- tone mapping ---------- */

static double srgb_oetf(double x) {
    return x <= 0.0031308 ? 12.92 * x : 1.055 * pow(x, 1 / 2.4) - 0.055;
}

/* Display mapping: an asinh stretch (Lupton et al. 2004), the standard for
   astronomical color images. It is linear near black and logarithmic toward
   the top, so the faint outer disk and the Doppler-boosted inner edge both
   stay visible. Scaling the three channels together keeps their colors. */
#define STRETCH 8.0
static void tonemap(const float *linear, double exposure, double out[3]) {
    double c[3] = {linear[0] * exposure, linear[1] * exposure, linear[2] * exposure};
    double L = 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2];
    if (L <= 0) {
        out[0] = out[1] = out[2] = 0;
        return;
    }
    double k = asinh(STRETCH * L) / (asinh(STRETCH) * L), m = 0;
    for (int i = 0; i < 3; i++) {
        c[i] = fmax(0, c[i] * k);
        m = fmax(m, c[i]);
    }
    for (int i = 0; i < 3; i++) out[i] = srgb_oetf(m > 1 ? c[i] / m : c[i]);
}

/* Disk brightness is normalized so the peak temperature at g = 1 has Y = 1, so
   one fixed exposure fits every view and keeps stars and disk in proportion. */
#define EXPOSURE 0.6

#endif
