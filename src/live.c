/*
 * live.c: interactive Kerr black hole viewer. Same physics as the still renderer
 * (physics.h), drawn with SDL2. Text uses CoreText, which ships with macOS.
 *
 * While the camera moves it renders a coarse frame; once it stops it refines
 * in passes up to full resolution with 4 samples per pixel. Hovering the
 * image traces the ray under the cursor and plots it.
 *
 *   drag        orbit            scroll    distance
 *   1 2 3       view mode        g         coordinate grid
 *   [ ]         disk temperature - =       field of view
 *   space       automatic orbit  h         hide panels
 *   s           save PNG         r         reset
 *   , .         spin             o         Orion view
 *   ?           key list
 *   q, esc      quit
 *
 * Build: make live
 */
#include <CoreText/CoreText.h>
#include <SDL.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include "physics.h"
#include "png.h"

/* ---------- look ---------- */

typedef struct { uint8_t r, g, b, a; } Color;

static const Color C_TEXT = {255, 255, 255, 255};
static const Color C_LABEL = {138, 138, 138, 255};
static const Color C_AXIS = {70, 70, 70, 255};
static const Color C_RULE = {34, 34, 34, 255};
static const Color C_ACCENT = {92, 200, 255, 255};
static const Color C_WARN = {255, 176, 72, 255};

#define FONT_PATH "/System/Library/Fonts/SFNSMono.ttf"
#define PANEL_W 344 /* points */
#define PLOTS_H 228
#define PAD 18
#define ROW_H 17

static const char *MODE_NAMES[] = {"TRUE COLOR", "REDSHIFT g", "IMAGE ORDER"};
enum { MODE_COLOR, MODE_REDSHIFT, MODE_ORDER, MODE_COUNT };

/* ---------- passes ---------- */

typedef struct {
    int div, ss;           /* resolution divisor, ss x ss samples per pixel */
    double step, step_max; /* integrator step, see trace() */
} Pass;

static const Pass PASSES[] = {
    {4, 1, 0.06, 0.5},
    {2, 1, 0.04, 0.5},
    {1, 1, 0.02, 0.5},
    {1, 2, 0.02, 0.5},
};
#define NPASSES (int)(sizeof PASSES / sizeof PASSES[0])
#define BANDS 24

static const View DEFAULT_VIEW = {.r_cam = 80, .incl_deg = 84, .azim_deg = 0, .fov_deg = 20,
                                  .t_peak = 9000, .r_out = 30, .spin = 0.9};
static const double SPINS[] = {-0.998, -0.9, -0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.7, 0.9, 0.998};
#define NSPINS (int)(sizeof SPINS / sizeof SPINS[0])
#define EDGE_N 180

/* ---------- small helpers ---------- */

static double clampd(double x, double lo, double hi) { return x < lo ? lo : x > hi ? hi : x; }
static double lerp(double a, double b, double t) { return a + (b - a) * t; }
static vec3 lerp3(vec3 a, vec3 b, double t) { return add(a, mul(sub(b, a), t)); }
static double now_s(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}
static Color rgba(vec3 c, uint8_t a) {
    return (Color){(uint8_t)lround(255 * clampd(c.x, 0, 1)), (uint8_t)lround(255 * clampd(c.y, 0, 1)),
                   (uint8_t)lround(255 * clampd(c.z, 0, 1)), a};
}
static Color alpha(Color c, double a) { return (Color){c.r, c.g, c.b, (uint8_t)lround(255 * clampd(a, 0, 1))}; }

/* Display color of a blackbody at temperature T, brightest channel scaled to 1. */
static vec3 bb_display(double T) {
    vec3 c = xyz_to_srgb(bb_xyz(T));
    c = v3(fmax(c.x, 0), fmax(c.y, 0), fmax(c.z, 0));
    double m = fmax(c.x, fmax(c.y, c.z));
    if (m <= 0) return v3(0, 0, 0);
    c = mul(c, 1 / m);
    return v3(srgb_oetf(c.x), srgb_oetf(c.y), srgb_oetf(c.z));
}

/* Diverging map for g: red for redshift (g < 1), blue for blueshift (g > 1). */
static vec3 redshift_color(double g) {
    static const vec3 deep_red = {0.50, 0.07, 0.06}, red = {0.97, 0.40, 0.20}, mid = {0.90, 0.90, 0.88},
                      blue = {0.36, 0.64, 1.00}, deep_blue = {0.10, 0.20, 0.66};
    double t = clampd((g - 1) / 0.5, -1, 1);
    if (t < -0.5) return lerp3(deep_red, red, (t + 1) * 2);
    if (t < 0) return lerp3(red, mid, (t + 0.5) * 2);
    if (t < 0.5) return lerp3(mid, blue, t * 2);
    return lerp3(blue, deep_blue, (t - 0.5) * 2);
}

static const vec3 ORDER_COLORS[3] = {{0.82, 0.82, 0.80}, {0.36, 0.78, 1.00}, {1.00, 0.69, 0.28}};

/* ---------- text (CoreText glyph atlas) ---------- */

typedef struct {
    SDL_Texture *tex;
    int cell_w, cell_h, cols;
    float adv; /* advance per glyph, pixels */
} Font;

static SDL_Renderer *R;
static float S = 1; /* pixels per point */
static Font F_MAIN, F_SMALL, F_TITLE;

static const char *EXTRA_GLYPHS = "°θφαΔλνμΩπ≈√∞·×−→←↑↓²³⁴¹⁻≥≤±⌀τ∝¼ω";
static int16_t glyph_index[0x10000];
static uint32_t glyph_cps[256];
static int glyph_count;

static uint32_t utf8_next(const char **s) {
    const unsigned char *p = (const unsigned char *)*s;
    uint32_t c = *p++;
    if (c >= 0xF0 && p[0] && p[1] && p[2]) {
        c = (c & 7) << 18 | (p[0] & 63u) << 12 | (p[1] & 63u) << 6 | (p[2] & 63u);
        p += 3;
    } else if (c >= 0xE0 && p[0] && p[1]) {
        c = (c & 15) << 12 | (p[0] & 63u) << 6 | (p[1] & 63u);
        p += 2;
    } else if (c >= 0xC0 && p[0]) {
        c = (c & 31) << 6 | (p[0] & 63u);
        p += 1;
    }
    *s = (const char *)p;
    return c;
}

static void glyphs_init(void) {
    memset(glyph_index, 0xFF, sizeof glyph_index);
    for (uint32_t c = 32; c < 127; c++) {
        glyph_index[c] = (int16_t)glyph_count;
        glyph_cps[glyph_count++] = c;
    }
    for (const char *s = EXTRA_GLYPHS; *s;) {
        uint32_t c = utf8_next(&s);
        if (c < 0x10000 && glyph_index[c] < 0) {
            glyph_index[c] = (int16_t)glyph_count;
            glyph_cps[glyph_count++] = c;
        }
    }
}

/* Rasterizes every glyph once into a white-on-transparent atlas texture.
   CTLineDraw falls back to other system fonts for glyphs SF Mono lacks. */
static bool font_load(Font *f, double px) {
    CGDataProviderRef dp = CGDataProviderCreateWithFilename(FONT_PATH);
    if (!dp) return false;
    CGFontRef cg = CGFontCreateWithDataProvider(dp);
    CGDataProviderRelease(dp);
    if (!cg) return false;
    CTFontRef ct = CTFontCreateWithGraphicsFont(cg, px, NULL, NULL);
    CGFontRelease(cg);

    double ascent = ceil(CTFontGetAscent(ct)), descent = ceil(CTFontGetDescent(ct));
    UniChar m = 'M';
    CGGlyph gm;
    CGSize adv;
    CTFontGetGlyphsForCharacters(ct, &m, &gm, 1);
    CTFontGetAdvancesForGlyphs(ct, kCTFontOrientationHorizontal, &gm, &adv, 1);
    f->adv = (float)round(adv.width);
    f->cell_w = (int)ceil(adv.width) + 6;
    f->cell_h = (int)(ascent + descent) + 2;
    f->cols = 16;
    int rows = (glyph_count + f->cols - 1) / f->cols;
    int W = f->cols * f->cell_w, H = rows * f->cell_h;

    uint8_t *bits = calloc((size_t)W * H, 1);
    CGContextRef ctx = CGBitmapContextCreate(bits, W, H, 8, W, NULL, (CGBitmapInfo)kCGImageAlphaOnly);
    CGContextSetShouldAntialias(ctx, true);
    CFStringRef keys[] = {kCTFontAttributeName};
    CFTypeRef vals[] = {ct};
    CFDictionaryRef attrs = CFDictionaryCreate(NULL, (const void **)keys, (const void **)vals, 1,
                                               &kCFTypeDictionaryKeyCallBacks, &kCFTypeDictionaryValueCallBacks);
    for (int i = 0; i < glyph_count; i++) {
        UniChar u = (UniChar)glyph_cps[i];
        CFStringRef str = CFStringCreateWithCharacters(NULL, &u, 1);
        CFAttributedStringRef as = CFAttributedStringCreate(NULL, str, attrs);
        CTLineRef line = CTLineCreateWithAttributedString(as);
        int col = i % f->cols, row = i / f->cols;
        CGContextSetTextPosition(ctx, col * f->cell_w + 3, H - (row * f->cell_h + 1 + ascent));
        CTLineDraw(line, ctx);
        CFRelease(line);
        CFRelease(as);
        CFRelease(str);
    }
    CFRelease(attrs);
    CGContextRelease(ctx);
    CFRelease(ct);

    uint32_t *argb = malloc(sizeof(uint32_t) * W * H);
    for (int i = 0; i < W * H; i++) argb[i] = (uint32_t)bits[i] << 24 | 0xFFFFFFu;
    f->tex = SDL_CreateTexture(R, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STATIC, W, H);
    SDL_UpdateTexture(f->tex, NULL, argb, 4 * W);
    SDL_SetTextureBlendMode(f->tex, SDL_BLENDMODE_BLEND);
    SDL_SetTextureScaleMode(f->tex, SDL_ScaleModeNearest);
    free(argb);
    free(bits);
    return f->tex != NULL;
}

static int utf8_len(const char *s) {
    int n = 0;
    while (*s) { utf8_next(&s); n++; }
    return n;
}
static float text_w(const Font *f, const char *s) { return utf8_len(s) * f->adv; }
static float line_h(const Font *f) { return (float)f->cell_h; }

/* Draws s with its top-left at (x, y), in pixels. Returns the pen position. */
static float text(const Font *f, float x, float y, Color c, const char *s) {
    SDL_SetTextureColorMod(f->tex, c.r, c.g, c.b);
    SDL_SetTextureAlphaMod(f->tex, c.a);
    float pen = roundf(x);
    int top = (int)roundf(y);
    while (*s) {
        uint32_t cp = utf8_next(&s);
        int gi = cp < 0x10000 ? glyph_index[cp] : -1;
        if (gi >= 0 && cp != ' ') {
            SDL_Rect src = {(gi % f->cols) * f->cell_w, (gi / f->cols) * f->cell_h, f->cell_w, f->cell_h};
            SDL_Rect dst = {(int)pen - 3, top, f->cell_w, f->cell_h};
            SDL_RenderCopy(R, f->tex, &src, &dst);
        }
        pen += f->adv;
    }
    return pen;
}

static float textf(const Font *f, float x, float y, Color c, const char *fmt, ...) {
    char buf[256];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    return text(f, x, y, c, buf);
}
static void text_right(const Font *f, float x_right, float y, Color c, const char *s) {
    text(f, x_right - text_w(f, s), y, c, s);
}
static void text_center(const Font *f, float x_mid, float y, Color c, const char *s) {
    text(f, x_mid - text_w(f, s) / 2, y, c, s);
}

/* ---------- drawing primitives (pixels) ---------- */

static void set_color(Color c) { SDL_SetRenderDrawColor(R, c.r, c.g, c.b, c.a); }
static void line(float x0, float y0, float x1, float y1, Color c) {
    set_color(c);
    SDL_RenderDrawLineF(R, x0, y0, x1, y1);
}
static void thick_line(float x0, float y0, float x1, float y1, float w, Color c) {
    float dx = x1 - x0, dy = y1 - y0, l = sqrtf(dx * dx + dy * dy);
    if (l <= 0) return;
    float nx = -dy / l, ny = dx / l;
    set_color(c);
    for (float o = -(w - 1) / 2; o <= (w - 1) / 2 + 1e-3f; o += 1)
        SDL_RenderDrawLineF(R, x0 + nx * o, y0 + ny * o, x1 + nx * o, y1 + ny * o);
}
static void fill(SDL_FRect r, Color c) {
    set_color(c);
    SDL_RenderFillRectF(R, &r);
}
static void circle(float cx, float cy, float rad, Color c, float dash) {
    int n = (int)clampd(rad * 0.8, 48, 720);
    set_color(c);
    for (int i = 0; i < n; i++) {
        if (dash > 0 && (int)(i * (2 * M_PI * rad / n) / dash) % 2) continue;
        double a0 = 2 * M_PI * i / n, a1 = 2 * M_PI * (i + 1) / n;
        SDL_RenderDrawLineF(R, cx + rad * cos(a0), cy + rad * sin(a0), cx + rad * cos(a1), cy + rad * sin(a1));
    }
}
static void disc(float cx, float cy, float rad, Color c) {
    set_color(c);
    for (float y = -rad; y <= rad; y += 1) {
        float hw = sqrtf(fmaxf(rad * rad - y * y, 0));
        SDL_RenderDrawLineF(R, cx - hw, cy + y, cx + hw, cy + y);
    }
}
static SDL_FRect inset(SDL_FRect r, float l, float t, float rt, float b) {
    return (SDL_FRect){r.x + l, r.y + t, r.w - l - rt, r.h - t - b};
}

/* ---------- app state ---------- */

#define HIST_N 72
#define THIST_N 256 /* log-spaced observed temperature bins, 300 K to 100000 K */
#define HIST_G0 0.3
#define HIST_G1 1.7
#define PATH_CAP 40000

typedef struct {
    float x, y; /* pixel in the composed image */
    double mag; /* magnification of the image there; 0 if the star has no image in view */
} Label;

typedef struct {
    Label labels[32];
    int hist[HIST_N], hist_max, thist[THIST_N];
    int disk_px, total_px, order_px[3];
    float gmin, gmax;
    float tmin, tmax; /* 1st and 99th percentile of the observed temperature */
} Stats;

typedef struct {
    double u, v; /* position in the view, 0..1 */
    Ray ray;
    Hit hit;
    vec3 n;
    double r, phi, g, t_emit, t_obs, ra, dec, deflection, r_min;
    vec3 e1, e2; /* plot plane: e1 toward the camera, e2 along the initial direction */
    Path path;
} Probe;

/* The shadow's edge as the camera sees it, in tan(angle) units from the center. */
typedef struct {
    double x[EDGE_N], y[EDGE_N];
    double width_deg, height_deg, shift_deg;
} Edge;

typedef struct {
    View view;
    int mode, pass;
    bool grid, hud, orbit, dragging, quit, keys;
    bool view_dirty, compose_dirty, probe_dirty, ui_dirty;

    /* the last finished pass ("front") and the one in progress ("work") */
    float *front_img, *work_img;
    Aux *front_aux, *work_aux;
    uint8_t *rgb;
    size_t cap;
    int img_w, img_h, img_pass;
    double pass_time, pass_started, last_draw;
    float progress;
    SDL_Texture *tex;
    int tex_w, tex_h;

    Stats stats;
    Probe probe;
    Edge edge;
    char saved[128];
    int threads;
    SDL_Window *win;
} App;

typedef struct { SDL_FRect view, panel, plots; } Layout;

static Layout layout(const App *a) {
    int pw, ph;
    SDL_GetRendererOutputSize(R, &pw, &ph);
    Layout L;
    if (!a->hud) {
        L.view = (SDL_FRect){0, 0, (float)pw, (float)ph};
        L.panel = L.plots = (SDL_FRect){0, 0, 0, 0};
        return L;
    }
    float panel_w = PANEL_W * S, plots_h = PLOTS_H * S;
    L.view = (SDL_FRect){0, 0, pw - panel_w, ph - plots_h};
    L.plots = (SDL_FRect){0, ph - plots_h, pw - panel_w, plots_h};
    L.panel = (SDL_FRect){pw - panel_w, 0, panel_w, (float)ph};
    return L;
}

/* Render resolution of the view for a pass: one sample per point at full res. */
static void pass_size(const App *a, int pass, int *w, int *h) {
    Layout L = layout(a);
    int div = PASSES[pass].div;
    *w = (int)fmax(1, floor(L.view.w / S / div));
    *h = (int)fmax(1, floor(L.view.h / S / div));
}

static void ensure_capacity(App *a, size_t px) {
    if (px <= a->cap) return;
    a->cap = px;
    a->front_img = realloc(a->front_img, sizeof(float) * 3 * px);
    a->work_img = realloc(a->work_img, sizeof(float) * 3 * px);
    a->front_aux = realloc(a->front_aux, sizeof(Aux) * px);
    a->work_aux = realloc(a->work_aux, sizeof(Aux) * px);
    a->rgb = realloc(a->rgb, 3 * px);
}

/* ---------- probe: one precise ray under the cursor ---------- */

/* Boyer-Lindquist r of a scene point (oblate spheroidal radius). */
static double scene_r(vec3 q) {
    double a2 = K.a * K.a, h = 0.5 * (dot(q, q) - a2);
    return sqrt(h + sqrt(h * h + a2 * q.y * q.y));
}

static Camera view_camera(const App *a) {
    Layout L = layout(a);
    return camera_make(a->view, (int)fmax(1, L.view.w), (int)fmax(1, L.view.h), 0.02, 0.02);
}

static void probe_update(App *a) {
    Probe *p = &a->probe;
    physics_sync(&a->view);
    Camera cam = view_camera(a);
    p->ray = camera_ray(&cam, 2 * p->u - 1, 1 - 2 * p->v, &p->n);
    p->path.n = 0;
    p->hit = trace_path(&p->ray, a->view.r_out, R_ESCAPE, 0.02, 0.02, &p->path);
    p->r_min = 1e30;
    for (int i = 0; i < p->path.n; i++) p->r_min = fmin(p->r_min, scene_r(p->path.pts[i]));

    p->e1 = normalize(cam.pos);
    vec3 t = sub(p->n, mul(p->e1, dot(p->n, p->e1)));
    p->e2 = len(t) > 1e-9 ? normalize(t) : cam.up;

    if (p->hit.kind == HIT_DISK) {
        p->r = 1 / p->hit.end.u;
        p->phi = atan2(p->hit.pos.z, p->hit.pos.x) * 180 / M_PI;
        p->g = disk_redshift(&p->ray, p->r);
        p->t_emit = disk_temp(p->r);
        p->t_obs = p->g * p->t_emit;
    } else if (p->hit.kind == HIT_SKY) {
        scene_to_radec(p->hit.dir, &p->ra, &p->dec);
        p->deflection = acos(clampd(dot(p->n, p->hit.dir), -1, 1)) * 180 / M_PI;
    }
}

static void edge_update(App *a) {
    physics_sync(&a->view);
    Camera cam = view_camera(a);
    double t[EDGE_N], xmin = 1e9, xmax = -1e9, ymin = 1e9, ymax = -1e9;
    shadow_edge(&cam, t, EDGE_N, a->threads);
    Edge *e = &a->edge;
    for (int i = 0; i < EDGE_N; i++) {
        double psi = 2 * M_PI * i / EDGE_N;
        e->x[i] = t[i] * cos(psi);
        e->y[i] = t[i] * sin(psi);
        xmin = fmin(xmin, e->x[i]);
        xmax = fmax(xmax, e->x[i]);
        ymin = fmin(ymin, e->y[i]);
        ymax = fmax(ymax, e->y[i]);
    }
    e->width_deg = (atan(xmax) - atan(xmin)) * 180 / M_PI;
    e->height_deg = (atan(ymax) - atan(ymin)) * 180 / M_PI;
    e->shift_deg = (atan(xmax) + atan(xmin)) / 2 * 180 / M_PI;
}

/* Speed of the disk gas at r measured by a local ZAMO, in units of c. */
static double disk_speed(double r) {
    double a = K.a, delta = r * r - 2 * r + a * a, A = (r * r + a * a) * (r * r + a * a) - a * a * delta;
    return A * (orbit_omega(r) - 2 * a * r / A) / (r * r * sqrt(delta));
}

/* ---------- named stars ---------- */

static const struct { int hr; const char *name; } NAMED[] = {
    {2061, "Betelgeuse"}, {1713, "Rigel"},   {1790, "Bellatrix"}, {2004, "Saiph"},   {1948, "Alnitak"},
    {1903, "Alnilam"},    {1852, "Mintaka"}, {2491, "Sirius"},    {1457, "Aldebaran"}, {2943, "Procyon"},
    {1708, "Capella"},    {2990, "Pollux"},  {2891, "Castor"},    {2326, "Canopus"}, {424, "Polaris"},
    {7001, "Vega"},       {5340, "Arcturus"}, {6134, "Antares"},  {5056, "Spica"},   {7924, "Deneb"},
    {7557, "Altair"},     {472, "Achernar"}, {8728, "Fomalhaut"}, {3982, "Regulus"},
};
#define N_NAMED (int)(sizeof NAMED / sizeof NAMED[0])

static vec3 named_eq[N_NAMED]; /* equatorial unit vectors */

static void named_init(void) {
    for (int k = 0; k < N_NAMED; k++)
        for (int i = 0; i < BSC5_COUNT; i++)
            if (BSC5[i].hr == NAMED[k].hr) {
                double ra = BSC5[i].ra, dec = BSC5[i].dec;
                named_eq[k] = v3(cos(dec) * cos(ra), cos(dec) * sin(ra), sin(dec));
            }
}

/* ---------- compose: modes, contour grid, statistics ---------- */

enum { F_R, F_PHI, F_G, F_DEC, F_RA };

static double field(const Aux *x, int f) {
    switch (f) {
    case F_R: return x->r;
    case F_PHI: return x->phi_deg;
    case F_G: return x->g;
    case F_DEC: return x->dec_deg;
    default: return x->ra_deg;
    }
}

/* Screen-space gradient magnitude of a field, using only neighbors that belong
   to the same image (same hit kind and order). Angles wrap at 360. */
static double gradient(const Aux *aux, int w, int h, int x, int y, int f, bool wrap) {
    const Aux *c = &aux[(size_t)y * w + x];
    double v = field(c, f), g[2] = {0, 0};
    for (int axis = 0; axis < 2; axis++) {
        double sum = 0;
        int n = 0;
        for (int s = -1; s <= 1; s += 2) {
            int nx = x + (axis == 0 ? s : 0), ny = y + (axis == 1 ? s : 0);
            if (nx < 0 || ny < 0 || nx >= w || ny >= h) continue;
            const Aux *o = &aux[(size_t)ny * w + nx];
            if (o->kind != c->kind || o->order != c->order) continue;
            double d = field(o, f) - v;
            if (wrap) d -= 360 * round(d / 360);
            sum += d * s;
            n++;
        }
        if (n) g[axis] = sum / n;
    }
    return sqrt(g[0] * g[0] + g[1] * g[1]);
}

/* Anti-aliased coverage of contour lines every `spacing` units. width is in
   full-resolution pixels and div is the pass's resolution divisor, so lines keep
   their width in coarse passes. Lines fade out where they would crowd. */
static double contour(double v, double grad, double spacing, double width, int div) {
    if (grad <= 1e-12) return 0;
    double d = fabs(remainder(v, spacing)) / grad;
    double cover = clampd(width / div / 2 + 0.5 - d, 0, 1) * fmin(1, width / div);
    double crowd = clampd((spacing / grad * div - 4) / 6, 0, 1);
    return cover * crowd;
}

static void blend(double *px, vec3 c, double a) {
    px[0] = lerp(px[0], c.x, a);
    px[1] = lerp(px[1], c.y, a);
    px[2] = lerp(px[2], c.z, a);
}

static void compose(App *a) {
    int w = a->img_w, h = a->img_h;
    const float *img = a->front_img;
    const Aux *aux = a->front_aux;
    Stats *st = &a->stats;
    memset(st, 0, sizeof *st);
    st->gmin = 1e30f;
    st->gmax = -1e30f;
    st->total_px = w * h;
    const vec3 accent = {C_ACCENT.r / 255.0, C_ACCENT.g / 255.0, C_ACCENT.b / 255.0};
    const vec3 white = {1, 1, 1}, black = {0, 0, 0};
    const int div = PASSES[a->img_pass].div;
    /* a pixel holds an image of a star when its ray lands within two pixel widths of it */
    double pix = 2 * tan(a->view.fov_deg * M_PI / 360) / h, near = 2 * pix * pix;

    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++) {
            size_t i = (size_t)y * w + x;
            const float *p = img + 3 * i;
            const Aux *q = &aux[i];
            double tc[3];
            tonemap(p, EXPOSURE, tc);
            double lum = 0.2126 * tc[0] + 0.7152 * tc[1] + 0.0722 * tc[2];
            double px[3] = {tc[0], tc[1], tc[2]};

            if (q->kind == HIT_DISK) {
                st->disk_px++;
                st->order_px[q->order > 2 ? 2 : q->order]++;
                st->gmin = fminf(st->gmin, q->g);
                st->gmax = fmaxf(st->gmax, q->g);
                if (q->t_obs > 300) {
                    int tb = (int)(log(q->t_obs / 300) / log(100000 / 300.0) * THIST_N);
                    st->thist[tb < THIST_N ? tb : THIST_N - 1]++;
                }
                int bin = (int)((q->g - HIST_G0) / (HIST_G1 - HIST_G0) * HIST_N);
                if (bin >= 0 && bin < HIST_N) st->hist[bin]++;

                if (a->mode == MODE_REDSHIFT) {
                    vec3 c = redshift_color(q->g);
                    px[0] = c.x; px[1] = c.y; px[2] = c.z;
                } else if (a->mode == MODE_ORDER) {
                    vec3 c = mul(ORDER_COLORS[q->order > 2 ? 2 : q->order], 0.22 + 0.78 * lum);
                    px[0] = c.x; px[1] = c.y; px[2] = c.z;
                }
                if (a->grid) {
                    double gr = gradient(aux, w, h, x, y, F_R, false);
                    double ring = q->r > 7.5 && q->r < a->view.r_out - 2.5 ? contour(q->r, gr, 5, 1.1, div) : 0;
                    double spoke = contour(q->phi_deg, gradient(aux, w, h, x, y, F_PHI, true), 30, 1.1, div);
                    vec3 lc = a->mode == MODE_COLOR ? accent : black;
                    blend(px, lc, 0.55 * fmax(ring, spoke));
                }
                if (a->mode == MODE_REDSHIFT) {
                    double gg = gradient(aux, w, h, x, y, F_G, false);
                    blend(px, black, 0.6 * contour(q->g - 1, gg, 0.1, 1.2, div));
                }
            } else if (q->kind == HIT_SKY) {
                if (a->mode != MODE_COLOR) for (int k = 0; k < 3; k++) px[k] *= 0.45;
                double ra = q->ra_deg * M_PI / 180, dec = q->dec_deg * M_PI / 180;
                vec3 eq = v3(cos(dec) * cos(ra), cos(dec) * sin(ra), sin(dec));
                for (int k = 0; k < N_NAMED; k++) {
                    if (1 - dot(eq, named_eq[k]) > near) continue;
                    /* Several images of one star can be in view. Label the most magnified:
                       there the sky changes least from one pixel to the next. */
                    double gd = gradient(aux, w, h, x, y, F_DEC, false) * M_PI / 180;
                    double gr = gradient(aux, w, h, x, y, F_RA, true) * M_PI / 180 * cos(dec);
                    double mag = pix * pix / (gd * gd + gr * gr + 1e-30);
                    if (mag > st->labels[k].mag) st->labels[k] = (Label){(float)x, (float)y, mag};
                }
                if (a->grid) {
                    /* equatorial grid: declination every 10°, right ascension every hour */
                    double la = contour(q->dec_deg, gradient(aux, w, h, x, y, F_DEC, false), 10, 1.0, div);
                    double lo = fabs(q->dec_deg) < 80
                                    ? contour(q->ra_deg, gradient(aux, w, h, x, y, F_RA, true), 15, 1.0, div)
                                    : 0;
                    blend(px, a->mode == MODE_COLOR ? accent : white, 0.16 * fmax(la, lo));
                }
            }
            uint8_t *o = a->rgb + 3 * i;
            for (int k = 0; k < 3; k++) o[k] = (uint8_t)lround(255 * clampd(px[k], 0, 1));
        }
    for (int b = 0; b < HIST_N; b++) st->hist_max = st->hist[b] > st->hist_max ? st->hist[b] : st->hist_max;
    int tcount = 0, acc = 0;
    for (int b = 0; b < THIST_N; b++) tcount += st->thist[b];
    for (int b = 0; b < THIST_N && tcount; b++) {
        double t_mid = 300 * pow(100000 / 300.0, (b + 0.5) / THIST_N);
        if (acc < 0.01 * tcount && acc + st->thist[b] >= 0.01 * tcount) st->tmin = (float)t_mid;
        if (acc < 0.99 * tcount && acc + st->thist[b] >= 0.99 * tcount) st->tmax = (float)t_mid;
        acc += st->thist[b];
    }

    if (!a->tex || a->tex_w != w || a->tex_h != h) {
        if (a->tex) SDL_DestroyTexture(a->tex);
        a->tex = SDL_CreateTexture(R, SDL_PIXELFORMAT_RGB24, SDL_TEXTUREACCESS_STREAMING, w, h);
        SDL_SetTextureScaleMode(a->tex, SDL_ScaleModeLinear);
        a->tex_w = w;
        a->tex_h = h;
    }
    SDL_UpdateTexture(a->tex, NULL, a->rgb, 3 * w);
}

/* ---------- view overlay ---------- */

static double nice_step(double approx) {
    static const double steps[] = {0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 15, 30, 45};
    for (int i = 0; i < (int)(sizeof steps / sizeof steps[0]); i++)
        if (steps[i] >= approx) return steps[i];
    return 90;
}

static void draw_view_overlay(const App *a, SDL_FRect v) {
    float cx = v.x + v.w / 2, cy = v.y + v.h / 2;
    double tan_half = tan(a->view.fov_deg * M_PI / 360), px_per_tan = (v.h / 2) / tan_half;
    SDL_Rect clip = {(int)v.x, (int)v.y, (int)v.w, (int)v.h};
    SDL_RenderSetClipRect(R, &clip);

    /* critical curve, traced numerically for this camera */
    const Edge *e = &a->edge;
    int top_i = 0;
    for (int i = 0; i < EDGE_N; i++) {
        int j = (i + 1) % EDGE_N;
        if (e->y[i] > e->y[top_i]) top_i = i;
        if (i % 2) continue;
        line(cx + (float)(e->x[i] * px_per_tan), cy - (float)(e->y[i] * px_per_tan),
             cx + (float)(e->x[j] * px_per_tan), cy - (float)(e->y[j] * px_per_tan), alpha(C_ACCENT, 0.8));
    }
    float lx = cx + (float)(e->x[top_i] * px_per_tan), ly = cy - (float)(e->y[top_i] * px_per_tan);
    float top = v.y + 44 * S;
    if (ly > top + 20 * S) {
        line(lx, ly, lx, top, alpha(C_ACCENT, 0.55));
        line(lx, top, lx + 10 * S, top, alpha(C_ACCENT, 0.55));
        disc(lx, ly, 2 * S, C_ACCENT);
        textf(&F_SMALL, lx + 14 * S, top - line_h(&F_SMALL) / 2, C_ACCENT, "CRITICAL CURVE  %.3f° × %.3f°",
              e->width_deg, e->height_deg);
    }

    /* reticle */
    float g = 6 * S, l = 16 * S;
    Color rc = alpha(C_TEXT, 0.35);
    line(cx - g - l, cy, cx - g, cy, rc);
    line(cx + g, cy, cx + g + l, cy, rc);
    line(cx, cy - g - l, cx, cy - g, rc);
    line(cx, cy + g, cx, cy + g + l, rc);

    /* angular scales along the left and bottom edges */
    double deg_per_px = a->view.fov_deg / v.h;
    double step = nice_step(deg_per_px * 90 * S);
    Color tc = alpha(C_TEXT, 0.45);
    for (int k = -200; k <= 200; k++) {
        for (int minor = 0; minor < 5; minor++) {
            double ang = (k + minor / 5.0) * step;
            if (fabs(ang) >= 89) continue;
            float off = (float)(tan(ang * M_PI / 180) * px_per_tan);
            float len_ = minor == 0 ? 9 * S : 4 * S;
            float yy = cy - off, xx = cx + off;
            if (yy > v.y && yy < v.y + v.h) line(v.x, yy, v.x + len_, yy, tc);
            if (xx > v.x && xx < v.x + v.w) line(xx, v.y + v.h - len_, xx, v.y + v.h, tc);
            if (minor == 0 && k != 0) {
                char lab[32];
                snprintf(lab, sizeof lab, step < 1 ? "%+.1f°" : "%+.0f°", ang);
                if (yy > v.y + 20 * S && yy < v.y + v.h - 30 * S)
                    text(&F_SMALL, v.x + 13 * S, yy - line_h(&F_SMALL) / 2, tc, lab);
                if (xx > v.x + 40 * S && xx < v.x + v.w - 40 * S)
                    text_center(&F_SMALL, xx, v.y + v.h - 13 * S - line_h(&F_SMALL), tc, lab);
            }
        }
    }

    /* named stars, at their most magnified image in view */
    if (a->grid && a->img_w) {
        for (int k = 0; k < N_NAMED; k++) {
            const Label *lb = &a->stats.labels[k];
            if (lb->mag <= 0) continue;
            float lx2 = v.x + (lb->x + 0.5f) / a->img_w * v.w, ly2 = v.y + (lb->y + 0.5f) / a->img_h * v.h;
            if (ly2 < v.y + 34 * S || ly2 > v.y + v.h - 34 * S || lx2 > v.x + v.w - 110 * S) continue;
            circle(lx2, ly2, 7 * S, alpha(C_TEXT, 0.55), 0);
            text(&F_SMALL, lx2 + 11 * S, ly2 - line_h(&F_SMALL) / 2, alpha(C_TEXT, 0.8), NAMED[k].name);
        }
    }

    /* probe marker */
    float px = v.x + (float)a->probe.u * v.w, py = v.y + (float)a->probe.v * v.h;
    float pg = 4 * S, pl = 7 * S;
    line(px - pg - pl, py, px - pg, py, C_ACCENT);
    line(px + pg, py, px + pg + pl, py, C_ACCENT);
    line(px, py - pg - pl, px, py - pg, C_ACCENT);
    line(px, py + pg, px, py + pg + pl, C_ACCENT);

    /* corner tags */
    float tx = text(&F_SMALL, v.x + 14 * S, v.y + 12 * S, C_TEXT, MODE_NAMES[a->mode]) + 14 * S;
    tx = textf(&F_SMALL, tx, v.y + 12 * S, alpha(C_TEXT, 0.7), "a = %+.3f", a->view.spin) + 14 * S;
    if (a->grid) text(&F_SMALL, tx, v.y + 12 * S, alpha(C_ACCENT, 0.9), "GRID");
    text_right(&F_SMALL, v.x + v.w - 14 * S, v.y + 12 * S, alpha(C_TEXT, a->keys ? 1 : 0.5), "? KEYS");
    if (a->orbit) text_right(&F_SMALL, v.x + v.w - 14 * S - text_w(&F_SMALL, "? KEYS") - 16 * S, v.y + 12 * S, C_WARN,
                             "AUTO ORBIT");
    if (a->keys) {
        static const char *keys[][2] = {{"drag", "orbit"},       {"scroll", "distance"}, {"arrows", "orbit 5°"},
                                        {", .", "spin"},         {"o", "Orion view"},
                                        {"1 2 3", "view mode"},  {"g", "grid"},          {"[ ]", "disk temperature"},
                                        {"- =", "field of view"}, {"space", "auto orbit"}, {"h", "panels"},
                                        {"s", "save png"},        {"r", "reset"},          {"q", "quit"}};
        int nk = (int)(sizeof keys / sizeof keys[0]);
        float kw = 30 * F_SMALL.adv, kh = nk * 15 * S + 20 * S;
        SDL_FRect box = {v.x + v.w - kw - 14 * S, v.y + 32 * S, kw, kh};
        fill(box, (Color){0, 0, 0, 215});
        for (int k = 0; k < nk; k++) {
            float yy = box.y + 10 * S + k * 15 * S;
            text(&F_SMALL, box.x + 12 * S, yy, C_TEXT, keys[k][0]);
            text(&F_SMALL, box.x + 12 * S + 9 * F_SMALL.adv, yy, C_LABEL, keys[k][1]);
        }
    }
    SDL_RenderSetClipRect(R, NULL);
}

/* ---------- plots ---------- */

typedef struct { SDL_FRect r; double x0, x1, y0, y1; } Axes;

static float ax_x(const Axes *a, double x) { return a->r.x + (float)((x - a->x0) / (a->x1 - a->x0)) * a->r.w; }
static float ax_y(const Axes *a, double y) {
    return a->r.y + a->r.h - (float)((y - a->y0) / (a->y1 - a->y0)) * a->r.h;
}

static void axes_frame(const Axes *a) {
    line(a->r.x, a->r.y, a->r.x, a->r.y + a->r.h, C_AXIS);
    line(a->r.x, a->r.y + a->r.h, a->r.x + a->r.w, a->r.y + a->r.h, C_AXIS);
}
static void xtick(const Axes *a, double x, const char *lab) {
    float px = ax_x(a, x), y = a->r.y + a->r.h;
    line(px, y, px, y + 4 * S, C_AXIS);
    text_center(&F_SMALL, px, y + 6 * S, C_LABEL, lab);
}
static void ytick(const Axes *a, double y, const char *lab) {
    float py = ax_y(a, y);
    line(a->r.x - 4 * S, py, a->r.x, py, C_AXIS);
    text_right(&F_SMALL, a->r.x - 7 * S, py - line_h(&F_SMALL) / 2, C_LABEL, lab);
}

/* Title row and plotting area for one plot cell. */
static Axes plot_cell(SDL_FRect cell, const char *title, const char *note, double x0, double x1, double y0,
                      double y1, float left) {
    text(&F_SMALL, cell.x, cell.y, C_ACCENT, title);
    if (note) text_right(&F_SMALL, cell.x + cell.w, cell.y, C_LABEL, note);
    return (Axes){inset(cell, left, line_h(&F_SMALL) + 12 * S, 0, line_h(&F_SMALL) + 8 * S), x0, x1, y0, y1};
}

static void plot_ray(const App *a, SDL_FRect cell) {
    const Probe *p = &a->probe;
    text(&F_SMALL, cell.x, cell.y, C_ACCENT, "RAY  projected");
    SDL_FRect area = inset(cell, 0, line_h(&F_SMALL) + 8 * S, 0, 0);
    double R_view = fmin(a->view.r_cam * 1.06, fmax(a->view.r_out * 0.62, 11));
    float scale = (float)(fmin(area.w, area.h) / 2 / R_view);
    float cx = area.x + area.w / 2, cy = area.y + area.h / 2;
#define MAPX(q) (cx - (float)dot((q), p->e1) * scale)
#define MAPY(q) (cy - (float)dot((q), p->e2) * scale)
    SDL_Rect clip = {(int)area.x, (int)area.y, (int)area.w, (int)area.h};
    SDL_RenderSetClipRect(R, &clip);

    /* where the disk cuts this plane, colored by its temperature */
    vec3 N = cross(p->e1, p->e2), d = cross(N, v3(0, 1, 0));
    if (len(d) > 1e-6) {
        d = normalize(d);
        for (int side = -1; side <= 1; side += 2)
            for (double r = K.r_isco; r < a->view.r_out; r += 0.25) {
                vec3 q0 = mul(d, side * r), q1 = mul(d, side * fmin(r + 0.25, a->view.r_out));
                thick_line(MAPX(q0), MAPY(q0), MAPX(q1), MAPY(q1), 3 * S, rgba(bb_display(disk_temp(r + 0.12)), 255));
            }
    }
    circle(cx, cy, (float)K.r_ph_min * scale, alpha(C_TEXT, 0.35), 3 * S);
    circle(cx, cy, (float)K.r_ph_max * scale, alpha(C_TEXT, 0.35), 3 * S);
    circle(cx, cy, (float)K.r_isco * scale, alpha(C_ACCENT, 0.3), 2 * S);
    disc(cx, cy, (float)K.r_plus * scale, (Color){0, 0, 0, 255});
    circle(cx, cy, (float)K.r_plus * scale, alpha(C_TEXT, 0.85), 0);

    /* the geodesic */
    for (int i = 1; i < p->path.n; i++) {
        vec3 q0 = p->path.pts[i - 1], q1 = p->path.pts[i];
        thick_line(MAPX(q0), MAPY(q0), MAPX(q1), MAPY(q1), 2 * S, C_ACCENT);
    }
    if (p->path.n) {
        vec3 end = p->path.pts[p->path.n - 1];
        if (p->hit.kind == HIT_DISK) disc(MAPX(end), MAPY(end), 3.5f * S, C_TEXT);
    }
    SDL_RenderSetClipRect(R, NULL);

    /* camera: on the left edge if it lies outside the frame */
    float camx = cx - (float)a->view.r_cam * scale;
    if (camx < area.x) {
        camx = area.x + 4 * S;
        text(&F_SMALL, area.x, area.y + area.h - line_h(&F_SMALL), C_LABEL, "← observer");
    }
    disc(camx, cy, 3 * S, C_TEXT);

    /* scale bar */
    float bar = 10 * scale, by = area.y + area.h - 4 * S;
    line(area.x + area.w - bar, by, area.x + area.w, by, C_LABEL);
    text_right(&F_SMALL, area.x + area.w, by - line_h(&F_SMALL) - 2 * S, C_LABEL, "10 M");
#undef MAPX
#undef MAPY
}

static void plot_temperature(const App *a, SDL_FRect cell) {
    double r1 = a->view.r_out, tmax = a->view.t_peak * 1.12;
    char note[48];
    snprintf(note, sizeof note, "T ∝ F¼");
    Axes ax = plot_cell(cell, "EMISSION  T(r)", note, K.r_isco, r1, 0, tmax, 40 * S);
    for (float x = ax.r.x; x < ax.r.x + ax.r.w; x += 1) {
        double r = lerp(K.r_isco, r1, (x - ax.r.x) / ax.r.w), t = disk_temp(r);
        line(x, ax_y(&ax, t), x, ax_y(&ax, 0), rgba(mul(bb_display(t), 0.55), 255));
    }
    float px = 0, py = 0;
    for (float x = ax.r.x; x <= ax.r.x + ax.r.w; x += 1) {
        double r = lerp(K.r_isco, r1, (x - ax.r.x) / ax.r.w);
        float y = ax_y(&ax, disk_temp(r));
        if (x > ax.r.x) thick_line(px, py, x, y, 2 * S, C_TEXT);
        px = x;
        py = y;
    }
    if (nt_rpeak < r1) {
        float xp = ax_x(&ax, nt_rpeak);
        for (float y = ax_y(&ax, a->view.t_peak); y < ax.r.y + ax.r.h; y += 6 * S)
            line(xp, y, xp, fminf(y + 3 * S, ax.r.y + ax.r.h), alpha(C_TEXT, 0.5));
    }
    if (a->probe.hit.kind == HIT_DISK) {
        float xp = ax_x(&ax, a->probe.r);
        line(xp, ax.r.y, xp, ax.r.y + ax.r.h, C_ACCENT);
        disc(xp, ax_y(&ax, a->probe.t_emit), 3 * S, C_ACCENT);
    }
    axes_frame(&ax);
    double step = r1 > 40 ? 20 : 5;
    for (double r = 10; r <= r1 + 1e-9; r += step) {
        char lab[16];
        snprintf(lab, sizeof lab, r + step > r1 + 1e-9 ? "%.0f M" : "%.0f", r);
        xtick(&ax, r, lab);
    }
    char isco[16];
    snprintf(isco, sizeof isco, "%.1f", K.r_isco);
    if (K.r_isco < 8.5) xtick(&ax, K.r_isco, isco);
    char lab[16];
    snprintf(lab, sizeof lab, "%.1fk", a->view.t_peak / 1000);
    ytick(&ax, a->view.t_peak, lab);
    snprintf(lab, sizeof lab, "%.1fk", a->view.t_peak / 2000);
    ytick(&ax, a->view.t_peak / 2, lab);
    ytick(&ax, 0, "0");
}

static void plot_histogram(const App *a, SDL_FRect cell) {
    const Stats *s = &a->stats;
    char note[48] = "";
    if (s->disk_px) snprintf(note, sizeof note, "%.1f%% of view", 100.0 * s->disk_px / s->total_px);
    Axes ax = plot_cell(cell, "REDSHIFT  g on disk", note, HIST_G0, HIST_G1, 0, 1.08, 14 * S);
    if (s->hist_max > 0) {
        for (int b = 0; b < HIST_N; b++) {
            double g0 = HIST_G0 + (HIST_G1 - HIST_G0) * b / HIST_N, g1 = g0 + (HIST_G1 - HIST_G0) / HIST_N;
            float x0 = ax_x(&ax, g0) + S * 0.5f, x1 = ax_x(&ax, g1) - S * 0.5f;
            float y = ax_y(&ax, (double)s->hist[b] / s->hist_max);
            fill((SDL_FRect){x0, y, fmaxf(x1 - x0, 1), ax.r.y + ax.r.h - y}, rgba(redshift_color((g0 + g1) / 2), 255));
        }
    }
    float x1 = ax_x(&ax, 1.0);
    for (float y = ax.r.y; y < ax.r.y + ax.r.h; y += 6 * S) line(x1, y, x1, y + 3 * S, alpha(C_TEXT, 0.5));
    if (a->probe.hit.kind == HIT_DISK) {
        float xp = ax_x(&ax, clampd(a->probe.g, HIST_G0, HIST_G1));
        line(xp, ax.r.y, xp, ax.r.y + ax.r.h, C_ACCENT);
    }
    axes_frame(&ax);
    for (double g = 0.5; g <= 1.51; g += 0.25) {
        char lab[16];
        snprintf(lab, sizeof lab, "%.2f", g);
        xtick(&ax, g, lab);
    }
}

/* R(r) / r^4 for the probe ray: it moves only where this is positive, and
   turns around where it touches zero. Tends to 1 far away. */
static void plot_potential(const App *a, SDL_FRect cell) {
    const Probe *p = &a->probe;
    const double r0 = floor(K.r_plus * 2) / 2, r1 = 14, y0 = -1.0, y1 = 1.15;
    Axes ax = plot_cell(cell, "RADIAL POTENTIAL", "R(r) / r⁴", r0, r1, y0, y1, 14 * S);
    SDL_Rect clip = {(int)ax.r.x, (int)ax.r.y, (int)ax.r.w + 1, (int)ax.r.h + 1};

    /* photon shell and horizon */
    float xa = ax_x(&ax, K.r_ph_min), xb = ax_x(&ax, K.r_ph_max), xh = ax_x(&ax, K.r_plus);
    fill((SDL_FRect){xa, ax.r.y, fmaxf(xb - xa, S), ax.r.h}, alpha(C_TEXT, 0.06));
    text(&F_SMALL, xa + 4 * S, ax.r.y, C_LABEL, "photon shell");
    line(xh, ax.r.y, xh, ax.r.y + ax.r.h, alpha(C_TEXT, 0.6));
    float yz = ax_y(&ax, 0);
    for (float x = ax.r.x; x < ax.r.x + ax.r.w; x += 6 * S) line(x, yz, x + 3 * S, yz, alpha(C_TEXT, 0.4));

    SDL_RenderSetClipRect(R, &clip);
    float px = 0, py = 0;
    for (float x = ax.r.x; x <= ax.r.x + ax.r.w; x += 1) {
        double r = lerp(r0, r1, (x - ax.r.x) / ax.r.w);
        double v = radial_potential(r, p->ray.l, p->ray.eta) / (r * r * r * r);
        float y = ax_y(&ax, clampd(v, y0 - 0.1, y1 + 0.1));
        if (v < 0 && r > K.r_plus) line(x, yz, x, y, alpha(C_WARN, 0.14));
        if (x > ax.r.x) thick_line(px, py, x, y, 2 * S, C_TEXT);
        px = x;
        py = y;
    }
    SDL_RenderSetClipRect(R, NULL);
    if (p->hit.kind != HIT_HORIZON && p->r_min < r1 && p->r_min > r0)
        disc(ax_x(&ax, p->r_min), yz, 3.5f * S, C_ACCENT);
    axes_frame(&ax);
    char lab[16];
    snprintf(lab, sizeof lab, "%.1f", K.r_plus);
    xtick(&ax, K.r_plus, lab);
    xtick(&ax, 6, "6");
    xtick(&ax, 10, "10");
    xtick(&ax, 14, "14 M");
}

static void draw_plots(const App *a, SDL_FRect r) {
    SDL_FRect inner = inset(r, PAD * S, 14 * S, PAD * S, 12 * S);
    float gap = 34 * S, sq = inner.h;
    SDL_FRect ray = {inner.x, inner.y, sq, inner.h};
    float rest = inner.w - sq - gap, w3 = (rest - 2 * gap) / 3;
    plot_ray(a, ray);
    float x = inner.x + sq + gap;
    plot_temperature(a, (SDL_FRect){x + 28 * S, inner.y, w3 - 28 * S, inner.h});
    x += w3 + gap;
    plot_histogram(a, (SDL_FRect){x, inner.y, w3, inner.h});
    x += w3 + gap;
    plot_potential(a, (SDL_FRect){x, inner.y, w3, inner.h});
}

/* ---------- panel ---------- */

typedef struct { float x, y, w; } Pen;

static void section(Pen *p, const char *title) {
    p->y += 16 * S;
    text(&F_SMALL, p->x, p->y, C_ACCENT, title);
    p->y += line_h(&F_SMALL) + 5 * S;
}

static void row(Pen *p, const char *label, Color vc, const char *unit, const char *fmt, ...) {
    char val[128];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(val, sizeof val, fmt, ap);
    va_end(ap);
    float right = p->x + p->w;
    text(&F_MAIN, p->x, p->y, C_LABEL, label);
    if (unit && *unit) {
        float uw = text_w(&F_MAIN, unit);
        text(&F_MAIN, right - uw, p->y, C_LABEL, unit);
        right -= uw + F_MAIN.adv;
    }
    text_right(&F_MAIN, right, p->y, vc, val);
    p->y += ROW_H * S;
}

static void colorbar(Pen *p, const App *a) {
    float x0 = p->x, w = p->w, h = 8 * S, y = p->y + 2 * S;
    if (a->mode == MODE_ORDER) {
        static const char *names[3] = {"n=0 direct", "n=1 secondary", "n≥2 ring"};
        float cx = x0;
        for (int k = 0; k < 3; k++) {
            fill((SDL_FRect){cx, y, 10 * S, h}, rgba(ORDER_COLORS[k], 255));
            float end = text(&F_SMALL, cx + 15 * S, y - 2 * S, C_TEXT, names[k]);
            cx = end + 14 * S;
        }
        p->y += 22 * S;
        if (a->stats.disk_px) {
            const Stats *s = &a->stats;
            row(p, "share of disk", C_TEXT, "%", "%.1f / %.2f / %.3f", 100.0 * s->order_px[0] / s->disk_px,
                100.0 * s->order_px[1] / s->disk_px, 100.0 * s->order_px[2] / s->disk_px);
        }
        return;
    }
    bool temp = a->mode == MODE_COLOR;
    double t0 = log(1500), t1 = log(15000);
    for (int i = 0; i < (int)w; i++) {
        double f = i / (w - 1);
        vec3 c = temp ? bb_display(exp(lerp(t0, t1, f))) : redshift_color(lerp(0.5, 1.5, f));
        line(x0 + i, y, x0 + i, y + h, rgba(c, 255));
    }
    /* markers for the range seen in this frame */
    const Stats *s = &a->stats;
    if (s->disk_px) {
        double lo = temp ? (log(s->tmin) - t0) / (t1 - t0) : (s->gmin - 0.5);
        double hi = temp ? (log(s->tmax) - t0) / (t1 - t0) : (s->gmax - 0.5);
        float xl = x0 + (float)clampd(lo, 0, 1) * w, xh = x0 + (float)clampd(hi, 0, 1) * w;
        line(xl, y + h + 2 * S, xh, y + h + 2 * S, C_TEXT);
        line(xl, y + h + S, xl, y + h + 5 * S, C_TEXT);
        line(xh, y + h + S, xh, y + h + 5 * S, C_TEXT);
    }
    float ty = y + h + 7 * S;
    if (temp) {
        static const double ticks[] = {2000, 3000, 5000, 7000, 10000, 15000};
        for (int k = 0; k < 6; k++) {
            float tx = x0 + (float)((log(ticks[k]) - t0) / (t1 - t0)) * w;
            char lab[16];
            snprintf(lab, sizeof lab, "%gk", ticks[k] / 1000);
            if (k == 5) text_right(&F_SMALL, tx, ty, C_LABEL, "15k K");
            else text_center(&F_SMALL, fmaxf(tx, x0 + text_w(&F_SMALL, lab) / 2), ty, C_LABEL, lab);
        }
    } else {
        for (int k = 0; k <= 4; k++) {
            float tx = x0 + k / 4.0f * w;
            char lab[16];
            snprintf(lab, sizeof lab, "%.2f", 0.5 + k * 0.25);
            if (k == 0) text(&F_SMALL, tx, ty, C_LABEL, lab);
            else if (k == 4) text_right(&F_SMALL, tx, ty, C_LABEL, lab);
            else text_center(&F_SMALL, tx, ty, C_LABEL, lab);
        }
    }
    p->y += 30 * S;
}

/* Planck B_lambda up to a constant factor. */
static double planck_lambda(double nm, double T) {
    double lam = nm * 1e-9;
    return 1.0 / (pow(lam * 1e6, 5) * expm1(1.438776877e-2 / (lam * T)));
}

/* Emitted and observed spectra of the probed disk point over the visible band. */
static void spectrum(const App *a, SDL_FRect r) {
    const Probe *pr = &a->probe;
    const double l0 = 200, l1 = 1600;
    Axes ax = {inset(r, 0, 0, 0, line_h(&F_SMALL) + 16 * S), l0, l1, 0, 1.08};
    /* visible band strip */
    float sy = ax.r.y + ax.r.h + 3 * S;
    for (float x = ax_x(&ax, 380); x < ax_x(&ax, 750); x += 1) {
        double nm = lerp(l0, l1, (x - ax.r.x) / ax.r.w);
        vec3 c = xyz_to_srgb(cmf(nm));
        c = v3(fmax(c.x, 0), fmax(c.y, 0), fmax(c.z, 0));
        double m = fmax(c.x, fmax(c.y, c.z));
        if (m > 0) c = mul(c, 1 / m);
        double edge = fmin(1, fmin(nm - 380, 750 - nm) / 40);
        line(x, sy, x, sy + 5 * S, rgba(mul(v3(srgb_oetf(c.x), srgb_oetf(c.y), srgb_oetf(c.z)), edge), 255));
    }
    if (pr->hit.kind != HIT_DISK || pr->t_emit < 300) {
        axes_frame(&ax);
        text_center(&F_SMALL, ax.r.x + ax.r.w / 2, ax.r.y + ax.r.h / 2 - line_h(&F_SMALL) / 2, C_LABEL,
                    "no disk emission on this ray");
    } else {
        double peak = 0;
        for (double nm = l0; nm <= l1; nm += 5)
            peak = fmax(peak, fmax(planck_lambda(nm, pr->t_emit), planck_lambda(nm, pr->t_obs)));
        /* emitted (dashed) and observed (solid) */
        for (int k = 0; k < 2; k++) {
            double T = k == 0 ? pr->t_emit : pr->t_obs;
            float px = 0, py = 0;
            for (float x = ax.r.x; x <= ax.r.x + ax.r.w; x += 1) {
                double nm = lerp(l0, l1, (x - ax.r.x) / ax.r.w);
                float y = ax_y(&ax, planck_lambda(nm, T) / peak);
                if (x > ax.r.x) {
                    if (k == 1) thick_line(px, py, x, y, 2 * S, C_TEXT);
                    else if ((int)((x - ax.r.x) / (4 * S)) % 2 == 0) line(px, py, x, y, C_LABEL);
                }
                px = x;
                py = y;
            }
        }
        double wien = 2.897771955e6 / pr->t_obs;
        if (wien > l0 && wien < l1) {
            float xw = ax_x(&ax, wien);
            for (float y = ax.r.y; y < ax.r.y + ax.r.h; y += 6 * S) line(xw, y, xw, y + 3 * S, C_ACCENT);
        }
        axes_frame(&ax);
        char note[64];
        snprintf(note, sizeof note, "λ peak %.0f nm", wien);
        text_right(&F_SMALL, ax.r.x + ax.r.w, ax.r.y, C_ACCENT, note);
        text_right(&F_SMALL, ax.r.x + ax.r.w, ax.r.y + line_h(&F_SMALL), C_LABEL, "- - emitted");
    }
    float ty = sy + 9 * S;
    static const double ticks[] = {400, 700, 1000, 1500};
    for (int k = 0; k < 4; k++) {
        char lab[16];
        snprintf(lab, sizeof lab, "%.0f", ticks[k]);
        text_center(&F_SMALL, ax_x(&ax, ticks[k]), ty, C_LABEL, lab);
    }
    text(&F_SMALL, ax.r.x, ty, C_LABEL, "nm");
}

static void draw_panel(const App *a, SDL_FRect r) {
    const View *v = &a->view;
    const Probe *pr = &a->probe;
    const Stats *st = &a->stats;
    Pen p = {r.x + PAD * S, r.y + 16 * S, r.w - 2 * PAD * S};
    Camera cam = view_camera(a);
    char label[48];

    text(&F_TITLE, p.x, p.y, C_TEXT, v->spin == 0 ? "SCHWARZSCHILD" : "KERR");
    text_right(&F_SMALL, p.x + p.w, p.y + (line_h(&F_TITLE) - line_h(&F_SMALL)) / 2 + S, C_LABEL,
               "G = c = M = 1");
    p.y += line_h(&F_TITLE);

    section(&p, "OBSERVER  ZAMO");
    row(&p, "distance r", C_TEXT, "M", "%.2f", v->r_cam);
    row(&p, "inclination / azimuth", C_TEXT, "°", "%.1f / %.1f", v->incl_deg,
        fmod(fmod(v->azim_deg, 360) + 360, 360));
    row(&p, "field of view", C_TEXT, "°", "%.1f", v->fov_deg);
    row(&p, "frame dragging ω", C_TEXT, "/M", "%.2e", cam.omega);
    row(&p, "clock rate dτ/dt", C_TEXT, "", "%.4f", cam.lapse);

    section(&p, "SPACETIME");
    row(&p, "spin a", C_TEXT, "M", "%+.3f", v->spin);
    row(&p, "horizon", C_TEXT, "M", "%.3f", K.r_plus);
    row(&p, "photon orbits", C_TEXT, "M", "%.3f / %.3f", K.r_ph_min, K.r_ph_max);
    row(&p, "ISCO", C_TEXT, "M", "%.3f", K.r_isco);
    row(&p, "shadow w × h", C_TEXT, "°", "%.3f × %.3f", a->edge.width_deg, a->edge.height_deg);
    row(&p, "shadow shift", C_TEXT, "°", "%+.3f", a->edge.shift_deg);

    section(&p, "DISK  Novikov-Thorne");
    row(&p, "inner / outer", C_TEXT, "M", "%.2f / %.1f", K.r_isco, v->r_out);
    snprintf(label, sizeof label, "peak T @ %.2f M", nt_rpeak);
    row(&p, label, C_TEXT, "K", "%.0f", v->t_peak);
    row(&p, "v orbit at ISCO", C_TEXT, "c", "%.3f", disk_speed(K.r_isco));
    row(&p, "efficiency 1 - E_isco", C_TEXT, "%", "%.2f", 100 * disk_efficiency());
    if (st->disk_px) {
        row(&p, "g in view", C_TEXT, "", "%.3f / %.3f", st->gmin, st->gmax);
        row(&p, "T seen, 1-99%", C_TEXT, "K", "%.0f / %.0f", st->tmin, st->tmax);
    } else {
        row(&p, "g in view", C_LABEL, "", "none");
        row(&p, "T seen, 1-99%", C_LABEL, "", "none");
    }

    section(&p, "PROBE  cursor ray");
    row(&p, "L/E, Q/E²", C_TEXT, "", "%.4f, %.4f", pr->ray.l, pr->ray.eta);
    if (pr->hit.kind == HIT_DISK) {
        const char *img = pr->hit.order == 0 ? "direct" : pr->hit.order == 1 ? "secondary" : "photon ring";
        row(&p, "hit", C_TEXT, "", "disk, %s n=%d", img, pr->hit.order);
        row(&p, "radius / azimuth", C_TEXT, "", "%.3f M / %.1f°", pr->r, pr->phi);
        row(&p, "redshift g", pr->g >= 1 ? rgba(redshift_color(1.35), 255) : rgba(redshift_color(0.65), 255),
            "", "%.4f", pr->g);
        row(&p, "T emit → seen", C_TEXT, "K", "%.0f → %.0f", pr->t_emit, pr->t_obs);
    } else if (pr->hit.kind == HIT_SKY) {
        row(&p, "hit", C_TEXT, "", "sky, crossed disk plane %d×", pr->hit.order);
        row(&p, "deflection", C_TEXT, "°", "%.3f", pr->deflection);
        double ra_h = pr->ra * 12 / M_PI;
        row(&p, "sky RA / Dec", C_TEXT, "", "%02dh%04.1fm / %+.2f°", (int)ra_h, fmod(ra_h, 1) * 60,
            pr->dec * 180 / M_PI);
        row(&p, "closest approach", C_TEXT, "M", "%.3f", pr->r_min);
        row(&p, "", C_TEXT, "", "");
    } else {
        row(&p, "hit", C_TEXT, "", "horizon, captured");
        row(&p, "", C_TEXT, "", "");
        row(&p, "", C_TEXT, "", "");
        row(&p, "", C_TEXT, "", "");
    }

    section(&p, "RENDER");
    const Pass *ps = &PASSES[a->img_pass];
    if (a->pass < NPASSES)
        row(&p, "pass", C_WARN, "", "%d / %d  %3.0f%%", a->pass + 1, NPASSES, 100 * a->progress);
    else
        row(&p, "pass", C_TEXT, "", "%d / %d  done", NPASSES, NPASSES);
    row(&p, "resolution", C_TEXT, "", "%d × %d  %d spp", a->img_w, a->img_h, ps->ss * ps->ss);
    row(&p, "time, geodesics/s", C_TEXT, "", "%.2f s, %.2f M", a->pass_time,
        a->pass_time > 0 ? a->img_w * a->img_h * ps->ss * ps->ss / a->pass_time / 1e6 : 0);
    row(&p, "integrator", C_TEXT, "", "RK4, Mino time, h %.2f", ps->step);

    section(&p, a->mode == MODE_COLOR ? "COLOR  observed blackbody" : a->mode == MODE_REDSHIFT
                                                                      ? "COLOR  g = ν obs / ν emit"
                                                                      : "COLOR  image order");
    colorbar(&p, a);

    section(&p, "SPECTRUM  probe ray");
    spectrum(a, (SDL_FRect){p.x, p.y, p.w, 118 * S});

    if (a->saved[0]) text(&F_SMALL, p.x, r.y + r.h - PAD * S - line_h(&F_SMALL), C_ACCENT, a->saved);
}

/* ---------- frame ---------- */

static void draw(App *a) {
    SDL_SetRenderDrawColor(R, 0, 0, 0, 255);
    SDL_RenderClear(R);
    SDL_SetRenderDrawBlendMode(R, SDL_BLENDMODE_BLEND);
    Layout L = layout(a);
    if (a->tex) SDL_RenderCopyF(R, a->tex, NULL, &L.view);
    draw_view_overlay(a, L.view);
    if (a->hud) {
        line(L.panel.x, 0, L.panel.x, L.panel.h, C_RULE);
        line(0, L.plots.y, L.plots.w, L.plots.y, C_RULE);
        draw_panel(a, L.panel);
        draw_plots(a, L.plots);
    }
    SDL_RenderPresent(R);
    a->last_draw = now_s();
}

/* ---------- input ---------- */

static void handle(App *a, const SDL_Event *e) {
    View *v = &a->view;
    switch (e->type) {
    case SDL_QUIT: a->quit = true; break;
    case SDL_WINDOWEVENT:
        if (e->window.event == SDL_WINDOWEVENT_SIZE_CHANGED) a->view_dirty = true;
        if (e->window.event == SDL_WINDOWEVENT_EXPOSED) a->ui_dirty = true;
        break;
    case SDL_MOUSEBUTTONDOWN: {
        Layout L = layout(a);
        float mx = e->button.x * S, my = e->button.y * S;
        if (e->button.button == SDL_BUTTON_LEFT && mx < L.view.x + L.view.w && my < L.view.y + L.view.h)
            a->dragging = true;
        break;
    }
    case SDL_MOUSEBUTTONUP:
        if (e->button.button == SDL_BUTTON_LEFT) a->dragging = false;
        break;
    case SDL_MOUSEMOTION:
        if (a->dragging) {
            v->azim_deg -= 0.25 * e->motion.xrel;
            v->incl_deg = clampd(v->incl_deg - 0.25 * e->motion.yrel, 0.5, 179.5);
            a->view_dirty = true;
        } else {
            Layout L = layout(a);
            float mx = e->motion.x * S, my = e->motion.y * S;
            if (mx >= L.view.x && my >= L.view.y && mx < L.view.x + L.view.w && my < L.view.y + L.view.h) {
                a->probe.u = (mx - L.view.x) / L.view.w;
                a->probe.v = (my - L.view.y) / L.view.h;
                a->probe_dirty = true;
            }
        }
        break;
    case SDL_MOUSEWHEEL: {
        double dy = e->wheel.preciseY;
        if (e->wheel.direction == SDL_MOUSEWHEEL_FLIPPED) dy = -dy;
        if (dy != 0) {
            v->r_cam = clampd(v->r_cam * pow(0.93, dy), 3.5, 2000);
            a->view_dirty = true;
        }
        break;
    }
    case SDL_KEYDOWN:
        switch (e->key.keysym.sym) {
        case SDLK_q: case SDLK_ESCAPE: a->quit = true; break;
        case SDLK_1: case SDLK_2: case SDLK_3:
            a->mode = e->key.keysym.sym - SDLK_1;
            a->compose_dirty = true;
            break;
        case SDLK_g: a->grid = !a->grid; a->compose_dirty = true; break;
        case SDLK_o: /* back off until Orion fits around the Einstein ring */
            *v = DEFAULT_VIEW;
            v->r_cam = 300;
            v->fov_deg = 40;
            a->orbit = false;
            a->view_dirty = true;
            break;
        case SDLK_SLASH: case SDLK_k: a->keys = !a->keys; a->ui_dirty = true; break;
        case SDLK_h: a->hud = !a->hud; a->view_dirty = true; break;
        case SDLK_SPACE: a->orbit = !a->orbit; a->view_dirty = true; break;
        case SDLK_r:
            *v = DEFAULT_VIEW;
            a->orbit = false;
            a->view_dirty = true;
            break;
        case SDLK_COMMA: case SDLK_PERIOD: {
            int k = 0;
            for (int i = 1; i < NSPINS; i++)
                if (fabs(SPINS[i] - v->spin) < fabs(SPINS[k] - v->spin)) k = i;
            k += e->key.keysym.sym == SDLK_PERIOD ? 1 : -1;
            v->spin = SPINS[k < 0 ? 0 : k >= NSPINS ? NSPINS - 1 : k];
            a->view_dirty = true;
            break;
        }
        case SDLK_LEFT: v->azim_deg += 5; a->view_dirty = true; break;
        case SDLK_RIGHT: v->azim_deg -= 5; a->view_dirty = true; break;
        case SDLK_UP: v->incl_deg = clampd(v->incl_deg - 5, 0.5, 179.5); a->view_dirty = true; break;
        case SDLK_DOWN: v->incl_deg = clampd(v->incl_deg + 5, 0.5, 179.5); a->view_dirty = true; break;
        case SDLK_LEFTBRACKET: case SDLK_RIGHTBRACKET:
            v->t_peak = clampd(v->t_peak * (e->key.keysym.sym == SDLK_RIGHTBRACKET ? 1.1 : 1 / 1.1), 1500, 100000);
            a->view_dirty = true;
            break;
        case SDLK_MINUS: v->fov_deg = clampd(v->fov_deg * 1.1, 2, 150); a->view_dirty = true; break;
        case SDLK_EQUALS: v->fov_deg = clampd(v->fov_deg / 1.1, 2, 150); a->view_dirty = true; break;
        case SDLK_s: {
            char name[96];
            time_t t = time(NULL);
            strftime(name, sizeof name, "blackhole-%Y%m%d-%H%M%S.png", localtime(&t));
            if (a->img_w && write_png(name, a->rgb, a->img_w, a->img_h))
                snprintf(a->saved, sizeof a->saved, "saved %s", name);
            else
                snprintf(a->saved, sizeof a->saved, "could not save %s", name);
            a->ui_dirty = true;
            break;
        }
        default: break;
        }
        break;
    default: break;
    }
}

/* ---------- rendering passes ---------- */

static void finish_pass(App *a, int w, int h) {
    const Pass *p = &PASSES[a->pass];
    float *ti = a->front_img;
    a->front_img = a->work_img;
    a->work_img = ti;
    if (p->ss == 1) { /* a supersampled pass reuses the previous pass's aux at the same size */
        Aux *ta = a->front_aux;
        a->front_aux = a->work_aux;
        a->work_aux = ta;
    }
    a->img_w = w;
    a->img_h = h;
    a->img_pass = a->pass;
    a->pass_time = now_s() - a->pass_started;
    compose(a);
    a->pass = a->orbit ? 0 : a->pass + 1;
    a->ui_dirty = true;
}

static void run_pass(App *a) {
    int w, h;
    pass_size(a, a->pass, &w, &h);
    ensure_capacity(a, (size_t)w * h);
    const Pass *p = &PASSES[a->pass];
    physics_sync(&a->view);
    Camera cam = camera_make(a->view, w, h, p->step, p->step_max);
    a->pass_started = now_s();
    a->progress = 0;
    for (int b = 0; b < BANDS; b++) {
        render_rows(&cam, a->work_img, p->ss == 1 ? a->work_aux : NULL, w, h, p->ss, h * b / BANDS,
                    h * (b + 1) / BANDS, a->threads);
        a->progress = (b + 1) / (float)BANDS;
        if (a->pass == 0) continue;
        SDL_Event e;
        while (SDL_PollEvent(&e)) handle(a, &e);
        if (a->view_dirty || a->quit) return;
        if (a->compose_dirty) {
            compose(a);
            a->compose_dirty = false;
            a->ui_dirty = true;
        }
        if (a->probe_dirty) {
            probe_update(a);
            a->probe_dirty = false;
            a->ui_dirty = true;
        }
        if (a->ui_dirty || now_s() - a->last_draw > 0.15) {
            draw(a);
            a->ui_dirty = false;
        }
    }
    finish_pass(a, w, h);
}

int main(void) {
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL_Init: %s\n", SDL_GetError());
        return 1;
    }
    SDL_Window *win = SDL_CreateWindow("Black hole", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1640, 940,
                                       SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
    R = win ? SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC) : NULL;
    if (!win || !R) {
        fprintf(stderr, "SDL window: %s\n", SDL_GetError());
        return 1;
    }
    SDL_SetWindowMinimumSize(win, 1100, 720);
    int ww, wh, pw, ph;
    SDL_GetWindowSize(win, &ww, &wh);
    SDL_GetRendererOutputSize(R, &pw, &ph);
    S = (float)pw / ww;

    glyphs_init();
    if (!font_load(&F_MAIN, 11.5 * S) || !font_load(&F_SMALL, 9.5 * S) || !font_load(&F_TITLE, 14 * S)) {
        fprintf(stderr, "could not load %s\n", FONT_PATH);
        return 1;
    }

    App a = {.view = DEFAULT_VIEW, .grid = true, .hud = true, .view_dirty = true, .win = win};
    a.threads = (int)sysconf(_SC_NPROCESSORS_ONLN);
    if (a.threads < 1) a.threads = 1;
    a.probe.u = 0.5;
    a.probe.v = 0.3;
    a.probe.path.cap = PATH_CAP;
    a.probe.path.pts = malloc(sizeof(vec3) * PATH_CAP);
    bb_init();
    stars_init();
    named_init();
    physics_sync(&a.view);

    double last_tick = now_s();
    while (!a.quit) {
        SDL_Event e;
        bool idle = a.pass >= NPASSES && !a.orbit && !a.ui_dirty && !a.view_dirty && !a.compose_dirty &&
                    !a.probe_dirty;
        if (idle && SDL_WaitEvent(&e)) handle(&a, &e);
        while (SDL_PollEvent(&e)) handle(&a, &e);

        double t = now_s();
        if (a.orbit) {
            a.view.azim_deg -= 12.0 * (t - last_tick);
            a.view_dirty = true;
        }
        last_tick = t;

        if (a.view_dirty) {
            a.view_dirty = false;
            a.pass = 0;
            a.probe_dirty = true;
            edge_update(&a);
        }
        if (a.probe_dirty) {
            probe_update(&a);
            a.probe_dirty = false;
            a.ui_dirty = true;
        }
        if (a.compose_dirty && a.img_w) {
            compose(&a);
            a.compose_dirty = false;
            a.ui_dirty = true;
        }
        if (a.pass < NPASSES && !a.quit) run_pass(&a);
        if (a.ui_dirty && !a.quit) {
            draw(&a);
            a.ui_dirty = false;
        }
    }

    if (a.tex) SDL_DestroyTexture(a.tex);
    SDL_DestroyRenderer(R);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
