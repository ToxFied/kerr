/* png.h: minimal PNG writer (stored deflate blocks, no dependencies). */
#ifndef PNG_H
#define PNG_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint32_t crc_table[256];
static void crc_init(void) {
    for (uint32_t n = 0; n < 256; n++) {
        uint32_t c = n;
        for (int k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320u ^ (c >> 1) : c >> 1;
        crc_table[n] = c;
    }
}
static uint32_t crc(uint32_t c, const uint8_t *buf, size_t n) {
    for (size_t i = 0; i < n; i++) c = crc_table[(c ^ buf[i]) & 0xFF] ^ (c >> 8);
    return c;
}
static void put32(uint8_t *p, uint32_t v) { p[0] = v >> 24; p[1] = v >> 16; p[2] = v >> 8; p[3] = v; }
static void chunk(FILE *f, const char *type, const uint8_t *data, uint32_t n) {
    uint8_t hdr[8];
    put32(hdr, n);
    memcpy(hdr + 4, type, 4);
    fwrite(hdr, 1, 8, f);
    if (n) fwrite(data, 1, n, f);
    uint32_t c = crc(0xFFFFFFFFu, hdr + 4, 4);
    c = crc(c, data, n) ^ 0xFFFFFFFFu;
    uint8_t tail[4];
    put32(tail, c);
    fwrite(tail, 1, 4, f);
}

/* Writes 8-bit RGB pixels. Returns 0 on failure. */
static int write_png(const char *path, const uint8_t *rgb, int w, int h) {
    crc_init();
    size_t raw_n = (size_t)h * (1 + 3 * (size_t)w);
    uint8_t *raw = malloc(raw_n);
    for (int y = 0; y < h; y++) {
        raw[y * (1 + 3 * (size_t)w)] = 0;
        memcpy(raw + y * (1 + 3 * (size_t)w) + 1, rgb + (size_t)y * 3 * w, 3 * (size_t)w);
    }
    size_t blocks = (raw_n + 65534) / 65535;
    size_t z_n = 2 + raw_n + 5 * blocks + 4;
    uint8_t *z = malloc(z_n), *p = z;
    *p++ = 0x78; *p++ = 0x01;
    uint32_t a = 1, b = 0;
    for (size_t off = 0; off < raw_n; off += 65535) {
        uint16_t n = (uint16_t)(raw_n - off < 65535 ? raw_n - off : 65535);
        *p++ = off + n >= raw_n;
        *p++ = n & 0xFF; *p++ = n >> 8; *p++ = ~n & 0xFF; *p++ = (uint16_t)~n >> 8;
        memcpy(p, raw + off, n);
        p += n;
    }
    for (size_t i = 0; i < raw_n; i++) { a = (a + raw[i]) % 65521; b = (b + a) % 65521; }
    put32(p, b << 16 | a);

    FILE *f = fopen(path, "wb");
    if (!f) { free(raw); free(z); return 0; }
    fwrite("\x89PNG\r\n\x1a\n", 1, 8, f);
    uint8_t ihdr[13];
    put32(ihdr, w); put32(ihdr + 4, h);
    ihdr[8] = 8; ihdr[9] = 2; ihdr[10] = ihdr[11] = ihdr[12] = 0;
    chunk(f, "IHDR", ihdr, 13);
    chunk(f, "IDAT", z, (uint32_t)z_n);
    chunk(f, "IEND", NULL, 0);
    fclose(f);
    free(raw);
    free(z);
    return 1;
}

#endif
