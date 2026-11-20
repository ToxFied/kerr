CC ?= cc
CFLAGS ?= -O3 -Wall -Wextra

SDL_CFLAGS := $(shell pkg-config --cflags sdl2 2>/dev/null)
SDL_LIBS := $(shell pkg-config --libs sdl2 2>/dev/null)
MAC_FRAMEWORKS := -framework CoreText -framework CoreGraphics -framework CoreFoundation

HEADERS := src/physics.h src/png.h src/bsc5.h

all: blackhole live

# Still renderer and test suite: needs only libc, libm and pthreads.
blackhole: src/blackhole.c $(HEADERS)
	$(CC) $(CFLAGS) -o $@ src/blackhole.c -lm -lpthread

# Interactive viewer: SDL2 for the window, CoreText (macOS) for text.
live: src/live.c $(HEADERS)
	$(CC) $(CFLAGS) $(SDL_CFLAGS) -o $@ src/live.c $(SDL_LIBS) $(MAC_FRAMEWORKS) -lm -lpthread

check: blackhole
	./blackhole --check

# Rebuilds src/bsc5.h from the Yale Bright Star Catalog at CDS.
catalog: tools/bsc5_to_h.c
	$(CC) -O2 -o tools/bsc5_to_h tools/bsc5_to_h.c -lm
	curl -sSfL https://cdsarc.cds.unistra.fr/ftp/V/50/catalog.gz | gunzip | tools/bsc5_to_h > src/bsc5.h

clean:
	rm -f blackhole live tools/bsc5_to_h

.PHONY: all check catalog clean
