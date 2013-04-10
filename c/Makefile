.PHONY: clean
LDFLAGS=-L/usr/lib/x86_64-linux-gnu -lSDL -lm
CFLAGS=-std=c99 -I/usr/include/SDL -D_GNU_SOURCE=1 -D_REENTRANT
EMCC=~/emscripten/emscripten/emcc
EMCCFLAGS=-O2 -s ASM_JS=1

sdltest: sdltest.c render.c render.h teapot.c teapot.h
	${CC} ${CFLAGS} sdltest.c render.c teapot.c ${LDFLAGS} -o sdltest

sdltest.html: sdltest.c render.c render.h teapot.c teapot.h
	${EMCC} ${EMCCFLAGS} sdltest.c render.c teapot.c -o sdltest.html

clean:
	rm -f sdltest sdltest.html
