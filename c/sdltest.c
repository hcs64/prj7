#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "SDL.h"

#include "render.h"
#include "teapot.h"

#ifdef EMSCRIPTEN
#include "emscripten.h"
#endif

#define FPS 60

void one_iter();

const int W = 1024;
const int H = 768;

// globals, the scene
SDL_Surface *screen;
double *zbuf;
Space * world;

Geometry * cube;
Geometry * cylinder;
Geometry * cylinder_orbiter;
Geometry * cylinder_orbiter2;
Geometry * sphere;
Geometry * teapot;

int main(int argc, char *argv[])
{
    if ( SDL_Init(SDL_INIT_VIDEO) < 0 ) {
        exit(1);
    }

    atexit(SDL_Quit);

    screen = SDL_SetVideoMode(W, H, 32, SDL_HWSURFACE|SDL_DOUBLEBUF);
    if (screen == NULL) {
        exit(1);
    }

    printf("ok!\n");

    zbuf = malloc(W*H*sizeof(double));
    world = makeSpace(4);

    sphere = makeGeometry(0);
    makeSphere(sphere);

    cube = makeGeometry(0);
    makeCube(cube);

    cylinder = makeGeometry(1);
    makeCylinder(cylinder);

    cylinder_orbiter = makeGeometry(1);
    makeSphere(cylinder_orbiter);

    cylinder_orbiter2 = makeGeometry(0);
    makeSphere(cylinder_orbiter2);

    teapot = makeGeometry(0);
    makeTeapot(teapot);

    addChild(&cylinder->space, &cylinder_orbiter->space);
    addChild(&cylinder_orbiter->space, &cylinder_orbiter2->space);

    addChild(world, &sphere->space);
    addChild(world, &cube->space);
    addChild(world, &cylinder->space);
    addChild(world, &teapot->space);

#ifdef EMSCRIPTEN
    emscripten_set_main_loop(one_iter, FPS, 1);
#else
    while (1) {
        one_iter();

        SDL_Delay(1000/FPS);
    }
#endif

    printf("done\n");

}

void one_iter() {
    SDL_Event event;

    while ( SDL_PollEvent(&event) ) {
        switch (event.type) {
            case SDL_KEYDOWN:
                break;
            case SDL_QUIT:
                exit(0);
                break;
            default:
                break;
        }
    }

    int start_t = SDL_GetTicks();
    double t = start_t/1000.;
    void * m;
    
    m = sphere->space.m;
    makeIdentity(m);
    translate(m,0,-2+sin(t*1.4),0);
    scaleMatrix(m,1,sin(t*1.4)+1.1,1);
    rotateY(m,t);

    m = cylinder->space.m;
    makeIdentity(m);
    rotateY(m,1);
    rotateX(m,t);
    translate(m,6,-2,-3);

    m = cube->space.m;
    makeIdentity(m);
    scaleMatrix(m,10,1,1);
    rotateX(m,t);
    rotateY(m,1);
    translate(m,0,-2,-10);

    // a little sphere orbiting the cylinder
    m = cylinder_orbiter->space.m;
    makeIdentity(m);
    scaleMatrix(m,.5,.5,.5);
    translate(m,1.5,0,0);
    rotateZ(m,t);
    
    // an orbiter for the orbiter
    m = cylinder_orbiter2->space.m;
    makeIdentity(m);
    scaleMatrix(m,.5,.5,.5);
    translate(m,0,0,-1.5);
    rotateX(m,t*5);
    
    m = teapot->space.m;
    makeIdentity(m);
    rotateY(m, -M_PI/2);
    translate(m, 3,-3,0);
    rotateY(m,-t);
    scaleMatrix(m,2,2,2);

    makeIdentity(world->m);
    world->update(world);

    // clear before locking!
    cls(screen, zbuf, 10);

    if ( SDL_MUSTLOCK(screen) ) {
        if ( SDL_LockSurface(screen) < 0) {
            return;
        }
    }

    fill(cube, screen, zbuf, 10);
    fill(cylinder, screen, zbuf, 10);
    fill(cylinder_orbiter, screen, zbuf, 10);
    fill(cylinder_orbiter2, screen, zbuf, 10);
    fill(sphere, screen, zbuf, 10);
    fill(teapot, screen, zbuf, 10);

    if ( SDL_MUSTLOCK(screen) ) {
        SDL_UnlockSurface(screen);
    }

    SDL_Flip(screen);

    int end_t = SDL_GetTicks();

    printf("%d\n", end_t-start_t);
}
