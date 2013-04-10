#ifndef _RENDER_H
#define _RENDER_H

#include <stdint.h>
#include "SDL.h"

// 4x4 transformation matrix
// m[row][col]

typedef double Matrix[4][4];
typedef double Vector[3];
typedef int Face[4];

void clearMatrix(Matrix m);
void makeIdentity(Matrix m);
void copyMatrix(Matrix dst, Matrix src);
void invertMatrixFrom(Matrix dst, Matrix src);
void transposeMatrix(Matrix m);
void translate(Matrix m, double x, double y, double z);
void rotateX(Matrix m, double radians);
void rotateY(Matrix m, double radians);
void rotateZ(Matrix m, double radians);
void scaleMatrix(Matrix m, double x, double y, double z);
void leftMultiply(Matrix m1, Matrix m2);
void rightMultiply(Matrix m1, Matrix m2);
void transform(Matrix m, Vector src, Vector dst);
void transformNorm(Matrix m, Vector src, Vector dst);

// a space with children and a transformation matrix

typedef struct Space {
    Matrix m, globalM, globalNM;
    struct Space ** children;
    int nChildren;
    int maxChildren;
    void (*update)(struct Space *);
    void (*updateWithParent)(struct Space *, Matrix);
} Space;

Space * makeSpace(int maxChildren);
void addChild(Space * s, Space * child);
Space * getChild(Space * s, int childIdx);
int getNumChildren(Space * s);

typedef struct Geometry {
    struct Space space;

    int nv;
    Vector * v; // 
    Vector * tv;// transformed vectors
    Vector * n; // vertex normals
    Vector * tn;// transformed normals

    int nf;
    Face * f;   // faces
} Geometry;

Geometry * makeGeometry(int maxChildren);
void setupInitialGeometry(Geometry * g, int nv, Vector * vertices, Vector * normals, int nf, Face * faces);

void makeSphere(Geometry * g);
void makeCylinder(Geometry * g);
void makeCube(Geometry * g);

// rendering
void cls(SDL_Surface * screen, double * zbuf, double FL);
void fill(Geometry * geo, SDL_Surface * screen, double * zbuf, double FL);

#endif
