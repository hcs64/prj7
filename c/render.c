#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "SDL.h"

#include "render.h"

#if 0
void print_matrix(Matrix m) {
    printf("%f %f %f %f\n", s->globalM[0][0],s->globalM[0][1],s->globalM[0][2],s->globalM[0][3]);
    printf("%f %f %f %f\n", s->globalM[1][0],s->globalM[1][1],s->globalM[1][2],s->globalM[1][3]);
    printf("%f %f %f %f\n", s->globalM[2][0],s->globalM[2][1],s->globalM[2][2],s->globalM[2][3]);
    printf("%f %f %f %f\n", s->globalM[3][0],s->globalM[3][1],s->globalM[3][2],s->globalM[3][3]);
}
#endif

#ifdef EMSCRIPTEN
#define set_pixel(pixel,r,g,b) (((uint8_t*)(pixel))[0]=(r),((uint8_t*)(pixel))[1]=(g),((uint8_t*)(pixel))[2]=(b))
#else
#define set_pixel(pixel,r,g,b) (*(pixel)=SDL_MapRGB(format,(r),(g),(b)))
#endif

void clearMatrix(Matrix m) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            m[i][j] = 0;
        }
    }
}

void makeIdentity(Matrix m) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            m[i][j] = i==j?1:0;
        }
    }
}

void copyMatrix(Matrix dst, Matrix src) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            dst[i][j] = src[i][j];
        }
    }
}

void invertMatrixFrom(Matrix dst, Matrix src) {
       // ported from Ken Perlin's Java implementation
       Matrix tmp, srcT;
       double max[4];
       int index[4];

      int i, j, k, m = 0, n = 4;

      for (i = 0 ; i < n ; i++) {
         index[i] = i;
         max  [i] = 0;
         for (j = 0 ; j < n ; j++) {
            max[i] = fmax(max[i], fabs(src[i][j]));
            tmp[i][j] = i == j ? 1 : 0;
            srcT[i][j] = src[i][j];
         }
      }

      for (j = 0 ; j < n - 1 ; j++) {
         double t = 0;
         for (i = j ; i < n ; i++) {
            double s = fabs(srcT[index[i]][j]) / max[index[i]];
            if (s > t) {
               t = s;
               m = i;
            }
         }
         int swap = index[j];
         index[j] = index[m];
         index[m] = swap;

         for (i = j + 1 ; i < n ; i++) {
            double p = srcT[index[i]][j] / srcT[index[j]][j];
            srcT[index[i]][j] = p;
            for (k = j + 1 ; k < n ; k++)
               srcT[index[i]][k] -= p * srcT[index[j]][k];
         }
      }

      for (i = 0     ; i < n - 1 ; i++)
      for (j = i + 1 ; j < n     ; j++)
      for (k = 0     ; k < n     ; k++)
         tmp[index[j]][k] -= srcT[index[j]][i] * tmp[index[i]][k];

      for (i = 0 ; i < n ; i++) {
         dst[n-1][i] = tmp[index[n-1]][i] / srcT[index[n-1]][n-1];
         for (j = n - 2 ; j >= 0 ; j--) {
            dst[j][i] = tmp[index[j]][i];
            for (k = j + 1 ; k < n ; k++)
               dst[j][i] -= srcT[index[j]][k] * dst[k][i];
            dst[j][i] /= srcT[index[j]][j];
         }
      }
}

void transposeMatrix(Matrix m) {
    // swap across the diagonal
    double swap;

    swap = m[0][1];
    m[0][1] = m[1][0];
    m[1][0] = swap;

    swap = m[0][2];
    m[0][2] = m[2][0];
    m[2][0] = swap;

    swap = m[0][3];
    m[0][3] = m[3][0];
    m[3][0] = swap;

    swap = m[1][2];
    m[1][2] = m[2][1];
    m[2][1] = swap;

    swap = m[1][3];
    m[1][3] = m[3][1];
    m[3][1] = swap;

    swap = m[2][3];
    m[2][3] = m[3][2];
    m[3][2] = swap;
}

void translate(Matrix m, double x, double y, double z) {
    m[0][3] += x;
    m[1][3] += y;
    m[2][3] += z;
}

void rotateX(Matrix m, double radians) {
    double c = cos(radians);
    double s = sin(radians);

    for (int i = 0; i < 4; i++) {
        // from Y (1) into  Z (2)
        // from Z (2) into -Y (0)
        double u = m[1][i];
        double d = m[2][i];

        m[1][i] = u*c + -d*s;
        m[2][i] = u*s +  d*c;
    }
}

void rotateY(Matrix m, double radians) {
    double c = cos(radians);
    double s = sin(radians);

    for (int i = 0; i < 4; i++) {
        // from X (0) into -Z (2)
        // from Z (2) into  X (0)
        double u = m[0][i];
        double d = m[2][i];

        m[0][i] =  u*c + d*s;
        m[2][i] = -u*s + d*c;
    }
}

void rotateZ(Matrix m, double radians) {
    double c = cos(radians);
    double s = sin(radians);

    for (int i = 0; i < 4; i++) {
        // from X (0) into  Y (1)
        // from Y (1) into -X (0)
        double u = m[0][i];
        double d = m[1][i];

        m[0][i] = u*c + -d*s;
        m[1][i] = u*s +  d*c;
    }
}

void scaleMatrix(Matrix m, double x, double y, double z) {
    for (int i = 0; i < 4; i++) {
        m[i][0] *= x;
        m[i][1] *= y;
        m[i][2] *= z;
    }
}

void leftMultiply(Matrix m1, Matrix m2) {
    Matrix t;

    // leftMultiply  means:  m1 <= m2 times m1
    for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 4; i++) {
            double sum = 0;
            for (int k = 0; k < 4; k++) {
                // walk across row of other, down col of this
                sum += m2[i][k] * m1[k][j];
            }

            t[i][j] = sum;
        }
    }

    copyMatrix(m1, t);
}

void rightMultiply(Matrix m1, Matrix m2) {
    Matrix t;

    // rightMultiply means:  m1 <= m1 times m2
    for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 4; i++) {
            double sum = 0;
            for (int k = 0; k < 4; k++) {
                // walk across row of this, down col of other
                sum += m1[i][k] * m2[k][j];
            }
            t[i][j] = sum;
        }
    }

    copyMatrix(m1, t);
}


void transform(Matrix m, Vector src, Vector dst) {
    for (int i = 0; i < 3; i++) {
        double sum = 0;
        for (int j = 0; j < 4; j++) {
            if (j < 3) {
                sum += m[i][j]*src[j];
            } else {
                sum += m[i][j];
            }
        }
        dst[i] = sum;
    }
}

void transformNorm(Matrix m, Vector src, Vector dst) {
    for (int i = 0; i < 3; i++) {
        double sum = 0;
        for (int j = 0; j < 4; j++) {
            if (j < 3) {
                sum += m[i][j]*src[j];
            }
        }
        dst[i] = sum;
    }

    // normalize
    double mag = sqrt(dst[0]*dst[0] + dst[1]*dst[1] + dst[2]*dst[2]);
    dst[0] /= mag;
    dst[1] /= mag;
    dst[2] /= mag;
}

//// Space
void spaceUpdate(Space * s);
void spaceUpdateWithParent(Space * s, Matrix parentM);

Space * makeSpace(int maxChildren) {
    Space * s = malloc(sizeof(Space));
    makeIdentity(s->m);
    clearMatrix(s->globalM);
    clearMatrix(s->globalNM);

    if (maxChildren > 0) {
        s->children = calloc(maxChildren, sizeof(Space *));
        s->maxChildren = maxChildren;
    } else {
        s->children = NULL;
        s->maxChildren = 0;
    }

    s->nChildren = 0;

    s->update = spaceUpdate;
    s->updateWithParent = spaceUpdateWithParent;

    return s;
}

void addChild(Space * s, Space * child) {
    assert(s->nChildren < s->maxChildren);
    s->children[s->nChildren++] = child;
}

Space * getChild(Space * s, int childIdx) {
    return s->children[childIdx];
}

void spaceUpdate(Space * s) {
    copyMatrix(s->globalM, s->m);

    for (int idx = 0; idx < s->nChildren; idx ++) {
        Space * child = s->children[idx];
        child->updateWithParent(child, s->globalM);
    }
}

void spaceUpdateWithParent(Space * s, Matrix parentM) {
    copyMatrix(s->globalM, parentM);
    rightMultiply(s->globalM, s->m);

    for (int idx = 0; idx < s->nChildren; idx ++) {
        Space * child = s->children[idx];
        child->updateWithParent(child, s->globalM);
    }
};

//// Geometry 
void geoUpdate(Space * s);
void geoUpdateWithParent(Space * s, Matrix parentM);
void computeNormals(Geometry * g);
void applyTransform(Geometry * g);

Geometry * makeGeometry(int maxChildren) {
    Space * s = makeSpace(maxChildren);
    s->update = geoUpdate;
    s->updateWithParent = geoUpdateWithParent;

    Geometry * g = malloc(sizeof(Geometry));
    memcpy(&g->space, s, sizeof(*s));

    return g;
}

void setupInitialGeometry(Geometry * g, int nv, Vector * vertices, Vector * normals, int nf, Face * faces) {

    g->nv = nv;
    g->v = vertices;
    g->tv = calloc(nv, sizeof(Vector));
    g->n = normals;
    g->tn = calloc(nv, sizeof(Vector));

    g->nf = nf;
    g->f = faces;

    if (normals == NULL) {
        computeNormals(g);
    }
}

void geoUpdate(Space * s) {
    Geometry * g = (Geometry *)s;
    spaceUpdate(s);
    applyTransform(g);
}

void geoUpdateWithParent(Space * s, Matrix parentM) {
    Geometry * g = (Geometry *)s;
    spaceUpdateWithParent(s, parentM);
    applyTransform(g);
}

void applyTransform(Geometry * g) {
    Space * s = &g->space;
    invertMatrixFrom(g->space.globalNM, g->space.globalM);
    transposeMatrix(g->space.globalNM);

    for (int i = 0; i < g->nv; i++) {
        transform(g->space.globalM, g->v[i], g->tv[i]);
        transformNorm(g->space.globalNM, g->n[i], g->tn[i]);
    }
}

void computeNormals(Geometry * g) {
    assert(g->n == NULL);

    Vector * v = g->v;
    Face * f = g->f;

    Vector * n = malloc(g->nv * sizeof(Vector));
    for (int i = 0; i < g->nv; i++) {
        n[i][0] = n[i][1] = n[i][2] = 0.;
    }

    for (int i = 0; i < g->nf; i++) {
        Vector fn = {0,0,0};

        int *tf = f[i];

        // weighted face normal
        // sum of cross products going counter-clockwise as viewed from outside
        for (int j = 0; j < 4; j++) {

            int pj = j-1;
            if (pj < 0) {
                pj = 3;
            }

            int nj = j+1;
            if (nj > 3) {
                nj = 0;
            }

            double *vp = v[tf[pj]];
            double *v0 = v[tf[j]];
            double *v1 = v[tf[nj]];

            Vector dp, d1;
            dp[0] = v0[0] - vp[0];
            dp[1] = v0[1] - vp[1];
            dp[2] = v0[2] - vp[2];

            d1[0] = v1[0] - v0[0];
            d1[1] = v1[1] - v0[1];
            d1[2] = v1[2] - v0[2];

            fn[0] += dp[1]*d1[2] - dp[2]*d1[1];
            fn[1] += dp[2]*d1[0] - dp[0]*d1[2];
            fn[2] += dp[0]*d1[1] - dp[1]*d1[0];
        }


        // contribute this to the vertex normals
        for (int j = 0; j < 4; j++) {
            double *tn = n[f[i][j]];
            tn[0] += fn[0];
            tn[1] += fn[1];
            tn[2] += fn[2];
        }
    }
 
    // normalize
    for (int i = 0; i < g->nv; i++) {
        double *tn = n[i];
        double mag = sqrt(tn[0]*tn[0]+tn[1]*tn[1]+tn[2]*tn[2]);
        tn[0] /= mag;
        tn[1] /= mag;
        tn[2] /= mag;
    }

    g->n = n;
};

void makeSphere(Geometry * g) {
    const int U = 20;
    const int V = 20;
    int i;
    const int num_verts = U*(V-1)+2;
    const int num_faces = U*V;

    Vector * v;
    Vector * n;
    Face * f;

    v = malloc(num_verts*sizeof(Vector));
    n = malloc(num_verts*sizeof(Vector));
    f = malloc(num_faces*sizeof(Face));

    // -z pole (vertex 0)
    v[0][0] = 0;
    v[0][1] = 0;
    v[0][2] = -1;
    // z pole (vertex 1)
    v[1][0] = 0;
    v[1][1] = 0;
    v[1][2] = 1;

    i = 2;
    for (int vi = 1; vi <= V; vi++) {
        // v ranges from -Pi/2 to Pi/2 (the poles)
        double fv = ((vi*1./U)-.5)*M_PI;
        double sv = sin(fv);
        double cv = cos(fv);

        for (int ui = 0; ui < U; ui++) {
            // u ranges from 0 to 2Pi
            double u = (ui*1./U)*M_PI*2;

            // 'previous' ui (with wrapping)
            int pui = ui-1;
            if (pui == -1) {
                pui = U-1;
            }

            int v0 = i+pui-U;
            int v1 = i+ ui-U;
            int v2 = i+ ui;
            int v3 = i+pui;
            int fi = (vi-1)*U+ui;

            // only adding vertices for non-polar rings
            if (vi != V) {
                v[i+ui][0] = cv*cos(u);
                v[i+ui][1] = cv*sin(u);
                v[i+ui][2] = sv;
            }

            if (vi == 1) {
                // all faces in this ring have a 0-length edge at the -z pole
                v0 = 0;
                v1 = 0;
            } else if (vi == V) {
                // all faces in this ring have a 0-length edge at the z pole
                v2 = 1;
                v3 = 1;
            }

            f[fi][0] = v0;
            f[fi][1] = v1;
            f[fi][2] = v2;
            f[fi][3] = v3;

        }

        i += U;
    }

    // trivially fill in normals
    for (i = 0; i < num_verts; i++) {
        n[i][0] = v[i][0];
        n[i][1] = v[i][1];
        n[i][2] = v[i][2];
    }
    
    setupInitialGeometry(g, num_verts, v, n, num_faces, f);
}

void makeCylinder(Geometry * g) {
    const int U = 20;
    const int num_verts = U*4+2;
    const int num_faces = U*3;
    int i;

    Vector * v;
    Vector * n;
    Face * f;

    v = malloc(num_verts*sizeof(Vector));
    n = malloc(num_verts*sizeof(Vector));
    f = malloc(num_faces*sizeof(Face));

    int fi = 0;

    // -z pole (vertex 0)
    v[0][0] = 0; v[0][1] = 0; v[0][2] = -1;
    n[0][0] = 0; n[0][0] = 0; n[0][2] = -1;

    // z pole (vertex 1)
    v[1][0] = 0; v[1][1] = 1; v[1][2] = 1;
    n[1][0] = 0; n[1][1] = 0; n[1][2] = 1;

    int vi = 2;
    int nvi = 2+4;

    for (int ui = 0; ui < U; ui++) {
        double u = (ui*1./U)*M_PI*2;
        double cu = cos(u);
        double su = sin(u);

        if (ui == U-1) {
            // wrap the "next" index
            nvi = 2;
        }

        // 
        v[vi+0][0] = cu; v[vi+0][1] = su; v[vi+0][2] = -1;
        n[vi+0][0] = cu; n[vi+0][1] = su; n[vi+0][2] = 0;
        v[vi+1][0] = cu; v[vi+1][1] = su; v[vi+1][2] = 1;
        n[vi+1][0] = cu; n[vi+1][1] = su; n[vi+1][2] = 0;

        // for endcaps
        v[vi+2][0] = cu; v[vi+2][1] = su; v[vi+2][2] = -1;
        n[vi+2][0] = 0;  n[vi+2][1] = 0;  n[vi+2][2] = -1;
        v[vi+3][0] = cu; v[vi+3][1] = su; v[vi+3][2] = 1;
        n[vi+3][0] = 0;  n[vi+3][1] = 0;  n[vi+3][2] = 1;

        // side
        f[fi+0][0] = vi + 0;
        f[fi+0][1] = nvi + 0;
        f[fi+0][2] = nvi + 1;
        f[fi+0][3] = vi + 1;

        // -z end cap (cwise)
        f[fi+1][0] = nvi + 2;
        f[fi+1][1] = vi + 2;
        f[fi+1][2] = 0;
        f[fi+1][3] = 0;

        // +z end cap (ccwise)
        f[fi+2][0] = vi + 3;
        f[fi+2][1] = nvi + 3;
        f[fi+2][2] = 1;
        f[fi+2][3] = 1;

        fi += 3;

        vi += 4;
        nvi += 4;
    }

    setupInitialGeometry(g, num_verts, v, n, num_faces, f);
}

void makeCube(Geometry * g) {
    const int num_verts = 24;
    const int num_faces = 6;

    Vector * v = malloc(num_verts*sizeof(Vector));
    Vector * n = malloc(num_verts*sizeof(Vector));
    Face * f = malloc(num_faces*sizeof(Face));

    for (int x = -1; x <= 1; x += 2) {
        for (int y = -1; y <= 1; y += 2) {
            for (int z = -1; z <= 1; z += 2) {
                int i = (x+1)*2+(y+1)+(z+1)/2;

                // normal along x
                v[i][0] = x; v[i][1] = y; v[i][2] = z;
                n[i][0] = x; n[i][1] = 0; n[i][2] = 0;
                i += 8;

                // normal along y
                v[i][0] = x; v[i][1] = y; v[i][2] = z;
                n[i][0] = 0; n[i][1] = y; n[i][2] = 0;
                i += 8;

                // normal along z
                v[i][0] = x; v[i][1] = y; v[i][2] = z;
                n[i][0] = 0; n[i][1] = 0; n[i][2] = z;
            }
        }
    }

    // -x
    f[0][0] = 2+1;      // -x,+y,+z
    f[0][1] = 2;        // -x,+y,-z
    f[0][2] = 0;        // -x,-y,-z
    f[0][3] = 1;        // -x,-y,+z

    // +x
    f[1][0] = 4+2;      // +x,+y,-z
    f[1][1] = 4+2+1;    // +x,+y,+z
    f[1][2] = 4+1;      // +x,-y,+z
    f[1][3] = 4;        // +x,-y,-z

    // -y
    f[2][0] = 8+4+1;    // +x,-y,+z
    f[2][1] = 8+1;      // -x,-y,+z
    f[2][2] = 8;        // -x,-y,-z
    f[2][3] = 8+4;      // +x,-y,-z

    // +y
    f[3][0] = 8+4+2;    // +x,+y,-z
    f[3][1] = 8+2;      // -x,+y,-z
    f[3][2] = 8+2+1;    // -x,+y,+z
    f[3][3] = 8+4+2+1;  // +x,+y,+z

    // -z
    f[4][0] = 16+4+2;   // +x,+y,-z
    f[4][1] = 16+4;     // +x,-y,-z
    f[4][2] = 16;       // -x,-y,-z
    f[4][3] = 16+2;     // -x,+y,-z

    // +z
    f[5][0] = 16+4+2+1; // +x,+y,+z
    f[5][1] = 16+2+1;   // -x,+y,+z
    f[5][2] = 16+1;     // -x,-y,+z
    f[5][3] = 16+4+1;   // +x,-y,+z

    setupInitialGeometry(g, num_verts, v, n, num_faces, f);
}

#if 0
var drawLine = function (imgData, sx, sy, ex, ey, c0, c1, c2, c3) {
    // an easy line plotter, everything's a double anyway

    // abandon if points are offscreen
    sx = Math.round(sx);
    sy = Math.round(sy);
    ex = Math.round(ex);
    ey = Math.round(ey);

    var W = imgData.width;
    var H = imgData.height;
    var data = imgData.data;

    if (sx < 0 || sx >= W || sy < 0 || sy >= H) return;
    if (ex < 0 || ex >= W || ey < 0 || ey >= H) return;

    var w = Math.abs(ex-sx);
    var h = Math.abs(ey-sy);
    var temp, dx, dy, cx, cy;
    var idx;

    if (w > h || (w > 0 && w == h)) {
        if (ex < sx) {
            temp = sx; sx = ex; ex = temp;
            temp = sy; sy = ey; ey = temp;
        }
        dy = ey-sy;
        for (cx = sx; cx <= ex; cx++) {
            cy = Math.round(sy+(dy * (cx-sx) / w));
            idx = (cx + cy*W)*4;
            data[idx+0] = c0;
            data[idx+1] = c1;
            data[idx+2] = c2;
            data[idx+3] = c3;
        }
    } else if (h > w) {
        if (ey < sy) {
            temp = sx; sx = ex; ex = temp;
            temp = sy; sy = ey; ey = temp;
        }
        dx = ex-sx;
        for (cy = sy; cy <= ey; cy++) {
            cx = Math.round(sx+(dx * (cy-sy) / h));
            idx = (cx + cy*W)*4;
            data[idx+0] = c0;
            data[idx+1] = c1;
            data[idx+2] = c2;
            data[idx+3] = c3;
        }
    }
};
#endif

void fillTrap(
    SDL_Surface * screen, double * zbuf,
    double minX0, double maxX0,
    int minY,
    double leftZ0, double rightZ0,
    int leftR0, int leftG0, int leftB0,
    int rightR0, int rightG0, int rightB0,
    double minX1, double maxX1,
    int maxY,
    double leftZ1, double rightZ1,
    int leftR1, int leftG1, int leftB1,
    int rightR1, int rightG1, int rightB1) {

    const int W = screen->w;
    const int pitch = screen->pitch;
    const int H = screen->h;
    uint32_t * pixels = screen->pixels;
    uint32_t * linepix;
    SDL_PixelFormat * format = screen->format;

    int dy = maxY - minY;
    double dMinX = minX1 - minX0;
    double dMaxX = maxX1 - maxX0;
    double dLeftZ =  leftZ1   -  leftZ0;
    double dRightZ = rightZ1   - rightZ0;
    int  dLeftR =  leftR1   -  leftR0;
    int dRightR = rightR1   - rightR0;
    int  dLeftG =  leftG1   -  leftG0;
    int dRightG = rightG1   - rightG0;
    int  dLeftB =  leftB1   -  leftB0;
    int dRightB = rightB1   - rightB0;

    int idxZ;
    int cMinX, cMaxX, cdx;
    double cLeftZ, cRightZ, cdz;
    int cLR, cRR, cLG, cRG, cLB, cRB, cdR, cdG, cdB;

    if (dy == 0) dy = 1;

    for (int y = minY, ay = 0; y <= maxY; y++, ay++) {
        double t = ay*1./dy;
        cMinX = ceil(minX0 + dMinX * t);
        cMaxX = floor(maxX0 + dMaxX * t);

        cLeftZ  = leftZ0    + dLeftZ*t;
        cRightZ = rightZ0   + dRightZ*t;
        
        cdx = cMaxX     - cMinX;
        cdz = cRightZ   - cLeftZ;
        if (cdx == 0) cdx = 1;

        cLR = leftR0    + dLeftR*t;
        cRR = rightR0   + dRightR*t;
        cLG = leftG0    + dLeftG*t;
        cRG = rightG0   + dRightG*t;
        cLB = leftB0    + dLeftB*t;
        cRB = rightB0   + dRightB*t;

        cdR = cRR - cLR;
        cdG = cRG - cLG;
        cdB = cRB - cLB;

        linepix = (uint32_t*)((char*)pixels + pitch*y);
        idxZ = cMinX+W*y;
        for (int x = cMinX, ax = 0; x <= cMaxX; x++, ax++, idxZ++) {
            double lt = ax*1./cdx;
            
            double z = cLeftZ + cdz * lt;

            if (z > zbuf[idxZ]) {
                zbuf[idxZ] = z;

                set_pixel(&linepix[x], 
                    cLR + cdR * lt,
                    cLG + cdG * lt,
                    cLB + cdB * lt);
            }
        }
    }
}

void fillTri(
    SDL_Surface * screen, double * zbuf,
    int minX, int minY, double minZ, Vector minN,
    int midX, int midY, double midZ, Vector midN,
    int maxX, int maxY, double maxZ, Vector maxN) {

    double midT = (midY-minY)*1./(maxY-minY);
    if (maxY == minY) midT = 0;
    double iMidX = minX + (maxX-minX)*midT;
    double iMidZ = minZ + (maxZ-minZ)*midT;

    // fake colors with normals
    int minR = fmax(0,minN[0] * 255);
    int minG = fmax(0,minN[1] * 255);
    int minB = fmax(0,minN[2] * 255);

    int midR = fmax(0,midN[0] * 255);
    int midG = fmax(0,midN[1] * 255);
    int midB = fmax(0,midN[2] * 255);

    int maxR = fmax(0,maxN[0] * 255);
    int maxG = fmax(0,maxN[1] * 255);
    int maxB = fmax(0,maxN[2] * 255);

    int iMidR = minR + (maxR - minR) * midT;
    int iMidG = minG + (maxG - minG) * midT;
    int iMidB = minB + (maxB - minB) * midT;

    if (midX > iMidX) {
        // knee is on the right (1)

        fillTrap(screen, zbuf,
            minX, minX,
            minY,
            minZ, minZ,
            minR, minG, minB,
            minR, minG, minB,

            iMidX, midX,
            midY,
            iMidZ, midZ,
            iMidR, iMidG, iMidB,
             midR,  midG,  midB);

        fillTrap(screen, zbuf,
            iMidX, midX,
            midY,
            iMidZ, midZ, 
            iMidR, iMidG, iMidB,
             midR,  midG,  midB,

            maxX, maxX,
            maxY, 
            maxZ, maxZ,
            maxR, maxG, maxB,
            maxR, maxG, maxB);

    } else {
        // knee is on the left (0)
        fillTrap(screen, zbuf,
            minX, minX,
            minY,
            minZ, minZ,
            minR, minG, minB,
            minR, minG, minB,

            midX, iMidX,
            midY,
            midZ, iMidZ,
             midR,  midG,  midB,
            iMidR, iMidG, iMidB);

        fillTrap(screen, zbuf,
            midX, iMidX,
            midY,
            midZ, iMidZ, 
             midR,  midG,  midB,
            iMidR, iMidG, iMidB,

            maxX, maxX,
            maxY, 
            maxZ, maxZ,
            maxR, maxG, maxB,
            maxR, maxG, maxB);
    }
}

int projectX(Vector p, int W, int H, double FL) {
    return floor((W / 2) + (H * p[0] / (FL - p[2])));
}

int projectY(Vector p, int W, int H, double FL) {
    return floor((H / 2) - (H * p[1] / (FL - p[2])));
}

double projectZ(Vector p, double FL) {
    return p[2] / (FL - p[2]);
}

#if 0
var drawProjLine = function (sp, ep, r, g, b) {
    var spx = projectX(sp);
    var spy = projectY(sp);

    var epx = projectX(ep);
    var epy = projectY(ep);

    return drawLine(imgData, spx, spy, epx, epy, r, g, b, 255);
};
#endif

void fillProjTri (SDL_Surface * screen, double * zbuf, double FL,
    Vector v0, Vector v1, Vector v2,
    Vector n0, Vector n1, Vector n2) {

    const int W = screen->w;
    const int H = screen->h;

    int minY, minX;
    double minZ;
    double *minN;
    int midY, midX;
    double midZ;
    double *midN;
    int maxY, maxX;
    double maxZ;
    double *maxN;

    int v0x = projectX(v0,W,H,FL);
    int v0y = projectY(v0,W,H,FL);
    double v0z = projectZ(v0,FL);

    int v1x = projectX(v1,W,H,FL);
    int v1y = projectY(v1,W,H,FL);
    double v1z = projectZ(v1,FL);

    int v2x = projectX(v2,W,H,FL);
    int v2y = projectY(v2,W,H,FL);
    double v2z = projectZ(v2,FL);

    // check that the triangle goes counterclockwise
    if ((v0y+v2y)*(v0x-v2x) + (v1y+v0y)*(v1x-v0x) + (v2y+v1y)*(v2x-v1x) <= 0) {
        return;
    }

    // sort vertices by Y
    if (v0y < v1y) {
        if (v1y < v2y) {
            // 0, 1, 2
            minY = v0y; minX = v0x; minZ = v0z; minN = n0;
            midY = v1y; midX = v1x; midZ = v1z; midN = n1;
            maxY = v2y; maxX = v2x; maxZ = v2z; maxN = n2;
        } else if (v0y < v2y) {
            // 0, 2, 1
            minY = v0y; minX = v0x; minZ = v0z; minN = n0;
            midY = v2y; midX = v2x; midZ = v2z; midN = n2;
            maxY = v1y; maxX = v1x; maxZ = v1z; maxN = n1;
        } else {
            // 2, 0, 1
            minY = v2y; minX = v2x; minZ = v2z; minN = n2;
            midY = v0y; midX = v0x; midZ = v0z; midN = n0;
            maxY = v1y; maxX = v1x; maxZ = v1z; maxN = n1;
        }
    } else {
        // 1 < 0
        if (v0y < v2y) {
            // 1, 0, 2
            minY = v1y; minX = v1x; minZ = v1z; minN = n1;
            midY = v0y; midX = v0x; midZ = v0z; midN = n0;
            maxY = v2y; maxX = v2x; maxZ = v2z; maxN = n2;
        } else if (v1y < v2y) {
            // 1, 2, 0
            minY = v1y; minX = v1x; minZ = v1z; minN = n1;
            midY = v2y; midX = v2x; midZ = v2z; midN = n2;
            maxY = v0y; maxX = v0x; maxZ = v0z; maxN = n0;
        } else {
            // 2, 1, 0
            minY = v2y; minX = v2x; minZ = v2z; minN = n2;
            midY = v1y; midX = v1x; midZ = v1z; midN = n1;
            maxY = v0y; maxX = v0x; maxZ = v0z; maxN = n0;
        }
    }

    fillTri(screen, zbuf, minX, minY, minZ, minN, midX, midY, midZ, midN, maxX, maxY, maxZ, maxN);
};

void cls(SDL_Surface *screen, double * zbuf, double FL) {
    const int H = screen->h;
    const int W = screen->w;

    SDL_Rect r = {.x=0,.y=0,.w=W,.h=H};
    SDL_PixelFormat * format = screen->format;
    uint32_t black = SDL_MapRGB(format,0,0,0);

    int x,y,z;
    z = 0;
    SDL_FillRect(screen, &r, black);
    for (y = 0; y < H; y++) {
        for (x = 0; x < W; x++) {
            zbuf[z] = -FL;
            z ++;
        }
    }
}

#if 0

    drawLine : function (sx, sy, ex, ey, r, g, b) {
        return drawLine(imgData, sx, sy, ex, ey, r, g, b, 255);
    },

    wireframe : function (geo, r, g, b) {
        var f = geo.f;
        var v = geo.tv;
        var i;

        for (i = 0; i < f.length; i++) {
            drawProjLine(v[f[i][0]], v[f[i][1]], r, g, b);
            drawProjLine(v[f[i][1]], v[f[i][2]], r, g, b);
            drawProjLine(v[f[i][2]], v[f[i][3]], r, g, b);
            drawProjLine(v[f[i][3]], v[f[i][0]], r, g, b);
        }
    },
#endif

void fill(Geometry * geo, SDL_Surface * screen, double * zbuf, double FL) {
    Face * f = geo->f;
    Vector * v = geo->tv;
    Vector * n = geo->tn;

    for (int i = 0; i < geo->nf; i++) {
        /*NOTE: assume the quad faces are convex so any triangleization works*/
        fillProjTri(screen, zbuf, FL,
                v[f[i][0]],
                v[f[i][1]],
                v[f[i][2]],
                n[f[i][0]],
                n[f[i][1]],
                n[f[i][2]]);
        fillProjTri(screen, zbuf, FL,
                v[f[i][0]],
                v[f[i][2]],
                v[f[i][3]],
                n[f[i][0]],
                n[f[i][2]],
                n[f[i][3]]);
    }
}
