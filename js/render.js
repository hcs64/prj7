"use strict";

////
var Matrix = function () {
    this.m = [0,0,0,0,
              0,0,0,0,
              0,0,0,0,
              0,0,0,0];

    return this;
};

Matrix.prototype = {
    copyFrom : function (m2) {
        var i;
        for (i = 0; i < 16; i++) {
            this.m[i] = m2.m[i];
        }
    },

    // these class-wide temporaries make these calls non-reentrant,
    // but they don't need to be reentrant anyway
    t : [0,0,0,0,
         0,0,0,0,
         0,0,0,0,
         0,0,0,0],

    invSrcT : [0,0,0,0,
          0,0,0,0,
          0,0,0,0,
          0,0,0,0],

    invertIndex : [0,0,0,0],
    invertMax : [0,0,0,0],

    invertFrom : function (m2) {
        // ported from Ken Perlin's Java implementation

        // INVERT A 4x4 MATRIX BY GAUSSIAN ELIMINATION

        var i, j, k, m = 0, n = 4;
        var t, s, p, swap;

        var src = m2.m;
        var srcT = this.invSrcT;
        var dst = this.m;
        var tmp = this.t;
        var max = this.invertMax;
        var index = this.invertIndex;

        for (i = 0 ; i < n ; i++) {
            index[i] = i;
            max  [i] = 0;
            for (j = 0 ; j < n ; j++) {
                max[i] = Math.max(max[i], Math.abs(src[i*n+j]));
                tmp[i*n+j] = ((i == j) ? 1 : 0);
                srcT[i*n+j] = src[i*n+j];
            }
        }

        src = srcT;

        for (j = 0 ; j < n - 1 ; j++) {
            t = 0;

            for (i = j ; i < n ; i++) {
                s = Math.abs(src[index[i]*n+j]) / max[index[i]];
                if (s > t) {
                    t = s;
                    m = i;
                }
            }
            swap = index[j];
            index[j] = index[m];
            index[m] = swap;

            for (i = j + 1 ; i < n ; i++) {
                p = src[index[i]*n+j] / src[index[j]*n+j];
                src[index[i]*n+j] = p;
                for (k = j + 1 ; k < n ; k++) {
                    src[index[i]*n+k] -= p * src[index[j]*n+k];
                }
            }
        }

        for (i = 0      ; i < n - 1 ; i++) {
        for (j = i + 1  ; j < n     ; j++) {
        for (k = 0      ; k < n     ; k++) {
            tmp[index[j]*n+k] -= src[index[j]*n+i] * tmp[index[i]*n+k];
        }}}

        for (i = 0 ; i < n ; i++) {
            dst[(n-1)*n+i] = tmp[index[n-1]*n+i] / src[index[n-1]*n+n-1];
            for (j = n - 2 ; j >= 0 ; j--) {
                dst[j*n+i] = tmp[index[j]*n+i];
                for (k = j + 1 ; k < n ; k++) {
                    dst[j*n+i] -= src[index[j]*n+k] * dst[k*n+i];
                }
                dst[j*n+i] /= src[index[j]*n+j];
            }
        }

    },

    set : function (c,r,v) {
        this.m[r*4+c] = v;
    },

    get : function (c,r) {
        return this.m[r*4+c];
    },

    inc : function (c,r,v) {
        this.m[r*4+c] += v;
    },

    mult : function (c,r,v) {
        this.m[r*4+c] *= v;
    },

    makeIdentity : function () {
        this.copyFrom(Matrix.identity);
    },

    transpose : function () {
        // swap across the diagonal
        var swap;
        var m = this.m;

        swap = m[0*4+1];
        m[0*4+1] = m[1*4+0];
        m[1*4+0] = swap;

        swap = m[0*4+2];
        m[0*4+2] = m[2*4+0];
        m[2*4+0] = swap;

        swap = m[0*4+3];
        m[0*4+3] = m[3*4+0];
        m[3*4+0] = swap;

        swap = m[1*4+2];
        m[1*4+2] = m[2*4+1];
        m[2*4+1] = swap;

        swap = m[1*4+3];
        m[1*4+3] = m[3*4+1];
        m[3*4+1] = swap;

        swap = m[2*4+3];
        m[2*4+3] = m[3*4+2];
        m[3*4+2] = swap;
    },

    translate : function (x, y, z) {
        this.inc(3,0,x);
        this.inc(3,1,y);
        this.inc(3,2,z);
    },

    rotateX : function (radians) {
        var cos = Math.cos(radians);
        var sin = Math.sin(radians);
        var i, u, d;

        for (i = 0; i < 4; i++) {
            // from Y (1) into  Z (2)
            // from Z (2) into -Y (0)
            u = this.get(i,1);
            d = this.get(i,2);

            this.set(i,1, u*cos + -d*sin);
            this.set(i,2, u*sin +  d*cos);
        }
    },

    rotateY : function (radians) {
        var cos = Math.cos(radians);
        var sin = Math.sin(radians);
        var i, u, d;

        for (i = 0; i < 4; i++) {
            // from X (0) into -Z (2)
            // from Z (2) into  X (0)
            u = this.get(i,0);
            d = this.get(i,2);

            this.set(i,0,  u*cos + d*sin);
            this.set(i,2, -u*sin + d*cos);
        }
    },

    rotateZ : function (radians) {
        var cos = Math.cos(radians);
        var sin = Math.sin(radians);
        var i, u, d;

        for (i = 0; i < 4; i++) {
            // from X (0) into  Y (1)
            // from Y (1) into -X (0)
            u = this.get(i,0);
            d = this.get(i,1);

            this.set(i,0, u*cos + -d*sin);
            this.set(i,1, u*sin +  d*cos);
        }
    },

    scale : function (x, y, z) {
        var i;
        for (i = 0; i < 4; i++) {
            this.mult(0,i,x);
            this.mult(1,i,y);
            this.mult(2,i,z);
        }
    },

    // leftMultiply  means:  this <= other times this
    leftMultiply : function (other) {
        var i,j,k,sum;
        var idx;

        idx = 0;
        for (j = 0; j < 4; j++) {
            for (i = 0; i < 4; i++) {
                sum = 0;
                for (k = 0; k < 4; k++) {
                    // walk across row of other, down col of this
                    sum += other.get(k,j) * this.get(i,k);
                }
                this.t[idx] = sum;
                idx ++;
            }
        }

        // update this Matrix
        for (idx = 0; idx < 16; idx++) {
            this.m[idx] = this.t[idx];
        }
    },

    // rightMultiply means:  this <= this times other
    rightMultiply : function (other) {
        var i,j,k,sum;
        var idx;

        idx = 0;
        for (j = 0; j < 4; j++) {
            for (i = 0; i < 4; i++) {
                sum = 0;
                for (k = 0; k < 4; k++) {
                    // walk across row of this, down col of other
                    sum += this.get(k,j) * other.get(i,k);
                }
                this.t[idx] = sum;
                idx ++;
            }
        }

        // update this Matrix
        for (idx = 0; idx < 16; idx++) {
            this.m[idx] = this.t[idx];
        }
    },

    transform : function (src, dst) {
        var i, j, sum, idx = 0;
        for (i = 0; i < 3; i++) {
            sum = 0;
            for (j = 0; j < 4; j++) {
                if (j < 3) {
                    sum += this.m[idx]*src[j];
                } else {
                    sum += this.m[idx];
                }
                idx ++;
            }
            dst[i] = sum;
        }
    },

    transformNorm : function (src, dst) {
        var i, j, sum, idx = 0, m;
        for (i = 0; i < 3; i++) {
            sum = 0;
            for (j = 0; j < 4; j++) {
                if (j < 3) {
                    sum += this.m[idx]*src[j];
                }
                idx ++;
            }
            dst[i] = sum;
        }

        // normalize
        m = Math.sqrt(dst[0]*dst[0] + dst[1]*dst[1] + dst[2]*dst[2]);
        dst[0] /= m;
        dst[1] /= m;
        dst[2] /= m;
    },

};

Matrix.identity = new Matrix;
Matrix.identity.m =
        [1,0,0,0,
         0,1,0,0,
         0,0,1,0,
         0,0,0,1];

////
var Space = function () {
    this.m = new Matrix;
    this.m.makeIdentity();
    this.globalM = new Matrix;
    this.globalNM = new Matrix;

    this.children = [];
};

Space.prototype = {
    add : function (child) {
        this.children.push(child);
    },

    remove : function (child) {
        var idx = this.children.indexOf(child);
        this.children.splice(idx, 1);
    },

    getChild : function (idx) {
        return this.children[idx];
    },

    getNumChildren : function () {
        return this.children.length;
    },

    update : function () {
        var idx;
        this.globalM.copyFrom(this.m);

        for (idx = 0; idx < this.children.length; idx ++) {
            this.children[idx].updateWithParent(this.globalM);
        }
    },

    updateWithParent : function (parentM) {
        var idx;
        this.globalM.copyFrom(parentM);
        this.globalM.rightMultiply(this.m);

        for (idx = 0; idx < this.children.length; idx ++) {
            this.children[idx].updateWithParent(this.globalM);
        }
    },
};

////
var Geometry = function () {
    Space.apply(this);
};

Geometry.prototype = new Space();
Geometry.prototype.setupInitialGeometry = function (vertices, faces, normals) {
    var i;

    this.v = vertices;
    this.tv = [];
    this.tn = [];
    for (i = 0; i < vertices.length; i++) {
        this.tv[i] = [0,0,0];
        this.tn[i] = [0,0,0];
    }

    this.f = faces;

    if (normals === undefined) {
        this.computeNormals();
    } else {
        this.n = normals;
    }
};

Geometry.prototype.update = function () {
    Space.prototype.update.apply(this);
    this.applyTransform();
};

Geometry.prototype.updateWithParent = function (parentM) {
    Space.prototype.updateWithParent.apply(this,[parentM]);
    this.applyTransform();
};

Geometry.prototype.applyTransform = function () {
    var i;
    var m = this.globalM;
    var nm = this.globalNM;
    var v = this.v, tv = this.tv, n = this.n, tn = this.tn;

    nm.invertFrom(m);
    nm.transpose();

    for (i = 0; i < v.length; i++) {
        m.transform(v[i], tv[i]);
        nm.transformNorm(n[i], tn[i]);
    }
};

Geometry.prototype.computeNormals = function () {
    var i, j, pj, nj;
    
    var fn, m;
    var f = this.f, tf;
    var v = this.v;
    var n, tn;
    var vp, v0, v1;
    var dp = [0,0,0];
    var d1 = [0,0,0];

    n = [];

    for (i = 0; i < f.length; i++) {
        fn = [0,0,0];

        tf = f[i];

        // weighted face normal
        // sum of cross products going counter-clockwise as viewed from outside
        for (j = 0; j < 4; j++) {

            pj = j-1;
            if (pj < 0) {
                pj = 3;
            }

            nj = j+1;
            if (nj > 3) {
                nj = 0;
            }

            vp = v[tf[pj]];
            v0 = v[tf[j]];
            v1 = v[tf[nj]];

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
        for (j = 0; j < 4; j++) {
            tn = n[f[i][j]];
            if (tn !== undefined) {
                tn[0] += fn[0];
                tn[1] += fn[1];
                tn[2] += fn[2];
            } else {
                n[f[i][j]] = fn;
            }
        }
    }
 
    // normalize
    for (i = 0; i < v.length; i++) {
        tn = n[i];
        m = Math.sqrt(tn[0]*tn[0]+tn[1]*tn[1]+tn[2]*tn[2]);
        tn[0] /= m;
        tn[1] /= m;
        tn[2] /= m;
    }

    this.n = n;
};

//// Some unit-sized parametric objects
var Sphere = function () {
    Geometry.apply(this);
    this.setupSphere();
    this.setupInitialGeometry(this.v, this.f, this.n);
};

Sphere.prototype = new Geometry;

Sphere.prototype.setupSphere = function () {
    var U = 20;
    var V = 20;
    var i;

    var u;
    var vi, ui, pui, fi;
    var fv, sv, cv;
    var v0, v1, v2, v3;

    this.v = [];
    this.f = [];

    // -z pole (vertex 0)
    this.v[0] = [0, 0, -1];
    // z pole (vertex 1)
    this.v[1] = [0, 0, 1];

    i = 2;
    for (vi = 1; vi <= V; vi++) {
        // v ranges from -Pi/2 to Pi/2 (the poles)
        fv = ((vi/U)-.5)*Math.PI;
        sv = Math.sin(fv);
        cv = Math.cos(fv);

        for (ui = 0; ui < U; ui++) {
            // u ranges from 0 to 2Pi
            u = (ui/U)*Math.PI*2;

            // 'previous' ui (with wrapping)
            pui = ui-1;
            if (pui == -1) {
                pui = U-1;
            }

            v0 = i+pui-U;
            v1 = i+ ui-U;
            v2 = i+ ui;
            v3 = i+pui;
            fi = (vi-1)*U+ui;

            // only adding vertices for non-polar rings
            if (vi != V) {
                this.v[i+ui] = [ cv*Math.cos(u), cv*Math.sin(u), sv ];
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

            this.f[fi] = [v0,v1,v2,v3];
        }

        i += U;
    }

    // trivially fill in normals
    this.n = [];
    for (i = 0; i < this.v.length; i++) {
        this.n[i] = this.v[i];
    }
};

var Cylinder = function () {
    Geometry.apply(this);
    this.setupCylinder();
    this.setupInitialGeometry(this.v, this.f, this.n);
};

Cylinder.prototype = new Geometry;

Cylinder.prototype.setupCylinder = function () {
    // sides
    var U = 20;

    this.v = [];
    this.f = [];
    this.n = [];

    var i, u, sin, cos;
    var fi, vi, nvi, ui;
    
    fi = 0;

    // -z pole (vertex 0)
    this.v[0] = [0, 0, -1];
    this.n[0] = [0, 0, -1];
    // z pole (vertex 1)
    this.v[1] = [0, 0, 1];
    this.n[1] = [0, 0, 1];

    vi = 2;
    nvi = 2+4;

    for (ui = 0; ui < U; ui++) {
        u = (ui/U)*Math.PI*2;
        cos = Math.cos(u);
        sin = Math.sin(u);

        if (ui == U-1) {
            // wrap the "next" index
            nvi = 2;
        }

        // 
        this.v[vi+0] = [cos, sin, -1];
        this.n[vi+0] = [cos, sin, 0];
        this.v[vi+1] = [cos, sin, 1];
        this.n[vi+1] = [cos, sin, 0];

        // for endcaps
        this.v[vi+2] = [cos, sin, -1];
        this.n[vi+2] = [0, 0, -1];
        this.v[vi+3] = [cos, sin, 1];
        this.n[vi+3] = [0, 0, 1];

        // side
        this.f[fi+0] = [vi + 0, nvi + 0, nvi + 1, vi + 1];

        // -z end cap (cwise)
        this.f[fi+1] = [nvi + 2, vi + 2, 0, 0];

        // +z end cap (ccwise)
        this.f[fi+2] = [vi + 3, nvi + 3, 1, 1];

        fi += 3;

        vi += 4;
        nvi += 4;
    }
};

var Cube = function () {
    Geometry.apply(this);
    this.setupCube();
    this.setupInitialGeometry(this.v, this.f, this.n);
};

Cube.prototype = new Geometry;

Cube.prototype.setupCube = function () {
    // copy a face and its vertices, applying some transformation
    var that = this;
    var copyTransFace = function (src_f, dst_v, dst_f, tm) {
        var i, v0;
        that.f[dst_f] = [0,0,0,0];
        for (i = 0; i < 4; i++) {
            v0 = that.f[src_f][i];
            that.v[dst_v+i] = [0,0,0];
            tm.transform(that.v[v0], that.v[dst_v+i]);
            that.f[dst_f][i] = dst_v+i;
        }
    };

    var rm;
    
    this.v = [];
    this.f = [];
    this.n = [];

    // facing +z
    this.v[0] = [ 1, 1, 1];
    this.v[1] = [-1, 1, 1];
    this.v[2] = [-1,-1, 1];
    this.v[3] = [ 1,-1, 1];

    this.n[0] = [ 0, 0, 1];
    this.n[1] = [ 0, 0, 1];
    this.n[2] = [ 0, 0, 1];
    this.n[3] = [ 0, 0, 1];

    this.f[0] = [0, 1, 2, 3];

    // facing x
    rm = new Matrix;
    rm.makeIdentity();
    rm.rotateY(Math.PI/2);
    copyTransFace(0, 4, 1, rm);

    this.n[4] = [ 1, 0, 0];
    this.n[5] = [ 1, 0, 0];
    this.n[6] = [ 1, 0, 0];
    this.n[7] = [ 1, 0, 0];

    // facing -z
    rm.rotateY(Math.PI/2);
    copyTransFace(0, 8, 2, rm);

    this.n[8] = [ 0, 0,-1];
    this.n[9] = [ 0, 0,-1];
    this.n[10]= [ 0, 0,-1];
    this.n[11]= [ 0, 0,-1];

    // facing -x
    rm.rotateY(Math.PI/2);
    copyTransFace(0, 12, 3, rm);

    this.n[12]= [-1, 0, 0];
    this.n[13]= [-1, 0, 0];
    this.n[14]= [-1, 0, 0];
    this.n[15]= [-1, 0, 0];

    // facing y
    rm.makeIdentity();
    rm.rotateX(-Math.PI/2);
    copyTransFace(0, 16, 4, rm);

    this.n[16] = [ 0, 1, 0];
    this.n[17] = [ 0, 1, 0];
    this.n[18] = [ 0, 1, 0];
    this.n[19] = [ 0, 1, 0];

    // facing -y
    rm.makeIdentity();
    rm.rotateX(Math.PI/2);
    copyTransFace(0, 20, 5, rm);

    this.n[20] = [ 0,-1, 0];
    this.n[21] = [ 0,-1, 0];
    this.n[22] = [ 0,-1, 0];
    this.n[23] = [ 0,-1, 0];

};

//// putting pixels onscreen
var graphics = (function(){
var o = {};

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

var fillTrap = function (imgData,
    minX0, maxX0,
    minY,
    leftZ0, rightZ0,
    leftR0, leftG0, leftB0,
    rightR0, rightG0, rightB0,
    minX1, maxX1,
    maxY,
    leftZ1, rightZ1,
    leftR1, leftG1, leftB1,
    rightR1, rightG1, rightB1) {

    var data = imgData.data;
    var zbuf = imgData.zbuf;
    var W = imgData.width;

    var dy = maxY - minY;
    var dMinX = minX1 - minX0;
    var dMaxX = maxX1 - maxX0;
    var  dLeftZ =  leftZ1   -  leftZ0;
    var dRightZ = rightZ1   - rightZ0;
    var  dLeftR =  leftR1   -  leftR0;
    var dRightR = rightR1   - rightR0
    var  dLeftG =  leftG1   -  leftG0;
    var dRightG = rightG1   - rightG0;
    var  dLeftB =  leftB1   -  leftB0;
    var dRightB = rightB1   - rightB0;

    var x, y, z, idx, idxZ;
    var ax, ay;
    var cMinX, cMaxX, cdx;
    var cLeftZ, cRightZ, cdz;
    var cLR, cRR, cLG, cRG, cLB, cRB, cdR, cdG, cdB;

    var t, lt;

    for (y = minY, ay = 0; y <= maxY; y++, ay++) {
        t = ay/dy;
        cMinX = Math.ceil(minX0 + dMinX * t);
        cMaxX = Math.floor(maxX0 + dMaxX * t);

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

        idx = (cMinX+W*y)*4;
        idxZ = cMinX+W*y;
        for (x = cMinX, ax = 0; x <= cMaxX; x++, ax++, idx+=4, idxZ++) {
            lt = ax/cdx;
            
            z = cLeftZ + cdz * lt;

            if (z > zbuf[idxZ]) {
                zbuf[idxZ] = z;

                data[idx+0] = cLR + cdR * lt;
                data[idx+1] = cLG + cdG * lt;
                data[idx+2] = cLB + cdB * lt;
                data[idx+3] = 255;
            }
        }
    }
};

var fillTri = function (imgData, minX, minY, minZ, minN, midX, midY, midZ, midN, maxX, maxY, maxZ, maxN, r, g, b) {
    var midT = (midY-minY)/(maxY-minY);
    var iMidX = minX + (maxX-minX)*midT;
    var iMidZ = minZ + (maxZ-minZ)*midT;
    
    // fake colors with normals
    var minR = Math.max(0,minN[0] * 255);
    var minG = Math.max(0,minN[1] * 255);
    var minB = Math.max(0,minN[2] * 255);

    var midR = Math.max(0,midN[0] * 255);
    var midG = Math.max(0,midN[1] * 255);
    var midB = Math.max(0,midN[2] * 255);

    var maxR = Math.max(0,maxN[0] * 255);
    var maxG = Math.max(0,maxN[1] * 255);
    var maxB = Math.max(0,maxN[2] * 255);

    var iMidR = minR + (maxR - minR) * midT;
    var iMidG = minG + (maxG - minG) * midT;
    var iMidB = minB + (maxB - minB) * midT;

    if (midX > iMidX) {
        // knee is on the right (1)

        fillTrap(imgData,
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

        fillTrap(imgData,
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
        fillTrap(imgData,
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

        fillTrap(imgData,
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
};

o.getGraphics = function (elem, W, H) {
var ctx, imgData, data, zbuf;

elem.width = W;
elem.height = H;

ctx = elem.getContext('2d');
imgData = ctx.getImageData(0,0,W,H);
zbuf = [];
imgData.zbuf = zbuf;
data = imgData.data;

var FL = 10;

var projectX = function (p) {
    return Math.floor((W / 2) + (H * p[0] / (FL - p[2])));
};

var projectY = function (p) {
    return Math.floor((H / 2) - (H * p[1] / (FL - p[2])));
};

var projectZ = function (p) {
    return p[2] / (FL - p[2]);
};

var drawProjLine = function (sp, ep, r, g, b) {
    var spx = projectX(sp);
    var spy = projectY(sp);

    var epx = projectX(ep);
    var epy = projectY(ep);

    return drawLine(imgData, spx, spy, epx, epy, r, g, b, 255);
};

var fillProjTri = function (v0, v1, v2, n0, n1, n2, r, g, b) {
    var minY, minX, minZ, minN;
    var midY, midX, midZ, midN;
    var maxY, maxX, maxZ, maxN;

    var v0x = projectX(v0);
    var v0y = projectY(v0);
    var v0z = projectZ(v0);

    var v1x = projectX(v1);
    var v1y = projectY(v1);
    var v1z = projectZ(v1);

    var v2x = projectX(v2);
    var v2y = projectY(v2);
    var v2z = projectZ(v2);

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

    fillTri(imgData, minX, minY, minZ, minN, midX, midY, midZ, midN, maxX, maxY, maxZ, maxN, r, g, b);

    //drawLine(imgData, minX, minY, midX, midY, 255, 0, 0, 255);
    //drawLine(imgData, midX, midY, maxX, maxY, 255, 255, 0, 255);
    //drawLine(imgData, maxX, maxY, minX, minY, 0, 0, 255, 255);

    //drawLine(imgData, v0x, v0y, v1x, v1y, 255, 0, 0, 255);
    //drawLine(imgData, v1x, v1y, v2x, v2y, 255, 255, 0, 255);
    //drawLine(imgData, v2x, v2y, v0x, v0y, 0, 0, 255, 255);

};

return {
    cls : function() {
        var x,y,i,z;
        i = 0, z = 0;
        for (y = 0; y < H; y++) {
            for (x = 0; x < W; x++) {
                data[i+0] = 0;
                data[i+1] = 0;
                data[i+2] = 0;
                data[i+3] = 255;
                i += 4;
                zbuf[z] = -FL;
                z ++;
            }
        }
    },

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

    fill : function (geo, r, g, b) {
        var f = geo.f;
        var v = geo.tv;
        var n = geo.tn;
        var i;

        for (i = 0; i < f.length; i++) {
            fillProjTri(
                v[f[i][0]],
                v[f[i][1]],
                v[f[i][2]],
                n[f[i][0]],
                n[f[i][1]],
                n[f[i][2]],
                r, g, b);
            fillProjTri(
                v[f[i][0]],
                v[f[i][2]],
                v[f[i][3]],
                n[f[i][0]],
                n[f[i][2]],
                n[f[i][3]],
                255 - r, 255 - g, 255 - b);
        }
    },
    
    render : function (t) {
        var idx = 0;
        var b = Math.sin(t/500)*127+127;

        for (var y = 0; y < H; y++) {
            for (var x = 0; x < W; x++) {
                data[idx+0]=x/W*255;
                data[idx+1]=y/H*255;
                data[idx+2]=b;
                data[idx+3]=255;
                idx += 4;
            }
        }
    },

    display : function (elem) {
        ctx.putImageData(imgData,0,0);
    }
};

}; // end of o.getGraphics


return o;

})();
