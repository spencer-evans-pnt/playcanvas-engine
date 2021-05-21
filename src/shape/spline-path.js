import { Vec3 } from '../math/vec3.js';
import { SPLINE_CONTROL_POINT_MODE_FREE } from './constants';

class Spline {
    static getPoint(p0, p1, p2, p3, t) {
        t = Math.clamp01(t);
        const oneMinusT = 1 - t;
        const v1 = p0.clone().mulScalar(oneMinusT * oneMinusT * oneMinusT);
        const v2 = p1.clone().mulScalar(3 * oneMinusT * oneMinusT * t);
        const v3 = p2.clone().mulScalar(3 * oneMinusT * t * t);
        const v4 = p3.clone().mulScalar(t * t * t);
        return v1.add(v2).add(v3).add(v4);
    }

    static getFirstDerivative(p0, p1, p2, p3, t) {
        t = Math.clamp01(t);
        const oneMinusT = 1 - t;
        const v1 = (new Vec3().sub2(p1, p0)).mulScalar(3 * oneMinusT * oneMinusT);
        const v2 = (new Vec3().sub2(p2, p1)).mulScalar(6 * oneMinusT * t);
        const v3 = (new Vec3().sub2(p3, p2)).mulScalar(3 * t * t);
        return v1.add(v2).add(v3);
    }
}

class SplinePoint {
    constructor(c1 = new Vec3(),
        p = new Vec3(),
        c2 = new Vec3(),
        m = SPLINE_CONTROL_POINT_MODE_FREE) {
        this.control1 = c1;
        this.point = p;
        this.control2 = c2;
        this.mode = m;
    }
}

class SplinePath {
    constructor(points = [], loop = false) {
        this.points = points;
        this.loop = loop;
    }

    get pointCount() {
        return this.points.length;
    }

    getPoint(t) {
        let i1, i2;
        if (t >= 1) {
            if (this.loop) {
                t  *= (this.pointCount - 1);
                i1 = this.pointCount - 1;
                i2 = 0;
                t -= i1;
            } else {
                t = 1;
                i1 = this.pointCount - 2;
                i2 = i1 + 1;
            }
        } else {
            t  = Math.clamp01(t) * (this.pointCount - 1);
            i1 = Math.floor(t);
            i2 = i1 + 1;
            t -= i1;
        }
        const position = Spline.getPoint(
            this.points[i1].point,
            this.points[i1].control2,
            this.points[i2].control1,
            this.points[i2].point, t);
        return position;
    }

    getDirection(t) {
        let i1, i2;
        if (t >= 1) {
            if (this.loop) {
                t *= (this.pointCount - 1);
                i1 = this.pointCount - 1;
                i2 = 0;
                t -= i1;
            } else {
                t = 1;
                i1 = this.pointCount - 2;
                i2 = i1 + 1;
            }
        } else {
            t = Math.clamp01(t) * (this.pointCount - 1);
            i1 = Math.floor(t);
            i2 = i1 + 1;
            t -= i1;
        }
        const direction = Spline.getFirstDerivative(
            this.points[i1].point,
            this.points[i1].control2,
            this.points[i2].control1,
            this.points[i2].point, t);
        return direction.normalize();
    }

    getControlPointMode(index) {
        return this.points[index].mode;
    }

    getControlPoint(index, control = -1) {
        switch (control) {
            case 1: return this.points[index].control1.clone();
            case 2: return this.points[index].control2.clone();
            default: return this.points[index].point.clone();
        }
    }

    // TODO

    addPoint(index) {}

    deletePoint(index = -1) {}

    setControlPoint(index, point, control = -1) {}

    alignPointControls(index, control = -1) {}

    reset() {}
}

export { SplinePath, SplinePoint };
