import { Vec3 } from '../math/vec3.js';
import { SPLINE_CONTROL_POINT_MODE_FREE,
    SPLINE_CONTROL_POINT_MODE_ALIGNED,
    SPLINE_CONTROL_POINT_MODE_MIRRORED } from './constants';

class Spline {
    static getPoint(p0, p1, p2, p3, t) {
        t = Math.clamp01(t);
        const oneMinusT = 1 - t;
        /*
        // De Casteljau's algorithm implementation
        const q0 = new Vec3().add2(p0.clone().mulScalar(oneMinusT), p1.clone().mulScalar(t));
        const q1 = new Vec3().add2(p1.clone().mulScalar(oneMinusT), p2.clone().mulScalar(t));
        const q2 = new Vec3().add2(p2.clone().mulScalar(oneMinusT), p3.clone().mulScalar(t));
        const r0 = new Vec3().add2(q0.mulScalar(oneMinusT), q1.clone().mulScalar(t));
        const r1 = new Vec3().add2(q1.clone().mulScalar(oneMinusT), q2.mulScalar(t));
        const b = new Vec3().add2(r0.mulScalar(oneMinusT), r1.mulScalar(t));
        return b;

        // De Casteljau's algorithm simplifies down to this for cubic bezier curve:
        */
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
    constructor(points = []) {
        this.points = points;
        this._chordLengths = [];
        this._reevaluateLength();
    }

    get pointCount() {
        return this.points.length;
    }

    get length() {
        return this._chordLengths.reduce((sum, len) => sum + len);
    }

    _reevaluateLength(method = 3) {
        this._chordLengths = [];
        var distance = new Vec3();
        for (var i = 1; i < this.points.length; i++) {
            var p1 = this.points[i - 1];
            var p2 = this.points[i];

            var length;
            switch (method) {
                case 0: {
                    // linear
                    length = distance.sub2(p1.point, p2.point).length();
                    break;
                }
                case 1: {
                    // break into linear segments
                    const iterations = 32;
                    let pLast = Spline.getPoint(
                        p1.point,
                        p1.control2,
                        p2.control1,
                        p2.point, 0);
                    let pCurrent;
                    let totalLen = 0;
                    for (let ii = 1; ii <= iterations; ii++) {
                        const t = ii / iterations;
                        pCurrent = Spline.getPoint(
                            p1.point,
                            p1.control2,
                            p2.control1,
                            p2.point, t);
                        const len = distance.sub2(pCurrent, pLast).length();
                        totalLen += len;
                        pLast = pCurrent.clone();
                    }
                    length = totalLen;
                    break;
                }
                case 2: {
                    // average chord and sum of net lengths
                    var chordLen = distance.sub2(p1.point, p2.point).length();
                    var s1Len = distance.sub2(p1.point, p1.control2).length();
                    var s2Len = distance.sub2(p1.control2, p2.control1).length();
                    var s3Len = distance.sub2(p2.control1, p2.point).length();
                    var netLen = s1Len + s2Len + s3Len;
                    length = (chordLen + netLen) / 2;
                    break;
                }
                case 3: {
                    // simpson's rule - use second degree polynomial
                    length = this._simpsonsArcLength(p1.point, p1.control2, p2.control1, p2.point);
                    break;
                }
            }

            this._chordLengths.push(length);
        }
    }

    _pointsForT(t0) {
        let i1, i2;
        t0 = Math.clamp01(t0);

        const totalLen = this.length;
        const lenT = totalLen * t0;
        let sum = 0;

        for (var i = 0; i < this._chordLengths.length; i++) {
            const cLen = this._chordLengths[i];
            if (lenT <= sum + cLen) {
                i1 = i;
                i2 = i + 1;
                t0 = (lenT - sum) / cLen;
                break;
            }
            sum += cLen;
        }
        return { i1, i2, t: t0 };
    }

    getPoint(t0) {
        const { i1, i2, t } = this._pointsForT(t0);
        const position = Spline.getPoint(
            this.points[i1].point,
            this.points[i1].control2,
            this.points[i2].control1,
            this.points[i2].point, t);
        return position;
    }

    getPointNewtons(d, totalLen, debug) {
        let tLast = d / totalLen;
        let tNext = tLast;

        let iterations = 0;
        const maxIterations = 16;

        if (debug) console.group(`searching for constant speed t with desired distance ${d.toFixed(4)} and total length ${totalLen.toFixed(4)}; initial t: ${tLast.toFixed(4)}`);

        // use bisection to get close
        let tMax = 1;
        let tMin = 0;
        let maxLen = this.getApproxArcLength(0, tMax);
        let minLen = this.getApproxArcLength(0, tMin);
        while (maxLen - minLen > 0.1) {
            const tMid = (tMax + tMin) / 2;
            const midLen = this.getApproxArcLength(0, tMid, debug);
            if (d - midLen < 0) {
                tMax = tMid;
                maxLen = this.getApproxArcLength(0, tMax);
            } else {
                tMin = tMid;
                minLen = this.getApproxArcLength(0, tMin);
            }
        }
        tNext = (tMax + tMin) / 2;
        if (debug) console.log(`found tNext ${tNext.toFixed(4)} using bisection`);

        // use newton's method to find a more accurate t
        do {
            tLast = tNext;
            tNext = tLast -
                ((this.getApproxArcLength(0, tLast, debug) - d) /
                 this.getArcLengthIntegrand(tLast, debug));
            iterations++;
            if (debug) console.log(`iteration ${iterations}: ${tLast.toFixed(4)}, ${tNext.toFixed(4)}`);
        } while (Math.abs(tLast - tNext) > 0.01 && iterations < maxIterations);

        if (debug) {
            if (iterations >= maxIterations) {
                console.warn(`failed to find suitable t for ${d.toFixed(4)}; last guess: ${tLast.toFixed(4)}, ${tNext.toFixed(4)}`);
            } else {
                console.log(`found constant speed t ${tLast.toFixed(4)} using newton's method after ${iterations} iterations`);
            }
            console.groupEnd();
        }

        // revert to our less accurate bisection result
        if (iterations === maxIterations) {
            tLast = (tMax + tMin) / 2;
        }

        const { i1, i2, t } = this._pointsForT(tLast);
        const position = Spline.getPoint(
            this.points[i1].point,
            this.points[i1].control2,
            this.points[i2].control1,
            this.points[i2].point, t);
        return position;
    }

    getDirection(t0) {
        const { i1, i2, t } = this._pointsForT(t0);
        const direction = Spline.getFirstDerivative(
            this.points[i1].point,
            this.points[i1].control2,
            this.points[i2].control1,
            this.points[i2].point, t);
        return direction.normalize();
    }

    getArcLengthIntegrand(t0, debug) {
        const { i1, i2, t } = this._pointsForT(t0);
        const len = this._arcLengthIntegrand(
            this.points[i1].point,
            this.points[i1].control2,
            this.points[i2].control1,
            this.points[i2].point, t);
        if (debug) console.log(`\tarc length integrand for ${t0.toFixed(4)} uses chords ${i1} and ${i2}, yields length ${len.toFixed(4)}`);
        return len;
    }

    getApproxArcLength(tStart, tEnd, debug = false) {
        const { i1: starti1, i2: starti2, t: startt } = this._pointsForT(tStart);
        const { i1: endi1, i2: endi2, t: endt } = this._pointsForT(tEnd);

        if (starti1 === endi1) {
            // both points are within the same chord
            const len = this._simpsonsArcLength(
                this.points[starti1].point,
                this.points[starti1].control2,
                this.points[starti2].control1,
                this.points[starti2].point,
                startt, endt
            );
            if (debug) console.log(`\t\ttEnd ${tEnd.toFixed(4)} requires a single chord\n` +
                                   `\t\tchord uses points ${starti1} and ${starti2}, with t values [${startt.toFixed(4)}, ${endt.toFixed(4)}] yields length ${len.toFixed(4)}`);
            return len;
        } else if (starti2 === endi1) {
            // points are on separate chords; combine
            const c1Len = this._simpsonsArcLength(
                this.points[starti1].point,
                this.points[starti1].control2,
                this.points[starti2].control1,
                this.points[starti2].point,
                startt, 1
            );
            const c2Len = this._simpsonsArcLength(
                this.points[endi1].point,
                this.points[endi1].control2,
                this.points[endi2].control1,
                this.points[endi2].point,
                0, endt
            );
            if (debug) {
                console.log(`\t\ttEnd ${tEnd.toFixed(4)} requires multiple chords\n` +
                            `\t\tchord 1 points ${starti1} and ${starti2}, with t values [${startt.toFixed(4)}, 1] yields length ${c1Len.toFixed(4)}\n` +
                            `\t\tchord 2 points ${endi1} and ${endi2}, with t values [0, ${endt.toFixed(4)}] yields length ${c2Len.toFixed(4)}\n` +
                            `\t\ttotal length: ${(c1Len + c2Len).toFixed(4)}`);
            }
            return c1Len + c2Len;
        }

        throw new Error(`Cannot reconcile arc length between ${tStart.toFixed(4)} and ${tEnd.toFixed(4)}`);
    }

    _segmentedArcLength(p0, p1, p2, p3, tStart = 0, tEnd = 1, resolution = 16) {
        tStart = Math.clamp01(tStart);
        tEnd = Math.clamp01(tEnd);

        var segments = Math.ceil((tEnd - tStart) * resolution);
        const distance = new Vec3();
        let pLast = Spline.getPoint(p0, p1, p2, p3, tStart);
        let pCurrent;
        let totalLen = 0;

        for (let ii = 1; ii <= segments; ii++) {
            const t = Math.min(tEnd, tStart + ii / segments);
            pCurrent = Spline.getPoint(p0, p1, p2, p3, t);
            const len = distance.sub2(pCurrent, pLast).length();
            totalLen += len;
            pLast = pCurrent.clone();
        }
        return totalLen;
    }

    _simpsonsArcLength(p0, p1, p2, p3, tStart = 0, tEnd = 1, resolution = 16) {
        // start dividing into sections
        const delta = (tEnd - tStart) / resolution;

        // length of first derivative is arc length integrand
        const endPoints = this._arcLengthIntegrand(p0, p1, p2, p3, tStart) +
            this._arcLengthIntegrand(p0, p1, p2, p3, tEnd);

        // Everything multiplied by 4
        let x4 = 0;
        for (let i = 1; i < resolution; i += 2) {
            const t = tStart + delta * i;
            x4 += this._arcLengthIntegrand(p0, p1, p2, p3, t);
        }

        // Everything multiplied by 2
        let x2 = 0;
        for (let i = 2; i < resolution; i += 2) {
            const t = tStart + delta * i;
            x2 += this._arcLengthIntegrand(p0, p1, p2, p3, t);
        }

        // The final length
        const length = (delta / 3) * (endPoints + 4 * x4 + 2 * x2);

        return length;
    }

    _arcLengthIntegrand(p0, p1, p2, p3, t) {
        return Spline.getFirstDerivative(p0, p1, p2, p3, t).length();
    }

    // helpers

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

    deletePoint(index = -1) {
        if (this.points.length <= 2) {
            return;
        }

        this.points.splice(index, 1);
    }

    setControlPoint(index, point, control = -1) {
        switch (control) {
            case 1:
                this.points[index].control1 = point;
                break;

            case 2:
                this.points[index].control2 = point;
                break;

            default: {
                // move the controls along with the point
                const delta = new Vec3().sub2(point, this.points[index].point);
                this.points[index].control1.add(delta);
                this.points[index].control2.add(delta);
                this.points[index].point = point;
                break;
            }
        }

        this.alignPointControls(index, control);
    }

    alignPointControls(index, control = -1) {
        // invalid index
        if (index == -1 || index >= this.points.length) {
            return;
        }

        // don't align free mode control points
        const mode = this.points[index].Mode;
        if (mode === SPLINE_CONTROL_POINT_MODE_FREE) {
            return;
        }

        // determine which control point to move
        var fixedPoint, enforcedPoint;
        if (control === 1) {
            fixedPoint = this.points[index].control1;
            enforcedPoint = this.points[index].control2;
        } else {
            fixedPoint = this.points[index].control2;
            enforcedPoint = this.points[index].control1;
        }

        // the point position, and the point control position/direction
        const middlePoint = this.points[index].point;
        const enforcedTangent = new Vec3().sub2(middlePoint, fixedPoint);

        // keep the distance of the control point if Aligned control point mode
        if (mode === SPLINE_CONTROL_POINT_MODE_ALIGNED) {
            enforcedTangent.normalize().mulScalar(new Vec3().sub2(middlePoint, enforcedPoint).length());
        }

        // now move the control position
        if (control == 1) {
            this.points[index].control2 = new Vec3().add2(middlePoint, enforcedTangent);
        } else {
            this.points[index].control1 = new Vec3().add2(middlePoint, enforcedTangent);
        }
    }

    reset() {
        this.points = [];

        this.points.push(new SplinePoint(
            Vec3.ZERO,
            new Vec3(-1, 0, 0),
            new Vec3(-0.5, 0, 1),
            SPLINE_CONTROL_POINT_MODE_MIRRORED,
        ));
        this.points.push(new SplinePoint(
            new Vec3(0.5, 0, 1),
            new Vec3(1, 0, 0),
            Vec3.ZERO,
            SPLINE_CONTROL_POINT_MODE_MIRRORED
        ));

        this.alignPointControls(0, 2);
        this.alignPointControls(1, 1);
    }

    clone() {
        const points = this.points.map((p) =>
            new SplinePoint(
                p.control1.clone(),
                p.point.clone(),
                p.control2.clone(),
                p.mode
            )
        );
        return new SplinePath(points);
    }
}

export { SplinePath, SplinePoint };
