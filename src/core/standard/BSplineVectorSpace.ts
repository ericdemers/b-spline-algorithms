

/*
const splineCurve2D = new BSpline(
    new Vector2DSpace(),
    new ControlPolygonCurve(controlPoints),
    new KnotsCurve([0, 1, 2, 3, 4, 5, 6]),
    new Vector2DSpace()
)
*/

// Type definitions
/*
export type Point1D = number;
export type Point2D = [number, number];
export type Point3D = [number, number, number];
export type Point4D = [number, number, number, number];
export type Complex = [number, number];
export type Point = Point1D | Point2D | Point3D | Point4D | Complex;


export class Index {}

export class IndexCurve extends Index {
    private _i: Readonly<number>

    constructor(i: Readonly<number>) {
        super()
        this._i = i
    }

    public get i(): Readonly<number> {
        return this._i
    }

}

export class IndexSurface extends Index {
    private _i: Readonly<number>
    private _j: Readonly<number>

    constructor(i: Readonly<number>, j: Readonly<number>) {
        super()
        this._i = i
        this._j = j
    }

    public get i(): Readonly<number> {
        return this._i
    }

    public get j(): Readonly<number> {
        return this._j
    }

    
}

export class IndexVolume extends Index {
    private _i: Readonly<number>
    private _j: Readonly<number>
    private _k: Readonly<number>

    constructor(i: Readonly<number>, j: Readonly<number>, k: Readonly<number>) {
        super()
        this._i = i
        this._j = j
        this._k = k
    }

    public get i(): Readonly<number> {
        return this._i
    }

    public get j(): Readonly<number> {
        return this._j
    }

    public get k(): Readonly<number> {
        return this._k
    }

}

export class ParametricPosition {}

export class ParametricPositionCurve extends ParametricPosition {
    private u: Readonly<number>

    constructor (u: Readonly<number>) {
        super()
        this.u = u
    }
}

export class ParametricPositionSurface extends ParametricPosition {
    private u: Readonly<number>
    private v: Readonly<number>

    constructor (u: Readonly<number>, v: Readonly<number>) {
        super()
        this.u = u
        this.v = v
    }
}

export class ParametricPositionVolume extends ParametricPosition {
    private u: Readonly<number>
    private v: Readonly<number>
    private w: Readonly<number>

    constructor (u: Readonly<number>, v: Readonly<number>, w: Readonly<number>) {
        super()
        this.u = u
        this.v = v
        this.w = w
    }
}

export class Knots {}

export class KnotsCurve extends Knots {
    private knots: ReadonlyArray<number>

    constructor(knots: ReadonlyArray<number>) {
        super()
        this.knots = knots
    }
}

export class KnotsSurface extends Knots {
    private knotsU: ReadonlyArray<number>
    private knotsV: ReadonlyArray<number>

    constructor(knotsU: number[], knotsV: number[]) {
        super()
        this.knotsU = knotsU
        this.knotsV = knotsV
    }
}

export class KnotsVolume extends Knots {
    private knotsU: ReadonlyArray<number>
    private knotsV: ReadonlyArray<number>
    private knotsW: ReadonlyArray<number>

    constructor(knotsU: number[], knotsV: number[], knotsW: number[]) {
        super()
        this.knotsU = knotsU
        this.knotsV = knotsV
        this.knotsW = knotsW
    }
}

export class BasisFunctions {

}

export interface ControlPolytope<P extends Point> {

    getPoint(index: Index): P 

}

export class ControlPolygonCurve<P extends Point> implements ControlPolytope<P> {
    private controlPoints: ReadonlyArray<P>
    private weights?: ReadonlyArray<P>

    constructor(controlPoints: ReadonlyArray<P>, weights?: ReadonlyArray<number> ) {
        this.controlPoints = controlPoints
    }

    getPoint(index: IndexCurve): P {
        return this.controlPoints[index.i]
    }
}

export class ControlPolygonSurface<P extends Point> extends ControlPolytope<P> {
    private controlPoints: ReadonlyArray<ReadonlyArray<P>>
    private weights?: ReadonlyArray<ReadonlyArray<P>>

    constructor(controlPoints: ReadonlyArray<ReadonlyArray<P>>, weights?: ReadonlyArray<ReadonlyArray<P>> ) {
        super()
        this.controlPoints = controlPoints
    }

    getPoint(index: IndexSurface): P {
        return this.controlPoints[index.i][index.j]
    }
}

export class ControlPolygonVolume<P extends Point> extends ControlPolytope<P> {
    private controlPoints: ReadonlyArray<ReadonlyArray<ReadonlyArray<P>>>
    private weights?: Readonly<ReadonlyArray<ReadonlyArray<P>>>

    constructor(controlPoints: ReadonlyArray<ReadonlyArray<ReadonlyArray<P>>>, weights?: ReadonlyArray<ReadonlyArray<ReadonlyArray<P>>>) {
        super()
        this.controlPoints = controlPoints
    }

    getPoint(index: IndexVolume): P {
        return this.controlPoints[index.i][index.j][index.k]
    }
}
*/

/*
export class BSpline<P extends Point> {

    private readonly controlPolygon: ControlPolytope<P>
    private readonly knots: Knots

    constructor(controlPolygon: ControlPolytope<P>, knots: Knots) {
        this.controlPolygon = controlPolygon
        this.knots = knots
    }

    controlPolygonPoint(index: Index): P {
        return this.controlPolygon.getPoint(index)
    }

    evaluate(p: ParametricPosition): P {
        // Evaluate the position using a pyramid algorithms from Ron Goldman book and Category theory concept

    }

}
*/

/*
// First, let's define a type class for vector spaces
interface VectorSpace<V> {
    zero: () => V;
    add: (a: V, b: V) => V;
    scale: (scalar: number, v: V) => V;
}

// A general affine combination type that works in any dimension
interface AffineCombinator<V> {
    combine: (points: V[], weights: number[]) => V;
}

// Implementation of the B-spline evaluation using these abstractions
class BSpline<V> {
    private vectorSpace: VectorSpace<V>;
    private controlPoints: V[];
    private knots: number[];
    private degree: number;

    constructor(
        vectorSpace: VectorSpace<V>,
        controlPoints: V[],
        knots: number[],
        degree: number
    ) {
        this.vectorSpace = vectorSpace;
        this.controlPoints = controlPoints;
        this.knots = knots;
        this.degree = degree;
    }

    // Affine combinator that works in any dimension
    private makeAffineCombinator(): AffineCombinator<V> {
        return {
            combine: (points: V[], weights: number[]): V => {
                return points.reduce((acc, point, i) => 
                    this.vectorSpace.add(
                        acc,
                        this.vectorSpace.scale(weights[i], point)
                    ),
                    this.vectorSpace.zero()
                );
            }
        };
    }

    evaluate(t: number): V {
        const combinator = this.makeAffineCombinator();
        const span = this.findSpanIndex(t);
        
        // Create the pyramid structure
        const pyramid: V[][] = Array(this.degree + 1)
            .fill(null)
            .map(() => Array(this.degree + 1));

        // Initialize with control points
        for (let i = 0; i <= this.degree; i++) {
            pyramid[0][i] = this.controlPoints[span - this.degree + i];
        }

        // Compute the pyramid using the affine combinator
        for (let r = 1; r <= this.degree; r++) {
            for (let i = 0; i <= this.degree - r; i++) {
                const alpha = (t - this.knots[span - this.degree + i + r]) /
                            (this.knots[span + i + 1] - this.knots[span - this.degree + i + r]);
                
                pyramid[r][i] = combinator.combine(
                    [pyramid[r-1][i], pyramid[r-1][i+1]],
                    [1 - alpha, alpha]
                );
            }
        }

        return pyramid[this.degree][0];
    }

    // ... rest of the implementation
}

// Example implementations for different dimensions:

// For 2D points
interface Point2D {
    x: number;
    y: number;
}

const point2DVectorSpace: VectorSpace<Point2D> = {
    zero: () => ({ x: 0, y: 0 }),
    add: (a, b) => ({ x: a.x + b.x, y: a.y + b.y }),
    scale: (s, v) => ({ x: s * v.x, y: s * v.y })
};

// For 3D points
interface Point3D {
    x: number;
    y: number;
    z: number;
}

const point3DVectorSpace: VectorSpace<Point3D> = {
    zero: () => ({ x: 0, y: 0, z: 0 }),
    add: (a, b) => ({ x: a.x + b.x, y: a.y + b.y, z: a.z + b.z }),
    scale: (s, v) => ({ x: s * v.x, y: s * v.y, z: s * v.z })
};

// Usage example:
const bspline2D = new BSpline<Point2D>(point2DVectorSpace, controlPoints2D, knots, 3);
const bspline3D = new BSpline<Point3D>(point3DVectorSpace, controlPoints3D, knots, 3);
*/

/*
// Functor type class
interface Functor<F> {
    map: <A, B>(fa: F<A>, f: (a: A) => B) => F<B>;
}

// Applicative type class
interface Applicative<F> extends Functor<F> {
    pure: <A>(a: A) => F<A>;
    ap: <A, B>(ff: F<(a: A) => B>, fa: F<A>) => F<B>;
}
*/

/*
// Vector space operations (only for vectors, not points)
interface VectorSpace<V> {
    zero: () => V;
    add: (a: V, b: V) => V;
    scale: (scalar: number, v: V) => V;
}

// Affine space operations
interface AffineSpace<P, V> {
    // Subtract two points to get a vector
    diff: (a: P, b: P) => V;
    // Translate a point by a vector
    translate: (p: P, v: V) => P;
    // Affine combination of points
    combine: (points: P[], weights: number[]) => P;
}

// Implementation using Affine Space
class BSpline<P, V> {
    private affineSpace: AffineSpace<P, V>;
    private vectorSpace: VectorSpace<V>;
    private controlPoints: P[];
    private knots: number[];
    private degree: number;

    constructor(
        affineSpace: AffineSpace<P, V>,
        vectorSpace: VectorSpace<V>,
        controlPoints: P[],
        knots: number[],
        degree: number
    ) {
        this.affineSpace = affineSpace;
        this.vectorSpace = vectorSpace;
        this.controlPoints = controlPoints;
        this.knots = knots;
        this.degree = degree;
    }

    evaluate(t: number): P {
        const span = this.findSpanIndex(t);
        
        // Create the pyramid structure
        const pyramid: P[][] = Array(this.degree + 1)
            .fill(null)
            .map(() => Array(this.degree + 1));

        // Initialize with control points
        for (let i = 0; i <= this.degree; i++) {
            pyramid[0][i] = this.controlPoints[span - this.degree + i];
        }

        // Compute the pyramid using affine combinations
        for (let r = 1; r <= this.degree; r++) {
            for (let i = 0; i <= this.degree - r; i++) {
                const alpha = (t - this.knots[span - this.degree + i + r]) /
                            (this.knots[span + i + 1] - this.knots[span - this.degree + i + r]);
                
                pyramid[r][i] = this.affineSpace.combine(
                    [pyramid[r-1][i], pyramid[r-1][i+1]],
                    [1 - alpha, alpha]
                );
            }
        }

        return pyramid[this.degree][0];
    }
}

// Example implementation for 2D
interface Point2D {
    x: number;
    y: number;
}

interface Vector2D {
    dx: number;
    dy: number;
}

const vector2DSpace: VectorSpace<Vector2D> = {
    zero: () => ({ dx: 0, dy: 0 }),
    add: (a, b) => ({ dx: a.dx + b.dx, dy: a.dy + b.dy }),
    scale: (s, v) => ({ dx: s * v.dx, dy: s * v.dy })
};

const affine2DSpace: AffineSpace<Point2D, Vector2D> = {
    diff: (a, b) => ({
        dx: a.x - b.x,
        dy: a.y - b.y
    }),
    translate: (p, v) => ({
        x: p.x + v.dx,
        y: p.y + v.dy
    }),
    combine: (points, weights) => {
        if (points.length === 0 || weights.length !== points.length) {
            throw new Error("Invalid points or weights");
        }

        // Check if weights sum to 1 (within numerical precision)
        const weightSum = weights.reduce((sum, w) => sum + w, 0);
        if (Math.abs(weightSum - 1) > 1e-10) {
            throw new Error("Weights must sum to 1");
        }

        // Compute weighted sum
        return points.reduce((acc, p, i) => ({
            x: acc.x + p.x * weights[i],
            y: acc.y + p.y * weights[i]
        }), { x: 0, y: 0 });
    }
};

// Example usage
const bspline2D = new BSpline<Point2D, Vector2D>(
    affine2DSpace,
    vector2DSpace,
    controlPoints2D,
    knots,
    3
);

*/