// ------------ Basic Examples ------------

import { BSpline} from "./bspline"
import { ControlNet, ControlPolygonCurve } from "./control-net"
import { Knots, KnotStructure } from "./knot-structure"
import { Complex, ComplexVector2D, Real, Vector, Vector1D, Vector2D, Vector3D, VectorSpace } from "./vector-space"



class ProductKnots implements KnotStructure {
    constructor(private readonly knots: ReadonlyArray<ReadonlyArray<number>>) {}

    getDimension(): number {
        return this.knots.length;
    }

    getKnotSequence(direction: number): ReadonlyArray<number> {
        return this.knots[direction];
    }
}

export class ControlNetSurface<V extends Vector> implements ControlNet<V> {

    constructor(    
        private readonly controlPoints: ReadonlyArray<ReadonlyArray<V>>,
        private readonly weights?: ReadonlyArray<ReadonlyArray<V>>) {}

    getPoint(index: number[]): V {
        return this.controlPoints[index[0]][index[1]]
    }

    getDimension() {
        return 2
    }

    /** Get size in each parametric direction 
     * Return the length of the outer array if direction === 0
     * and return the length of the inner array if the direction === 1
    */
    getSize(direction: number) : number {
        if (direction === 0) {
            return this.controlPoints.length;
        } else if (direction === 1) {
            return this.controlPoints[0].length;
        } else {
            throw new Error("Invalid direction");
        }
    }
}

class ControlNetVolume<V extends Vector> implements ControlNet<V> {

    constructor(
        private readonly controlPoints: ReadonlyArray<ReadonlyArray<ReadonlyArray<V>>>,
        private readonly weights?: ReadonlyArray<ReadonlyArray<V>>) {}

    getPoint(index: number[]): V {
        return this.controlPoints[index[0]][index[1]][index[2]]
    }

    getDimension() {
        return 3
    }

    getSize(direction: number) : number {
        if (direction === 0) {
            return this.controlPoints.length;
        } else if (direction === 1) {
            return this.controlPoints[0].length;
        } else if (direction === 2) {
            return this.controlPoints[0][0].length;
        } else {
            throw new Error("Invalid direction");
        }
    }
}


class RealVectorSpace1D implements VectorSpace<Real, Vector1D> {
    zero(): Vector1D {
        return [0];
    }

    add(a: Vector1D, b: Vector1D): Vector1D {
        return [a[0] + b[0]];
    }

    scale(scalar: Real, v: Vector1D): Vector1D {
        return [scalar * v[0]];
    }

    dimension() {
        return 1
    }
}

class RealVectorSpace2D implements VectorSpace<Real, Vector2D> {
    zero(): Vector2D {
        return [0, 0];
    }

    add(a: Vector2D, b: Vector2D): Vector2D {
        return [a[0] + b[0], a[1] + b[1]];
    }

    scale(scalar: Real, v: Vector2D): Vector2D {
        return [scalar * v[0], scalar * v[1]];
    }

    dimension() {
        return 2
    }
}

class RealVectorSpace3D implements VectorSpace<Real, Vector3D> {
    zero(): Vector3D {
        return [0, 0, 0];
    }

    add(a: Vector3D, b: Vector3D): Vector3D {
        return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
    }

    scale(scalar: Real, v: Vector3D): Vector3D {
        return [scalar * v[0], scalar * v[1], scalar * v[2]];
    }

    dimension() {
        return 3
    }
}

class ComplexVectorSpace2D implements VectorSpace<Complex, ComplexVector2D> {
    zero(): ComplexVector2D {
        return [[0, 0], [0, 0]];
    }

    add(a: ComplexVector2D, b: ComplexVector2D): ComplexVector2D {
        return [
            [a[0][0] + b[0][0], a[0][1] + b[0][1]],
            [a[1][0] + b[1][0], a[1][1] + b[1][1]]
        ];
    }

    scale(scalar: Complex, v: ComplexVector2D): ComplexVector2D {
        return [
            [scalar[0] * v[0][0] - scalar[1]*v[0][1], scalar[0] * v[0][1] + scalar[1] * v[0][0]],
            [scalar[0] * v[1][0] - scalar[1]*v[1][1], scalar[0] * v[1][1] + scalar[1] * v[1][0]]
        ];
    }

    dimension() {
        return 2
    }
}

// 1. Simple 2D curve with real coordinates
const curve2DReal = () => {
    // Define control points in R²
    const controlPoints: Vector2D[] = [
        [0, 0],
        [1, 2],
        [3, 1],
        [4, 3]
    ];

    // Define uniform knot vector for cubic curve
    const knots = [0, 0, 0, 0, 1, 1, 1, 1];

    // Create control net
    const controlNet = new ControlPolygonCurve(controlPoints);

    // Create knot structure
    const knotStructure = new Knots(knots);

    // Create vector space
    const vectorSpace = new RealVectorSpace2D();

    // Create B-spline (cubic degree)
    const bspline = new BSpline(
        vectorSpace,
        controlNet,
        knotStructure,
        [3] // cubic degree
    );

    // Evaluate at different parameters
    const point1 = bspline.evaluate([0.5]);
    console.log('Point at t=0.5:', point1);

    // Evaluate multiple points for plotting
    const points = Array.from({length: 101}, (_, i) => i / 100)
        .map(t => bspline.evaluate([t]));
    
    return points;
};

// 2. Complex 2D curve
const curve2DComplex = () => {
    // Define control points in C²
    const controlPoints: ComplexVector2D[] = [
        [[0, 0], [0, 0]],     // (0 + 0i, 0 + 0i)
        [[1, 1], [0, 1]],     // (1 + i, 0 + i)
        [[2, -1], [1, 0]],    // (2 - i, 1 + 0i)
        [[3, 0], [1, -1]]     // (3 + 0i, 1 - i)
    ];

    // Create complex vector space
    const vectorSpace = new ComplexVectorSpace2D();

    // Create control net and knot structure
    const controlNet = new ControlPolygonCurve(controlPoints);
    const knotStructure = new Knots([0, 0, 0, 0, 1, 1, 1, 1]);

    // Create B-spline
    const bspline = new BSpline(
        vectorSpace,
        controlNet,
        knotStructure,
        [3]
    );

    // Evaluate curve
    const point = bspline.evaluate([0.5]);
    return point; // Returns complex coordinates
};

// 3. 3D Surface with real coordinates
const surface3DReal = () => {
    // Define control points grid in R³
    const controlPoints: Vector3D[][] = [
        [[0,0,0], [0,1,1], [0,2,0]],
        [[1,0,1], [1,1,2], [1,2,1]],
        [[2,0,0], [2,1,1], [2,2,0]]
    ];

    // Create control net for surface
    const controlNet = new ControlNetSurface(controlPoints);

    // Create knot structure for both directions
    const knotStructure = new ProductKnots([
        [0, 0, 0, 1, 1, 1], // u direction
        [0, 0, 0, 1, 1, 1]  // v direction
    ]);

    // Create vector space
    const vectorSpace = new RealVectorSpace3D();

    // Create B-spline surface
    const bspline = new BSpline(
        vectorSpace,
        controlNet,
        knotStructure,
        [2, 2] // quadratic in both directions
    );

    // Evaluate at point
    const point = bspline.evaluate([0.5, 0.5]);

    // Generate surface points for visualization
    const surfacePoints: Vector3D[][] = [];
    for (let u = 0; u <= 1; u += 0.1) {
        const row: Vector3D[] = [];
        for (let v = 0; v <= 1; v += 0.1) {
            row.push(bspline.evaluate([u, v]));
        }
        surfacePoints.push(row);
    }

    return surfacePoints;
};

// ------------ Advanced Examples ------------

// 4. Volumetric B-spline with real coordinates
const volume3DReal = () => {
    // Define 3D grid of control points
    const controlPoints: Vector3D[][][] = Array(3).fill(0).map((_, i) =>
        Array(3).fill(0).map((_, j) =>
            Array(3).fill(0).map((_, k) => [i, j, k])
        )
    );

    // Create control net for volume
    const controlNet = new ControlNetVolume(controlPoints);

    // Create knot structure for three directions
    const knotStructure = new ProductKnots([
        [0, 0, 0, 1, 1, 1],
        [0, 0, 0, 1, 1, 1],
        [0, 0, 0, 1, 1, 1]
    ]);

    // Create B-spline volume
    const bspline = new BSpline(
        new RealVectorSpace3D(),
        controlNet,
        knotStructure,
        [2, 2, 2] // quadratic in all directions
    );

    // Evaluate at point
    return bspline.evaluate([0.5, 0.5, 0.5]);
};

/*

// 5. Example with derivatives
const curveWithDerivatives = () => {
    // ... similar setup as curve2DReal
    const bspline = // ... create B-spline

    // Evaluate point and derivatives
    const derivatives = bspline.evaluateWithDerivatives([0.5], 2);
    console.log('Point:', derivatives[0]);
    console.log('First derivative:', derivatives[1]);
    console.log('Second derivative:', derivatives[2]);
};

// 6. Example with error handling
const errorHandlingExample = () => {
    try {
        const bspline = // ... create B-spline

        // Try to evaluate outside parameter range
        const point = bspline.evaluate([1.5]); // Should throw error
    } catch (error) {
        console.error('Evaluation error:', error.message);
    }

    try {
        // Try to create invalid B-spline
        const invalidBspline = new BSpline(
            vectorSpace,
            controlNet,
            knotStructure,
            [3, 2] // Wrong number of degrees
        );
    } catch (error) {
        console.error('Construction error:', error.message);
    }
};

// 7. Performance testing example
const performanceTest = () => {
    const bspline = // ... create B-spline
    
    console.time('B-spline evaluation');
    for (let i = 0; i < 1000; i++) {
        bspline.evaluate([i / 1000]);
    }
    console.timeEnd('B-spline evaluation');
};

*/
