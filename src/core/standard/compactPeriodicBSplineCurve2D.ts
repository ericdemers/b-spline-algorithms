import { BSpline } from './bspline';
import { RealVectorSpace, Vector2D, VectorSpace } from "./vector-space";
import { PeriodicKnots } from './periodicKnots';
import { ControlNet } from './control-net';
import { KnotValue } from './knot-structure';

/**
 * Vector space implementation for 2D vectors
 */
class Vector2DSpace implements VectorSpace<number, Vector2D> {
    add(a: Vector2D, b: Vector2D): Vector2D {
        return [a[0] + b[0], a[1] + b[1]];
    }

    multiply(a: number, v: Vector2D): Vector2D {
        return [a * v[0], a * v[1]];
    }

    zero(): Vector2D {
        return [0, 0];
    }

    scale(scalar: number, v: Vector2D): Vector2D {
        return (v as number[]).map(val => scalar * val) as Vector2D;
    }

    subtract(a: Vector2D, b: Vector2D): Vector2D {
        return (a as number[]).map((val, i) => val - (b as number[])[i]) as Vector2D;
    }

    dimension(): number {
        return 2;
    }
}

/**
 * Represents a control net for periodic B-spline curves
 */
class PeriodicControlNet implements ControlNet<Vector2D> {
    constructor(private readonly basePoints: ReadonlyArray<Vector2D>) {
        if (!basePoints || basePoints.length < 2) {
            throw new Error('At least 2 control points are required');
        }
    }

    getBasePoints(): ReadonlyArray<Vector2D> {
        return this.basePoints;
    }

    getDimension(): number {
        return 1; // 1D curve
    }

    getSize(direction: number): number {
        if (direction !== 0) {
            throw new Error('Invalid direction');
        }
        return this.basePoints.length;
    }

    getPoint(indices: ReadonlyArray<number>): Vector2D {
        if (indices.length !== 1) {
            throw new Error('Invalid number of indices');
        }
        const index = indices[0];
        // Handle periodic indexing
        const normalizedIndex = ((index % this.basePoints.length) + this.basePoints.length) % this.basePoints.length;
        return this.basePoints[normalizedIndex];
    }


}

/**
 * Represents a compact periodic B-spline curve in 2D
 * Extends the base BSpline class with periodic behavior
 */
export class CompactPeriodicBSplineCurve2D extends BSpline<number, Vector2D> {
    private readonly patternLength: number;
    private readonly patternInterval: number;
    private readonly periodicKnots: PeriodicKnots;

    /**
     * Creates a new compact periodic B-spline curve
     * 
     * @param knotPattern - The basic knot pattern that repeats
     * @param controlPoints - Array of control points
     * @param degree - Degree of the B-spline curve
     */
    constructor(
        knotPattern: ReadonlyArray<number>,
        controlPoints: ReadonlyArray<Vector2D>,
        degree: number
    ) {
        // Create vector space for 2D points
        const vectorSpace = new Vector2DSpace();

        // Create periodic control net
        const periodicControlNet = new PeriodicControlNet(controlPoints);

        // Calculate period from knot pattern
        const period = knotPattern[knotPattern.length - 1] - knotPattern[0];

        // Create periodic knot structure
        const periodicKnots = new PeriodicKnots(knotPattern, period);

        // Initialize base class
        super(vectorSpace, periodicControlNet, periodicKnots, [degree]);

        this.patternLength = knotPattern.length;
        this.patternInterval = period;
        this.periodicKnots = periodicKnots;
    }

    /**
     * Evaluates the curve at parameter u
     * @param parameters - Array containing single parameter value
     */
    override evaluate(parameters: ReadonlyArray<number>): Vector2D {
        if (parameters.length !== 1) {
            throw new Error('Expected exactly one parameter');
        }

        const u = parameters[0];
        const normalizedU = this.normalizeParameter(u);
        return super.evaluate([normalizedU]);
    }

    /**
     * Gets the domain of the curve
     */
    getDomain(): { min: number; max: number } {
        return this.periodicKnots.getDomain(0);
    }

    /**
     * Gets the period of the curve
     */
    getPeriod(): number {
        return this.patternInterval;
    }

    /**
     * Gets distinct knots with their multiplicities
     */
    getDistinctKnots(): ReadonlyArray<KnotValue> {
        return this.periodicKnots.getDistinctKnots(0);
    }

    /**
     * Creates a new curve with an inserted knot
     * @param u - Parameter value where to insert the knot
     */
    withInsertedKnot(u: number): CompactPeriodicBSplineCurve2D {
        const newKnotStructure = this.periodicKnots.withInsertedKnot(0, u) as PeriodicKnots;
        return new CompactPeriodicBSplineCurve2D(
            newKnotStructure.getPattern(),
            (this.controlNet as PeriodicControlNet).getBasePoints(),
            this.degrees[0]
        );
    }

    /**
     * Creates a new curve with removed knots
     */
    withRemovedKnots(): CompactPeriodicBSplineCurve2D {
        const newKnotStructure = this.periodicKnots.withRemovedKnots(0) as PeriodicKnots;
        return new CompactPeriodicBSplineCurve2D(
            newKnotStructure.getPattern(),
            (this.controlNet as PeriodicControlNet).getBasePoints(),
            this.degrees[0]
        );
    }

    /**
     * Normalizes a parameter value to the fundamental domain
     */
    private normalizeParameter(u: number): number {
        const { min, max } = this.getDomain();
        return min + ((u - min) % this.patternInterval + this.patternInterval) % this.patternInterval;
    }
}

// Example usage:
/*
const knotPattern = [0, 0, 0, 1, 2, 3, 3, 3];
const controlPoints = [
    new Vector2D(0, 0),
    new Vector2D(1, 1),
    new Vector2D(2, 0),
    new Vector2D(1, -1),
    new Vector2D(0, 0)
];
const degree = 2;

const periodicCurve = new CompactPeriodicBSplineCurve2D(knotPattern, controlPoints, degree);

// Evaluate at some parameter
const point = periodicCurve.evaluate([1.5]);

// Insert a knot
const refinedCurve = periodicCurve.withInsertedKnot(1.5);

// Get domain information
const domain = periodicCurve.getDomain();
const period = periodicCurve.getPeriod();
*/
