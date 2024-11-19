import { Vector, VectorSpace } from "./vector-space";
import { BSpline } from "./bspline";
import { ControlPolygon } from "./controlNet";
import { Knots } from "./knotStructure";

/**
 * Implementation of a B-spline curve
 * Represents a parametric curve defined by control points, knot vector, and degree
 */
export class BSplineCurve<K extends number, V extends Vector> extends BSpline<K, V> {
    //private readonly vectorSpace: VectorSpace<K, V>;
    //private readonly controlPoints: V[];
    //private readonly knots: number[];
    //private readonly degree: number;

    /**
     * Constructs a new B-spline curve
     * 
     * @param vectorSpace - Vector space for control points
     * @param controlPoints - Array of control points
     * @param knots - Knot vector
     * @param degree - Degree of the B-spline
     */
    constructor(
        vectorSpace: VectorSpace<K, V>,
        controlPoints: V[],
        knots: number[],
        degree: number
    ) {
        /*
        this.vectorSpace = vectorSpace;
        this.controlPoints = [...controlPoints];
        this.knots = [...knots];
        this.degree = degree;
        this.validateConstructorParams();
        */
        // Convert to general form
        super(
            vectorSpace,
            new ControlPolygon(controlPoints),
            new Knots(knots ?? generateUniformKnots(controlPoints.length, degree)),
            [degree]
        );
    }

    /**
     * Evaluates the B-spline curve at parameter t
     * 
     * @param t - Parameter value
     * @returns Point on the curve
     */
    public evaluate(t: number): V {
        // Clamp parameter to domain
        t = Math.max(this.knots[this.degree], 
            Math.min(t, this.knots[this.knots.length - this.degree - 1]));

        // Find knot span
        const span = this.findSpanIndex(t);

        // Compute basis functions
        const basis = this.computeBasisFunctions(t, span);

        // Compute point
        let point = this.vectorSpace.zero();
        for (let i = 0; i <= this.degree; i++) {
            point = this.vectorSpace.add(
                point,
                this.vectorSpace.scale(basis[i] as K, this.controlPoints[span - this.degree + i])
            );
        }

        return point;
    }

    /**
     * Inserts a new knot into the B-spline curve
     * 
     * @param t - Parameter value for new knot
     * @returns New B-spline curve with inserted knot
     */
    public insertKnot(t: number): BSplineCurve<K, V> {
        // Validate knot insertion
        if (t < this.knots[0] || t > this.knots[this.knots.length - 1]) {
            throw new Error('Knot value outside of domain');
        }

        // Find span containing new knot
        const span = this.findSpanIndex(t);

        // Calculate alpha factors
        const alphas = this.calculateAlphaFactors(t, span);

        // Create new control points
        const newControlPoints: V[] = [];
        
        // Copy points before affected region
        for (let i = 0; i <= span - this.degree; i++) {
            newControlPoints.push(this.controlPoints[i]);
        }

        // Calculate new points in affected region
        for (let i = span - this.degree + 1; i <= span; i++) {
            const alpha = alphas[i - (span - this.degree)];
            const newPoint = this.vectorSpace.add(
                this.vectorSpace.scale(1 - alpha as K, this.controlPoints[i - 1]),
                this.vectorSpace.scale(alpha as K, this.controlPoints[i])
            );
            newControlPoints.push(newPoint);
        }

        // Copy remaining points
        for (let i = span + 1; i < this.controlPoints.length; i++) {
            newControlPoints.push(this.controlPoints[i]);
        }

        // Create new knot vector
        const newKnots = [...this.knots];
        newKnots.splice(span + 1, 0, t);

        return new BSplineCurve(
            this.vectorSpace,
            newControlPoints,
            newKnots,
            this.degree
        );
    }

    /**
     * Gets the derivative curve
     * 
     * @returns B-spline curve representing the derivative
     */
    public getDerivative(): BSplineCurve<K, V> {
        if (this.degree < 1) {
            throw new Error('Cannot derive curve of degree 0');
        }

        // Calculate derivative control points
        const derivativePoints: V[] = [];
        for (let i = 0; i < this.controlPoints.length - 1; i++) {
            const coefficient = this.degree / (this.knots[i + this.degree + 1] - this.knots[i + 1]);
            const diff = this.vectorSpace.subtract(
                this.controlPoints[i + 1],
                this.controlPoints[i]
            );
            derivativePoints.push(this.vectorSpace.scale(coefficient as K, diff));
        }

        // Create new knot vector
        const derivativeKnots = this.knots.slice(1, -1);

        return new BSplineCurve(
            this.vectorSpace,
            derivativePoints,
            derivativeKnots,
            this.degree - 1
        );
    }

    /**
     * Gets the control points
     */
    public getControlPoints(): V[] {
        return [...this.controlPoints];
    }

    /**
     * Gets the knot vector
     */
    public getKnots(): number[] {
        return [...this.knots];
    }

    /**
     * Gets the degree
     */
    public getDegree(): number {
        return this.degree;
    }

    /**
     * Validates constructor parameters
     */
    private validateConstructorParams(): void {
        // Check degree
        if (this.degree < 0) {
            throw new Error('Degree must be non-negative');
        }

        // Check number of control points
        if (this.controlPoints.length <= this.degree) {
            throw new Error(
                'Number of control points must be greater than degree'
            );
        }

        // Check knot vector length
        const requiredKnots = this.controlPoints.length + this.degree + 1;
        if (this.knots.length !== requiredKnots) {
            throw new Error(
                `Invalid number of knots. Expected ${requiredKnots}, got ${this.knots.length}`
            );
        }

        // Check knot vector is non-decreasing
        for (let i = 1; i < this.knots.length; i++) {
            if (this.knots[i] < this.knots[i - 1]) {
                throw new Error('Knot vector must be non-decreasing');
            }
        }

        // Check end knot multiplicities
        const startMultiplicity = this.countMultiplicity(this.knots[0]);
        const endMultiplicity = this.countMultiplicity(
            this.knots[this.knots.length - 1]
        );

        if (startMultiplicity < this.degree + 1 || endMultiplicity < this.degree + 1) {
            throw new Error(
                'End knots must have multiplicity at least degree + 1'
            );
        }
    }

    /**
     * Finds the knot span index containing parameter t
     */
    private findSpanIndex(t: number): number {
        // Special case for end of curve
        if (t >= this.knots[this.knots.length - this.degree - 1]) {
            return this.controlPoints.length - 1;
        }

        // Binary search
        let low = this.degree;
        let high = this.knots.length - this.degree - 1;
        let mid = Math.floor((low + high) / 2);

        while (t < this.knots[mid] || t >= this.knots[mid + 1]) {
            if (t < this.knots[mid]) {
                high = mid;
            } else {
                low = mid;
            }
            mid = Math.floor((low + high) / 2);
        }

        return mid;
    }

    /**
     * Computes basis functions for parameter t
     */
    private computeBasisFunctions(t: number, span: number): number[] {
        const N: number[] = new Array(this.degree + 1).fill(0);
        const left: number[] = new Array(this.degree + 1).fill(0);
        const right: number[] = new Array(this.degree + 1).fill(0);

        // Initialize degree 0
        N[0] = 1.0;

        // Compute basis functions
        for (let j = 1; j <= this.degree; j++) {
            left[j] = t - this.knots[span + 1 - j];
            right[j] = this.knots[span + j] - t;
            let saved = 0.0;

            for (let r = 0; r < j; r++) {
                const temp = N[r] / (right[r + 1] + left[j - r]);
            }
        }
        return N;
    }

    /**
     * Calculates alpha factors for knot insertion
     */
    private calculateAlphaFactors(t: number, span: number): number[] {
        const alphas: number[] = new Array(this.degree + 1);
        
        for (let i = 0; i <= this.degree; i++) {
            const numerator = t - this.knots[span - this.degree + i];
            const denominator = this.knots[span + 1 + i - this.degree] - 
                              this.knots[span - this.degree + i];
            alphas[i] = denominator === 0 ? 0 : numerator / denominator;
        }
        
        return alphas;
    }

    /**
     * Counts multiplicity of a knot value
     */
    private countMultiplicity(value: number): number {
        return this.knots.filter(k => Math.abs(k - value) < Number.EPSILON).length;
    }
}

