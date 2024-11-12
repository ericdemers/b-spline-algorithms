import { ControlNet, ControlPolygonCurve } from "./control-net";
import { KnotStructure, ProductKnots } from "./knot-structure";
import { Complex, Scalar, Vector, VectorSpace } from "./vector-space";

// ------------ B-spline Implementation ------------

/**
 * B-spline class implementing the core algorithm
 * Supports arbitrary dimension, real and complex vector spaces
 * 
 * @typeParam K - Scalar type (Real or Complex)
 * @typeParam V - Vector type
 */
export class BSpline<K extends Scalar, V extends Vector> {
    /**
     * Creates a new B-spline instance
     * 
     * @param vectorSpace - The vector space implementation
     * @param controlNet - The control net defining the B-spline
     * @param knots - The knot structure
     * @param degrees - Array of polynomial degrees for each parameter
     */
    constructor(
        private readonly vectorSpace: VectorSpace<K, V>,
        private readonly controlNet: ControlNet<V>,
        private readonly knots: KnotStructure,
        private readonly degrees: ReadonlyArray<number>
    ) {
        this.validateConstructorParams();
    }

    public getKnots() {
        return this.knots
    }

    public getDegrees() {
        return this.degrees
    }

    public getControlNet() {
        return this.controlNet
    }



    /**
     * Evaluates the B-spline at given parameters using de Boor's algorithm
     * 
     * @param parameters - Parameter values for evaluation
     * @returns Evaluated point on the B-spline
     */
    evaluate(parameters: number[]): V {
        this.validateParameters(parameters);
        return this.evaluateRecursive(parameters, 0);
    }

    /**
     * Recursive implementation of de Boor's algorithm
     * Handles arbitrary dimensional B-splines through dimension reduction
     * 
     * @param parameters - Array of parameter values for each dimension
     * @param depth - Current recursion depth (dimension being processed)
     * @returns Evaluated point on the B-spline
     */
    private evaluateRecursive(parameters: number[], depth: number): V {
        // Base case: univariate evaluation
        if (depth === this.controlNet.getDimension() - 1) {
            return this.evaluateUnivariateCase(
                parameters[depth],
                this.getCurrentControlPoints(parameters, depth),
                this.degrees[depth],
                this.knots.getKnotSequence(depth)
            );
        }
    
        // Get the number of control points needed for the current dimension
        const degree = this.degrees[depth];
        const span = this.findSpanIndex(
            parameters[depth],
            degree,
            this.knots.getKnotSequence(depth)
        );
    
        // Create intermediate control points through recursive evaluation
        const intermediatePoints: V[] = new Array(degree + 1);
        
        // For each point in the current span
        for (let i = 0; i <= degree; i++) {
            // Create parameter set for the next dimension
            const nextParameters = this.createParameterSetForNextDimension(
                parameters,
                depth,
                span - degree + i
            );
    
            // Recursively evaluate in the next dimension
            intermediatePoints[i] = this.evaluateRecursive(nextParameters, depth + 1);
        }
    
        // Compute basis functions for current parameter
        const basisFunctions = computeBasisFunctions(
            parameters[depth],
            span,
            degree,
            this.knots.getKnotSequence(depth)
        );
    
        // Compute final point as linear combination of intermediate points
        let result = this.vectorSpace.zero();
        for (let i = 0; i <= degree; i++) {
            const scaledPoint = this.vectorSpace.scale(
                basisFunctions[i] as K,
                intermediatePoints[i]
            );
            result = this.vectorSpace.add(result, scaledPoint);
        }
    
        return result;
    }
    
    /**
     * Helper method to get current control points for a specific parameter set
     * 
     * @param parameters - Full parameter set
     * @param depth - Current dimension
     * @returns Array of control points for the current evaluation
     */
    private getCurrentControlPoints(parameters: number[], depth: number): V[] {
        const indices = new Array(this.controlNet.getDimension()).fill(0);
        const result: V[] = [];
        
        // Get the size of the control net in the current dimension
        const size = this.controlNet.getSize(depth);
        
        // Fill in known parameters
        for (let i = 0; i < depth; i++) {
            const span = this.findSpanIndex(
                parameters[i],
                this.degrees[i],
                this.knots.getKnotSequence(i)
            );
            indices[i] = span;
        }
    
        // Collect control points
        for (let i = 0; i < size; i++) {
            indices[depth] = i;
            result.push(this.controlNet.getPoint(indices));
        }
    
        return result;
    }
    
    /**
     * Creates parameter set for the next dimension evaluation
     * 
     * @param parameters - Current parameter set
     * @param depth - Current dimension
     * @param index - Index in current dimension
     * @returns New parameter set for next dimension
     */
    private createParameterSetForNextDimension(
        parameters: number[],
        depth: number,
        index: number
    ): number[] {
        const nextParameters = [...parameters];
        nextParameters[depth] = index;
        return nextParameters;
    }
    
    
    /**
     * Validates that the parameters array matches the dimension of the B-spline
     * 
     * @param parameters - Array of parameter values
     * @throws Error if parameters array is invalid
     */
    private validateParameters(parameters: number[]): void {
        if (parameters.length !== this.controlNet.getDimension()) {
            throw new Error(
                `Invalid number of parameters. Expected ${
                    this.controlNet.getDimension()
                }, got ${parameters.length}`
            );
        }
    
        // Check each parameter is within its knot vector range
        parameters.forEach((param, i) => {
            const knots = this.knots.getKnotSequence(i);
            if (param < knots[0] || param > knots[knots.length - 1]) {
                throw new Error(
                    `Parameter ${param} at dimension ${i} is outside valid range [${
                        knots[0]
                    }, ${knots[knots.length - 1]}]`
                );
            }
        });
    }
    

    /**
     * Implements de Boor's algorithm for the univariate case
     * This is the base case for the recursive evaluation
     */
    private evaluateUnivariateCase(
        t: number,
        controlPoints: V[],
        degree: number,
        knots: ReadonlyArray<number>
    ): V {
        const span = this.findSpanIndex(t, degree, knots);

        // Basis functions are always real-valued
        const realBasisFunctions = computeBasisFunctions(t, span, degree, knots);
        
        // Convert real basis functions to scalar field K if needed
        const scalarBasisFunctions = this.convertToScalarField(realBasisFunctions);
        
        let point = this.vectorSpace.zero();
        for (let i = 0; i <= degree; i++) {
            const scaledPoint = this.vectorSpace.scale(
                scalarBasisFunctions[i],
                controlPoints[span - degree + i]
            );
            point = this.vectorSpace.add(point, scaledPoint);
        }
        
        return point;
    }

    /**
     * Converts real-valued basis functions to scalar field K if needed
     * 
     * @param realBasisFunctions - Array of real-valued basis functions
     * @returns Array of values in scalar field K
     */
    private convertToScalarField(realBasisFunctions: number[]): K[] {
        if (this.isRealField()) {
            return realBasisFunctions as K[];
        } else {
            // For complex field, convert to complex numbers [real, 0]
            return realBasisFunctions.map(x => [x, 0] as Complex) as K[];
        }
    }

    /**
     * Determines if we're working over the real field
     */
    private isRealField(): boolean {
        return !this.isComplexNumber(this.vectorSpace.zero()[0]);
    }

    private isComplexNumber(value: unknown): value is Complex {
        return Array.isArray(value) && value.length === 2;
    }

    /**
     * Validates constructor parameters for mathematical consistency
     * Ensures the B-spline is well-defined according to spline theory
     * 
     * @throws Error if any parameter is invalid or inconsistent
     */
    private validateConstructorParams(): void {
        // 1. Validate dimensions consistency
        this.validateDimensions();

        // 2. Validate degrees
        this.validateDegrees();

        // 3. Validate knot structure
        this.validateKnotStructure();

        // 4. Validate control net
        this.validateControlNet();
    }

    /**
     * Validates dimensional consistency across all components
     */
    private validateDimensions(): void {
        const netDimension = this.controlNet.getDimension();
        const knotDimension = this.knots.getDimension();
        const degreeDimension = this.degrees.length;

        if (netDimension !== knotDimension || netDimension !== degreeDimension) {
            throw new Error(
                `Dimensional mismatch: Control net (${netDimension}D), ` +
                `Knot structure (${knotDimension}D), ` +
                `Degrees (${degreeDimension}D)`
            );
        }
    }

    /**
     * Validates degrees against control points and knot vectors
     */
    private validateDegrees(): void {
        this.degrees.forEach((degree, i) => {
            // Degree must be non-negative
            if (degree < 0) {
                throw new Error(
                    `Degree must be non-negative in dimension ${i}: ${degree}`
                );
            }

            // Get number of control points in this dimension
            const numControlPoints = this.controlNet.getSize(i);

            // Get number of knots in this dimension
            const numKnots = this.knots.getKnotSequence(i).length;

            // Validate degree against number of control points
            if (numControlPoints <= degree) {
                throw new Error(
                    `Dimension ${i}: Number of control points (${numControlPoints}) ` +
                    `must be greater than degree (${degree})`
                );
            }

            // Validate degree against knot vector length
            // For a B-spline of degree p with n+1 control points,
            // we need m+1 knots where m = n + p + 1
            const requiredKnots = numControlPoints + degree + 1;
            if (numKnots !== requiredKnots) {
                throw new Error(
                    `Dimension ${i}: Invalid number of knots. ` +
                    `Expected ${requiredKnots}, got ${numKnots}`
                );
            }
        });
    }

    /**
     * Validates knot vector properties
     */
    private validateKnotStructure(): void {
        for (let dim = 0; dim < this.knots.getDimension(); dim++) {
            const knots = this.knots.getKnotSequence(dim);

            // Check knot vector is non-decreasing
            for (let i = 1; i < knots.length; i++) {
                if (knots[i] < knots[i - 1]) {
                    throw new Error(
                        `Dimension ${dim}: Knot sequence must be non-decreasing. ` +
                        `Invalid at index ${i}: ${knots[i - 1]} > ${knots[i]}`
                    );
                }
            }

            // Check for correct multiplicity at ends
            const degree = this.degrees[dim];
            const startMultiplicity = this.countMultiplicity(knots, knots[0]);
            const endMultiplicity = this.countMultiplicity(
                knots,
                knots[knots.length - 1]
            );

            if (startMultiplicity < degree + 1) {
                throw new Error(
                    `Dimension ${dim}: Start knot multiplicity (${startMultiplicity}) ` +
                    `must be at least degree + 1 (${degree + 1})`
                );
            }

            if (endMultiplicity < degree + 1) {
                throw new Error(
                    `Dimension ${dim}: End knot multiplicity (${endMultiplicity}) ` +
                    `must be at least degree + 1 (${degree + 1})`
                );
            }
        }
    }

    /**
     * Validates control net structure
     */
    private validateControlNet(): void {
        // Verify control points are valid vectors in the vector space
        try {
            const testPoint = this.controlNet.getPoint(
                Array(this.controlNet.getDimension()).fill(0)
            );
            
            // Verify we can perform vector space operations
            this.vectorSpace.add(testPoint, this.vectorSpace.zero());
            this.vectorSpace.scale(1 as K, testPoint);
        } catch (error) {
            throw new Error(
                `Invalid control points: Control points must be valid vectors in ` +
                `the specified vector space. ${error}`
            );
        }

        // Verify control net size in each dimension
        for (let dim = 0; dim < this.controlNet.getDimension(); dim++) {
            const size = this.controlNet.getSize(dim);
            if (size < 2) {
                throw new Error(
                    `Dimension ${dim}: Control net must have at least 2 points ` +
                    `in each dimension. Got ${size}`
                );
            }
        }
    }

    /**
     * Counts multiplicity of a value in a knot sequence
     */
    private countMultiplicity(knots: ReadonlyArray<number>, value: number): number {
        return knots.filter(k => Math.abs(k - value) < Number.EPSILON).length;
    }


    private findSpanIndex(
        t: number,
        degree: number,
        knots: ReadonlyArray<number>
    ): number {
        return findSpanIndex(t, degree, knots);
    }

    private computeBasisFunctions(
        t: number,
        i: number,
        p: number,
        knots: ReadonlyArray<number>
    ): number[] {
        return computeBasisFunctions(t, i, p, knots);
    }

    
    /**
     * Helper method to get all control points in current (first) dimension
     * Returns a flattened array of control points for the current parametric dimension
     * 
     * @returns Array of control points
     */
    public getAllControlPoints(): V[] {
        const indices = new Array(this.controlNet.getDimension()).fill(0);
        const points: V[] = [];
        
        const size = this.controlNet.getSize(0);
        for (let i = 0; i < size; i++) {
            indices[0] = i;
            points.push(this.controlNet.getPoint(indices));
        }
        
        return points;
    }
    



}

// ------------ Utility Functions ------------

/**
 * Computes basis functions using the Cox-de Boor recursion formula
 * 
 * @param t - Parameter value
 * @param span - Knot span index
 * @param degree - Degree of the B-spline
 * @param knots - Knot vector
 * @returns Array of basis function values N_{i,p}(t)
 */
function computeBasisFunctions(
    t: number,
    span: number,
    degree: number,
    knots: ReadonlyArray<number>
): number[] {
    // Initialize the basis functions array
    const N: number[] = new Array(degree + 1).fill(0);
    
    // Initialize the zeroth degree basis functions
    N[0] = 1.0;
    
    // Temporary arrays to store divided differences
    const left: number[] = new Array(degree + 1).fill(0);
    const right: number[] = new Array(degree + 1).fill(0);

    // Compute basis functions of increasing degree
    for (let j = 1; j <= degree; j++) {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        
        let saved = 0.0;

        // Compute basis functions of degree j
        for (let r = 0; r < j; r++) {
            const temp = N[r] / (right[r + 1] + left[j - r]);
            N[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }

        N[j] = saved;
    }

    return N;
}

/**
 * Alternative implementation using explicit recursion
 * This version is more readable but less efficient
 * 
 * @param t - Parameter value
 * @param i - Basis function index
 * @param p - Degree
 * @param knots - Knot vector
 * @returns Value of basis function N_{i,p}(t)
 */
function computeBasisFunctionRecursive(
    t: number,
    i: number,
    p: number,
    knots: ReadonlyArray<number>
): number {
    // Base case: degree 0
    if (p === 0) {
        return (t >= knots[i] && t < knots[i + 1]) ? 1.0 : 0.0;
    }

    // Compute the two terms of the recursion formula
    let term1 = 0.0;
    let term2 = 0.0;

    // First term: (t - t_i)/(t_{i+p} - t_i) * N_{i,p-1}(t)
    const denom1 = knots[i + p] - knots[i];
    if (denom1 !== 0) {
        term1 = ((t - knots[i]) / denom1) * 
                computeBasisFunctionRecursive(t, i, p - 1, knots);
    }

    // Second term: (t_{i+p+1} - t)/(t_{i+p+1} - t_{i+1}) * N_{i+1,p-1}(t)
    const denom2 = knots[i + p + 1] - knots[i + 1];
    if (denom2 !== 0) {
        term2 = ((knots[i + p + 1] - t) / denom2) * 
                computeBasisFunctionRecursive(t, i + 1, p - 1, knots);
    }

    return term1 + term2;
}

/**
 * Computes all basis functions and their derivatives up to a given order
 * 
 * @param t - Parameter value
 * @param span - Knot span index
 * @param degree - Degree of the B-spline
 * @param knots - Knot vector
 * @param derivOrder - Maximum derivative order to compute
 * @returns Array of arrays containing basis functions and their derivatives
 */
function computeBasisFunctionsWithDerivatives(
    t: number,
    span: number,
    degree: number,
    knots: ReadonlyArray<number>,
    derivOrder: number
): number[][] {
    const ndu: number[][] = Array(degree + 1)
        .fill(null)
        .map(() => Array(degree + 1).fill(0));
    
    const left: number[] = new Array(degree + 1).fill(0);
    const right: number[] = new Array(degree + 1).fill(0);
    
    // Array of derivatives
    const derivatives: number[][] = Array(derivOrder + 1)
        .fill(null)
        .map(() => Array(degree + 1).fill(0));

    // Initialize the zeroth degree basis function
    ndu[0][0] = 1.0;

    // Compute triangular table of basis functions
    for (let j = 1; j <= degree; j++) {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        
        let saved = 0.0;

        for (let r = 0; r < j; r++) {
            // Lower triangle
            ndu[j][r] = right[r + 1] + left[j - r];
            const temp = ndu[r][j - 1] / ndu[j][r];
            
            // Upper triangle
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }

        ndu[j][j] = saved;
    }

    // Load the basis functions
    for (let j = 0; j <= degree; j++) {
        derivatives[0][j] = ndu[j][degree];
    }

    // Compute the derivatives
    for (let r = 0; r <= degree; r++) {
        let s1 = 0;
        let s2 = 1;
        
        const a: number[][] = Array(2)
            .fill(null)
            .map(() => Array(degree + 1).fill(0));

        // Loop to compute kth derivative
        for (let k = 1; k <= derivOrder; k++) {
            let d = 0.0;
            const rk = r - k;
            const pk = degree - k;

            if (r >= k) {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }

            const j1 = rk >= -1 ? 1 : -rk;
            const j2 = (r - 1 <= pk) ? k - 1 : degree - r;

            for (let j = j1; j <= j2; j++) {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][j] * ndu[rk + j][pk];
            }

            if (r <= pk) {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }

            derivatives[k][r] = d;

            // Swap rows
            [s1, s2] = [s2, s1];
        }
    }

    // Multiply through by the correct factors
    let r = degree;
    for (let k = 1; k <= derivOrder; k++) {
        for (let j = 0; j <= degree; j++) {
            derivatives[k][j] *= r;
        }
        r *= (degree - k);
    }

    return derivatives;
}


/**
 * Finds the knot span index using binary search
 */
function findSpanIndex(
    t: number,
    degree: number,
    knots: ReadonlyArray<number>
): number {
    const n = knots.length - degree - 2 // number of control points
    if (t >= knots[n + 1]) return n;
    if (t <= knots[degree]) return degree;
    let low = degree;
    let high = n + 1;
    let mid = Math.floor((low + high) / 2);
    while (t < knots[mid] || t >= knots[mid + 1]) {
        if (t < knots[mid]) {
            high = mid;
        } else {
            low = mid;
        }
        mid = Math.floor((low + high) / 2);
    }
    return mid;
}


