import { ControlNet } from "./controlNet";
import { KnotStructure } from "./knotStructure";
import { Scalar, Vector, VectorSpace } from "./vector-space";

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
        protected readonly vectorSpace: VectorSpace<K, V>,
        protected readonly controlNet: ControlNet<V>,
        protected readonly knots: KnotStructure,
        protected readonly degrees: ReadonlyArray<number>
    ) {
        this.validateConstructorParams();
    }

    // Getters
    public getKnots(): KnotStructure {
        return this.knots;
    }

    public getDegrees(): ReadonlyArray<number> {
        return this.degrees;
    }

    public getControlNet(): ControlNet<V> {
        return this.controlNet;
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
    protected evaluateRecursive(parameters: number[], depth: number): V {
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
            const nextParameters = this.createParameterSetForNextDimension(
                parameters,
                depth,
                span - degree + i
            );
            intermediatePoints[i] = this.evaluateRecursive(nextParameters, depth + 1);
        }
    
        // Compute basis functions for current parameter
        const basisFunctions = this.computeBasisFunctions(
            parameters[depth],
            span,
            degree,
            this.knots.getKnotSequence(depth)
        );
    
        // Compute final point as linear combination of intermediate points
        return this.computeLinearCombination(intermediatePoints, basisFunctions);
    }

    /**
     * Computes linear combination of points with given coefficients
     */
    protected computeLinearCombination(points: V[], coefficients: number[]): V {
        let result = this.vectorSpace.zero();
        for (let i = 0; i < points.length; i++) {
            const scaledPoint = this.vectorSpace.scale(
                coefficients[i] as K,
                points[i]
            );
            result = this.vectorSpace.add(result, scaledPoint);
        }
        return result;
    }

    /**
     * Implements de Boor's algorithm for the univariate case
     */
    protected evaluateUnivariateCase(
        t: number,
        controlPoints: V[],
        degree: number,
        knots: ReadonlyArray<number>
    ): V {
        const span = this.findSpanIndex(t, degree, knots);
        const basisFunctions = this.computeBasisFunctions(t, span, degree, knots);
        return this.computeLinearCombination(
            controlPoints.slice(span - degree, span + 1),
            basisFunctions
        );
    }

    /**
     * Computes B-spline basis functions
     */
    protected computeBasisFunctions(
        u: number,
        span: number,
        degree: number,
        knots: ReadonlyArray<number>
    ): number[] {
        const basis: number[] = new Array(degree + 1).fill(0);
        const left: number[] = new Array(degree + 1).fill(0);
        const right: number[] = new Array(degree + 1).fill(0);

        // Initialize zeroth-degree basis functions
        basis[0] = 1.0;

        // Compute basis functions
        for (let j = 1; j <= degree; j++) {
            left[j] = u - knots[span + 1 - j];
            right[j] = knots[span + j] - u;
            let saved = 0.0;

            for (let r = 0; r < j; r++) {
                const temp = basis[r] / (right[r + 1] + left[j - r]);
                basis[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }

            basis[j] = saved;
        }

        return basis;
    }

    /**
     * Finds the knot span index for a given parameter
     */
    protected findSpanIndex(
        u: number,
        degree: number,
        knots: ReadonlyArray<number>
    ): number {
        const n = knots.length - degree - 1;

        // Handle boundary cases
        if (u >= knots[n]) return n - 1;
        if (u <= knots[degree]) return degree;

        // Binary search
        let low = degree;
        let high = n;
        let mid = Math.floor((low + high) / 2);

        while (u < knots[mid] || u >= knots[mid + 1]) {
            if (u < knots[mid]) {
                high = mid;
            } else {
                low = mid;
            }
            mid = Math.floor((low + high) / 2);
        }

        return mid;
    }

    /**
     * Helper method to get current control points for a specific parameter set
     */
    protected getCurrentControlPoints(parameters: number[], depth: number): V[] {
        const indices = new Array(this.controlNet.getDimension()).fill(0);
        const result: V[] = [];
        
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
        const size = this.controlNet.getSize(depth);
        for (let i = 0; i < size; i++) {
            indices[depth] = i;
            result.push(this.controlNet.getPoint(indices));
        }
    
        return result;
    }

    // ... (continued in next part)

    /**
     * Creates parameter set for next dimension evaluation
     * 
     * @param currentParams - Current parameter values
     * @param currentDepth - Current dimension being processed
     * @param nextIndex - Index for the next dimension
     * @returns New parameter set for next dimension evaluation
     */
    protected createParameterSetForNextDimension(
        currentParams: number[],
        currentDepth: number,
        nextIndex: number
    ): number[] {
        const newParams = [...currentParams];
        
        // Map the index to a parameter value in the next dimension
        const knotSequence = this.knots.getKnotSequence(currentDepth + 1);
        const paramValue = knotSequence[nextIndex];
        newParams[currentDepth + 1] = paramValue;
        
        return newParams;
    }
     
 

    /**
     * Inserts a knot into the B-spline
     * 
     * @param dimension - The dimension in which to insert the knot
     * @param u - Parameter value where to insert the knot
     * @returns New B-spline with inserted knot
     */
    insertKnot(dimension: number, u: number): BSpline<K, V> {
        this.validateDimension(dimension);
        
        const degree = this.degrees[dimension];
        const knots = this.knots.getKnotSequence(dimension);
        const span = this.findSpanIndex(u, degree, knots);
        
        // Create new knot vector
        const newKnots = [...knots];
        newKnots.splice(span + 1, 0, u);
        
        // Calculate new control points
        const newControlNet = this.computeNewControlPoints(dimension, span, u);
        
        // Create new knot structure
        const newKnotStructure = this.knots.withInsertedKnot(dimension, u);
        
        return new BSpline(
            this.vectorSpace,
            newControlNet,
            newKnotStructure,
            this.degrees
        );
    }

    /**
     * Computes derivative of the B-spline
     * 
     * @param dimension - The dimension along which to take the derivative
     * @returns New B-spline representing the derivative
     */
    getDerivative(dimension: number): BSpline<K, V> {
        this.validateDimension(dimension);
        
        const degree = this.degrees[dimension];
        if (degree < 1) {
            throw new Error('Cannot take derivative of degree 0 B-spline');
        }
        
        // Compute new degrees
        const newDegrees = [...this.degrees];
        newDegrees[dimension]--;
        
        // Compute new control points
        const newControlNet = this.computeDerivativeControlPoints(dimension);
        
        // Create new knot structure
        const newKnotStructure = this.knots.withRemovedKnots(dimension);
        
        return new BSpline(
            this.vectorSpace,
            newControlNet,
            newKnotStructure,
            newDegrees
        );
    }

    /**
     * Computes new control points for knot insertion
     */
    protected computeNewControlPoints(
        dimension: number,
        span: number,
        u: number
    ): ControlNet<V> {
        // Implementation of control point computation for knot insertion
        // This would involve creating a new control net with updated points
        throw new Error('Not implemented');
    }

    /**
     * Computes control points for derivative
     */
    protected computeDerivativeControlPoints(
        dimension: number
    ): ControlNet<V> {
        // Implementation of control point computation for derivative
        // This would involve creating a new control net with derivative points
        throw new Error('Not implemented');
    }

    /**
     * Validates constructor parameters
     * Checks degrees, control net, and knot structure consistency
     * @throws Error if parameters are invalid
     */
    protected validateConstructorParams(): void {
        // Check degrees array
        if (!this.degrees || this.degrees.length === 0) {
            throw new Error('Degrees array must not be empty');
        }
        
        if (this.degrees.some(d => d < 0)) {
            throw new Error('Degrees must be non-negative');
        }

        // Check control net dimensions
        if (!this.controlNet || this.controlNet.getDimension() !== this.degrees.length) {
            throw new Error(
                `Control net dimension mismatch: ` +
                `expected ${this.degrees.length}, ` +
                `got ${this.controlNet?.getDimension()}`
            );
        }

        // Check knot structure dimensions
        if (!this.knots || this.knots.getDimension() !== this.degrees.length) {
            throw new Error(
                `Knot structure dimension mismatch: ` +
                `expected ${this.degrees.length}, ` +
                `got ${this.knots?.getDimension()}`
            );
        }

        // Validate each dimension's knot sequence
        for (let i = 0; i < this.degrees.length; i++) {
            const degree = this.degrees[i];
            const knotSequence = this.knots.getKnotSequence(i);
            const controlPointCount = this.controlNet.getSize(i);

            // Check knot vector length
            if (knotSequence.length !== controlPointCount + degree + 1) {
                throw new Error(
                    `Invalid knot sequence length in dimension ${i}: ` +
                    `expected ${controlPointCount + degree + 1}, ` +
                    `got ${knotSequence.length}`
                );
            }

            // Check knot sequence is non-decreasing
            for (let j = 1; j < knotSequence.length; j++) {
                if (knotSequence[j] < knotSequence[j - 1]) {
                    throw new Error(
                        `Knot sequence must be non-decreasing in dimension ${i}`
                    );
                }
            }
        }
    }
    


    /**
     * Validates evaluation parameters
     */
    protected validateParameters(parameters: number[]): void {
        if (!parameters || parameters.length !== this.degrees.length) {
            throw new Error('Number of parameters must match B-spline dimension');
        }
        
        parameters.forEach((param, i) => {
            const knots = this.knots.getKnotSequence(i);
            if (param < knots[0] || param > knots[knots.length - 1]) {
                throw new Error(`Parameter ${param} out of bounds for dimension ${i}`);
            }
        });
    }

    /**
     * Validates dimension index
     */
    protected validateDimension(dimension: number): void {
        if (dimension < 0 || dimension >= this.degrees.length) {
            throw new Error('Invalid dimension index');
        }
    }
}
    
