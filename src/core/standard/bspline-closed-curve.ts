import { Vector, VectorSpace } from "./vector-space";

/**
 * Implementation of a closed B-spline curve
 * The curve is periodic and the first and last points are connected
 */
export class ClosedBSplineCurve<K extends number, V extends Vector> {
    private readonly vectorSpace: VectorSpace<K, V>;
    private readonly controlPoints: V[];
    private readonly knots: number[];
    private readonly degree: number;
    private readonly wrappedPoints: V[]; // Internal wrapped control points

    /**
     * Constructs a new closed B-spline curve
     * 
     * @param vectorSpace - Vector space for control points
     * @param controlPoints - Array of control points
     * @param degree - Degree of the B-spline
     */
    constructor(
        vectorSpace: VectorSpace<K, V>,
        controlPoints: V[],
        degree: number
    ) {
        this.vectorSpace = vectorSpace;
        this.controlPoints = [...controlPoints];
        this.degree = degree;

        // Create wrapped control points by repeating first 'degree' points
        this.wrappedPoints = [...controlPoints];
        for (let i = 0; i < degree; i++) {
            this.wrappedPoints.push(controlPoints[i]);
        }

        // Create uniform periodic knot vector
        this.knots = this.createPeriodicKnots();

        this.validateConstructorParams();
    }

    /**
     * Evaluates the closed B-spline curve at parameter t
     * 
     * @param t - Parameter value [0,1]
     * @returns Point on the curve
     */
    public evaluate(t: number): V {
        // Map t to [0,1] by taking modulo 1
        t = t - Math.floor(t);
        
        // Map t to knot domain
        const startKnot = this.knots[this.degree];
        const endKnot = this.knots[this.knots.length - this.degree - 1];
        t = startKnot + t * (endKnot - startKnot);

        // Find knot span
        const span = this.findSpanIndex(t);

        // Compute basis functions
        const basis = this.computeBasisFunctions(t, span);

        // Compute point using wrapped control points
        let point = this.vectorSpace.zero();
        for (let i = 0; i <= this.degree; i++) {
            const idx = (span - this.degree + i) % this.controlPoints.length;
            point = this.vectorSpace.add(
                point,
                this.vectorSpace.scale(
                    basis[i] as K, 
                    this.controlPoints[idx]
                )
            );
        }
        return point;
    }

    /**
     * Gets the original control points (unwrapped)
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
     * Creates a periodic knot vector for the closed curve
     */
    private createPeriodicKnots(): number[] {
        const n = this.wrappedPoints.length - 1; // Number of wrapped control points - 1
        const total = n + this.degree + 1;
        const knots = new Array(total);

        // Create uniform knot vector
        for (let i = 0; i < total; i++) {
            knots[i] = i / (n - this.degree + 1);
        }

        return knots;
    }

    /**
     * Validates constructor parameters
     */
    private validateConstructorParams(): void {
        // Check degree
        if (this.degree < 1) {
            throw new Error('Degree must be positive');
        }

        // Check number of control points
        if (this.controlPoints.length <= this.degree) {
            throw new Error(
                'Number of control points must be greater than degree'
            );
        }
    }

    /**
     * Finds the knot span index containing parameter t
     */
    private findSpanIndex(t: number): number {
        // Special case for end of curve
        if (t >= this.knots[this.knots.length - this.degree - 1]) {
            return this.wrappedPoints.length - this.degree - 1;
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
                N[r] = saved + right[r + 1] * temp
                saved = left[j - r] * temp;
            }
            N[j] = saved;
        }

        return N;
    }

    /**
     * Gets the derivative curve
     */
    public getDerivative(): ClosedBSplineCurve<K, V> {
        if (this.degree < 2) {
            throw new Error('Degree must be at least 2 for derivative of closed curve');
        }

        // Calculate derivative control points
        const derivativePoints: V[] = [];
        const n = this.controlPoints.length;

        for (let i = 0; i < n; i++) {
            const nextIndex = (i + 1) % n;
            const coefficient = this.degree / (this.knots[i + this.degree + 1] - this.knots[i + 1]);
            const diff = this.vectorSpace.subtract(
                this.controlPoints[nextIndex],
                this.controlPoints[i]
            );
            derivativePoints.push(this.vectorSpace.scale(coefficient as K, diff));
        }

        return new ClosedBSplineCurve(
            this.vectorSpace,
            derivativePoints,
            this.degree - 1
        );
    }

    /**
     * Inserts a new knot into the closed B-spline curve
     */
    public insertKnot(t: number): ClosedBSplineCurve<K, V> {
        // Map t to [0,1]
        t = t - Math.floor(t);
        
        // Map t to knot domain
        t = this.knots[this.degree] + 
            t * (this.knots[this.knots.length - this.degree - 1] - this.knots[this.degree]);

        // Create new control points with knot insertion
        const newPoints = this.insertKnotHelper(t);

        // Create new closed curve with updated control points
        return new ClosedBSplineCurve(
            this.vectorSpace,
            newPoints,
            this.degree
        );
    }

    /**
     * Helper method for knot insertion
     */
    private insertKnotHelper(t: number): V[] {
        const span = this.findSpanIndex(t);
        const alphas = this.calculateAlphaFactors(t, span);
        const newPoints: V[] = [...this.controlPoints];

        // Update affected control points
        for (let i = span - this.degree + 1; i <= span; i++) {
            const alpha = alphas[i - (span - this.degree)];
            const index = i % this.controlPoints.length;
            const prevIndex = (i - 1) % this.controlPoints.length;
            
            newPoints[index] = this.vectorSpace.add(
                this.vectorSpace.scale(1 - alpha as K, this.controlPoints[prevIndex]),
                this.vectorSpace.scale(alpha as K, this.controlPoints[index])
            );
        }

        return newPoints;
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
}


/**
 * Implementation of a non-uniform closed B-spline curve
 * Supports custom knot vectors while maintaining closure
 */
export class NonUniformClosedBSplineCurve<K extends number, V extends Vector> {
    private readonly vectorSpace: VectorSpace<K, V>;
    private readonly controlPoints: V[];
    private readonly knots: number[];
    private readonly degree: number;
    private readonly wrappedPoints: V[];

    /**
     * Constructs a new non-uniform closed B-spline curve
     * 
     * @param vectorSpace - Vector space for control points
     * @param controlPoints - Array of control points
     * @param knots - Custom knot vector (will be adjusted for closure)
     * @param degree - Degree of the B-spline
     */
    constructor(
        vectorSpace: VectorSpace<K, V>,
        controlPoints: V[],
        knots: number[],
        degree: number
    ) {
        this.vectorSpace = vectorSpace;
        this.controlPoints = [...controlPoints];
        this.degree = degree;

        // Create wrapped control points
        this.wrappedPoints = this.createWrappedPoints();

        // Create closed knot vector
        this.knots = this.createClosedKnotVector(knots);

        this.validateConstructorParams();
    }

    /**
     * Creates wrapped control points for closure
     */
    private createWrappedPoints(): V[] {
        const wrapped = [...this.controlPoints];
        // Append first 'degree' points to end
        for (let i = 0; i < this.degree; i++) {
            wrapped.push(this.controlPoints[i]);
        }
        return wrapped;
    }

    /**
     * Creates a closed knot vector from input knots
     */
    private createClosedKnotVector(inputKnots: number[]): number[] {
        // Normalize input knots to [0,1]
        const normalizedKnots = this.normalizeKnots(inputKnots);
        
        const n = this.wrappedPoints.length - 1; // Number of wrapped control points - 1
        const m = n + this.degree + 1; // Number of knots needed
        const closedKnots = new Array(m);

        // Create periodic knot vector
        for (let i = 0; i < m; i++) {
            const baseIndex = i % (this.controlPoints.length - 1);
            const period = Math.floor(i / (this.controlPoints.length - 1));
            closedKnots[i] = normalizedKnots[baseIndex] + period;
        }

        return closedKnots;
    }

    /**
     * Normalizes knot vector to [0,1] range
     */
    private normalizeKnots(knots: number[]): number[] {
        const min = Math.min(...knots);
        const max = Math.max(...knots);
        const range = max - min;
        
        return knots.map(k => (k - min) / range);
    }

    /**
     * Evaluates the curve at parameter t
     */
    public evaluate(t: number): V {
        // Map t to [0,1] by taking modulo 1
        t = t - Math.floor(t);
        
        // Map t to knot domain
        const startKnot = this.knots[this.degree];
        const endKnot = this.knots[this.knots.length - this.degree - 1];
        t = startKnot + t * (endKnot - startKnot);

        // Find knot span
        const span = this.findSpanIndex(t);

        // Compute basis functions
        const basis = this.computeBasisFunctions(t, span);

        // Compute point
        let point = this.vectorSpace.zero();
        for (let i = 0; i <= this.degree; i++) {
            const idx = (span - this.degree + i) % this.controlPoints.length;
            point = this.vectorSpace.add(
                point,
                this.vectorSpace.scale(basis[i] as K, this.wrappedPoints[idx])
            );
        }

        return point;
    }

    /**
     * Inserts a new knot into the curve
     */
    public insertKnot(t: number): NonUniformClosedBSplineCurve<K, V> {
        // Map t to [0,1]
        t = t - Math.floor(t);
        
        // Map t to knot domain
        const startKnot = this.knots[this.degree];
        const endKnot = this.knots[this.knots.length - this.degree - 1];
        t = startKnot + t * (endKnot - startKnot);

        // Find span and compute alpha factors
        const span = this.findSpanIndex(t);
        const alphas = this.calculateAlphaFactors(t, span);

        // Create new control points
        const newPoints = this.insertKnotHelper(t, span, alphas);

        // Create new knot vector
        const newKnots = [...this.knots];
        newKnots.splice(span + 1, 0, t);

        return new NonUniformClosedBSplineCurve(
            this.vectorSpace,
            newPoints,
            newKnots,
            this.degree
        );
    }

    /**
     * Gets the derivative curve
     */
    public getDerivative(): NonUniformClosedBSplineCurve<K, V> {
        if (this.degree < 2) {
            throw new Error('Degree must be at least 2 for derivative of closed curve');
        }

        const derivativePoints: V[] = [];
        const n = this.controlPoints.length;

        for (let i = 0; i < n; i++) {
            const nextIndex = (i + 1) % n;
            const coefficient = this.degree / 
                (this.knots[i + this.degree + 1] - this.knots[i + 1]);
            const diff = this.vectorSpace.subtract(
                this.controlPoints[nextIndex],
                this.controlPoints[i]
            );
            derivativePoints.push(
                this.vectorSpace.scale(coefficient as K, diff)
            );
        }

        // Create new knot vector for derivative
        const derivativeKnots = this.knots.slice(1, -1);

        return new NonUniformClosedBSplineCurve(
            this.vectorSpace,
            derivativePoints,
            derivativeKnots,
            this.degree - 1
        );
    }

    private validateConstructorParams(): void {
        if (this.degree < 1) {
            throw new Error('Degree must be positive');
        }

        if (this.controlPoints.length <= this.degree) {
            throw new Error(
                'Number of control points must be greater than degree'
            );
        }

        // Validate knot vector
        if (!this.isValidKnotVector()) {
            throw new Error('Invalid knot vector for closed B-spline');
        }
    }

    private isValidKnotVector(): boolean {
        // Check knot vector is non-decreasing
        for (let i = 1; i < this.knots.length; i++) {
            if (this.knots[i] < this.knots[i - 1]) {
                return false;
            }
        }

        // Check periodicity conditions
        const n = this.controlPoints.length;
        for (let i = 0; i < this.degree; i++) {
            const diff1 = this.knots[i + n] - this.knots[i];
            const diff2 = this.knots[i + n + 1] - this.knots[i + 1];
            if (Math.abs(diff1 - diff2) > Number.EPSILON) {
                return false;
            }
        }

        return true;
    }

    private findSpanIndex(t: number): number {
        if (t >= this.knots[this.knots.length - this.degree - 1]) {
            return this.wrappedPoints.length - this.degree - 1;
        }

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

    private computeBasisFunctions(t: number, span: number): number[] {
        const N: number[] = new Array(this.degree + 1).fill(0);
        const left: number[] = new Array(this.degree + 1).fill(0);
        const right: number[] = new Array(this.degree + 1).fill(0);

        N[0] = 1.0;

        // Compute basis functions
        for (let j = 1; j <= this.degree; j++) {
            left[j] = t - this.knots[span + 1 - j];
            right[j] = this.knots[span + j] - t;
            let saved = 0.0;

            for (let r = 0; r < j; r++) {
                const temp = N[r] / (right[r + 1] + left[j - r]);
                N[r] = saved + right[r + 1] * temp
                saved = left[j - r] * temp;
            }
            N[j] = saved;
        }

        return N
    }

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

    private insertKnotHelper(t: number, span: number, alphas: number[]): V[] {
        const newPoints: V[] = [...this.controlPoints];

        for (let i = span - this.degree + 1; i <= span; i++) {
            const alpha = alphas[i - (span - this.degree)];
            const index = i % this.controlPoints.length;
            const prevIndex = (i - 1) % this.controlPoints.length;
            
            newPoints[index] = this.vectorSpace.add(
                this.vectorSpace.scale(1 - alpha as K, this.controlPoints[prevIndex]),
                this.vectorSpace.scale(alpha as K, this.controlPoints[index])
            );
        }

        return newPoints;
    }

    // Getters
    public getControlPoints(): V[] {
        return [...this.controlPoints];
    }

    public getKnots(): number[] {
        return [...this.knots];
    }

    public getDegree(): number {
        return this.degree;
    }
}

