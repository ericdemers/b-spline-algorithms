// Type for 3D point/vector
interface Point3D {
    x: number;
    y: number;
    z: number;
}

// Represents a triangular Bézier patch
class TriangularBezierPatch {
    private controlPoints: Point3D[][];
    private degree: number;

    constructor(controlPoints: Point3D[][]) {
        this.controlPoints = controlPoints;
        this.degree = controlPoints.length - 1;
    }

    // Evaluate the patch at barycentric coordinates (u, v, w)
    // where u + v + w = 1
    evaluate(u: number, v: number, w: number): Point3D {
        // Validate barycentric coordinates
        if (Math.abs(u + v + w - 1.0) > 1e-10) {
            throw new Error("Invalid barycentric coordinates: must sum to 1");
        }

        // Initialize result point
        const result: Point3D = { x: 0, y: 0, z: 0 };

        // Apply de Casteljau's algorithm for triangular patches
        for (let i = 0; i <= this.degree; i++) {
            for (let j = 0; j <= this.degree - i; j++) {
                const k = this.degree - i - j;
                
                // Calculate multinomial coefficient
                const coef = this.multinomial(this.degree, i, j, k) * 
                           Math.pow(u, i) * Math.pow(v, j) * Math.pow(w, k);
                
                // Add weighted contribution of control point
                result.x += coef * this.controlPoints[i][j].x;
                result.y += coef * this.controlPoints[i][j].y;
                result.z += coef * this.controlPoints[i][j].z;
            }
        }

        return result;
    }

    // Calculate multinomial coefficient (n! / (i!j!k!))
    private multinomial(n: number, i: number, j: number, k: number): number {
        return this.factorial(n) / 
               (this.factorial(i) * this.factorial(j) * this.factorial(k));
    }

    // Helper function to calculate factorial
    private factorial(n: number): number {
        if (n <= 1) return 1;
        let result = 1;
        for (let i = 2; i <= n; i++) {
            result *= i;
        }
        return result;
    }
}

// Example usage:
function example() {
    // Example for a quadratic (degree 2) triangular Bézier patch
    const controlPoints: Point3D[][] = [
        [   // First row (1 point)
            { x: 0, y: 0, z: 0 }
        ],
        [   // Second row (2 points)
            { x: 1, y: 0, z: 1 },
            { x: 1, y: 1, z: 0 }
        ],
        [   // Third row (3 points)
            { x: 2, y: 0, z: 0 },
            { x: 2, y: 1, z: 1 },
            { x: 2, y: 2, z: 0 }
        ]
    ];

    const patch = new TriangularBezierPatch(controlPoints);
    
    // Evaluate at a point using barycentric coordinates
    const point = patch.evaluate(1/3, 1/3, 1/3);
    console.log('Evaluated point:', point);
}

interface Point3D {
    x: number;
    y: number;
    z: number;
}

class TriangularBezierPatch {
    private controlPoints: Point3D[][];

    constructor(controlPoints: Point3D[][]) {
        this.controlPoints = controlPoints;
    }

    // Main evaluation function
    evaluate(u: number, v: number, w: number): Point3D {
        // Validate barycentric coordinates
        if (Math.abs(u + v + w - 1.0) > 1e-10) {
            throw new Error("Invalid barycentric coordinates: must sum to 1");
        }

        return this.deCasteljau(
            this.controlPoints,
            u,
            v,
            w
        );
    }

    // Recursive de Casteljau's algorithm
    private deCasteljau(
        points: Point3D[][],
        u: number,
        v: number,
        w: number
    ): Point3D {
        // Base case: if we have only one point, return it
        if (points.length === 1) {
            return points[0][0];
        }

        // Create a new array for the next level of control points
        const newPoints: Point3D[][] = [];
        
        // For each row except the last
        for (let i = 0; i < points.length - 1; i++) {
            const newRow: Point3D[] = [];
            
            // For each point in the row
            for (let j = 0; j < points[i].length - 1; j++) {
                // Compute new control point using barycentric combination
                const p1 = points[i][j];
                const p2 = points[i][j + 1];
                const p3 = points[i + 1][j];

                const newPoint = this.interpolate(p1, p2, p3, u, v, w);
                newRow.push(newPoint);
            }
            newPoints.push(newRow);
        }

        // Recursive call with reduced control net
        return this.deCasteljau(newPoints, u, v, w);
    }

    // Interpolate between three points using barycentric coordinates
    private interpolate(
        p1: Point3D,
        p2: Point3D,
        p3: Point3D,
        u: number,
        v: number,
        w: number
    ): Point3D {
        return {
            x: u * p1.x + v * p2.x + w * p3.x,
            y: u * p1.y + v * p2.y + w * p3.y,
            z: u * p1.z + v * p2.z + w * p3.z
        };
    }
}

// Example usage
function example() {
    // Example for a quadratic (degree 2) triangular Bézier patch
    const controlPoints: Point3D[][] = [
        [   // First row (3 points)
            { x: 0, y: 0, z: 0 },
            { x: 1, y: 0, z: 1 },
            { x: 2, y: 0, z: 0 }
        ],
        [   // Second row (2 points)
            { x: 0, y: 1, z: 1 },
            { x: 1, y: 1, z: 2 }
        ],
        [   // Third row (1 point)
            { x: 0, y: 2, z: 0 }
        ]
    ];

    const patch = new TriangularBezierPatch(controlPoints);
    
    // Evaluate at the center of the patch
    const point = patch.evaluate(1/3, 1/3, 1/3);
    console.log('Evaluated point:', point);
}

interface Point3D {
    x: number;
    y: number;
    z: number;
}

class TriangularBPatch {
    private controlPoints: Point3D[][];
    private degree: number;
    private knots: number[];

    constructor(controlPoints: Point3D[][], degree: number) {
        this.controlPoints = controlPoints;
        this.degree = degree;
        
        // Initialize uniform knot vector
        const numKnots = controlPoints.length + degree + 1;
        this.knots = [];
        for (let i = 0; i < numKnots; i++) {
            this.knots.push(i / (numKnots - 1));
        }
    }

    evaluate(u: number, v: number, w: number): Point3D {
        // Validate barycentric coordinates
        if (Math.abs(u + v + w - 1.0) > 1e-10) {
            throw new Error("Invalid barycentric coordinates: must sum to 1");
        }

        const result: Point3D = { x: 0, y: 0, z: 0 };
        const basisValues = this.computeBSplineBasis(u, v, w);

        let index = 0;
        for (let i = 0; i < this.controlPoints.length; i++) {
            for (let j = 0; j < this.controlPoints[i].length; j++) {
                const point = this.controlPoints[i][j];
                const basis = basisValues[index++];

                result.x += point.x * basis;
                result.y += point.y * basis;
                result.z += point.z * basis;
            }
        }

        return result;
    }

    private computeBSplineBasis(u: number, v: number, w: number): number[] {
        const basisValues: number[] = [];
        let totalPoints = 0;
        
        for (let i = 0; i < this.controlPoints.length; i++) {
            totalPoints += this.controlPoints[i].length;
        }

        // Initialize basis values
        for (let i = 0; i < totalPoints; i++) {
            basisValues[i] = 0;
        }

        // Compute basis functions using de Boor-Cox recursion
        this.computeTriangularBSplineBasis(u, v, w, this.degree, basisValues);

        return basisValues;
    }

    private computeTriangularBSplineBasis(
        u: number, 
        v: number, 
        w: number, 
        k: number, 
        values: number[]
    ): void {
        // Base case for degree 0
        if (k === 0) {
            let index = 0;
            for (let i = 0; i < this.controlPoints.length; i++) {
                for (let j = 0; j < this.controlPoints[i].length; j++) {
                    // Check if point is in the current triangle
                    const domain = this.isInDomain(u, v, w, i, j);
                    values[index++] = domain ? 1 : 0;
                }
            }
            return;
        }

        // Recursive case
        const lowerValues: number[] = new Array(values.length).fill(0);
        this.computeTriangularBSplineBasis(u, v, w, k - 1, lowerValues);

        // Compute higher degree basis functions
        let index = 0;
        for (let i = 0; i < this.controlPoints.length - k; i++) {
            for (let j = 0; j < this.controlPoints[i].length - k; j++) {
                const weights = this.computeWeights(u, v, w, i, j, k);
                
                values[index] = 0;
                for (let l = 0; l < weights.length; l++) {
                    values[index] += weights[l] * lowerValues[index + l];
                }
                index++;
            }
        }
    }

    private isInDomain(u: number, v: number, w: number, i: number, j: number): boolean {
        // Check if point (u,v,w) is in the domain triangle of control point (i,j)
        const domainU = this.knots[i];
        const domainV = this.knots[j];
        const domainW = 1 - domainU - domainV;

        const epsilon = 1e-10;
        return u >= domainU - epsilon && 
               v >= domainV - epsilon && 
               w >= domainW - epsilon;
    }

    private computeWeights(
        u: number, 
        v: number, 
        w: number, 
        i: number, 
        j: number, 
        k: number
    ): number[] {
        // Compute barycentric weights for the recursive formula
        const weights: number[] = [];
        
        // Get knot values
        const u0 = this.knots[i];
        const u1 = this.knots[i + k];
        const v0 = this.knots[j];
        const v1 = this.knots[j + k];
        
        // Compute denominators
        const du = u1 - u0;
        const dv = v1 - v0;
        
        // Compute weights based on barycentric coordinates
        if (du !== 0) weights.push((u - u0) / du);
        if (weights.length < 3) weights.push(w);

        return weights;
    }
}

// Example usage
function example() {
    const controlPoints: Point3D[][] = [
        [
            { x: 0, y: 0, z: 0 },
            { x: 1, y: 0, z: 1 },
            { x: 2, y: 0, z: 0 }
        ],
        [
            { x: 0, y: 1, z: 1 },
            { x: 1, y: 1, z: 2 }
        ],
        [
            { x: 0, y: 2, z: 0 }
        ]
    ];

    const degree = 2;
    const patch = new TriangularBPatch(controlPoints, degree);
    const point = patch.evaluate(1/3, 1/3, 1/3);
    console.log('Evaluated point:', point);
}



interface Point3D {
    x: number;
    y: number;
    z: number;
}

class TriangularBPatch {
    private controlPoints: Point3D[][];
    private degree: number;
    private knots: number[];

    constructor(controlPoints: Point3D[][], degree: number) {
        this.controlPoints = controlPoints;
        this.degree = degree;
        
        // Initialize knot vector with multiplicity at endpoints
        this.initializeKnots();
    }

    private initializeKnots(): void {
        // Create knot vector with multiplicity k+1 at endpoints
        this.knots = [];
        // Add degree+1 knots at start
        for (let i = 0; i <= this.degree; i++) {
            this.knots.push(0);
        }
        // Add interior knots
        const interior = this.controlPoints.length - this.degree;
        for (let i = 1; i < interior; i++) {
            this.knots.push(i / interior);
        }
        // Add degree+1 knots at end
        for (let i = 0; i <= this.degree; i++) {
            this.knots.push(1);
        }
    }

    evaluate(u: number, v: number, w: number): Point3D {
        // Validate barycentric coordinates
        if (Math.abs(u + v + w - 1.0) > 1e-10) {
            throw new Error("Invalid barycentric coordinates: must sum to 1");
        }

        return this.blossom(u, v, w);
    }

    private blossom(u: number, v: number, w: number): Point3D {
        // Create array of parameter values for blossoming
        const params: [number, number, number][] = [];
        for (let i = 0; i < this.degree; i++) {
            params.push([u, v, w]);
        }

        return this.recursiveBlossom(
            this.controlPoints,
            params,
            0,
            this.controlPoints.length - 1
        );
    }

    private recursiveBlossom(
        points: Point3D[][],
        params: [number, number, number][],
        level: number,
        depth: number
    ): Point3D {
        // Base case: if we've used all parameters
        if (level === this.degree) {
            return points[0][0];
        }

        // Get current parameter values
        const [u, v, w] = params[level];

        // Create new control net for next level
        const newPoints: Point3D[][] = [];
        for (let i = 0; i < depth; i++) {
            const newRow: Point3D[] = [];
            for (let j = 0; j < depth - i; j++) {
                // Compute affine combination using polar form
                const p1 = points[i][j];
                const p2 = points[i][j + 1];
                const p3 = points[i + 1][j];

                // Get knot intervals for barycentric coordinates
                const a1 = this.getKnotInterval(i, j, level);
                const a2 = this.getKnotInterval(i, j + 1, level);
                const a3 = this.getKnotInterval(i + 1, j, level);

                // Compute barycentric coordinates for this level
                const [alpha1, alpha2, alpha3] = this.computePolarWeights(
                    u, v, w,
                    a1, a2, a3
                );

                newRow.push(this.affineCombo(p1, p2, p3, alpha1, alpha2, alpha3));
            }
            newPoints.push(newRow);
        }

        // Recursive call with new control net
        return this.recursiveBlossom(newPoints, params, level + 1, depth - 1);
    }

    private getKnotInterval(i: number, j: number, level: number): [number, number, number] {
        // Get knot interval for control point (i,j) at given level
        const u = this.knots[i + level];
        const v = this.knots[j + level];
        const w = 1 - u - v;
        return [u, v, w];
    }

    private computePolarWeights(
        u: number, v: number, w: number,
        [u1, v1, w1]: [number, number, number],
        [u2, v2, w2]: [number, number, number],
        [u3, v3, w3]: [number, number, number]
    ): [number, number, number] {
        // Compute weights using polar form
        const det = (u2 - u1) * (w3 - w1) - (w2 - w1) * (u3 - u1);
        
        if (Math.abs(det) < 1e-10) {
            // Handle degenerate case
        }

        const alpha1 = ((u2 - u) * (w3 - w) - (w2 - w) * (u3 - u)) / det;
        const alpha2 = ((u3 - u) * (w1 - w) - (w3 - w) * (u1 - u)) / det;
        const alpha3 = 1 - alpha1 - alpha2;

        return [alpha1, alpha2, alpha3];
    }

    private affineCombo(
        p1: Point3D,
        p2: Point3D,
        p3: Point3D,
        w1: number,
        w2: number,
        w3: number
    ): Point3D {
        return {
            x: w1 * p1.x + w2 * p2.x + w3 * p3.x,
            y: w1 * p1.y + w2 * p2.y + w3 * p3.y,
            z: w1 * p1.z + w2 * p2.z + w3 * p3.z
        };
    }
}

// Example usage
function example() {
    const controlPoints: Point3D[][] = [
        [
            { x: 0, y: 0, z: 0 },
            { x: 1, y: 0, z: 1 },
            { x: 2, y: 0, z: 0 }
        ],
        [
            { x: 0, y: 1, z: 1 },
            { x: 1, y: 1, z: 2 }
        ],
        [
            { x: 0, y: 2, z: 0 }
        ]
    ];

    const degree = 2;
    const patch = new TriangularBPatch(controlPoints, degree);
    const point = patch.evaluate(1/3, 1/3, 1/3);
    console.log('Evaluated point:', point);
}

interface Point3D {
    x: number;
    y: number;
    z: number;
}

class TensorProductBSpline {
    private controlPoints: Point3D[][];
    private degreeU: number;
    private degreeV: number;
    private knotsU: number[];
    private knotsV: number[];

    constructor(
        controlPoints: Point3D[][],
        degreeU: number,
        degreeV: number,
        knotsU?: number[],
        knotsV?: number[]
    ) {
        this.controlPoints = controlPoints;
        this.degreeU = degreeU;
        this.degreeV = degreeV;

        // Initialize uniform knot vectors if not provided
        this.knotsU = knotsU || this.initializeKnots(controlPoints.length, degreeU);
        this.knotsV = knotsV || this.initializeKnots(controlPoints[0].length, degreeV);
    }

    private initializeKnots(numControlPoints: number, degree: number): number[] {
        const knots: number[] = [];
        const numKnots = numControlPoints + degree + 1;

        // Add degree+1 knots at start
        for (let i = 0; i <= degree; i++) {
            knots.push(0);
        }

        // Add interior knots
        const interior = numControlPoints - degree;
        for (let i = 1; i < interior; i++) {
            knots.push(i / interior);
        }

        // Add degree+1 knots at end
        for (let i = 0; i <= degree; i++) {
            knots.push(1);
        }

        return knots;
    }

    evaluate(u: number, v: number): Point3D {
        // Validate parameter range
        if (u < 0 || u > 1 || v < 0 || v > 1) {
            throw new Error("Parameters must be in range [0,1]");
        }

        return this.blossom(u, v);
    }

    private blossom(u: number, v: number): Point3D {
        // Create parameter sequences for blossoming
        const uParams: number[] = Array(this.degreeU).fill(u);
        const vParams: number[] = Array(this.degreeV).fill(v);

        // Start with original control points
        let points = this.controlPoints.map(row => [...row]);

        // Blossom in u direction first
        for (let i = 0; i < points.length; i++) {
            points[i] = this.blossomCurve(points[i], uParams, this.knotsU, this.degreeU);
        }

        // Convert result to single column
        const column = points.map(row => row[0]);

        // Blossom in v direction
        return this.blossomCurve(column, vParams, this.knotsV, this.degreeV)[0];
    }

    private blossomCurve(
        points: Point3D[],
        params: number[],
        knots: number[],
        degree: number
    ): Point3D[] {
        if (degree === 0) {
            return points;
        }

        const n = points.length - 1;
        const newPoints: Point3D[] = new Array(n);

        for (let i = 0; i < n; i++) {
            const t = params[degree - 1];
            const alpha = this.computeBlossomWeight(t, knots[span], knots[span + 1]);

            newPoints[i] = this.interpolate(points[i], points[i + 1], alpha);
        }

        return this.blossomCurve(newPoints, params, knots, degree - 1);
    }

    private findKnotSpan(t: number, knots: number[], degree: number): number {
        const n = knots.length - degree - 2;
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

    private computeBlossomWeight(t: number, t0: number, t1: number): number {
        if (Math.abs(t1 - t0) < 1e-10) return 0;
        return (t - t0) / (t1 - t0);
    }

    private interpolate(p0: Point3D, p1: Point3D, alpha: number): Point3D {
        return {
            x: (1 - alpha) * p0.x + alpha * p1.x,
            y: (1 - alpha) * p0.y + alpha * p1.y,
            z: (1 - alpha) * p0.z + alpha * p1.z
        };
    }
}

// Example usage
function example() {
    const controlPoints: Point3D[][] = [
        [
            { x: 0, y: 0, z: 0 },
            { x: 0, y: 1, z: 1 },
            { x: 0, y: 2, z: 0 }
        ],
        [
            { x: 1, y: 0, z: 1 },
            { x: 1, y: 1, z: 2 },
            { x: 1, y: 2, z: 1 }
        ],
        [
            { x: 2, y: 0, z: 0 },
            { x: 2, y: 1, z: 1 },
            { x: 2, y: 2, z: 0 }
        ]
    ];

    const degreeU = 2;
    const degreeV = 2;
    const surface = new TensorProductBSpline(controlPoints, degreeU, degreeV);

    // Evaluate at a point
    const point = surface.evaluate(0.5, 0.5);
    console.log('Evaluated point:', point);
}


/// T-splines are a generalization of NURBS that allow T-junctions in the control mesh, enabling local refinement. Here's a basic implementation focusing on the core concepts:

interface Point3D {
    x: number;
    y: number;
    z: number;
    w: number;  // Weight for rational surfaces
}

interface TJunction {
    u: number;
    v: number;
    direction: 'horizontal' | 'vertical';
}

class TSplineSurface {
    private controlPoints: Map<string, Point3D>;
    private tJunctions: TJunction[];
    private degreeU: number;
    private degreeV: number;
    private knotIntervals: Map<string, [number, number]>;

    constructor(
        controlPoints: Map<string, Point3D>,
        tJunctions: TJunction[],
        degreeU: number = 3,
        degreeV: number = 3
    ) {
        this.controlPoints = controlPoints;
        this.tJunctions = tJunctions;
        this.degreeU = degreeU;
        this.degreeV = degreeV;
        this.knotIntervals = new Map();
        
        this.initializeKnotIntervals();
    }

    evaluate(u: number, v: number): Point3D {
        // Find active blending functions at (u,v)
        const activeBlendingFunctions = this.findActiveBlendingFunctions(u, v);
        
        // Evaluate using blossoming
        return this.evaluatePoint(u, v, activeBlendingFunctions);
    }

    private findActiveBlendingFunctions(u: number, v: number): Map<string, number> {
        const active = new Map<string, number>();
        
        // Find all control points that influence this parameter value
        for (const [key, point] of this.controlPoints) {
            if (this.isActive(key, u, v)) {
                const blendingValue = this.computeBlendingFunction(key, u, v);
                active.set(key, blendingValue);
            }
        }

        return active;
    }

    private isActive(controlPointKey: string, u: number, v: number): boolean {
        const knotInterval = this.knotIntervals.get(controlPointKey);
        if (!knotInterval) return false;

        const [uMin, uMax] = knotInterval[0];
        const [vMin, vMax] = knotInterval[1];

        return u >= uMin && u <= uMax && v >= vMin && v <= vMax;
    }

    private computeBlendingFunction(key: string, u: number, v: number): number {
        // Get knot intervals for this control point
        const intervals = this.getLocalKnotIntervals(key);
        
        // Compute blending function using blossoming
        return this.blossom(u, v, intervals);
    }

    private blossom(
        u: number, 
        v: number, 
        intervals: { u: number[], v: number[] }
    ): number {
        // Create parameter sequences for blossoming
        const uParams = Array(this.degreeU).fill(u);
        const vParams = Array(this.degreeV).fill(v);

        // Compute univariate blossoming in both directions
        const uBlossom = this.univariateBlossom(uParams, intervals.u);
        const vBlossom = this.univariateBlossom(vParams, intervals.v);

        return uBlossom * vBlossom;
    }

    private univariateBlossom(
        params: number[],
        knots: number[]
    ): number {
        if (params.length === 0) return 1;

        const t = params[0];
        const remainingParams = params.slice(1);
        let result = 0;

        // Compute weights using polar form
        for (let i = 0; i < knots.length - 1; i++) {
            const alpha = this.computeWeight(t, knots[i], knots[i + 1]);
            result += alpha * this.univariateBlossom(remainingParams, knots.slice(1));
        }

        return result;
    }

    private evaluatePoint(
        u: number,
        v: number,
        blendingFunctions: Map<string, number>
    ): Point3D {
        let result = { x: 0, y: 0, z: 0, w: 0 };
        let weightSum = 0;

        // Compute weighted sum of control points
        for (const [key, blend] of blendingFunctions) {
            const point = this.controlPoints.get(key)!;
            const weight = point.w * blend;

            result.x += weight * point.x;
            result.y += weight * point.y;
            result.z += weight * point.z;
            weightSum += weight;
        }

        // Rationalize the result
        if (Math.abs(weightSum) > 1e-10) {
            result.x /= weightSum;
            result.y /= weightSum;
            result.z /= weightSum;
            result.w = 1;
        }

        return result;
    }

    private initializeKnotIntervals(): void {
        // Initialize knot intervals for each control point
        for (const [key, point] of this.controlPoints) {
            const intervals = this.computeLocalKnotIntervals(key);
            this.knotIntervals.set(key, intervals);
        }
    }

    private computeLocalKnotIntervals(key: string): [number, number][] {
        // Compute local knot intervals considering T-junctions
        const intervals: [number, number][] = [];
        
        // Parse control point index from key
        const [i, j] = key.split(',').map(Number);
        
        // Get base intervals
        const uInterval = this.getBaseInterval(i, 'u');
        const vInterval = this.getBaseInterval(j, 'v');
        
        // Modify intervals based on T-junctions
        const modifiedIntervals = this.modifyIntervalsForTJunctions(
            uInterval,
            vInterval,
            i,
            j
        );
        
        return modifiedIntervals;
    }

    private getBaseInterval(index: number, direction: 'u' | 'v'): [number, number] {
        // Compute base interval without considering T-junctions
        const degree = direction === 'u' ? this.degreeU : this.degreeV;
        const start = index / (degree + 1);
        const end = (index + 1) / (degree + 1);
        return [start, end];
    }

    private modifyIntervalsForTJunctions(
        uInterval: [number, number],
        vInterval: [number, number],
        i: number,
        j: number
    ): [number, number][] {
        // Modify intervals based on nearby T-junctions
        const relevantTJunctions = this.findRelevantTJunctions(i, j);
        
        // Apply T-junction modifications
        let modified = [uInterval, vInterval];
        for (const junction of relevantTJunctions) {
        }
        
        return modified;
    }

    private findRelevantTJunctions(i: number, j: number): TJunction[] {
        // Find T-junctions that affect this control point
        return this.tJunctions.filter(tj => 
            this.isTJunctionRelevant(tj, i, j)
        );
    }

    private isTJunctionRelevant(tj: TJunction, i: number, j: number): boolean {
        // Check if T-junction affects the given control point
        const distance = tj.direction === 'horizontal' 
            ? Math.abs(j - tj.v)
            : Math.abs(i - tj.u);
            
        return distance <= Math.max(this.degreeU, this.degreeV);
    }

    private applyTJunctionModification(
        intervals: [number, number][],
        junction: TJunction
    ): [number, number][] {
        // Modify intervals based on T-junction
        // This is a simplified version - real implementation would be more complex
        return intervals;
    }
}

// Analysis-suitable T-splines (ASTS) are crucial for ensuring the linear independence of blending functions and proper behavior in isogeometric analysis. Here's a detailed discussion with implementation:

interface TJunction {
    u: number;
    v: number;
    direction: 'horizontal' | 'vertical';
    index: number; // Index in the mesh
}

interface Face {
    vertices: number[];  // Indices of vertices
    tJunctions: TJunction[];
}

class AnalysisSuitableTSpline {
    private controlPoints: Map<string, Point3D>;
    private tJunctions: TJunction[];
    private faces: Face[];
    private degreeU: number;
    private degreeV: number;

    constructor(
        controlPoints: Map<string, Point3D>,
        tJunctions: TJunction[],
        faces: Face[],
        degreeU: number = 3,
        degreeV: number = 3
    ) {
        this.controlPoints = controlPoints;
        this.tJunctions = tJunctions;
        this.faces = faces;
        this.degreeU = degreeU;
        this.degreeV = degreeV;

        if (!this.checkAnalysisSuitability()) {
            throw new Error("T-spline mesh is not analysis-suitable");
        }
    }

    private checkAnalysisSuitability(): boolean {
        return (
            this.checkIntersectionRule() &&
            this.checkParallelRule() &&
            this.checkClosureRule() &&
            this.checkPartitionOfUnity()
        );
    }

    private checkIntersectionRule(): boolean {
        // Check if any T-junction extensions intersect
        for (let i = 0; i < this.tJunctions.length; i++) {
            for (let j = i + 1; j < this.tJunctions.length; j++) {
                if (this.doExtensionsIntersect(
                    this.tJunctions[i],
                    this.tJunctions[j]
                )) {
                    return false;
                }
            }
        }
        return true;
    }

    private doExtensionsIntersect(t1: TJunction, t2: TJunction): boolean {
        // Get extension spans
        const ext1 = this.computeExtensionSpan(t1);
        const ext2 = this.computeExtensionSpan(t2);

        // Check for intersection based on direction
        if (t1.direction === t2.direction) {
            // Parallel extensions shouldn't overlap
            return this.doIntervalsOverlap(ext1, ext2);
        } else {
            // Perpendicular extensions - check if they cross
            return this.doPerpendicularExtensionsIntersect(ext1, ext2, t1, t2);
        }
    }

    private checkParallelRule(): boolean {
        // Group T-junctions by their parameter lines
        const horizontalLines = new Map<number, TJunction[]>();
        const verticalLines = new Map<number, TJunction[]>();

        for (const tj of this.tJunctions) {
            if (tj.direction === 'horizontal') {
                const line = horizontalLines.get(tj.v) || [];
                line.push(tj);
                horizontalLines.set(tj.v, line);
            } else {
                const line = verticalLines.get(tj.u) || [];
                line.push(tj);
                verticalLines.set(tj.u, line);
            }
        }

        // Check each parameter line
        for (const [_, tjs] of horizontalLines) {
            if (!this.checkParallelTJunctions(tjs)) return false;
        }
        for (const [_, tjs] of verticalLines) {
            if (!this.checkParallelTJunctions(tjs)) return false;
        }

        return true;
    }

    private checkParallelTJunctions(tjs: TJunction[]): boolean {
        // Sort T-junctions along their parameter line
        tjs.sort((a, b) => 
            a.direction === 'horizontal' ? a.u - b.u : a.v - b.v
        );

        // Check for proper spacing
        for (let i = 1; i < tjs.length; i++) {
            const spacing = this.computeTJunctionSpacing(tjs[i-1], tjs[i]);
            if (spacing < this.degreeU + 1) return false;
        }

        return true;
    }

    private checkClosureRule(): boolean {
        // For each T-junction, verify its extension reaches valid edges
        for (const tj of this.tJunctions) {
            if (!this.checkExtensionClosure(tj)) return false;
        }
        return true;
    }

    private checkExtensionClosure(tj: TJunction): boolean {
        const extension = this.computeExtensionSpan(tj);
        
        // Find faces intersected by the extension
        const intersectedFaces = this.findIntersectedFaces(tj, extension);
        
        // Check if extension ends at valid edges
        return this.areExtensionEndpointsValid(tj, intersectedFaces);
    }

    private checkPartitionOfUnity(): boolean {
        // Sample points for checking partition of unity
        const samplePoints = this.generateSamplePoints();
        
        for (const [u, v] of samplePoints) {
            const sum = this.computeBlendingFunctionSum(u, v);
            if (Math.abs(sum - 1.0) > 1e-10) return false;
        }
        
        return true;
    }

    private computeExtensionSpan(tj: TJunction): [number, number] {
        // Compute the span of a T-junction extension
        const degree = tj.direction === 'horizontal' ? this.degreeU : this.degreeV;
        const span = degree + 1;

        if (tj.direction === 'horizontal') {
            return [tj.u - span/2, tj.u + span/2];
        } else {
            return [tj.v - span/2, tj.v + span/2];
        }
    }

    private doIntervalsOverlap(
        interval1: [number, number],
        interval2: [number, number]
    ): boolean {
        return !(
            interval1[1] < interval2[0] ||
            interval2[1] < interval1[0]
        );
    }

    private doPerpendicularExtensionsIntersect(
        ext1: [number, number],
        ext2: [number, number],
        t1: TJunction,
        t2: TJunction
    ): boolean {
        // Check if perpendicular extensions intersect
        if (t1.direction === 'horizontal') {
            return (
                ext1[0] <= t2.u && t2.u <= ext1[1] &&
                ext2[0] <= t1.v && t1.v <= ext2[1]
            );
        } else {
            return (
                ext1[0] <= t2.v && t2.v <= ext1[1] &&
                ext2[0] <= t1.u && t1.u <= ext2[1]
            );
        }
    }

    private computeTJunctionSpacing(t1: TJunction, t2: TJunction): number {
        return t1.direction === 'horizontal' 
            ? Math.abs(t2.u - t1.u)
            : Math.abs(t2.v - t1.v);
    }

    private findIntersectedFaces(
        tj: TJunction,
        extension: [number, number]
    ): Face[] {
        return this.faces.filter(face => 
            this.doesExtensionIntersectFace(tj, extension, face)
        );
    }

    private areExtensionEndpointsValid(
        tj: TJunction,
        faces: Face[]
    ): boolean {
        // Check if extension endpoints terminate at valid edges
        const endpoints = this.getExtensionEndpoints(tj);
        return endpoints.every(endpoint => 
            this.isValidEndpoint(endpoint, faces)
        );
    }

    private generateSamplePoints(): [number, number][] {
        // Generate sample points for partition of unity test
        const points: [number, number][] = [];
        const numSamples = 10;

        for (let i = 0; i <= numSamples; i++) {
            for (let j = 0; j <= numSamples; j++) {
                points.push([i/numSamples, j/numSamples]);
            }
        }

        return points;
    }

    private computeBlendingFunctionSum(u: number, v: number): number {
        let sum = 0;
        for (const [_, point] of this.controlPoints) {
            const blend = this.computeBlendingFunction(u, v, point);
            sum += blend;
        }
        return sum;
    }
}



