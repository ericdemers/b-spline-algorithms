// Core types and interfaces
type Dimension = number;
type ParameterValue = number[];
type Index = number[];

interface Point {
    coordinates: number[];
    weight: number;
}

interface ParametricDomain {
    dimension: number;
    evaluate(parameters: ParameterValue): boolean;
    project(point: Point): ParameterValue;
}

interface KnotVector {
    knots: number[];
    degree: number;
    domain: [number, number];
}

// Base classes and interfaces
abstract class BasisFunction {
    abstract evaluate(parameters: ParameterValue): number;
    abstract support: ParametricDomain;
}

abstract class ParametricSpace {
    abstract dimension: number;
    abstract domain: ParametricDomain;
    
    abstract evaluateBasis(parameters: ParameterValue, index: Index): number;
    abstract getActiveBasis(parameters: ParameterValue): Map<string, BasisFunction>;
}

abstract class SplineBase {
    protected space: ParametricSpace;
    protected controlPoints: Map<string, Point>;

    constructor(space: ParametricSpace) {
        this.space = space;
        this.controlPoints = new Map();
    }

    abstract evaluate(parameters: ParameterValue): Point;
}

// B-spline specific implementations
class BSplineBasis extends BasisFunction {
    private knotVector: KnotVector;
    private index: number;

    constructor(knotVector: KnotVector, index: number) {
        super();
        this.knotVector = knotVector;
        this.index = index;
    }

    evaluate(parameters: ParameterValue): number {
        return this.deBoorCox(parameters[0], this.index, this.knotVector.degree);
    }

    get support(): ParametricDomain {
        const start = this.knotVector.knots[this.index];
        const end = this.knotVector.knots[this.index + this.knotVector.degree + 1];
        return {
            dimension: 1,
            evaluate: (params: ParameterValue) => params[0] >= start && params[0] < end,
            project: (point: Point) => [Math.max(start, Math.min(end, point.coordinates[0]))]
        };
    }

    private deBoorCox(t: number, i: number, p: number): number {
        if (p === 0) {
            return (t >= this.knotVector.knots[i] && t < this.knotVector.knots[i+1]) ? 1 : 0;
        }
        
        let left = 0, right = 0;
        
        if (this.knotVector.knots[i+p] !== this.knotVector.knots[i]) {
            left = (t - this.knotVector.knots[i]) / (this.knotVector.knots[i+p] - this.knotVector.knots[i]);
        }
        
        if (this.knotVector.knots[i+p+1] !== this.knotVector.knots[i+1]) {
            right = (this.knotVector.knots[i+p+1] - t) / (this.knotVector.knots[i+p+1] - this.knotVector.knots[i+1]);
        }
        
        return left * this.deBoorCox(t, i, p-1) + right * this.deBoorCox(t, i+1, p-1);
    }
}

class BSplineSpace extends ParametricSpace {
    private knotVectors: KnotVector[];

    constructor(knotVectors: KnotVector[]) {
        super();
        this.knotVectors = knotVectors;
    }

    get dimension(): number {
        return this.knotVectors.length;
    }

    get domain(): ParametricDomain {
        return {
            dimension: this.dimension,
            evaluate: (params: ParameterValue) => 
                params.every((t, i) => t >= this.knotVectors[i].domain[0] && t <= this.knotVectors[i].domain[1]),
            project: (point: Point) => 
                point.coordinates.map((coord, i) => 
                    Math.max(this.knotVectors[i].domain[0], 
                             Math.min(this.knotVectors[i].domain[1], coord)))
        };
    }

        return this.knotVectors.reduce((acc, kv, i) => 
            acc * new BSplineBasis(kv, index[i]).evaluate([parameters[i]]), 1);
    }

    getActiveBasis(parameters: ParameterValue): Map<string, BasisFunction> {
        const activeBasis = new Map<string, BasisFunction>();
        // Implementation depends on how you want to handle multi-dimensional basis functions
        // This is a simplified version
        this.knotVectors.forEach((kv, dim) => {
            for (let i = 0; i < kv.knots.length - kv.degree - 1; i++) {
                const basis = new BSplineBasis(kv, i);
                if (basis.support.evaluate([parameters[dim]])) {
                    activeBasis.set(`${dim}_${i}`, basis);
                }
            }
        });
        return activeBasis;
    }
}

class BSpline extends SplineBase {
    constructor(knotVectors: KnotVector[]) {
        super(new BSplineSpace(knotVectors));
    }

    evaluate(parameters: ParameterValue): Point {
        const activeBasis = this.space.getActiveBasis(parameters);
        let result = { coordinates: new Array(this.space.dimension).fill(0), weight: 0 };
        
        for (const [key, basis] of activeBasis) {
            const controlPoint = this.controlPoints.get(key);
            if (controlPoint) {
                const value = basis.evaluate(parameters);
                result.coordinates = result.coordinates.map((coord, i) => 
                    coord + value * controlPoint.coordinates[i] * controlPoint.weight);
                result.weight += value * controlPoint.weight;
            }
        }
        
        result.coordinates = result.coordinates.map(coord => coord / result.weight);
        return result;
    }

    setControlPoint(index: Index, point: Point): void {
        this.controlPoints.set(index.join('_'), point);
    }
}

// Usage example
function createSimpleBSplineCurve(): BSpline {
    const knotVector: KnotVector = {
        knots: [0, 0, 0, 1, 2, 3, 4, 4, 4],
        degree: 2,
        domain: [0, 4]
    };
    
    const bspline = new BSpline([knotVector]);
    
    bspline.setControlPoint([0], { coordinates: [0, 0], weight: 1 });
    bspline.setControlPoint([1], { coordinates: [1, 1], weight: 1 });
    bspline.setControlPoint([2], { coordinates: [2, -1], weight: 1 });
    bspline.setControlPoint([3], { coordinates: [3, 0], weight: 1 });
    bspline.setControlPoint([4], { coordinates: [4, 1], weight: 1 });
    
    return bspline;
}

// Example usage
const curve = createSimpleBSplineCurve();
console.log(curve.evaluate([2])); // Evaluate at parameter 2

export default createSimpleBSplineCurve
