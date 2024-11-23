class BSplineTensorProduct<K extends Scalar, V extends Vector> extends BSpline<K, V> {
    constructor(
        vectorSpace: VectorSpace<K, V>,
        controlNet: ControlGrid<V>,
        knotStructures: KnotStructure[],
        degrees: number[]
    ) {
        super(vectorSpace, controlNet, 
              new ProductKnotStructure(knotStructures), 
              degrees);
    }

    evaluate(parameters: number[]): V {
        if (parameters.length !== this.getDimension()) {
            throw new Error('Invalid number of parameters');
        }

        // Evaluate using tensor product algorithm
        return this.evaluateTensorProduct(parameters);
    }

    private evaluateTensorProduct(parameters: number[]): V {
        // Implementation of tensor product evaluation
        const dimensions = this.getDimension();
        
        // Start with the control points
        let currentPoints = this.controlNet.getControlPoints();
        
        // For each parametric direction
        for (let dim = 0; dim < dimensions; dim++) {
            const u = parameters[dim];
            const degree = this.degrees[dim];
            const knots = this.knots.getKnotSequence(dim);
            
            // Apply de Boor algorithm in this direction
            currentPoints = this.evaluateInDirection(
                currentPoints, 
                u, 
                degree, 
                knots, 
                dim
            );
        }

        return currentPoints[0];
    }

    private evaluateInDirection(
        points: V[], 
        u: number, 
        degree: number,
        knots: number[],
        direction: number
    ): V[] {
        // De Boor algorithm implementation
        // ...
    }
}

// Helper class for tensor product knot structure
class ProductKnotStructure implements KnotStructure {
    constructor(
        private readonly knotStructures: KnotStructure[]
    ) {}

    getDimension(): number {
        return this.knotStructures.length;
    }

    getKnotSequence(direction: number): ReadonlyArray<number> {
        return this.knotStructures[direction].getKnotSequence(0);
    }

    // ... other KnotStructure methods
}

// Helper class for grid of control points
class ControlGrid<V> implements ControlNet<V> {
    constructor(
        private readonly points: V[][],
        private readonly dimensions: number[]
    ) {}

    getDimension(): number {
        return this.dimensions.length;
    }

    getControlPoints(): ReadonlyArray<V> {
        // Flatten multi-dimensional array
        return this.flattenGrid(this.points);
    }

    private flattenGrid(grid: any[]): V[] {
        return grid.flat(this.dimensions.length - 1);
    }

    // ... other ControlNet methods
}

// Usage example
const surface = new BSplineTensorProduct(
    vector3DSpace,
    new ControlGrid([
        [[0,0,0], [0,1,0], [0,2,0]],
        [[1,0,0], [1,1,1], [1,2,0]],
        [[2,0,0], [2,1,0], [2,2,0]]
    ], [3, 3]),
    [
        new UnivariateKnots([0,0,0,1,2,2,2]),
        new PeriodicKnots([0, 0.5, 1], 1.0)
    ],
    [2, 2] // degrees in each direction
);

// Evaluate at parameter values
const point = surface.evaluate([0.5, 0.3]);


class RecursiveTensorProduct<K extends Scalar, V extends Vector> {
    // Main evaluation entry point
    evaluate(parameters: number[]): V {
        // Start recursion with all control points and dimension 0
        return this.evaluateRecursive(
            this.controlNet.getControlPoints(),
            parameters,
            0
        );
    }

    private evaluateRecursive(
        points: V[],
        parameters: number[],
        currentDimension: number
    ): V {
        // Base case: when we've processed all dimensions
        if (currentDimension === parameters.length) {
            // Only one point should remain after all evaluations
            return points[0];
        }

        // Get current dimension's parameters
        const u = parameters[currentDimension];
        const degree = this.degrees[currentDimension];
        const knots = this.knots.getKnotSequence(currentDimension);

        // Evaluate along current dimension
        const newPoints = this.evaluateInDirection(
            points, 
            u, 
            degree, 
            knots, 
            currentDimension
        );

        // Recurse to next dimension with reduced point set
        return this.evaluateRecursive(
            newPoints,
            parameters,
            currentDimension + 1
        );
    }

    private evaluateInDirection(
        points: V[], 
        u: number, 
        degree: number,
        knots: number[],
        direction: number
    ): V[] {
        // Find the knot span containing u
        const span = this.findSpan(u, degree, knots);
        
        // Compute basis function values
        const basis = this.computeBasisValues(span, u, degree, knots);
        
        // Get the stride for this direction
        //stride represents the distance (number of elements) you need to skip to move to the next point in a particular direction within a flattened multi-dimensional array.
        const stride = this.getDirectionalStride(direction);
        
        // Number of points to process in this direction
        const numPoints = this.getNumPointsInDirection(direction);
        
        // Result array for this direction's evaluation
        const result: V[] = [];
        
        // Process each group of points
        for (let i = 0; i < points.length / stride; i++) {
            let sum = this.vectorSpace.zero();
            
            // Apply basis functions to points in span
            for (let j = 0; j <= degree; j++) {
                const pointIndex = i * stride + (span - degree + j) * this.getStride(direction);
                const weighted = this.vectorSpace.multiply(basis[j], points[pointIndex]);
                sum = this.vectorSpace.add(sum, weighted);
            }
            
            result.push(sum);
        }
        
        return result;
    }

    // Helper method to get stride for a direction
    private getDirectionalStride(direction: number): number {
        const sizes = this.controlNet.getSizes();
        let stride = 1;
        for (let i = direction + 1; i < sizes.length; i++) {
            stride *= sizes[i];
        }
        return stride;
    }

    // Helper method to get number of points in a direction
    private getNumPointsInDirection(direction: number): number {
        return this.controlNet.getSizes()[direction];
    }
}

