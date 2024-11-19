/**
 * De Boor algorithm implementation using blossoming
 */
class DeBoorBlossomEvaluation<T> implements EvaluationAlgorithm<T> {
    constructor(private readonly vectorSpace: VectorSpace<T>) {}

    evaluate(bspline: BSplineData<T>, parameter: number): T {
        const degree = bspline.degree;
        const knots = bspline.knots;
        
        // Create array of parameter values for blossoming
        // For de Boor, we use the parameter value repeated (degree + 1) times
        const parameters = Array(degree + 1).fill(parameter);
        
        // Create blossom context for evaluation
        const context = this.createBlossomContext(knots, parameters, degree);
        
        // Apply blossom to relevant control points
        return context.coefficients.apply(
            context.influencingIndices.map(i => bspline.controlPoints[i])
        );
    }

    private createBlossomContext(
        knots: ReadonlyArray<number>,
        parameters: ReadonlyArray<number>,
        degree: number
    ): BlossomContext<T> {
        // Find the knot span containing the parameter
        const spanIndex = this.findSpanIndex(knots, parameters[0], degree);
        
        // Compute blossom coefficients using recursive approach
        const coefficients = this.computeBlossomCoefficients(
            knots,
            parameters,
            spanIndex,
            degree
        );

        return {
            spans: this.getRelevantKnotSpans(knots, spanIndex, degree),
            coefficients: new BlossomCoefficients(this.vectorSpace, coefficients),
            influencingIndices: this.getInfluencingIndices(spanIndex, degree)
        };
    }

    /**
     * Computes blossom coefficients recursively
     */
    private computeBlossomCoefficients(
        knots: ReadonlyArray<number>,
        parameters: ReadonlyArray<number>,
        spanIndex: number,
        degree: number
    ): ReadonlyArray<number> {
        if (degree === 0) {
            return [1]; // Base case: constant function
        }

        const u = parameters[0];
        const remainingParams = parameters.slice(1);
        
        // Recursive computation using blossom properties
        const prevCoeffs = this.computeBlossomCoefficients(
            knots,
            remainingParams,
            spanIndex,
            degree - 1
        );

        // Compute linear interpolation coefficients
        const coefficients: number[] = [];
        for (let i = 0; i <= degree; i++) {
            const leftIndex = spanIndex - degree + i;
            const rightIndex = leftIndex + degree;
            
            // Handle boundary cases
            if (leftIndex < 0 || rightIndex >= knots.length) {
                coefficients.push(0);
                continue;
            }

            // Compute linear interpolation coefficient
            const denom = knots[rightIndex] - knots[leftIndex];
            const alpha = denom === 0 ? 0 : 
                (knots[rightIndex] - u) / denom;

            // Apply coefficient using blossom properties
            if (i > 0) {
                coefficients[i] = (1 - alpha) * prevCoeffs[i - 1];
            }
            if (i < degree) {
                coefficients[i] = (coefficients[i] || 0) + alpha * prevCoeffs[i];
            }
        }

        return coefficients;
    }

    /**
     * Implementation of BlossomCoefficients for de Boor algorithm
     */
    private class BlossomCoefficients implements BlossomCoefficients<T> {
        constructor(
            readonly vectorSpace: VectorSpace<T>,
            readonly values: ReadonlyArray<number>
        ) {}

        apply(points: ReadonlyArray<T>): T {
            return points.reduce(
                (result, point, i) => this.vectorSpace.add(
                    result,
                    this.vectorSpace.multiply(this.values[i], point)
                ),
                this.vectorSpace.zero()
            );
        }

        combine<U>(other: BlossomCoefficients<U>): BlossomCoefficients<[T, U]> {
            return new ProductBlossomCoefficients(this, other);
        }
    }

    /**
     * Gets indices of control points that influence the result
     */
    private getInfluencingIndices(spanIndex: number, degree: number): number[] {
        return Array.from(
            { length: degree + 1 },
            (_, i) => spanIndex - degree + i
        );
    }

    /**
     * Gets relevant knot spans for the evaluation
     */
    private getRelevantKnotSpans(
        knots: ReadonlyArray<number>,
        spanIndex: number,
        degree: number
    ): ReadonlyArray<KnotSpan> {
        const spans: KnotSpan[] = [];
        for (let i = spanIndex - degree; i <= spanIndex; i++) {
            if (i >= 0 && i < knots.length - 1) {
                spans.push({
                    start: knots[i],
                    end: knots[i + 1],
                    startIndex: i,
                    multiplicity: {
                        start: this.getMultiplicity(knots, knots[i]),
                        end: this.getMultiplicity(knots, knots[i + 1])
                    }
                });
            }
        }
        return spans;
    }

    /**
     * Finds the knot span containing the parameter
     */
    private findSpanIndex(
        knots: ReadonlyArray<number>,
        parameter: number,
        degree: number
    ): number {
        // Binary search implementation
        let low = degree;
        let high = knots.length - degree - 2;

        while (low <= high) {
            const mid = Math.floor((low + high) / 2);
            if (parameter >= knots[mid] && parameter < knots[mid + 1]) {
                return mid;
            }
            if (parameter < knots[mid]) {
                high = mid - 1;
            } else {
                low = mid + 1;
            }
        }

        return knots.length - degree - 2;
    }

    /**
     * Gets multiplicity of a knot value
     */
    private getMultiplicity(
        knots: ReadonlyArray<number>,
        value: number
    ): number {
        return knots.filter(k => Math.abs(k - value) < Number.EPSILON).length;
    }
}

// Usage example
const evaluator = new DeBoorBlossomEvaluation(vector2DSpace);

// Evaluate a curve
const curve = {
    controlPoints: [[0,0], [1,1], [2,0], [3,1]] as [number, number][],
    knots: [0, 0, 0, 0, 1, 1, 1, 1],
    degree: 3
};

const point = evaluator.evaluate(curve, 0.5);
