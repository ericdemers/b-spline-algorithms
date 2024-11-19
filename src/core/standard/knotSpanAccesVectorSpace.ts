/**
 * Represents coefficients in a specific vector space
 */
interface BlossomCoefficients<T> {
    /**
     * The vector space these coefficients operate in
     */
    readonly vectorSpace: VectorSpace<T>;

    /**
     * Raw coefficient values
     */
    readonly values: ReadonlyArray<number>;

    /**
     * Applies coefficients to vector space elements
     * @param points Elements to combine
     */
    apply(points: ReadonlyArray<T>): T;

    /**
     * Combines with other coefficients (useful for tensor product)
     * @param other Other coefficients in same or different vector space
     */
    combine<U>(other: BlossomCoefficients<U>): BlossomCoefficients<[T, U]>;
}

/**
 * Context for efficient polar form evaluation in a vector space
 */
interface BlossomContext<T> {
    /**
     * Gets the sequence of knot spans needed for evaluation
     */
    readonly spans: ReadonlyArray<KnotSpan>;
    
    /**
     * Gets precomputed coefficients in the vector space
     */
    readonly coefficients: BlossomCoefficients<T>;
    
    /**
     * Gets indices of control points that influence the result
     */
    readonly influencingIndices: ReadonlyArray<number>;
}

/**
 * Implementation for real-valued coefficients
 */
class RealBlossomCoefficients implements BlossomCoefficients<number> {
    readonly vectorSpace = realVectorSpace;

    constructor(readonly values: ReadonlyArray<number>) {}

    apply(points: ReadonlyArray<number>): number {
        return points.reduce(
            (sum, p, i) => sum + p * this.values[i],
            0
        );
    }

    combine<U>(other: BlossomCoefficients<U>): BlossomCoefficients<[number, U]> {
        return new ProductBlossomCoefficients(this, other);
    }
}

/**
 * Implementation for vector-valued coefficients
 */
class VectorBlossomCoefficients implements BlossomCoefficients<[number, number]> {
    readonly vectorSpace = vector2DSpace;

    constructor(readonly values: ReadonlyArray<number>) {}

    apply(points: ReadonlyArray<[number, number]>): [number, number] {
        return points.reduce(
            (sum, p, i) => this.vectorSpace.add(
                sum,
                this.vectorSpace.multiply(this.values[i], p)
            ),
            [0, 0]
        );
    }

    combine<U>(other: BlossomCoefficients<U>): BlossomCoefficients<[[number, number], U]> {
        return new ProductBlossomCoefficients(this, other);
    }
}

/**
 * Tensor product coefficients for surfaces
 */
class ProductBlossomCoefficients<T, U> implements BlossomCoefficients<[T, U]> {
    readonly vectorSpace: VectorSpace<[T, U]>;

    constructor(
        private coeff1: BlossomCoefficients<T>,
        private coeff2: BlossomCoefficients<U>
    ) {
        this.vectorSpace = new ProductVectorSpace(
            coeff1.vectorSpace,
            coeff2.vectorSpace
        );
    }

    get values(): ReadonlyArray<number> {
        // Tensor product of coefficient values
        return this.coeff1.values.flatMap(
            v1 => this.coeff2.values.map(v2 => v1 * v2)
        );
    }

    apply(points: ReadonlyArray<[T, U]>): [T, U] {
        // Apply coefficients in each direction
        const n1 = this.coeff1.values.length;
        const n2 = this.coeff2.values.length;

        // Separate points into T and U components
        const points1: T[] = points.map(p => p[0]);
        const points2: U[] = points.map(p => p[1]);

        // Apply coefficients in each direction
        return [
            this.coeff1.apply(points1.slice(0, n1)),
            this.coeff2.apply(points2.slice(0, n2))
        ];
    }

    combine<V>(other: BlossomCoefficients<V>): BlossomCoefficients<[[T, U], V]> {
        return new ProductBlossomCoefficients(this, other);
    }
}

/**
 * Extended knot structure with vector-space aware blossoming
 */
interface BlossomKnotStructure<T> extends KnotStructure {
    /**
     * Creates a polar form evaluation context
     */
    createBlossomContext(
        direction: number,
        parameters: ReadonlyArray<number>,
        vectorSpace: VectorSpace<T>
    ): BlossomContext<T>;
}

// Usage examples
class BlossomKnots<T> extends BaseKnotStructure implements BlossomKnotStructure<T> {
    createBlossomContext(
        direction: number,
        parameters: ReadonlyArray<number>,
        vectorSpace: VectorSpace<T>
    ): BlossomContext<T> {
        // ... implementation ...
    }
}

// Example usage for curves
const curve = new BlossomKnots<[number, number]>();
const context2D = curve.createBlossomContext(0, [0.5], vector2DSpace);
const point = context2D.coefficients.apply(controlPoints);

// Example usage for surfaces
const surface = new BlossomKnots<[[number, number], number]>();
const contextSurface = surface.createBlossomContext(0, [0.5], 
    new ProductVectorSpace(vector2DSpace, realVectorSpace)
);

// Combining coefficients for tensor product surfaces
const uCoeffs = new VectorBlossomCoefficients([0.2, 0.6, 0.2]);
const vCoeffs = new RealBlossomCoefficients([0.5, 0.5]);
const surfaceCoeffs = uCoeffs.combine(vCoeffs);

// Evaluating with combined coefficients
const surfacePoint = surfaceCoeffs.apply(controlPoints);
