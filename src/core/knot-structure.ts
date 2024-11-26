/**
 * Represents a knot value with its properties
 */
export interface KnotValue {
    readonly value: number;
    readonly multiplicity: number;
}

/**
 * Interface for distinct knots with multiplicities
 */

export interface DistinctKnotsWithMultiplicities {
    readonly knots: ReadonlyArray<number>;
    readonly multiplicities: ReadonlyArray<number>;
}

/**
 * Base interface for all knot structures
 */
export interface KnotStructure {
    /**
     * Gets the dimension of the knot structure
     */
    getDimension(): number;

    /**
     * Gets the knot sequence for a specific direction
     * @param direction - The parametric direction
     */
    getKnotSequence(direction: number): ReadonlyArray<number>;

    /**
     * Creates a new knot structure with an inserted knot
     * @param dimension - The dimension where to insert the knot
     * @param u - The knot value to insert
     */
    withInsertedKnot(u: number, dimension: number): KnotStructure;

    /**
     * Creates a new knot structure with removed knots
     * @param dimension - The dimension where to remove knots
     */
    withRemovedKnot(index: number, dimension: number): KnotStructure;

    /**
     * Gets the domain bounds for a specific direction
     * @param direction - The parametric direction
     */
    getDomain(direction: number): Domain;

    /**
     * Gets distinct knots with their multiplicities
     * @param direction - The parametric direction
     */
    getDistinctKnots(direction: number): ReadonlyArray<KnotValue>;
}


/**
 * Represents a domain interval
 */
export interface Domain {
    readonly min: number;
    readonly max: number;
}

/**
 * Abstract base class providing common functionality for knot structures
 */
export abstract class BaseKnotStructure implements KnotStructure {
    abstract getDimension(): number;
    abstract getKnotSequence(direction: number): ReadonlyArray<number>;
    abstract withInsertedKnot(u: number, dimension: number): KnotStructure;
    abstract withRemovedKnot(index: number, dimension: number): KnotStructure;
    abstract getDomain(direction: number): Domain;
    abstract getDistinctKnots(direction: number): ReadonlyArray<KnotValue>;

    /**
     * Validates a direction index
     */
    protected validateDirection(direction: number): void {
        if (direction < 0 || direction >= this.getDimension()) {
            throw new Error(`Invalid direction: ${direction}`);
        }
    }

    /**
     * Checks if a sequence is non-decreasing
     */
    protected isNonDecreasing(knots: ReadonlyArray<number>): boolean {
        return knots.every((val, i) => i === 0 || val >= knots[i - 1]);
    }

    /**
     * Computes distinct knots with their multiplicities
     */
    protected computeDistinctKnots(knots: ReadonlyArray<number>): ReadonlyArray<KnotValue> {
        const result: KnotValue[] = [];
        let currentValue = knots[0];
        let multiplicity = 1;

        for (let i = 1; i < knots.length; i++) {
            if (Math.abs(knots[i] - currentValue) < Number.EPSILON) {
                multiplicity++;
            } else {
                result.push({ value: currentValue, multiplicity });
                currentValue = knots[i];
                multiplicity = 1;
            }
        }
        result.push({ value: currentValue, multiplicity });

        return result;
    }

    /**
     * Finds insertion index for a new knot value
     */
    protected findInsertionIndex(knots: ReadonlyArray<number>, u: number): number {
        return knots.findIndex(k => k > u);
    }
}

/**
 * Basic implementation for single-dimensional knot vectors
 */
export class Knots extends BaseKnotStructure {
    constructor(private readonly knots: ReadonlyArray<number>) {
        super();
        this.validateKnots(knots);
    }

    getDimension(): number {
        return 1;
    }

    getKnotSequence(direction: number = 0): ReadonlyArray<number> {
        this.validateDirection(direction);
        return this.knots;
    }

    withInsertedKnot(u: number, dimension: number = 0): Knots {
        this.validateDirection(dimension);
        const insertIndex = this.findInsertionIndex(this.knots, u);
        const newKnots = [...this.knots];
        newKnots.splice(insertIndex >= 0 ? insertIndex : this.knots.length, 0, u);
        return new Knots(newKnots);
    }

    withRemovedKnot(index: number, dimension: number = 0): Knots {
        this.validateDirection(dimension);
        const newKnots = this.knots.filter((_, i) => i !== index);
        return new Knots(newKnots);
    }

    getDomain(direction: number): Domain {
        this.validateDirection(direction);
        return {
            min: this.knots[0],
            max: this.knots[this.knots.length - 1]
        };
    }

    getDistinctKnots(direction: number): ReadonlyArray<KnotValue> {
        this.validateDirection(direction);
        return this.computeDistinctKnots(this.knots);
    }

    private validateKnots(knots: ReadonlyArray<number>): void {
        if (!knots || knots.length < 2) {
            throw new Error('Knot vector must have at least 2 values');
        }
        if (!this.isNonDecreasing(knots)) {
            throw new Error('Knot vector must be non-decreasing');
        }
    }
}

/**
 * Product knot structure for multi-dimensional non-periodic B-splines
 */
export class ProductKnots extends BaseKnotStructure {
    constructor(private readonly knotVectors: ReadonlyArray<ReadonlyArray<number>>) {
        super();
        this.validateKnotVectors(knotVectors);
    }

    getDimension(): number {
        return this.knotVectors.length;
    }

    getKnotSequence(direction: number): ReadonlyArray<number> {
        this.validateDirection(direction);
        return this.knotVectors[direction];
    }

    withInsertedKnot(dimension: number, u: number): KnotStructure {
        this.validateDirection(dimension);
        const newKnotVectors = [...this.knotVectors];
        const knots = this.knotVectors[dimension];
        const insertIndex = this.findInsertionIndex(knots, u);
        const newKnots = [...knots];
        newKnots.splice(insertIndex >= 0 ? insertIndex : knots.length, 0, u);
        newKnotVectors[dimension] = newKnots;
        return new ProductKnots(newKnotVectors);
    }

    withRemovedKnot(index: number, dimension: number = 0): KnotStructure {
        this.validateDirection(dimension);
        const newKnotVectors = [...this.knotVectors];
        const knots = this.knotVectors[dimension];
        const newKnots = knots.filter((_, i) => i !== index);
        newKnotVectors[dimension] = newKnots;
        return new ProductKnots(newKnotVectors);
    }

    getDomain(direction: number): Domain {
        this.validateDirection(direction);
        const knots = this.knotVectors[direction];
        return {
            min: knots[0],
            max: knots[knots.length - 1]
        };
    }

    getDistinctKnots(direction: number): ReadonlyArray<KnotValue> {
        this.validateDirection(direction);
        return this.computeDistinctKnots(this.knotVectors[direction]);
    }

    private validateKnotVectors(knotVectors: ReadonlyArray<ReadonlyArray<number>>): void {
        if (!knotVectors || knotVectors.length === 0) {
            throw new Error('At least one knot vector must be provided');
        }
        knotVectors.forEach((knots, i) => {
            if (!knots || knots.length < 2) {
                throw new Error(`Knot vector ${i} must have at least 2 values`);
            }
            if (!this.isNonDecreasing(knots)) {
                throw new Error(`Knot vector ${i} must be non-decreasing`);
            }
        });
    }
}
