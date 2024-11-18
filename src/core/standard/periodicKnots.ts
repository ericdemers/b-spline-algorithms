import { BaseKnotStructure, KnotStructure, Domain, KnotValue } from './knot-structure';

/**
 * Interface for distinct knots with multiplicities
 */
export interface DistinctKnotsWithMultiplicities {
    readonly knots: ReadonlyArray<number>;
    readonly multiplicities: ReadonlyArray<number>;
}

/**
 * Represents a periodic knot vector for closed B-spline curves.
 * A periodic knot vector repeats with a specified period and is used
 * to define closed B-spline curves.
 */
export class PeriodicKnots extends BaseKnotStructure {
    private readonly domain: Domain;
    private readonly distinctKnotsCache: ReadonlyArray<KnotValue>;

    /**
     * Creates a new PeriodicKnots instance
     * @param pattern - Array of knot values that define the basic pattern
     * @param period - The period after which the pattern repeats
     * @throws {Error} If the knot vector is invalid
     */
    constructor(
        private readonly pattern: ReadonlyArray<number>,
        private readonly period: number,
    ) {
        super();
        this.validateConstructorParams(pattern, period);
        this.domain = this.computeDomain();
        this.distinctKnotsCache = this.computeDistinctKnots(pattern);
    }

    getDimension(): number {
        return 1;
    }

    getKnotSequence(direction: number): ReadonlyArray<number> {
        this.validateDirection(direction);
        return this.unrollKnots(3); // Default to 3 periods
    }

    withInsertedKnot(dimension: number, u: number): KnotStructure {
        this.validateDirection(dimension);
        const normalizedU = this.normalizeParameter(u);
        const insertIndex = this.findInsertionIndex(this.pattern, normalizedU);
        const newPattern = [...this.pattern];
        newPattern.splice(insertIndex, 0, normalizedU);
        return new PeriodicKnots(newPattern, this.period);
    }

    withRemovedKnots(dimension: number): KnotStructure {
        this.validateDirection(dimension);
        // For periodic knots, maintain pattern structure while removing interior knots
        const newPattern = [
            this.pattern[0],
            ...this.pattern.slice(1, -1).filter((_, index) => index % 2 === 0),
            this.pattern[this.pattern.length - 1]
        ];
        return new PeriodicKnots(newPattern, this.period);
    }

    getDomain(direction: number): Domain {
        this.validateDirection(direction);
        return { ...this.domain };
    }

    getDistinctKnots(direction: number): ReadonlyArray<KnotValue> {
        this.validateDirection(direction);
        return this.distinctKnotsCache;
    }

    /**
     * Gets the basic pattern of the periodic knot vector
     */
    getPattern(): ReadonlyArray<number> {
        return this.pattern;
    }

    /**
     * Gets the period of the knot vector
     */
    getPeriod(): number {
        return this.period;
    }

    /**
     * Creates a periodic knot structure from distinct knots with multiplicities
     */
    static fromDistinctKnots(
        values: DistinctKnotsWithMultiplicities,
        period: number
    ): PeriodicKnots {
        const pattern = values.knots.flatMap((knot, index) =>
            Array(values.multiplicities[index]).fill(knot)
        );
        return new PeriodicKnots(pattern, period);
    }

    /**
     * Unrolls the knot pattern for evaluation
     * @param numPeriods - Number of periods to unroll
     */
    unrollKnots(numPeriods: number): ReadonlyArray<number> {
        const result: number[] = [];
        for (let i = -numPeriods; i <= numPeriods; i++) {
            this.pattern.forEach(knot => {
                result.push(knot + i * this.period);
            });
        }
        return result;
    }

    /**
     * Unrolls knots for knot insertion operation
     * @param degree - Degree of the B-spline
     * @param u - Parameter value for insertion
     */
    unrollForKnotInsertion(degree: number, u: number): {
        knots: KnotStructure;
        knotsToBeInserted: ReadonlyArray<number>;
    } {
        const numUnroll = Math.ceil((degree + 1) / this.pattern.length);
        const startIndex = -this.pattern.length * numUnroll;
        const endIndex = this.pattern.length * (numUnroll + 1)
        
        const newPattern: number[] = [];
        for (let i = startIndex; i < endIndex; i++) {
            newPattern.push(this.getKnotValue(i));
        }

        const knotsToBeInserted = Array.from(
            { length: 2 * numUnroll + 1 },
            (_, i) => u + (i - numUnroll) * this.period
        );

        return {
            knots: new PeriodicKnots(newPattern, this.period * (2 * numUnroll + 1)),
            knotsToBeInserted
        };
    }

    /**
     * Gets a knot value at a specific index, handling periodic repetition
     */
    getKnotValue(index: number): number {
        const baseIndex = ((index % this.pattern.length) + this.pattern.length) % this.pattern.length;
        const numPeriods = Math.floor(index / this.pattern.length);
        return this.pattern[baseIndex] + numPeriods * this.period;
    }

    /**
     * Finds the knot span containing a parameter value
     */
    findKnotSpan(u: number): number {
        const normalizedU = this.normalizeParameter(u);
        const periodOffset = Math.floor(u / this.period) * this.pattern.length;
        
        // Binary search in the pattern
        let low = 0;
        let high = this.pattern.length - 1;
        
        if (normalizedU >= this.pattern[high]) {
            return high + periodOffset;
        }
        
        if (normalizedU <= this.pattern[0]) {
            return periodOffset;
        }
        
        while (low < high) {
            const mid = Math.floor((low + high) / 2);
            if (normalizedU >= this.pattern[mid] && normalizedU < this.pattern[mid + 1]) {
                return mid + periodOffset;
            }
            if (normalizedU < this.pattern[mid]) {
                high = mid;
            } else {
                low = mid + 1;
            }
        }
        
        return low + periodOffset;
    }

    private normalizeParameter(u: number): number {
        const { min, max } = this.domain;
        return min + ((u - min) % this.period + this.period) % this.period;
    }

    private computeDomain(): Domain {
        return {
            min: this.pattern[0],
            max: this.pattern[0] + this.period
        };
    }

    private validateConstructorParams(
        pattern: ReadonlyArray<number>,
        period: number
    ): void {
        if (!pattern || pattern.length < 2) {
            throw new Error('Pattern must have at least 2 values');
        }

        if (period <= 0) {
            throw new Error('Period must be positive');
        }

        if (!this.isNonDecreasing(pattern)) {
            throw new Error('Pattern must be non-decreasing');
        }

        if (pattern[pattern.length - 1] - pattern[0] >= period) {
            throw new Error(
                `Pattern range (${pattern[pattern.length - 1] - pattern[0]}) ` +
                `exceeds specified period (${period})`
            );
        }
    }

    /**
     * Gets relevant knots and indices for operations like knot insertion
     */
    getRelevantKnotsAndIndices(insertPosition: number, degree: number): {
        knots: ReadonlyArray<number>;
        startIndex: number;
        endIndex: number;
    } {
        const patternLength = this.pattern.length;
        const segment = Math.floor(insertPosition / patternLength);
        const unrolledKnots = this.unrollKnots(Math.ceil(degree / patternLength));
        const insertSpan = this.findKnotSpan(insertPosition);
        
        const startIndex = Math.max(0, insertSpan - degree);
        const endIndex = Math.min(unrolledKnots.length - 1, insertSpan + degree + 1);

        return {
            knots: unrolledKnots.slice(startIndex, endIndex + 1),
            startIndex,
            endIndex
        };
    }
}
