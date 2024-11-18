/**
 * Represents a periodic knot vector for closed B-spline curves.
 * A periodic knot vector repeats with a specified period and is used
 * to define closed B-spline curves.
 * Note : The number of periodic knot in the pattern is equal
 * the number of periodic control point
 */

import { BaseKnotStructure, KnotStructure } from "./knot-structure2";

export class PeriodicKnots extends BaseKnotStructure {

    private readonly domain: { min: number; max: number };


    /**
     * Creates a new PeriodicKnots instance
     * @param pattern - Array of knot values that define the basic pattern
     * @param period - The period after which the pattern repeats
     * @throws {Error} If the knot vector is invalid
     */
    constructor(
        private readonly pattern: number[],
        private readonly period: number,
    ) {
        super();
        if (period <= 0) {
            throw new Error('Period must be a positive number');
        }
        this.validateKnots(pattern);
        this.domain = {
            min: pattern[0],
            max: pattern[0] + period
        };
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
        const insertIndex = this.findInsertionIndex(normalizedU);
        const newPattern = [...this.pattern];
        newPattern.splice(insertIndex, 0, normalizedU);
        return new PeriodicKnots(newPattern, this.period);
    }

    withRemovedKnots(dimension: number): KnotStructure {
        this.validateDirection(dimension);
        // For periodic knots, we maintain the pattern structure
        const newPattern = [
            this.pattern[0],
            ...this.pattern.slice(1, -1).filter((_, index) => index % 2 === 0),
            this.pattern[this.pattern.length - 1]
        ];
        return new PeriodicKnots(newPattern, this.period);
    }

    // Factory method
    public static fromDistinctKnots(values: DistinctKnotsWithMultiplicities, period: number): PeriodicKnots {
        const pattern = values.knots.flatMap((knot, index) => 
            Array(values.multiplicities[index]).fill(knot)
        );
        return new PeriodicKnots(pattern, period);
    }

    // Query methods
    public getDomain(): KnotDomain {
        const {pattern, period} = this;
        return {minValue: pattern[0], maxValue: pattern[0] + period}
    }

    public getDistinctKnots() {
        return Array.from(new Set(this.pattern))
    }

    public getPatternLength() {
        return this.pattern.length;
    }

    public getKnotPattern() {
        return this.pattern;
    }

    public getMultiplicities() {
        const multiplicities = new Map<number, number>();
        this.pattern.forEach(knot => {
            multiplicities.set(knot, (multiplicities.get(knot) || 0) + 1);
        });
        return Array.from(multiplicities.values());
    }

    public getKnotValue(index: number): number {        
        const baseIndex = ((index % this.pattern.length) + this.pattern.length) % this.pattern.length;
        const numPeriods = Math.floor(index / this.pattern.length);
        
        return this.pattern[baseIndex] + numPeriods * this.period;
    }

    public getDistinctKnotsWithMultiplicities(): DistinctKnotsWithMultiplicities {
        const knots = this.getDistinctKnots();
        const multiplicities = this.getMultiplicities();
        return { knots, multiplicities };
    }

    // Modification methods

    public insertKnot(u: number) : PeriodicKnots {
        const uInsideDomain = this.shiftInsideDomain(u);
        const index = this.findKnotSpan(uInsideDomain);
        const newPattern = [...this.pattern];
        newPattern.splice(index + 1, 0, uInsideDomain);
        return new PeriodicKnots(newPattern, this.period);
    }

    public elevateDegree() : PeriodicKnots {
        const {pattern, period} = this;
        const { knots, multiplicities } = this.getDistinctKnotsWithMultiplicities()
        const newMultiplicities = multiplicities.map(m => m + 1);
        return PeriodicKnots.fromDistinctKnots({knots, multiplicities: newMultiplicities}, this.period)
    }


    /**
     * 
     * @param degree Degree of the b-spline curve
     * @param u Position to insert a new knot
     * @returns An unroll knots
     */
    public unrollForKnotInsertion(degree: number, u: number) {
        const numUnroll = Math.ceil((degree + 1) / this.pattern.length)
        const startIndex = - this.pattern.length * numUnroll
        const endIndex = this.pattern.length * (numUnroll + 1)
        let newPattern: number[] = []
        for (let i = startIndex; i < endIndex; i++) {
            newPattern.push(this.getKnotValue(i))
        }
        const knotsToBeInserted: number[] = []
        for (let i = -numUnroll; i <= numUnroll; i++) {
            knotsToBeInserted.push(u + i * this.period)
        }
        return {
            knots: new PeriodicKnots(newPattern, this.period * (1 + 2 * numUnroll)),
            knotsToBeInserted
        }
    }


    // Utility methods
    public findKnotSpan(u: number): number {
        // For periodic behavior, we need to handle values outside the base period
        const periodOffset = Math.floor(u / this.period) * this.pattern.length;
        
        // Normalize u to the first period
        const uInsideDomain = this.shiftInsideDomain(u);
        
        // Binary search to find the span index
        let low = 0;
        let high = this.pattern.length - 1;
        
        // If u is beyond the last knot in the pattern
        if (uInsideDomain >= this.pattern[high]) {
            return high + periodOffset;
        }
        
        // If u is at or before the first knot
        if (uInsideDomain <= this.pattern[0]) {
            return periodOffset;
        }
        
        // Binary search
        while (low < high) {
            const mid = Math.floor((low + high) / 2);
            
            if (uInsideDomain >= this.pattern[mid] && uInsideDomain < this.pattern[mid + 1]) {
                return mid + periodOffset;
            }
            
            if (uInsideDomain < this.pattern[mid]) {
                high = mid;
            } else {
                low = mid + 1;
            }
        }
        
        return low + periodOffset;
    }
    
    


    public shiftInsideDomain(u: number): number {
        const {minValue, maxValue} = this.getDomain();
        const width = maxValue - minValue;
        return minValue + ((u - minValue) % width + width) % width;
    }


    // Private helper methods

    /**
     * Validates the knot vector during construction.
     * 
     * @param knots - Knot vector to validate
     * @throws {Error} If knot vector is invalid
     */
    private validateKnots(knots: number[]): void {
        if (knots.length < 2) {
            throw new Error('Knot vector must have at least 2 values for a valid B-spline');
        }

        if (!this.isNonDecreasing(knots)) {
            throw new Error('Knot vector values must be in non-decreasing order');
        }

        if (this.exceedsPeriod(knots)) {
            throw new Error(
                `Pattern range (${knots[knots.length - 1] - knots[0]}) ` +
                `exceeds specified period (${this.period})`
            );
        }
    }

    protected isNonDecreasing(knots: readonly number[]): boolean {
        return knots.every((val, i) => i === 0 || val >= knots[i - 1]);
    }
    
    private exceedsPeriod(knots: number[]): boolean {
        return knots[knots.length - 1] - knots[0] >= this.period;
    }

    public getRelevantKnotsAndIndices(insertPosition: number, degree: number): {
        knots: number[],
        startIndex: number,
        endIndex: number
    } {
        // Get the pattern length
        const patternLength = this.getKnotPattern().length;
        
        // Find which periodic segment contains our insertion point
        const segment = Math.floor(insertPosition / patternLength);
        
        // Unroll enough knots for the insertion
        const unrolledKnots = this.unrollKnotsAroundPosition(insertPosition, degree);
        
        // Find the knot span containing the insertion position
        const insertSpan = this.findKnotSpan(insertPosition);
        
        // Calculate the relevant range
        const startIndex = Math.max(0, insertSpan - degree);
        const endIndex = Math.min(unrolledKnots.length - 1, insertSpan + degree + 1);

        return {
            knots: unrolledKnots,
            startIndex,
            endIndex
        }

    }

    private unrollKnotsAroundPosition(position: number, degree: number): number[] {
        // Similar to the previous implementation but as a method of PeriodicKnot
        // ...
    }

    private validatePattern(pattern: ReadonlyArray<number>, period: number): void {
        if (pattern.length < 2) {
            throw new Error('Pattern must have at least 2 values');
        }
        if (!this.isNonDecreasing(pattern)) {
            throw new Error('Pattern must be non-decreasing');
        }
        if (pattern[pattern.length - 1] - pattern[0] >= period) {
            throw new Error('Pattern range exceeds period');
        }
        if (period <= 0) {
            throw new Error('Period must be positive');
        }
    }

}

interface KnotDomain {
    minValue: number;
    maxValue: number;
}

interface DistinctKnotsWithMultiplicities {
    knots: number[];
    multiplicities: number[];
}



