/**
 * Represents a knot vector for B-spline curves.
 * Handles both periodic and non-periodic cases with proper domain mapping.
 */
export class Knots {
    private readonly values: number[];
    private readonly _isPeriodic: boolean;
    private readonly period: number;
    private readonly degree: number; 
    private readonly numControlPoints: number;
    private readonly EPSILON = 1e-10;

    /**
     * Creates a new knot vector.
     * For periodic B-splines, the knot sequence must satisfy:
     * u[i + numControlPoints] = u[i] + T, where T is the period
     *
     * @param knots - Array of knot values in non-decreasing order
     * @param degree - Degree of the B-spline curve
     * @param isPeriodic - Whether the knot vector represents a periodic B-spline
     * @param numControlPoints - Number of distinct control points (required for periodic case)
     * @throws {Error} If knot vector is not in non-decreasing order
     * @throws {Error} If periodic knot vector doesn't satisfy periodicity conditions
     */
    constructor(
        knots: number[], 
        degree: number,
        isPeriodic: boolean = false,
        numControlPoints?: number,
    ) {
        this.values = [...knots];
        this._isPeriodic = isPeriodic;
        this.degree = degree;

        if (isPeriodic && !numControlPoints) {
            throw new Error('Number of control points must be specified for periodic knot vector');
        }
        this.numControlPoints = numControlPoints || knots.length - degree - 1;
        this.validateKnots(knots);

    
        if (isPeriodic) {
            this.period = this.calculatePeriod();
            this.validatePeriodicity();
        } else {
            this.period = 0;
        }
    
    }

    /**
     * Gets the degree of the B-spline curve.
     * 
     * @returns The degree of the curve
     */
        public getDegree(): number {
            return this.degree;
        }


    /**
     * Gets the knot value at the specified index, handling periodic cases.
     * For periodic knots, implements the formula: u[i + k] = u[i] + T
     * where T is the period and k is the number of control points.
     * 
     * @param index - The index of the knot value to retrieve
     * @returns The knot value at the specified index
     */
    public getValue(index: number): number {
        if (this._isPeriodic) {
            return this.getPeriodicValue(index);
        }
        return this.values[index];
    }

    public getNumberOfControlPoints() : number {
        return this.numControlPoints
    }

    public getNumberOfWrappedControlPoints() : number {
        if (this._isPeriodic) {
            return this.numControlPoints + this.degree;
        }
        else {
            throw new Error('Cannot get number of wrapped control points for non-periodic knot vector');
        }
    }

    public getDomain() : {minValue: number, maxValue: number} {
        return { minValue: this.values[this.degree], maxValue: this.values[this.values.length - this.degree - 1] };
    }

    /**
     * Finds the knot span index for a given parameter value.
     * Uses binary search to efficiently locate the span.
     * 
     * @param u - Parameter value
     * @returns Index i where u ∈ [uᵢ, uᵢ₊₁)
     */
    private findSpanIndex(u: number): number {
        const n = this.values.length - 1;
        
        // Handle boundary cases
        if (u >= this.values[n]) return n - 1;
        if (u <= this.values[0]) return 0;

        // Binary search
        let low = 0;
        let high = n;
        let mid = Math.floor((low + high) / 2);

        while (u < this.values[mid] || u >= this.values[mid + 1]) {
            if (u < this.values[mid]) {
                high = mid;
            } else {
                low = mid;
            }
            mid = Math.floor((low + high) / 2);
        }

        return mid;
    }

    /**
     * Inserts a single knot value in a non-periodic knot vector.
     * 
     * @param u - Knot value to insert
     * @returns New Knots instance with inserted knot
     */
    private insertSingle(u: number): Knots {
        const newKnots = [...this.values];
        const insertIndex = this.findSpanIndex(u) + 1;
        newKnots.splice(insertIndex, 0, u);
        return new Knots(newKnots, this.degree);
    }

    /**
     * Inserts a knot value in a periodic knot vector.
     * Maintains periodicity by inserting corresponding knots in all periods.
     * 
     * @param u - Knot value to insert
     * @returns New Knots instance with inserted knots
     */
    private insertPeriodic(u: number): Knots {
        // Map u to fundamental domain
        const uMapped = this.mapToFundamentalDomain(u);
        
        // Find base insertion point
        const baseIndex = this.findSpanIndex(uMapped) + 1;
        
        // Create new knot array
        const newKnots = [...this.values];
        
        // Insert knot in fundamental domain
        newKnots.splice(baseIndex, 0, uMapped);
        
        // Insert corresponding periodic knots
        const numPeriods = Math.ceil(
            (this.values.length - baseIndex) / this.values.length
        );
        
        for (let i = 1; i <= numPeriods; i++) {
            const periodicValue = uMapped + i * this.period;
            const periodicIndex = baseIndex + i * (this.numControlPoints + 1); // after degree elevation the period is numControlPoints + 1
            newKnots.splice(periodicIndex, 0, periodicValue);
        }

        return new Knots(newKnots, this.degree, true, this.numControlPoints + 1);
    }

    /**
     * Elevates degree for standard (non-periodic) knot vector.
     * Increases multiplicities of all knots by one.
     * 
     * @returns New Knots instance with elevated degree
     */
    private elevateStandardKnots(): Knots {
        const newKnots: number[] = [];
        
        // Process each distinct knot value
        let i = 0;
        while (i < this.values.length) {
            const currentKnot = this.values[i];
            let multiplicity = 1;
            
            // Count multiplicity
            while (i + 1 < this.values.length && 
                   Math.abs(this.values[i + 1] - currentKnot) < this.EPSILON) {
                multiplicity++;
                i++;
            }
            
            // Add knot with increased multiplicity
            for (let j = 0; j <= multiplicity; j++) {
                newKnots.push(currentKnot);
            }
            
            i++;
        }
        
        return new Knots(newKnots, this.degree + 1);
    }

    /**
     * Elevates degree for periodic knot vector.
     * Maintains periodicity while increasing multiplicities.
     * 
     * @returns New Knots instance with elevated degree
     */
    private elevatePeriodicKnots(): Knots {

        function elevateKnots(distinctKnots: number[], multiplicities: number[]): number[] {
            const elevatedKnots: number[] = [];
            for (let i = 0; i < distinctKnots.length; i++) {
                const knot = distinctKnots[i];
                const multiplicity = multiplicities[i];
                for (let j = 0; j <= multiplicity; j++) {
                    elevatedKnots.push(knot);
                }
            }
            return elevatedKnots;
        }

        // Extract the knot in the domain
        const {minValue, maxValue} = this.getDomain()
        const domainKnots = this.values.filter(knot => knot >= minValue && knot < maxValue)
        const distinctDomainKnots = Array.from(new Set(domainKnots));
        const domainKnotsMultiplicities = getKnotMultiplicities(domainKnots)
        const period = this.getPeriod()

        const newDomainKnots = elevateKnots(distinctDomainKnots, domainKnotsMultiplicities)

        const newDegree = this.degree + 1;
        const newNumControlPoints = newDomainKnots.length
        const newKnots = extendPeriodicArray(newDomainKnots, period, newDegree, newDegree + 1)
        return new Knots(newKnots, newDegree, true, newNumControlPoints)
    }
    
    

    /**
     * Gets all distinct knot values in the vector.
     * 
     * @returns Array of distinct knot values
     */
    private getDistinctKnots(): number[] {
        return Array.from(new Set(this.values));
    }


    /**
     * Finds the knot span index for a given parameter value.
     * For periodic cases, first maps the parameter to the fundamental domain.
     * 
     * @param u - Parameter value
     * @returns Index i where u ∈ [uᵢ, uᵢ₊₁)
     */
    public getSpan(u: number): number {
        if (this._isPeriodic) {
            return this.findSpanIndex(this.mapToFundamentalDomain(u));
        }
        return this.findSpanIndex(u);
    }

    /**
     * Inserts a new knot value while maintaining the knot vector properties.
     * For periodic cases, ensures periodicity is preserved.
     * 
     * @param u - Knot value to insert
     * @returns New Knots instance with inserted knot
     */
    public insert(u: number): Knots {
        if (this._isPeriodic) {
            return this.insertPeriodic(u);
        }
        return this.insertSingle(u);
    }

    /**
     * Performs degree elevation on the knot vector.
     * Increases multiplicities appropriately for both periodic and non-periodic cases.
     * 
     * @returns New Knots instance with elevated degree
     */
    public elevate(): Knots {
        if (this._isPeriodic) {
            return this.elevatePeriodicKnots();
        }
        return this.elevateStandardKnots();
    }

    /**
     * Gets the multiplicity of a knot value at a given index.
     * 
     * @param index - Index to check multiplicity
     * @returns Number of times the knot value appears consecutively
     */
    public getMultiplicity(index: number): number {
        let mult = 1;
        const value = this.getValue(index);
        
        let i = index + 1;
        while (i < this.values.length && 
               Math.abs(this.getValue(i) - value) < this.EPSILON) {
            mult++;
            i++;
        }
        
        return mult;
    }
    

    /**
     * Checks if the knot vector is periodic.
     * 
     * @returns True if the knot vector is periodic
     */
    public isPeriodic(): boolean {
        return this._isPeriodic;
    }

    /**
     * Gets the period of the knot vector.
     * 
     * @returns Period value for periodic knot vectors, 0 for non-periodic
     */
    public getPeriod(): number {
        return this.period;
    }

    /**
     * Maps a parameter value to the fundamental domain [u₀, u₀ + T)
     * where T is the period.
     * 
     * @param u - Parameter value to map
     * @returns Mapped parameter value in fundamental domain
     */
    private mapToFundamentalDomain(u: number): number {
        return u - Math.floor((u - this.values[0]) / this.period) * this.period;
    }

    /**
     * Gets a knot value for periodic knot vectors.
     * Implements the periodic extension formula.
     * 
     * @param index - Index of the knot value
     * @returns Periodic knot value
     */
    private getPeriodicValue(index: number): number {
        // Helper function for proper mathematical modulo
        const mod = (n: number, m: number): number => {
            return ((n % m) + m) % m;
        };

        // Get the wrapped index within the knot vector
        const baseIndex = mod(index, this.numControlPoints);
        
        // Calculate how many full periods we've moved
        const period = Math.floor(index / this.numControlPoints);
        
        // Apply the base value plus the period shift
        return this.values[baseIndex] + period * this.period;
    }


    /**
     * Validates the knot vector during construction.
     * 
     * @param knots - Knot vector to validate
     * @throws {Error} If knot vector is invalid
     */
    private validateKnots(knots: number[]): void {
        if (knots.length < 2) {
            throw new Error('Knot vector must have at least 2 values');
        }

        for (let i = 1; i < knots.length; i++) {
            if (knots[i] < knots[i - 1]) {
                throw new Error('Knot vector must be non-decreasing');
            }
        }
    }

    /**
     * Validates periodicity conditions.
     * For periodic B-splines:
     * 1. The knot sequence must be non-decreasing
     * 2. u[i + numControlPoints] = u[i] + T for all applicable i
     */
    private validatePeriodicity(): void {
        // Check minimum length requirement
        const minLength = this.degree + this.numControlPoints + 1;
        if (this.values.length < minLength) {
            throw new Error(
                `Periodic knot vector must have at least ${minLength} knots`
            );
        }

        // Check non decreasing periodic case
        for (let i = 1; i < this.values.length; i++) {
            if (this.values[i] - this.values[i - 1] < 0) {
                throw new Error('Periodic knot vector must be non decreasing');
            }
        }

        // Check periodic extension
        for (let i = 0; i < this.values.length - this.numControlPoints; i++) {
            const diff = this.values[i + this.numControlPoints] - this.values[i];
            if (Math.abs(diff - this.period) > this.EPSILON) {
                //console.log(this.numControlPoints)
                //console.log(i)
                throw new Error(
                    'Periodic knot vector must satisfy u[i + n] = u[i] + T'
                );
            }
        }
    }

    /**
     * Calculates the period of a periodic knot vector.
     * Period T is defined as u[i + numControlPoints] - u[i]
     */
    private calculatePeriod(): number {
        // Take the difference between knots separated by numControlPoints
        const period = this.values[this.numControlPoints] - this.values[0];
        return period;
    }

    
}


/**
 * Gets the multiplicities of knots in ascending order.
 * 
 * @param knots Array of knot values in non-decreasing order
 * @param epsilon Optional tolerance for floating point comparison (default: 1e-10)
 * @returns Array of multiplicities for each distinct knot value
 * 
 * @example
 * getKnotMultiplicities([1,1,2,2,2,3,4,4]) // returns [2,3,1,2]
 * getKnotMultiplicities([0,1,1,2,2,3,3,4,4,5,5,6,6]) // returns [1,2,2,2,2,2,2]
 */
function getKnotMultiplicities(knots: number[], epsilon: number = 1e-10): number[] {
    if (knots.length === 0) return [];
    
    const multiplicities: number[] = [];
    
    let i = 0;
    while (i < knots.length) {
        const currentKnot = knots[i];
        let multiplicity = 1;
        
        // Count multiplicity of current knot value
        while (i + 1 < knots.length && 
               Math.abs(knots[i + 1] - currentKnot) < epsilon) {
            multiplicity++;
            i++;
        }
        
        multiplicities.push(multiplicity);
        i++;
    }
    
    return multiplicities;
}

/**
 * Extends a periodic array by adding periodic values before and after
 * 
 * @param values - Original array of numbers
 * @param period - The period value
 * @param numBefore - Number of values to add before
 * @param numAfter - Number of values to add after
 * @returns Extended array with periodic values
 */
export function extendPeriodicArray(
    values: number[], 
    period: number,
    numBefore: number, 
    numAfter: number
): number[] {
    if (values.length === 0) return [];
    
    const result: number[] = [...values];
    
    // Add values before
    for (let i = 0; i < numBefore; i++) {
        const baseIndex = ((values.length - 1) - (i % values.length));
        const periodShift = Math.ceil((i + 1) / values.length);
        result.unshift(values[baseIndex] - period * periodShift);
    }
    
    // Add values after
    for (let i = 0; i < numAfter; i++) {
        const baseIndex = i % values.length;
        const periodShift = Math.floor(i / values.length) + 1;
        result.push(values[baseIndex] + period * periodShift);
    }
    
    return result;
}


/**
 * Represents a periodic knot sequence defined by a repeating pattern and a period length.
 * 
 * @interface PeriodicKnots
 * @property {number[]} pattern - The base sequence of knot values that repeats
 * @property {number} period - The length of one complete period (the spacing between repetitions)
 * 
 * @example
 * // For sequence [-5, -3, -2, 0, 1, 3, 4, 6, 7, 9, ...]
 * const knots: PeriodicKnots = {
 *     pattern: [1, 3],  // Base pattern
 *     period: 3         // Each repetition adds 3
 * };
 */
export interface PeriodicKnots {
    pattern: number[];    
    period: number;      
}

/**
 * Calculates the knot value at a given index in a periodic knot sequence.
 * 
 * @param {PeriodicKnots} knots - The periodic knot sequence definition
 * @param {number} index - The index at which to calculate the knot value
 * @returns {number} The knot value at the specified index
 * @throws {Error} If the pattern array is empty
 */

export function getKnotValue(knots: PeriodicKnots, index: number): number {
    const { pattern, period } = knots;
    
    if (pattern.length === 0) {
        throw new Error('Pattern array cannot be empty');
    }

    const baseIndex = ((index % pattern.length) + pattern.length) % pattern.length;
    const numPeriods = Math.floor(index / pattern.length);
    
    return pattern[baseIndex] + numPeriods * period;
}

//note :
// The number of periodic knot inside a period should be the number of periodic control point





