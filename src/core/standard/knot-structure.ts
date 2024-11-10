
/**
 * Knot structure for B-spline definition
 * Represents the knot vectors that define the B-spline
 */
export interface KnotStructure {
    /** Get the dimension of the knot structure */
    getDimension(): number;
    
    /** Get knot sequence for a given parametric direction */
    getKnotSequence(direction: number): ReadonlyArray<number>;
    
    /** Get multiplicity of a knot value */
    getMultiplicity?(direction: number, value: number): number;
}

export class Knots implements KnotStructure {
    constructor(private readonly knots: ReadonlyArray<number>) {}

    getDimension(): number {
        return 1;
    }

    getKnotSequence(direction: number): ReadonlyArray<number> {
        return this.knots;
    }
}


export class ProductKnots implements KnotStructure {
    constructor(private readonly knots: ReadonlyArray<ReadonlyArray<number>>) {}

    getDimension(): number {
        return this.knots.length;
    }

    getKnotSequence(direction: number): ReadonlyArray<number> {
        return this.knots[direction];
    }
}