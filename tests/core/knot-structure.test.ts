import { Knots, ProductKnots, KnotValue } from '../../src/core/knot-structure';

describe('Knots', () => {
    let knots: Knots;

    beforeEach(() => {
        // Initialize a Knots instance with a valid knot vector
        knots = new Knots([0, 0, 0, 0, 1, 1, 1, 1]);
    });

    test('getDimension returns correct dimension', () => {
        expect(knots.getDimension()).toBe(1);
    });

    test('getKnotSequence returns correct sequence', () => {
        expect(knots.getKnotSequence()).toEqual([0, 0, 0, 0, 1, 1, 1, 1]);
    });

    test('withInsertedKnot inserts a knot correctly', () => {
        const newKnots = knots.withInsertedKnot(0.5);
        expect(newKnots.getKnotSequence()).toEqual([0, 0, 0, 0, 0.5, 1, 1, 1, 1]);
    });

    test('withRemovedKnot removes a knot correctly', () => {
        const newKnots = knots.withRemovedKnot(3); // Remove the knot at index 3
        expect(newKnots.getKnotSequence()).toEqual([0, 0, 0, 1, 1, 1, 1]);
    });

    test('getDomain returns correct domain', () => {
        const domain = knots.getDomain(0);
        expect(domain).toEqual({ min: 0, max: 1 });
    });

    test('getDistinctKnots returns correct distinct knots', () => {
        const distinctKnots = knots.getDistinctKnots(0);
        expect(distinctKnots).toEqual([
            { value: 0, multiplicity: 4 },
            { value: 1, multiplicity: 4 },
        ]);
    });

    test('constructor throws error for invalid knot vector', () => {
        expect(() => {
            new Knots([2, 1, 3]); // Not non-decreasing
        }).toThrow('Knot vector must be non-decreasing');
    });
});

describe('ProductKnots', () => {
    let productKnots: ProductKnots;

    beforeEach(() => {
        // Initialize a ProductKnots instance with valid knot vectors
        productKnots = new ProductKnots([
            [0, 0, 0, 1, 1],
            [0, 0, 1, 1, 1],
        ]);
    });

    test('getDimension returns correct dimension', () => {
        expect(productKnots.getDimension()).toBe(2);
    });

    test('getKnotSequence returns correct sequence for a specific direction', () => {
        expect(productKnots.getKnotSequence(0)).toEqual([0, 0, 0, 1, 1]);
        expect(productKnots.getKnotSequence(1)).toEqual([0, 0, 1, 1, 1]);
    });

    test('withInsertedKnot inserts a knot correctly', () => {
        const newProductKnots = productKnots.withInsertedKnot(0, 0.5);
        expect(newProductKnots.getKnotSequence(0)).toEqual([0, 0, 0, 0.5, 1, 1]);
    });

    test('withRemovedKnot removes knots correctly', () => {
        const newProductKnots = productKnots.withRemovedKnot(1, 0);
        expect(newProductKnots.getKnotSequence(0)).toEqual([0, 0, 1, 1]);
    });

    test('getDomain returns correct domain for a specific direction', () => {
        const domain = productKnots.getDomain(0);
        expect(domain).toEqual({ min: 0, max: 1 });
    });

    test('getDistinctKnots returns correct distinct knots for a specific direction', () => {
        const distinctKnots = productKnots.getDistinctKnots(0);
        expect(distinctKnots).toEqual([
            { value: 0, multiplicity: 3 },
            { value: 1, multiplicity: 2 },
        ]);
    });

    test('constructor throws error for invalid knot vectors', () => {
        expect(() => {
            new ProductKnots([[1, 0], [0, 1]]); // Invalid knot vector
        }).toThrow('Knot vector 0 must be non-decreasing');
    });
});
