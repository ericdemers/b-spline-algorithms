import { DistinctKnotsWithMultiplicities } from '../../src/core/knot-structure';
import { PeriodicKnots } from '../../src/core/periodic-knots';

describe('PeriodicKnots', () => {
  let periodicKnots: PeriodicKnots;

  beforeEach(() => {
    // Initialize a PeriodicKnots instance with a basic pattern and period
    periodicKnots = new PeriodicKnots([0, 1, 2], 3); // Example pattern and period
  });

  test('constructor initializes correctly', () => {
    expect(periodicKnots.getPattern()).toEqual([0, 1, 2]);
    expect(periodicKnots.getPeriod()).toBe(3);
  });

  test('getDimension returns correct dimension', () => {
    expect(periodicKnots.getDimension()).toBe(1);
  });

  test('getKnotSequence returns correct unrolled knots', () => {
    const sequence = periodicKnots.getKnotSequence(); 
    expect(sequence).toEqual([0, 1, 2, 3, 4, 5, 6, 7]); 
  });

  
  test('withInsertedKnot inserts a knot correctly', () => {
    const newKnots = periodicKnots.withInsertedKnot(1.5); // Insert a knot at 1.5
    expect(newKnots.getPattern()).toEqual([0, 1, 1.5, 2]); // New pattern should include the inserted knot
  });

  test('withRemovedKnots removes knots correctly', () => {
    const newKnots = periodicKnots.withRemovedKnot(1); // Remove knot at index 1
    expect(newKnots.getPattern()).toEqual([0, 2]); 
  });

  
  test('getDomain returns correct domain', () => {
    const domain = periodicKnots.getDomain(0);
    expect(domain).toEqual({ min: 0, max: 3 }); // Domain based on pattern and period
  });

  test('getDistinctKnots returns cached distinct knots', () => {
    const distinctKnots = periodicKnots.getDistinctKnots();
    expect(distinctKnots).toEqual([{value: 0, multiplicity: 1}, {value: 1, multiplicity: 1}, {value: 2, multiplicity: 1}]); // Assuming distinct knots are the same as the pattern
  });

  test('fromDistinctKnots creates PeriodicKnots correctly', () => {
    const distinctKnots: DistinctKnotsWithMultiplicities = {
      knots: [0, 1, 2],
      multiplicities: [1, 1, 1],
    };
    const newPeriodicKnots = PeriodicKnots.fromDistinctKnots(distinctKnots, 3);
    expect(newPeriodicKnots.getPattern()).toEqual([0, 1, 2]);
    expect(newPeriodicKnots.getPeriod()).toBe(3);
  });

  test('constructor throws error for invalid pattern', () => {
    expect(() => {
      new PeriodicKnots([1], 3); // Invalid pattern with less than 2 values
    }).toThrow('Pattern must have at least 2 values');
  });

  test('constructor throws error for negative period', () => {
    expect(() => {
      new PeriodicKnots([0, 1, 2], -1); // Invalid negative period
    }).toThrow('Period must be positive');
  });

  test('findKnotSpan returns correct span for a parameter', () => {
    const span = periodicKnots.findKnotSpan(1.5);
    expect(span).toBe(1); // Should return the index of the knot span containing 1.5
  });
  
  
  
});
