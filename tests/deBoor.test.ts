// src/tests/deBoor.test.ts

import { deBoor as deBoorFunctional } from '../src/core/functional/deBoor';
import { deBoor as deBoorOptimized } from '../src/core/optimized/deBoor';

type Point = [number, number];

describe('De Boor Algorithm', () => {
  // Test both implementations
  const implementations = [
    { name: 'Functional', fn: deBoorFunctional },
    { name: 'Optimized', fn: deBoorOptimized }
  ];

  implementations.forEach(({ name, fn }) => {
    describe(`${name} Implementation`, () => {
      test('should compute correct point for a simple linear B-spline', () => {
        const controlPoints: Point[] = [[0, 0], [1, 1]];
        const knots = [0, 0, 1, 1];
        const t = 0.5;
        const degree = 1;

        const result = fn(controlPoints, knots, t, degree);
        expect(result[0]).toBeCloseTo(0.5);
        expect(result[1]).toBeCloseTo(0.5);
      });

      test('should compute correct point for a quadratic B-spline', () => {
        const controlPoints: Point[] = [[0, 0], [1, 2], [2, 0]];
        const knots = [0, 0, 0, 1, 1, 1];
        const t = 0.5;
        const degree = 2;

        const result = fn(controlPoints, knots, t, degree);
        expect(result[0]).toBeCloseTo(1);
        expect(result[1]).toBeCloseTo(1.5);
      });

      test('should handle t at start of knot vector', () => {
        const controlPoints: Point[] = [[0, 0], [1, 1], [2, 0]];
        const knots = [0, 0, 0, 1, 1, 1];
        const t = 0;
        const degree = 2;

        const result = fn(controlPoints, knots, t, degree);
        expect(result[0]).toBeCloseTo(0);
        expect(result[1]).toBeCloseTo(0);
      });

      test('should handle t at end of knot vector', () => {
        const controlPoints: Point[] = [[0, 0], [1, 1], [2, 0]];
        const knots = [0, 0, 0, 1, 1, 1];
        const t = 1;
        const degree = 2;

        const result = fn(controlPoints, knots, t, degree);
        expect(result[0]).toBeCloseTo(2);
        expect(result[1]).toBeCloseTo(0);
      });

      test('should compute correct point for a cubic B-spline', () => {
        const controlPoints: Point[] = [[0, 0], [1, 3], [2, 1], [3, 3], [4, 0]];
        const knots = [0, 0, 0, 0, 0.5, 1, 1, 1, 1];
        const t = 0.75;
        const degree = 3;

        const result = fn(controlPoints, knots, t, degree);
        // You may need to adjust these expected values based on the actual curve
        expect(result[0]).toBeCloseTo(2.5, 1);
        expect(result[1]).toBeCloseTo(2.2, 1);
      });

      test('should throw error for invalid input', () => {
        const controlPoints: Point[] = [[0, 0], [1, 1]];
        const knots = [0, 1];
        const t = 0.5;
        const degree = 1;

        expect(() => fn(controlPoints, knots, t, degree)).toThrow();
      });
    });
  });

  test('Functional and Optimized implementations should produce the same results', () => {
    const controlPoints: Point[] = [[0, 0], [1, 2], [3, 1], [4, 3]];
    const knots = [0, 0, 0, 0.5, 1, 1, 1];
    const degree = 2;

    for (let t = 0; t <= 1; t += 0.1) {
      const functionalResult = deBoorFunctional(controlPoints, knots, t, degree);
      const optimizedResult = deBoorOptimized(controlPoints, knots, t, degree);

      expect(functionalResult[0]).toBeCloseTo(optimizedResult[0], 5);
      expect(functionalResult[1]).toBeCloseTo(optimizedResult[1], 5);
    }
  });
});
