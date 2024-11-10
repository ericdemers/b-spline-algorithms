import {dividedDifferences as dividedDifferencesFunctional} from '../../src/core/standard/dividedDifferences';
import {dividedDifferences as dividedDifferencesOptimized, dividedDifferencesOptimizedForLargeInput as dividedDifferencesOptimizedForLargeInput} from '../../src/core/optimized/dividedDifferences';

const implementations = [
    { name: 'Functional dividedDifferences', func: dividedDifferencesFunctional },
    //{ name: 'More Functional dividedDifferences', func: dividedDifferencesMoreFunctional },
    { name: 'Optimized dividedDifferences', func: dividedDifferencesOptimized },
    { name: 'Optimized for large input dividedDifferences', func: dividedDifferencesOptimizedForLargeInput },
  ];

implementations.forEach(({ name, func }) => {

    describe(`${name}`, () => {
    test('calculates divided differences for a linear function', () => {
        const points: [number, number][] = [[0, 0], [1, 1], [2, 2], [3, 3]];
        const result = func(points);
        expect(result).toHaveLength(4);
        expect(result[0]).toBeCloseTo(0);
        expect(result[1]).toBeCloseTo(1);
        expect(result[2]).toBeCloseTo(0);
        expect(result[3]).toBeCloseTo(0);
    });

    test('calculates divided differences for a quadratic function', () => {
        const points: [number, number][] = [[0, 0], [1, 1], [2, 4], [3, 9]];
        const result = func(points);
        expect(result).toHaveLength(4);
        expect(result[0]).toBeCloseTo(0);
        expect(result[1]).toBeCloseTo(1);
        expect(result[2]).toBeCloseTo(1);
        expect(result[3]).toBeCloseTo(0);
    });

    test('handles a single point correctly', () => {
        const points: [number, number][] = [[1, 5]];
        const result = func(points);
        expect(result).toHaveLength(1);
        expect(result[0]).toBeCloseTo(5);
    });

    test('handles two points correctly', () => {
        const points: [number, number][] = [[1, 2], [3, 4]];
        const result = func(points);
        expect(result).toHaveLength(2);
        expect(result[0]).toBeCloseTo(2);
        expect(result[1]).toBeCloseTo(1);
    });

    test('calculates divided differences for non-uniform x values', () => {
        const points: [number, number][] = [[1, 1], [2, 8], [5, 125]];
        const result = func(points);
        expect(result).toHaveLength(3);
        expect(result[0]).toBeCloseTo(1);
        expect(result[1]).toBeCloseTo(7);
        expect(result[2]).toBeCloseTo(8);
    });

    });

});


implementations.forEach(({ name, func }) => {
    describe(`${name}`, () => {  
      test('handles negative x and y values correctly', () => {
        const points: [number, number][] = [[-2, 4], [-1, 1], [0, 0], [1, 1], [2, 4]];
        const result = func(points);
        expect(result).toHaveLength(5);
        expect(result[0]).toBeCloseTo(4);
        expect(result[1]).toBeCloseTo(-3);
        expect(result[2]).toBeCloseTo(1);
        expect(result[3]).toBeCloseTo(0);
        expect(result[4]).toBeCloseTo(0);
      });
  
      test('calculates divided differences for exponential function', () => {
        const points: [number, number][] = [[0, 1], [1, Math.E], [2, Math.E ** 2]];
        const result = func(points);
        expect(result).toHaveLength(3);
        expect(result[0]).toBeCloseTo(1);
        expect(result[1]).toBeCloseTo(Math.E - 1);
        expect(result[2]).toBeCloseTo((Math.E ** 2 - 2 * Math.E + 1) / 2);
      });
  
      test('handles very large x differences', () => {
        const points: [number, number][] = [[0, 1], [1000000, 1000001]];
        const result = func(points);
        expect(result).toHaveLength(2);
        expect(result[0]).toBeCloseTo(1);
        expect(result[1]).toBeCloseTo(1, 6);
      });
  
      test('handles very small x differences', () => {
        const points: [number, number][] = [[0, 1], [1e-10, 1.000000001]];
        const result = func(points);
        expect(result).toHaveLength(2);
        expect(result[0]).toBeCloseTo(1);
        expect(result[1]).toBeCloseTo(10);
      });
  
      test('calculates divided differences for trigonometric function', () => {
        const points: [number, number][] = [
          [0, 0],
          [Math.PI / 6, 0.5],
          [Math.PI / 4, Math.sqrt(2) / 2],
          [Math.PI / 3, Math.sqrt(3) / 2]
        ];
        const result = func(points);
        expect(result).toHaveLength(4);
        expect(result[0]).toBeCloseTo(0);
        expect(result[1]).toBeCloseTo(3 / Math.PI);
        // The following expectations might need adjustment based on precision
        expect(result[2]).toBeCloseTo(-0.20860760161962227, 5);
        expect(result[3]).toBeCloseTo(-0.13648909830897232, 5);
      });
  
      test('handles repeated y values correctly', () => {
        const points: [number, number][] = [[1, 5], [2, 5], [3, 5], [4, 5]];
        const result = func(points);
        expect(result).toHaveLength(4);
        expect(result[0]).toBeCloseTo(5);
        expect(result[1]).toBeCloseTo(0);
        expect(result[2]).toBeCloseTo(0);
        expect(result[3]).toBeCloseTo(0);
      });
  
  
      test('handles very large number of points', () => {
        const points: [number, number][] = Array.from({ length: 1000 }, (_, i) => [i, i * i]);
        const result = func(points);
        expect(result).toHaveLength(1000);
        expect(result[0]).toBeCloseTo(0);
        expect(result[1]).toBeCloseTo(1);
        expect(result[2]).toBeCloseTo(1);
        // All higher-order differences should be very close to 0
        for (let i = 3; i < 1000; i++) {
          expect(result[i]).toBeCloseTo(0, 5);
        }
      });
    });
  });