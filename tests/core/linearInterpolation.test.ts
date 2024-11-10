import { linearInterpolation as linearInterpolationFunctional} from '../../src/core/standard/linearInterpolation';
import fc from 'fast-check';

describe('linearInterpolation', () => {
  test('interpolates correctly between two 2D points', () => {
    const p1 = [0, 0];
    const p2 = [10, 20];
    
    expect(linearInterpolationFunctional(p1, p2, 0)).toEqual([0, 0]);
    expect(linearInterpolationFunctional(p1, p2, 1)).toEqual([10, 20]);
    expect(linearInterpolationFunctional(p1, p2, 0.5)).toEqual([5, 10]);
    expect(linearInterpolationFunctional(p1, p2, 0.25)).toEqual([2.5, 5]);
    expect(linearInterpolationFunctional(p1, p2, 0.75)).toEqual([7.5, 15]);
  });

  test('interpolates correctly between two 3D points', () => {
    const p1 = [0, 0, 0];
    const p2 = [10, 20, 30];
    
    expect(linearInterpolationFunctional(p1, p2, 0)).toEqual([0, 0, 0]);
    expect(linearInterpolationFunctional(p1, p2, 1)).toEqual([10, 20, 30]);
    expect(linearInterpolationFunctional(p1, p2, 0.5)).toEqual([5, 10, 15]);
  });

  test('handles edge case of identical points', () => {
    const p = [5, 5];
    
    expect(linearInterpolationFunctional(p, p, 0)).toEqual(p);
    expect(linearInterpolationFunctional(p, p, 0.5)).toEqual(p);
    expect(linearInterpolationFunctional(p, p, 1)).toEqual(p);
  });

  test('throws error for points with different dimensions', () => {
    const p1 = [0, 0];
    const p2 = [10, 20, 30];
    
    expect(() => linearInterpolationFunctional(p1, p2, 0.5)).toThrow('Points must have the same number of dimensions');
  });

  test('handles t values outside [0, 1] range', () => {
    const p1 = [0, 0];
    const p2 = [10, 20];
    
    expect(linearInterpolationFunctional(p1, p2, -0.5)).toEqual([-5, -10]);
    expect(linearInterpolationFunctional(p1, p2, 1.5)).toEqual([15, 30]);
  });
});

import { linearInterpolation, linearInterpolationUnrolled } from '../../src/core/optimized/linearInterpolation';

describe('Linear Interpolation Functions', () => {
  // Helper function to compare arrays with a small epsilon for floating point precision
  const arraysEqual = (arr1: number[], arr2: number[], epsilon = 1e-6) => {
    if (arr1.length !== arr2.length) return false;
    return arr1.every((val, i) => Math.abs(val - arr2[i]) < epsilon);
  };

  // Test cases
  const testCases = [
    { p1: [0], p2: [10], t: 0.5, expected: [5] },
    { p1: [0, 0], p2: [10, 20], t: 0.5, expected: [5, 10] },
    { p1: [0, 0, 0], p2: [10, 20, 30], t: 0.5, expected: [5, 10, 15] },
    { p1: [0, 0, 0, 0], p2: [10, 20, 30, 40], t: 0.5, expected: [5, 10, 15, 20] },
    { p1: [0, 0, 0, 0, 0], p2: [10, 20, 30, 40, 50], t: 0.5, expected: [5, 10, 15, 20, 25] },
    { p1: [1, 2, 3], p2: [4, 5, 6], t: 0.25, expected: [1.75, 2.75, 3.75] },
    { p1: [0, 0], p2: [10, 10], t: 0, expected: [0, 0] },
    { p1: [0, 0], p2: [10, 10], t: 1, expected: [10, 10] },
  ];

  describe('linearInterpolation', () => {
    test.each(testCases)('interpolates correctly for p1=$p1, p2=$p2, t=$t', ({ p1, p2, t, expected }) => {
      const result = linearInterpolation(p1, p2, t);
      expect(arraysEqual(result, expected)).toBe(true);
    });

    test('throws error when dimensions do not match', () => {
      expect(() => linearInterpolation([0, 0], [10], 0.5)).toThrow('Points must have the same number of dimensions');
    });
  });

  describe('linearInterpolationUnrolled', () => {
    test.each(testCases)('interpolates correctly for p1=$p1, p2=$p2, t=$t', ({ p1, p2, t, expected }) => {
      const result = linearInterpolationUnrolled(p1, p2, t);
      expect(arraysEqual(result, expected)).toBe(true);
    });

    test('throws error when dimensions do not match', () => {
      expect(() => linearInterpolationUnrolled([0, 0], [10], 0.5)).toThrow('Points must have the same number of dimensions');
    });
  });

  describe('Comparison of linearInterpolation and linearInterpolationUnrolled', () => {
    test.each(testCases)('both functions produce the same result for p1=$p1, p2=$p2, t=$t', ({ p1, p2, t }) => {
      const result1 = linearInterpolation(p1, p2, t);
      const result2 = linearInterpolationUnrolled(p1, p2, t);
      expect(arraysEqual(result1, result2)).toBe(true);
    });
  });

  describe('Edge cases', () => {
    test('handles very large numbers', () => {
      const p1 = [Number.MAX_SAFE_INTEGER, Number.MAX_SAFE_INTEGER];
      const p2 = [Number.MAX_SAFE_INTEGER * 2, Number.MAX_SAFE_INTEGER * 2];
      const t = 0.5;
      const expected = [Number.MAX_SAFE_INTEGER * 1.5, Number.MAX_SAFE_INTEGER * 1.5];
      
      const result1 = linearInterpolation(p1, p2, t);
      const result2 = linearInterpolationUnrolled(p1, p2, t);
      
      expect(arraysEqual(result1, expected)).toBe(true);
      expect(arraysEqual(result2, expected)).toBe(true);
    });

    test('handles very small numbers', () => {
      const p1 = [Number.MIN_VALUE, Number.MIN_VALUE];
      const p2 = [Number.MIN_VALUE * 2, Number.MIN_VALUE * 2];
      const t = 0.5;
      const expected = [Number.MIN_VALUE * 1.5, Number.MIN_VALUE * 1.5];
      
      const result1 = linearInterpolation(p1, p2, t);
      const result2 = linearInterpolationUnrolled(p1, p2, t);
      
      expect(arraysEqual(result1, expected)).toBe(true);
      expect(arraysEqual(result2, expected)).toBe(true);
    });
  });
});

describe('Linear Interpolation Functions - Consistency Tests', () => {
    // Helper function to compare arrays with a small epsilon for floating point precision
    const arraysEqual = (arr1: number[], arr2: number[], epsilon = 1e-6) => {
      if (arr1.length !== arr2.length) return false;
      return arr1.every((val, i) => Math.abs(val - arr2[i]) < epsilon);
    };
  
    // Test cases
    const testCases = [
      { p1: [0], p2: [10], t: 0.5 },
      { p1: [0, 0], p2: [10, 20], t: 0.5 },
      { p1: [0, 0, 0], p2: [10, 20, 30], t: 0.5 },
      { p1: [1, 2, 3], p2: [4, 5, 6], t: 0.25 },
      { p1: [0, 0], p2: [10, 10], t: 0 },
      { p1: [0, 0], p2: [10, 10], t: 1 },
      { p1: [-5, -10], p2: [5, 10], t: 0.75 },
    ];
  
    test.each(testCases)('all functions produce the same result for p1=$p1, p2=$p2, t=$t', ({ p1, p2, t }) => {
      const result1 = linearInterpolationFunctional(p1, p2, t);
      const result2 = linearInterpolation(p1, p2, t);
      const result3 = linearInterpolationUnrolled(p1, p2, t);
  
      expect(arraysEqual(result1, result2)).toBe(true);
      expect(arraysEqual(result1, result3)).toBe(true);
      expect(arraysEqual(result2, result3)).toBe(true);
    });
  
    test('all functions produce the same result for random inputs', () => {
      fc.assert(
        fc.property(
          fc.array(fc.float({noNaN: true, noDefaultInfinity: true}), {minLength: 1, maxLength: 10}),
          fc.array(fc.float({noNaN: true, noDefaultInfinity: true}), {minLength: 1, maxLength: 10}),
          fc.float({noNaN: true, noDefaultInfinity: true, min: 0, max: 1}),
          (p1, p2, t) => {
            if (p1.length !== p2.length) return true; // Skip if lengths don't match
  
            const result1 = linearInterpolationFunctional(p1, p2, t);
            const result2 = linearInterpolation(p1, p2, t);
            const result3 = linearInterpolationUnrolled(p1, p2, t);
  
            return arraysEqual(result1, result2) && arraysEqual(result1, result3) && arraysEqual(result2, result3);
          }
        )
      );
    });
  
    test('all functions handle edge cases consistently', () => {
      const edgeCases = [
        { p1: [Number.MAX_SAFE_INTEGER], p2: [Number.MAX_SAFE_INTEGER * 2], t: 0.5 },
        { p1: [Number.MIN_VALUE], p2: [Number.MIN_VALUE * 2], t: 0.5 },
        { p1: [-Number.MAX_SAFE_INTEGER], p2: [-Number.MAX_SAFE_INTEGER * 2], t: 0.5 },
      ];
  
      edgeCases.forEach(({ p1, p2, t }) => {
        const result1 = linearInterpolationFunctional(p1, p2, t);
        const result2 = linearInterpolation(p1, p2, t);
        const result3 = linearInterpolationUnrolled(p1, p2, t);
  
        expect(arraysEqual(result1, result2)).toBe(true);
        expect(arraysEqual(result1, result3)).toBe(true);
        expect(arraysEqual(result2, result3)).toBe(true);
      });
    });
  });

