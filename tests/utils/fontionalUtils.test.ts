import { zipWith } from "../../src/utils/fonctionalUtils";

describe('zipWith', () => {
  test('combines two arrays of numbers with addition', () => {
    const result = zipWith((a, b) => a + b, [1, 2, 3], [4, 5, 6]);
    expect(result).toEqual([5, 7, 9]);
  });

  test('combines two arrays of strings with concatenation', () => {
    const result = zipWith((a, b) => a + b, ['a', 'b', 'c'], ['1', '2', '3']);
    expect(result).toEqual(['a1', 'b2', 'c3']);
  });

  test('combines arrays of different types', () => {
    const result = zipWith((a, b) => a.toString() + b, [1, 2, 3], ['a', 'b', 'c']);
    expect(result).toEqual(['1a', '2b', '3c']);
  });

  test('handles empty arrays', () => {
    const result = zipWith((a, b) => a + b, [], []);
    expect(result).toEqual([]);
  });

  test('handles arrays of length 1', () => {
    const result = zipWith((a, b) => a + b, [1], [2]);
    expect(result).toEqual([3]);
  });

  test('uses complex combining function', () => {
    const result = zipWith((a, b) => ({ sum: a + b, product: a * b }), [1, 2, 3], [4, 5, 6]);
    expect(result).toEqual([
      { sum: 5, product: 4 },
      { sum: 7, product: 10 },
      { sum: 9, product: 18 }
    ]);
  });

  test('handles arrays of different lengths (uses shorter length)', () => {
    const result = zipWith((a, b) => a + b, [1, 2, 3], [4, 5]);
    expect(result).toEqual([5, 7, NaN]);
  });

  test('works with arrays of functions', () => {
    const fns1 = [(x: number) => x + 1, (x: number) => x * 2];
    const fns2 = [(x: number) => x - 1, (x: number) => x / 2];
    const result = zipWith((f, g) => (x: number) => f(g(x)), fns1, fns2);
    
    expect(result[0](5)).toBe(5);  // (5 - 1) + 1 = 5
    expect(result[1](6)).toBe(6);  // (6 / 2) * 2 = 6
  });

test('produces NaN when combining function is called with undefined', () => {
  const result = zipWith((a, b) => a + b, [1, 2, 3], [4]);
  expect(result).toEqual([5, NaN, NaN]);
});

});
