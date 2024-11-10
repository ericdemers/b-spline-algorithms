/**
 * A TypeScript implementation of Haskell's zipWith function.
 * 
 * This function applies a given function to pairs of elements from two input arrays,
 * producing an array of results.
 * 
 * @template T The type of elements in the first array
 * @template U The type of elements in the second array
 * @template R The type of elements in the resulting array
 * 
 * @param f A function that takes an element from each input array and returns a value
 * @param as The first input array
 * @param bs The second input array
 * @returns An array containing the results of applying f to corresponding elements of as and bs
 * 
 * @precondition The input arrays 'as' and 'bs' must have the same length.
 *               If the arrays have different lengths, the behavior is undefined
 *               and may result in errors or unexpected results.
 * 
 * @example
 * zipWith((a, b) => a + b, [1, 2, 3], [1, 2, 3])
 * // Returns [2, 4, 6]
 * 
 * @example
 * 
 * zipWith((a, b) => a + b, [1, 2, 3], ['a', 'b', 'c'])
 * // Returns ["1a", "2b", "3c"]
 */
export function zipWith<T, U, R>(
    f: (a: T, b: U) => R,
    as: T[],
    bs: U[]
  ): R[] {
    // Use map to iterate over the first array (as)
    // For each element in 'as', we apply the function f to it and the corresponding element in 'bs'
    return as.map((a, i) => f(a, bs[i]));
  }
