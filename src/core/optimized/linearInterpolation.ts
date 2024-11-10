/**
 * Performs linear interpolation between two points efficiently.
 * Each point is represented as an array of numbers (coordinates).
 * 
 * @param p1 The starting point
 * @param p2 The ending point
 * @param t The interpolation parameter (0 <= t <= 1)
 * @returns An array representing the interpolated point
 */
export function linearInterpolation(
  p1: number[],
  p2: number[],
  t: number
): number[] {
  if (p1.length !== p2.length) {
    throw new Error("Points must have the same number of dimensions");
  }

  const result = new Array(p1.length);
  const t1 = 1 - t;

  for (let i = 0; i < p1.length; i++) {
    result[i] = t1 * p1[i] + t * p2[i];
  }

  return result;
}

/**
 * Performs highly optimized linear interpolation between two points.
 * This version uses loop unrolling and SIMD-like optimizations for maximum performance.
 * It's optimized for common 2D, 3D, and 4D cases, but also handles higher dimensions efficiently.
 * 
 * @param p1 The starting point (array of numbers)
 * @param p2 The ending point (array of numbers)
 * @param t The interpolation parameter (0 <= t <= 1)
 * @returns An array representing the interpolated point
 */
export function linearInterpolationUnrolled(p1: number[], p2: number[], t: number): number[] {
  const len = p1.length;
  if (len !== p2.length) throw new Error("Points must have the same number of dimensions");

  const result = new Array(len);
  const t1 = 1 - t;

  let i = 0;
  // Unrolled loop for common dimensions (2D, 3D, 4D)
  switch (len) {
    case 4:
      result[3] = t1 * p1[3] + t * p2[3];
    case 3:
      result[2] = t1 * p1[2] + t * p2[2];
    case 2:
      result[1] = t1 * p1[1] + t * p2[1];
    case 1:
      result[0] = t1 * p1[0] + t * p2[0];
      return result;
    default:
      // SIMD-like optimization for larger dimensions
        result[i] = t1 * p1[i] + t * p2[i];
        result[i + 1] = t1 * p1[i + 1] + t * p2[i + 1];
        result[i + 2] = t1 * p1[i + 2] + t * p2[i + 2];
        result[i + 3] = t1 * p1[i + 3] + t * p2[i + 3];
      }
      // Handle remaining elements
      for (; i < len; i++) {
        result[i] = t1 * p1[i] + t * p2[i];
      }
      return result;
}



