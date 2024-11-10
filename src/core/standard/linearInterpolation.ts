import { zipWith } from "../../utils/fonctionalUtils";

/**
 * Performs linear interpolation between two points.
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
  return zipWith((a, b) => interpolate(t, a, b), p1, p2);
}

/**
 * Interpolates between two numbers.
 * 
 * @param t The interpolation parameter (0 <= t <= 1)
 * @param a The starting number
 * @param b The ending number
 * @returns The interpolated value
 */
function interpolate(t: number, a: number, b: number): number {
  return (1 - t) * a + t * b;
}
