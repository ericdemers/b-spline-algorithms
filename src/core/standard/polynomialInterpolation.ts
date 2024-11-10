/**
 * Represents a point with x and y coordinates.
 */
type Point = [number, number];

/**
 * Calculates the Lagrange basis polynomial for a given index.
 * 
 * The Lagrange basis polynomial is defined as:
 * Lᵢ(x) = ∏ (x - xⱼ) / (xᵢ - xⱼ), for j ≠ i
 * 
 * Where:
 * - i is the index of the current point
 * - j iterates over all other points
 * - xᵢ is the x-coordinate of the i-th point
 * - xⱼ is the x-coordinate of the j-th point
 * 
 * This function has the property that:
 * Lᵢ(xᵢ) = 1 and Lᵢ(xⱼ) = 0 for all j ≠ i
 * 
 * @param points - Array of points used for interpolation
 * @param index - Index of the current point
 * @param x - The x-value at which to evaluate the basis polynomial
 * @returns The value of the Lagrange basis polynomial at x
 */
function lagrangeBasis(points: Point[], index: number, x: number): number {
  return points.reduce((basis, [xi, _], i) => {
    if (i === index) return basis;
    return basis * (x - xi) / (points[index][0] - xi);
  }, 1);
}

/**
 * Performs polynomial interpolation using the Lagrange method.
 * 
 * The Lagrange interpolation polynomial is defined as:
 * P(x) = ∑ yᵢ * Lᵢ(x)
 * 
 * Where:
 * - yᵢ is the y-coordinate of the i-th point
 * - Lᵢ(x) is the Lagrange basis polynomial for the i-th point
 * 
 * This polynomial has the property that P(xᵢ) = yᵢ for all input points (xᵢ, yᵢ).
 * 
 * @param points - Array of points to interpolate
 * @param x - The x-value at which to interpolate
 * @returns The interpolated y-value
 */
export function interpolate(points: Point[], x: number): number {
  return points.reduce((sum, [_, yi], i) => {
    return sum + yi * lagrangeBasis(points, i, x);
  }, 0);
}
