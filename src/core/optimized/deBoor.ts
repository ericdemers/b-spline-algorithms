// src/core/optimized/deBoor.ts
// Note: This is an optimized version prioritizing performance.
// For a more readable implementation, see src/core/functional/deBoor.ts

type Point = [number, number];

/**
 * Computes a point on a B-spline curve using De Boor's algorithm.
 * This is an optimized implementation prioritizing performance.
 *
 * @param controlPoints - Array of control points defining the B-spline curve
 * @param knots - Array of knot values
 * @param t - Parameter value at which to evaluate the B-spline
 * @param degree - Degree of the B-spline curve
 * @returns The point on the B-spline curve at parameter t
 */
export function deBoor(
  controlPoints: Point[],
  knots: number[],
  t: number,
  degree: number
): Point {
  // Find the knot span index
  let spanIndex = degree;
  const n = knots.length - degree - 1;
  
  if (t === knots[n]) {
    spanIndex = n - 1;
  } else {
    while (t >= knots[spanIndex + 1]) spanIndex++;
  }

  // Create a working array of points
  const points: Point[] = new Array(degree + 1);
  for (let i = 0; i <= degree; i++) {
    points[i] = [...controlPoints[spanIndex - degree + i]];
  }

  // Compute the point using optimized De Boor's algorithm
  for (let r = 1; r <= degree; r++) {
    for (let i = degree; i >= r; i--) {
      const alpha = (t - knots[spanIndex - degree + i]) / 
                    (knots[spanIndex + i + 1 - r] - knots[spanIndex - degree + i]);
      points[i][0] = (1 - alpha) * points[i - 1][0] + alpha * points[i][0];
      points[i][1] = (1 - alpha) * points[i - 1][1] + alpha * points[i][1];
    }
  }

  return points[degree];
}
