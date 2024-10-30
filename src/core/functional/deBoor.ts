// src/core/functional/deBoor.ts

type Point = [number, number];

/**
 * Computes a point on a B-spline curve using De Boor's algorithm.
 * This is a functional, easy-to-read implementation prioritizing clarity over performance.
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
  const spanIndex = findKnotSpanIndex(knots, t, degree);

  // Extract the relevant control points
  const relevantControlPoints = controlPoints.slice(spanIndex - degree, spanIndex + 1);

  // Compute the point using recursive De Boor's algorithm
  return deBoorRecursive(relevantControlPoints, knots, t, degree, spanIndex, 0);
}

/**
 * Recursive implementation of de Boor's algorithm.
 */
function deBoorRecursive(
  points: Point[],
  knots: number[],
  t: number,
  degree: number,
  spanIndex: number,
  r: number
): Point {
  if (r === degree) {
    return points[0];
  }

  const newPoints: Point[] = [];

  for (let i = 0; i < points.length - 1; i++) {
    const alpha = computeAlpha(t, knots[spanIndex + i - degree + r + 1], knots[spanIndex + i + 1]);
    newPoints.push(interpolate(points[i], points[i + 1], alpha));
  }

  return deBoorRecursive(newPoints, knots, t, degree, spanIndex, r + 1);
}

/**
 * Finds the knot span index for a given parameter t.
 */
function findKnotSpanIndex(knots: number[], t: number, degree: number): number {
  for (let i = degree; i < knots.length - 1; i++) {
    if (t >= knots[i] && t < knots[i + 1]) {
      return i;
    }
  }
  return knots.length - degree - 2;
}

/**
 * Computes the alpha value for De Boor's algorithm.
 */
function computeAlpha(t: number, knotI: number, knotIPlus1: number): number {
  return (t - knotI) / (knotIPlus1 - knotI);
}

/**
 * Linearly interpolates between two points.
 */
function interpolate(p1: Point, p2: Point, alpha: number): Point {
  return [
    (1 - alpha) * p1[0] + alpha * p2[0],
    (1 - alpha) * p1[1] + alpha * p2[1]
  ];
}
