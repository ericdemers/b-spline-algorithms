# De Boor's Algorithm

## Overview

De Boor's algorithm is a method for evaluating points on a B-spline or NURBS curve. It's a generalization of the de Casteljau algorithm used for Bézier curves. The algorithm recursively computes a point on the B-spline curve for a given parameter value.

## Implementation

Our implementation in `src/core/functional/deBoor.ts` prioritizes readability and understanding over performance. Here's a breakdown of the algorithm:

### Main Function: `deBoor`

```typescript
export function deBoor(
  controlPoints: Point[],
  knots: number[],
  t: number,
  degree: number
): Point {
  // ...
}
```

This function takes four parameters:

- controlPoints: An array of control points that define the shape of the B-spline curve.

- knots: An array of knot values that determine how the control points affect the curve.

- t: The parameter value at which to evaluate the curve.

- degree: The degree of the B-spline curve.

The function performs these steps:

1. Finds the knot span index for the given parameter t.

2. Extracts the relevant control points.

3. Calls the recursive De Boor function to compute the final point.

```typescript
Recursive Function: deBoorRecursive
function deBoorRecursive(
  points: Point[],
  knots: number[],
  t: number,
  degree: number,
  spanIndex: number,
  r: number
): Point {
  // ...
}
```

This function implements the core of de Boor's algorithm. It recursively computes new sets of points, each set getting closer to the final point on the curve. The recursion continues until we reach the degree of the curve.

Parameters:

- points: The current set of points being processed.

- knots, t, degree: Same as in the main function.

- spanIndex: The index of the knot span containing t.

- r: The current recursion depth.

## Helper Functions

1. findKnotSpanIndex: Determines which knot span contains the parameter t.

2. computeAlpha: Calculates the interpolation factor between two knots.

3. interpolate: Performs linear interpolation between two points.

## Mathematical Explanation

De Boor's algorithm is based on the recursive definition of B-spline basis functions. For a B-spline of degree p, we start with p+1 control points and recursively compute new points using linear interpolation.

The key formula is:

$$ d[i,r](t) = (1 - α) _ d[i-1,r-1](t) + α _ d[i,r-1](t) $$

$d[i,r]$ is the i-th point at recursion level r

$$α = (t - knot[i]) / (knot[i+p-r+1] - knot[i])$$

This process is repeated p times, resulting in the final point on the curve.

Advantages of This Implementation

- Readability : The recursive approach closely mirrors the mathematical definition, making it easier to understand the algorithm's logic.

- Modularity : Helper functions break down the algorithm into clear, manageable steps.

- Educational Value : This implementation serves as an excellent reference for learning and teaching the algorithm.

## Performance Considerations

While this implementation prioritizes clarity, it's not optimized for performance. For high-performance applications, consider using the optimized version in src/core/optimized/deBoor.ts, which uses techniques like in-place computation and loop unrolling to improve efficiency.

## Key optimizations in this version include:

1. In-place computation : Instead of creating new arrays in each iteration, we use a single working array ( points) and modify it in-place.

2. Loop unrolling : We've unrolled the inner loop for the 2D point calculation, avoiding the need for a separate interpolation function.

3. Efficient knot span search : The knot span index finding is optimized with a simple while loop, and it handles the edge case where t is at the end of the knot vector.

4. Minimized function calls : Helper functions have been inlined to reduce function call overhead.

5. Reduced memory allocation : We only allocate one array ( points) instead of creating new arrays in each recursive step.

6. Direct array access : We use direct array indexing instead of array methods like slice().
