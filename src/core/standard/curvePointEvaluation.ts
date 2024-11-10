/*
function evaluateBSplinePoint(
  t: number,
  degree: number,
  knots: number[],
  controlPoints: number[][]
): number[] {
  return deBoorRecursive(t, degree, 0, controlPoints.length - 1, knots, controlPoints);
}

function deBoorRecursive(
  t: number,
  degree: number,
  i: number,
  k: number,
  knots: number[],
  controlPoints: number[][]
): number[] {
  // Base case: degree 0
  if (degree === 0) {
    return controlPoints[i];
  }

  // Recursive case
  const alpha = (t - knots[i]) / (knots[i + degree] - knots[i]);
  
  const left = deBoorRecursive(t, degree - 1, i, k, knots, controlPoints);
  const right = deBoorRecursive(t, degree - 1, i + 1, k, knots, controlPoints);

  return left.map((leftValue, index) => 
    (1 - alpha) * leftValue + alpha * right[index]
  );
}
  */
