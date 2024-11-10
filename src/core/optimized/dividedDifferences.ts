type Point = [number, number]; // [x, y]

/**
 * Calculates the divided differences for a set of points.
 * 
 * Divided differences are used in numerical analysis to describe
 * the coefficients in the Newton form of the interpolation polynomial.
 * 
 * @param points An array of [x, y] coordinate pairs
 * @returns An array of coefficients representing the divided differences
 * @precondition The input array 'points' must not be empty.
 * @precondition The x-coordinates in the input points must be distinct and in ascending order.
 * @throws {Error} If the input array is empty.
 */
export function dividedDifferences(points: Point[]): number[] {
    const n = points.length;
    if (n === 0) {
        throw new Error("Input array must not be empty");
    }

    // Use a single array to store the coefficients
    const coefficients = new Array(n);
    
    // Initialize with y values
    for (let i = 0; i < n; i++) {
        coefficients[i] = points[i][1];
    }

    // Calculate divided differences in-place
    for (let j = 1; j < n; j++) {
        for (let i = n - 1; i >= j; i--) {
            coefficients[i] = (coefficients[i] - coefficients[i - 1]) / (points[i][0] - points[i - j][0]);
        }
    }

    return coefficients;
}

/**
 * Calculates the divided differences for a set of points.
 * 
 * @param points An array of [x, y] coordinate pairs
 * @returns An array of coefficients representing the divided differences
 * @precondition The input array 'points' must not be empty.
 * @precondition The x-coordinates in the input points must be distinct and in ascending order.
 * @throws {Error} If the input array is empty.
 */
export function dividedDifferencesOptimizedForLargeInput(points: Point[]): number[] {
    const n = points.length;
    if (n === 0) throw new Error("Input array must not be empty");

    const coeffs = new Float64Array(n);
    
    // Initialize with y values
    for (let i = 0; i < n; i++) {
        coeffs[i] = points[i][1];
    }

    // Calculate divided differences in-place
    for (let j = 1; j < n; j++) {
        for (let i = n - 1; i >= j; i--) {
            coeffs[i] = (coeffs[i] - coeffs[i - 1]) / (points[i][0] - points[i - j][0]);
        }
    }

    return Array.from(coeffs);
}
