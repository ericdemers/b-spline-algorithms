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
 */
export function dividedDifferences(points: Point[]): number[] {
    const n = points.length;
    // Initialize a 2D array to store divided differences
    const f: number[][] = new Array(n).fill(0).map(() => new Array(n).fill(0));

    // Initialize the first column with y values
    // f[i][0] represents f(x_i)
    for (let i = 0; i < n; i++) {
        f[i][0] = points[i][1];
    }

    // Calculate divided differences
    for (let j = 1; j < n; j++) {
        for (let i = 0; i < n - j; i++) {
            // Apply the divided difference formula:
            // f[x_i, ..., x_{i+j}] = (f[x_{i+1}, ..., x_{i+j}] - f[x_i, ..., x_{i+j-1}]) / (x_{i+j} - x_i)
            f[i][j] = (f[i + 1][j - 1] - f[i][j - 1]) / (points[i + j][0] - points[i][0]);
        }
    }

    // Return the first row, which contains the coefficients
    // These coefficients can be used in Newton's interpolation polynomial
    // The function returns an array because it represents the coefficients of the Newton form of the interpolation polynomial:
    // P(x) = a₀ + a₁(x - x₀) + a₂(x - x₀)(x - x₁) + ... + aₙ₋₁(x - x₀)(x - x₁)...(x - xₙ₋₂)
    // Where a₀, a₁, a₂, ..., aₙ₋₁ are the coefficients in the returned array.
    // This array format allows for efficient storage, easy access to individual coefficients, and direct use in constructing the polynomial.
    return f[0];
}

/*
// Example usage
const points: Point[] = [[0, 1], [1, 2], [2, 4], [3, 8]];
const coefficients = dividedDifferences(points);
console.log(coefficients);
*/







