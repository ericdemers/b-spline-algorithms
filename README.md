# B-spline Algorithms

A comprehensive TypeScript library for working with B-splines of arbitrary dimension in arbitrary vector spaces. Built with mathematical rigor and computational efficiency in mind, this library provides a flexible foundation for B-spline operations and algorithms.

## Core Features

### Fundamental B-spline Operations

Each operation supports multiple algorithm implementations that can be selected based on your needs:

#### Evaluation Algorithms

- De Boor's algorithm (stable and general purpose)
- Optimized evaluation for uniform B-splines
- Educational step-by-step implementation

```typescript
const curve = new BSplineCurve2D(controlPoints, degree, {
  evaluation: new DeBoorEvaluation(), // Default stable algorithm
  // OR
  evaluation: new OptimizedEvaluation(), // Performance-focused
  // OR
  evaluation: new EducationalEvaluation(), // Clear step-by-step implementation
});
```

### Knot Operations

- Boehm's knot insertion algorithm
- Oslo algorithm for knot refinement
- Knot removal with error control

```typescript
const refinedCurve = curve.insertKnot(0.5, {
    algorithm: new BoehmsInsertion()     // Classical algorithm
    // OR
    algorithm: new OsloAlgorithm()       // Better for multiple insertions
    // OR
    algorithm: new AdaptiveInsertion()   // Automatic error control
});
```

### Degree Manipulation

- Degree elevation with various bases
- Degree reduction with error control

```typescript
const elevatedCurve = curve.elevateDegree({
  method: new BernsteinElevation(), // Using Bernstein basis
  // OR
  method: new PowerBasisElevation(), // Using power basis
  // OR
  method: new AdaptiveElevation(), // Automatic error control
});
```

### Derivative Computation

- Symbolic differentiation
- Automatic differentiation
- Finite difference methods

```typescript
const derivative = curve.derivative({
  method: new SymbolicDerivative(), // Exact derivatives
  // OR
  method: new AutomaticDerivative(), // Efficient for higher orders
  // OR
  method: new FiniteDifference(), // Numerical approximation
});
```

### Multiplication Algorithms

- Direct multiplication
- FFT-based multiplication for uniform B-splines
- Adaptive multiplication with error control

```typescript
const product = f.multiply(g, {
  algorithm: new DirectMultiplication(), // Standard algorithm
  // OR
  algorithm: new FFTMultiplication(), // Fast for uniform B-splines
  // OR
  algorithm: new AdaptiveMultiplication(), // Automatic precision control
});
```

### Blossom (Polar Form) Algorithms

- Classical polar form computation
- Optimized evaluation for specific cases

```typescript
const blossomValue = curve.blossom([u1, u2, u3], {
  method: new ClassicalBlossom(), // Standard implementation
  // OR
  method: new OptimizedBlossom(), // Performance-focused
  // OR
  method: new SymbolicBlossom(), // Symbolic computation
});
```

### Flexible Architecture

- Support for B-splines of any dimension (functions, curves, surfaces, volumes)
- Works with arbitrary vector spaces (ℝⁿ, ℂⁿ, or custom spaces)
- Extensible design for adding new algorithms and spaces

### Implementation Choices

- Type-safe implementation in TypeScript
- Both readable "educational" and optimized implementations
- Comprehensive test coverage
- Well-documented mathematical foundations

## Design Philosophy

### Flexible Algorithm Selection

The library supports interchangeable algorithms for core operations:

```typescript
// Choose between different algorithm implementations
const curve = new BSplineCurve2D(controlPoints, degree, {
  evaluation: new OptimizedEvaluation(), // Fast evaluation
  knotInsertion: new BoehmsInsertion(), // Boehm's algorithm
  multiplication: new FastMultiplication(), // Optimized multiplication
});

// Or use convenient factory methods
const educationalCurve = BSplineFactory.createEducational(
  controlPoints,
  degree
);
const optimizedCurve = BSplineFactory.createOptimized(controlPoints, degree);
```

## Simple Example

### Easy-to-Use Curve Creation

```typescript
// Create a cubic B-spline curve in 2D with automatic knot vector generation
const curve = new BSplineCurve2D(
  [
    [0, 0], // Control points as [x, y] coordinates
    [1, 1],
    [2, 0],
    [1, -1],
    [0, 0],
  ],
  3
); // Degree 3 (cubic)

// Evaluate the curve at parameter value
const point = curve.evaluate(0.5); // Returns [x, y] point

// Get curve derivatives
const tangent = curve.derivative(1); // First derivative
const curvature = curve.derivative(2); // Second derivative

// Modify the curve
const refinedCurve = curve.insertKnot(0.3); // Knot insertion
const elevatedCurve = curve.elevateDegree(); // Degree elevation
```

### B-Spline Function Multiplication

```typescript
// Create two B-spline functions (1D B-splines)
const f = new BSplineFunction(
  [
    1, // Control points as scalar values
    2,
    0.5,
    1,
  ],
  2
); // Quadratic B-spline

const g = new BSplineFunction(
  [
    0.5, // Control points as scalar values
    1,
    1.5,
  ],
  1
); // Linear B-spline

// Multiply the functions
const product = f.multiply(g);

// Evaluate the product at a parameter value
const value = product.evaluate(0.5); // Returns scalar value

// The resulting B-spline maintains mathematical properties:
console.log(product.degree); // Sum of input degrees (2 + 1 = 3)
console.log(product.domain); // Intersection of input domains
console.log(product.controlPoints.length); // New control point count

// Functions can be manipulated before/after multiplication
const refined = f.insertKnot(0.5).multiply(g.elevateDegree());
```

### Create a cubic B-spline curve in 3D space

```typescript
// Create a cubic B-spline curve in 3D space
const curve = new BSpline(
  vectorSpace3D, // Vector space definition
  controlPoints, // Control points in 3D
  knotVector, // Knot vector
  3 // Degree
);

// Evaluate the curve
const point = curve.evaluate(0.5);

// Compute the derivative
const derivative = curve.derivative();

// Insert a knot
const refinedCurve = curve.insertKnot(0.3);
```

## Installation

```bash
npm install b-spline-algorithms
```

## Current Status

This library is under active development. Current focus:

- Core B-spline evaluation algorithms
- Basic operations (knot insertion, degree elevation)
- Foundation for arbitrary vector spaces
- Documentation of mathematical concepts

Coming soon:

- Additional vector space implementations
- Advanced algorithms (subdivision, intersection)
- Performance optimizations
- More examples and tutorials

## Contributing

Contributions are welcome! Please read our Contributing Guide for details on our code of conduct and the process for submitting pull requests.

### Development Setup

1. Fork the repository

2. Clone your fork

3. Install dependencies

4. Create a feature branch

5. Make your changes

6. Run tests

7. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

List any references, papers, or other libraries that inspired or informed your implementations

## Contact

- Create an issue on our GitHub repository for bug reports or feature requests

- Submit pull requests for contributions

## Version History

- 1.0.0
  - Initial release
  - Basic B-spline functionality
  - Periodic curve support
