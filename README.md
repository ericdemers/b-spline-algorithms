# b-spline-algorithms

A comprehensive TypeScript library for B-spline and NURBS algorithms, offering both functional and optimized implementations.

## Features

- Implements core B-spline algorithms including De Boor's algorithm, knot insertion, and more
- Provides both functional (readable) and optimized (performant) versions of each algorithm
- Comprehensive test suite ensuring correctness and performance
- Written in TypeScript with full type definitions
- Thoroughly documented with explanations of underlying mathematical concepts

## Installation

```bash
npm install b-spline-algorithms
```

## Usage

Basic usage example:

```typescript
import { deBoor } from "b-spline-algorithms";

const controlPoints = [
  [0, 0],
  [1, 1],
  [2, 0],
  [3, 1],
];
const knots = [0, 0, 0, 0, 1, 1, 1, 1];
const t = 0.5;
const degree = 3;

const point = deBoor(controlPoints, knots, t, degree);
console.log(point); // Output: [1.5, 0.75]
```

For debugging or educational purposes, you can use the functional version:

```typescript
import { deBoorDebug } from "b-spline-algorithms";

const point = deBoorDebug(controlPoints, knots, t, degree);
```

## API Reference

See our API documentation for detailed information on all available functions.

## Algorithm Explanations

For in-depth explanations of the implemented algorithms, check our algorithm guide.

## Performance

Our optimized implementations are designed for high performance. For benchmarks and performance tips, see our performance guide.

## Contributing

We welcome contributions! Please see our contributing guidelines for details on how to get started.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

List any references, papers, or other libraries that inspired or informed your implementations

## Contact

For questions and feedback, please open an issue on our GitHub repository.
