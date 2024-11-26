# Detailed API Documentation

# Detailed API Documentation

## Core Architecture

### Vector Spaces

The foundation of the library is built on abstract vector spaces, allowing operations on different types of geometric data:

```typescript
interface VectorSpace<T> {
  add(a: T, b: T): T;
  multiply(scalar: number, vector: T): T;
  zero(): T;
}

// Built-in vector spaces
const vector2DSpace: VectorSpace<[number, number]>;
const vector3DSpace: VectorSpace<[number, number, number]>;
const realVectorSpace: VectorSpace<number>;
```

## Base B-Spline Structure

The library uses a generic base class that handles multi-dimensional B-splines:

```typescript
export class BSpline<K extends Scalar, V extends Vector> {
  constructor(
    protected readonly vectorSpace: VectorSpace<K, V>,
    protected readonly controlNet: ControlNet<V>,
    protected readonly knots: KnotStructure,
    protected readonly degrees: ReadonlyArray<number>
  ) {
    // Core B-spline properties independent of dimension
  }

  evaluate(parameters: ReadonlyArray<number> | readonly number): V;
  derivative(orders: ReadonlyArray<number> | readonly number): BSpline<K, V>;
}
```

## Specialize Implementations

### B-Spline Curves

```typescript
function adaptParameter() {
  return function (
    target: any,
    propertyKey: string,
    descriptor: PropertyDescriptor
  ) {
    const originalMethod = descriptor.value;
    descriptor.value = function (position: number | number[]) {
      const input = Array.isArray(position) ? position : [position];
      return originalMethod.call(this, input);
    };
  };
}
```

```typescript
export class BSplineCurve<K extends Scalar, V extends Vector> extends BSpline<
  K,
  V
> {
  constructor(
    vectorSpace: VectorSpace<K, V>,
    controlPoints: ReadonlyArray<V>,
    degree: number,
    knots?: ReadonlyArray<number>
  ) {
    // Convert to general form
    super(
      vectorSpace,
      new ControlPointArray(controlPoints),
      new UnivariateKnots(
        knots ?? generateUniformKnots(controlPoints.length, degree)
      ),
      [degree]
    );
  }

  // Simplified curve-specific methods
  @adaptParameter()
  evaluate(position: number | number[]): number {
    return super.evaluate(position as number[]);
  }

  @adaptParameter()
  derivative(order: number | number[] = 1): BSplineCurve<K, V> {
    // Curve-specific implementation using general machinery
  }
}
```

### BSpline<T>

Base class for B-splines in arbitrary vector spaces.

#### Type Parameters

- `T`: The vector space element type (e.g., number for 1D, [number, number] for 2D)

#### Constructor

```typescript
constructor(
    vectorSpace: VectorSpace<T>,
    controlPoints: ReadonlyArray<T>,
    knots: ReadonlyArray<number>,
    degree: number,
    options?: {
        evaluation?: EvaluationAlgorithm<T>;
        knotInsertion?: KnotInsertionAlgorithm<T>;
        multiplication?: MultiplicationAlgorithm<T>;
    }
)
```

```typescript

readonly degree: number;
readonly domain: [number, number];
readonly controlPoints: ReadonlyArray<T>;
readonly knots: ReadonlyArray<number>;

```

### Utility Types

#### BSplineData&lt;T&gt;

Interface representing B-spline data structure.

```typescript
interface BSplineData<T> {
  controlPoints: ReadonlyArray<T>;
  knots: ReadonlyArray<number>;
  degree: number;
}
```

### Algorithm Interfaces

```typescript
interface EvaluationAlgorithm<T> {
  evaluate(bspline: BSplineData<T>, parameter: number): T;
}
```

#### KnotInsertionAlgorithm&lt;T&gt;

```typescript
interface KnotInsertionAlgorithm<T> {
  insertKnot(bspline: BSplineData<T>, parameter: number): BSplineData<T>;
}
```

### Algorithm Selection

#### Using Different Evaluation Algorithms

```typescript
// For educational purposes
const educational = new BSplineCurve2D(points, degree, {
  evaluation: new ClearEvaluation(),
});

// For performance
const optimized = new BSplineCurve2D(points, degree, {
  evaluation: new OptimizedEvaluation(),
});
```

### Error Handling

```typescript
try {
  const curve = new BSplineCurve2D(points, degree);
  const value = curve.evaluate(0.5);
} catch (error) {
  if (error instanceof InvalidParameterError) {
    // Handle parameter out of domain
  } else if (error instanceof DegenerateBSplineError) {
    // Handle invalid B-spline configuration
  }
}
```

### Memory Management

```typescript
// Immutable operations return new instances
const refined = curve.insertKnot(0.5); // Original curve unchanged
const elevated = curve.elevateDegree(); // Original curve unchanged
```

## Knot Structure Module

The knot structure module provides interfaces and implementations for handling knot vectors in B-spline representations.

### Interfaces

#### KnotValue

Represents a knot value with its multiplicity.

```typescript
interface KnotValue {
  readonly value: number; // The knot parameter value
  readonly multiplicity: number; // Number of times the value appears
}
```

#### KnotStructure

Base interface for all knot structures.

```typescript
interface KnotStructure {
  getDimension(): number;
  getKnotSequence(direction: number): ReadonlyArray<number>;
  withInsertedKnot(dimension: number, u: number): KnotStructure;
  withRemovedKnots(dimension: number): KnotStructure;
  getDomain(direction: number): Domain;
  getDistinctKnots(direction: number): ReadonlyArray<KnotValue>;
}
```

### Classes

#### BaseKnotStructure

Abstract base class providing common functionality for knot structures.

```typescript
abstract class BaseKnotStructure implements KnotStructure {
  protected validateDirection(direction: number): void;
  protected isNonDecreasing(knots: ReadonlyArray<number>): boolean;
  protected computeDistinctKnots(
    knots: ReadonlyArray<number>
  ): ReadonlyArray<KnotValue>;
  protected findInsertionIndex(knots: ReadonlyArray<number>, u: number): number;
}
```

#### Knots

Implementation for single-dimensional knot vectors.

```typescript
class Knots extends BaseKnotStructure {
  constructor(knots: ReadonlyArray<number>);
}
```

#### Example:

```typescript
// Create a knot vector for a cubic B-spline
const knots = new Knots([0, 0, 0, 0, 1, 2, 3, 3, 3, 3]);

// Get domain bounds
const domain = knots.getDomain(0); // { min: 0, max: 3 }

// Get distinct knots with multiplicities
const distinct = knots.getDistinctKnots(0);
// [
//   { value: 0, multiplicity: 4 },
//   { value: 1, multiplicity: 1 },
//   { value: 2, multiplicity: 1 },
//   { value: 3, multiplicity: 4 }
// ]
```

### ProductKnots

Implementation for multi-dimensional tensor-product knot structures.

```typescript
class ProductKnots extends BaseKnotStructure {
  constructor(knotVectors: ReadonlyArray<ReadonlyArray<number>>);
}
```

#### Example:

```typescript
// Create a knot structure for a B-spline surface
const surfaceKnots = new ProductKnots([
  [0, 0, 0, 1, 1, 1], // u-direction
  [0, 0, 1, 2, 2], // v-direction
]);

// Get knot sequence for v-direction
const vKnots = surfaceKnots.getKnotSequence(1);

// Insert a knot in u-direction
const refined = surfaceKnots.withInsertedKnot(0, 0.5);
```

### Error Handling

```typescript
try {
  const knots = new Knots([2, 1, 3]); // Will throw - not non-decreasing
} catch (error) {
  console.error("Invalid knot vector");
}

try {
  const knots = new Knots([0, 1, 2]);
  knots.getKnotSequence(1); // Will throw - invalid direction
} catch (error) {
  console.error("Invalid direction");
}
```

### Periodic Knot Structure

The library provides specialized support for periodic knot vectors, commonly used in closed B-spline curves:

```typescript
interface DistinctKnotsWithMultiplicities {
  readonly knots: ReadonlyArray<number>;
  readonly multiplicities: ReadonlyArray<number>;
}

class PeriodicKnots extends BaseKnotStructure {
  constructor(
    private readonly pattern: ReadonlyArray<number>,
    private readonly period: number
  );

  // Core KnotStructure interface
  getDimension(): number;
  getKnotSequence(direction: number): ReadonlyArray<number>;
  withInsertedKnot(dimension: number, u: number): KnotStructure;
  getDomain(direction: number): Domain;
  getDistinctKnots(direction: number): ReadonlyArray<KnotValue>;

  // Periodic-specific methods
  getPattern(): ReadonlyArray<number>;
  getPeriod(): number;

  // Static factory method
  static fromDistinctKnots(
    values: DistinctKnotsWithMultiplicities,
    period: number
  ): PeriodicKnots;
}
```

#### Key Features

1. Pattern-Based Construction

```typescript
// Create periodic knots from a repeating pattern
const knots = new PeriodicKnots(
  [0, 0.25, 0.5, 0.75], // pattern
  1.0 // period
);
```

2. Distinct Knots Factory

```typescript
// Create from distinct knots with multiplicities
const periodicKnots = PeriodicKnots.fromDistinctKnots(
  {
    knots: [0, 0.25, 0.5, 0.75],
    multiplicities: [2, 1, 1, 2],
  },
  1.0
);
```

3. Knot Operations

```typescript
// Get unrolled knot sequence
const sequence = knots.getKnotSequence(0); // direction = 0 for curves
```

#### Advanced Features

1. Parameter Normalization

```typescript
// Parameters are automatically normalized to the periodic domain
const span = knots.findKnotSpan(1.25); // Automatically wraps to [0, period)
```

2. Knot Insertion Support

```typescript
// Get relevant knots for insertion operation
const { knots: relevantKnots, knotsToBeInserted } =
  periodicKnots.unrollForKnotInsertion(3, 0.5); // degree 3, parameter 0.5
```

3. Domain Manageement

```typescript
// Get the valid parameter domain
const domain = knots.getDomain(0); // { min: 0, max: 1.0 }
```

#### Implementation Details

1. Validation

- Pattern must be non-decreasing
- Pattern range must not exceed period
- Pattern must have at least 2 values
- Period must be positive

2. Knot Span Finding

```typescript
// Efficient binary search implementation
const span = knots.findKnotSpan(u);
```

3. Knot Value Access

```typescript
// Get knot value at any index (handles periodic repetition)
const value = knots.getKnotValue(index);
```

### Best Practices For Periodic Knots

1. Period Selection

- Choose period based on geometric meaning
- Typically use 1.0 or 2π for normalized parameters
- Ensure period matches control point structure

2. Pattern Design

- Use uniform patterns for regular parameterization
- Consider multiplicity for continuity control
- Ensure pattern covers desired parameter range

3. Error Handling

```typescript
try {
  const knots = new PeriodicKnots(pattern, period);
} catch (error) {
  if (error.message.includes("Pattern must be non-decreasing")) {
    // Handle invalid pattern
  } else if (error.message.includes("Period must be positive")) {
    // Handle invalid period
  }
}
```

### Best Practices

#### Immutability

All knot structure operations return new instances:

```typescript
const original = new Knots([0, 0, 1, 1]);
const modified = original.withInsertedKnot(0, 0.5);
// original is unchanged
```

#### Validation

Always validate knot vectors:

```typescript
// Good: Non-decreasing knot vector
const valid = new Knots([0, 0, 0, 1, 2, 2, 2]);

// Bad: Will throw error
const invalid = new Knots([1, 0, 2]);
```

#### Memory Efficiency

Use shared knot vectors for multiple B-splines when possible:

```typescript
const sharedKnots = new Knots([0, 0, 1, 1]);
const curve1 = new BSpline(controlPoints1, sharedKnots, 1);
const curve2 = new BSpline(controlPoints2, sharedKnots, 1);
```

## Control Net Module

The control net module provides structures and operations for managing control points in B-spline representations.

### Interfaces

#### ControlPoint&lt;T&gt;

Represents a control point with optional weight for rational B-splines.

```typescript
interface ControlPoint<T> {
  readonly position: T; // Position in vector space
  readonly weight?: number; // Weight for rational B-splines
}
```

#### ControlNet&lt;T&gt;

Base interface for control point structures of any dimension.

```typescript
interface ControlNet<T> {
  getDimension(): number;
  getControlPoints(): ReadonlyArray<ControlPoint<T>>;
  getControlPointAt(...indices: number[]): ControlPoint<T>;
  withInsertedPoint(
    point: ControlPoint<T>,
    ...indices: number[]
  ): ControlNet<T>;
  withModifiedPoint(
    point: ControlPoint<T>,
    ...indices: number[]
  ): ControlNet<T>;
  map<U>(transform: (point: T) => U): ControlNet<U>;
}
```

### Classes

#### BaseControlNet&lt;T&gt;

Abstract base class providing common functionality for control nets.

```typescript
abstract class BaseControlNet<T> implements ControlNet<T> {
  protected validateIndices(...indices: number[]): void;
  protected computeBoundingBox(): BoundingBox<T>;
  protected isValid(): boolean;
}
```

#### ControlPointArray&lt;T&gt;

Implementation for one-dimensional control point arrays (curves).

```typescript
class ControlPointArray<T> extends BaseControlNet<T> {
  constructor(points: ReadonlyArray<ControlPoint<T>>);

  getSize(): number;
  getBoundingBox(): BoundingBox<T>;
}
```

#### Example:

```typescript
// Create a control point array for a curve
const points = new ControlPointArray([
  { position: [0, 0] },
  { position: [1, 1] },
  { position: [2, 0] },
]);

// Get a specific control point
const point = points.getControlPointAt(1); // { position: [1, 1] }

// Create new array with modified point
const modified = points.withModifiedPoint({ position: [1, 2] }, 1);
```

#### ControlPointGrid&lt;T&gt;

Implementation for two-dimensional control point grids (surfaces).

```typescript
class ControlPointGrid<T> extends BaseControlNet<T> {
  constructor(points: ReadonlyArray<ReadonlyArray<ControlPoint<T>>>);

  getSize(): [number, number];
  getBoundingBox(): BoundingBox<T>;
}
```

#### Example:

```typescript
// Create a control point grid for a surface
const grid = new ControlPointGrid([
  [{ position: [0, 0, 0] }, { position: [1, 0, 0] }],
  [{ position: [0, 1, 1] }, { position: [1, 1, 1] }],
]);

// Get a specific control point
const point = grid.getControlPointAt(1, 0); // { position: [0, 1, 1] }

// Create new grid with inserted row/column
const refined = grid.withInsertedPoint({ position: [0.5, 0.5, 0.5] }, 1, 1);
```

### Utility Types

#### BoundingBox&lt;T&gt;

Represents the bounding box of a control net.

```typescript
interface BoundingBox<T> {
  readonly min: T;
  readonly max: T;
}
```

### Operations

#### Transform Operations

```typescript
// Transform all control points
const scaled = controlNet.map((p) => vectorSpace.multiply(2, p));

// Project to 2D
const projected = controlNet3D.map(([x, y, z]) => [x, y]);
```

#### Rational Operations

```typescript
// Create rational control points
const rational = new ControlPointArray([
  { position: [0, 0], weight: 1 },
  { position: [1, 1], weight: 0.7071 },
  { position: [2, 0], weight: 1 },
]);

// Convert to homogeneous coordinates
const homogeneous = rational.map((p) => [
  p.position[0] * p.weight,
  p.position[1] * p.weight,
  p.weight,
]);
```

#### Error Handling

```typescript
try {
  const grid = new ControlPointGrid([
    [{ position: [0, 0] }],
    [{ position: [1, 0] }, { position: [1, 1] }],
  ]); // Will throw - inconsistent row lengths
} catch (error) {
  console.error("Invalid control point grid");
}

try {
  const points = new ControlPointArray([]);
  points.getControlPointAt(0); // Will throw - index out of bounds
} catch (error) {
  console.error("Invalid index");
}
```

### Best Practices

#### Immutability

All control net operations return new instances:

```typescript
const original = new ControlPointArray(points);
const modified = original.withModifiedPoint(newPoint, 1);
// original is unchanged
```

#### Type Safety

Use type parameters to ensure consistent vector spaces:

```typescript
// 2D control points
const points2D: ControlPointArray<[number, number]> = /*...*/;

// 3D control points
const points3D: ControlPointArray<[number, number, number]> = /*...*/;

```

#### Memory Efficiency

Share control point data when possible:

```typescript
const sharedPoints = new ControlPointArray(points);
const curve1 = new BSpline(sharedPoints, knots1, degree);
const curve2 = new BSpline(sharedPoints, knots2, degree);
```

#### Validation

Always validate control points structures:

```typescript
// Good: Consistent grid dimensions
const valid = new ControlPointGrid([
  [p00, p01, p02],
  [p10, p11, p12],
]);

// Bad: Inconsistent dimensions (will throw)
const invalid = new ControlPointGrid([
  [p00, p01],
  [p10, p11, p12],
]);
const sharedPoints = new ControlPointArray(points);
const curve1 = new BSpline(sharedPoints, knots1, degree);
const curve2 = new BSpline(sharedPoints, knots2, degree);
```
