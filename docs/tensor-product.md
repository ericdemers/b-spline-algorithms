# Different approaches for implementing tensor product evaluation in B-splines

## Multiple strategies with their pros and cons

### Each implementation has its advantages:

1. [Sequential : Simple to understand and implement](#1-sequential-direction-evaluation)

2. [Recursive : Clean code structure, good for debugging](#2-recursive-implementation)

3. [Basis Product : More mathematical approach, good for analysis](#3-basis-function-product-approach)

4. [Block Processing : Better cache utilization](#4-optimized-block-processing)

5. [Cache-Aware : Best performance for repeated evaluations](#5-cache-aware-implementation)

6. [Category Theory View](#6-category-theory-view)

## 1. Sequential Direction Evaluation

```typescript
class BSplineTensorProduct<K extends Scalar, V extends Vector> {
  private evaluateTensorProduct(parameters: number[]): V {
    const dimensions = this.getDimension();
    let currentPoints = this.controlNet.getControlPoints();

    // Evaluate one direction at a time
    for (let dim = 0; dim < dimensions; dim++) {
      const u = parameters[dim];
      const degree = this.degrees[dim];
      const knots = this.knots.getKnotSequence(dim);

      currentPoints = this.evaluateInDirection(
        currentPoints,
        u,
        degree,
        knots,
        dim
      );
    }

    return currentPoints[0];
  }

  private evaluateInDirection(
    points: V[],
    u: number,
    degree: number,
    knots: number[],
    direction: number
  ): V[] {
    const span = this.findSpan(u, degree, knots);
    const basis = this.computeBasisValues(span, u, degree, knots);

    // Compute new points using basis functions
    const newPoints: V[] = [];
    // stride represents the distance (number of elements) you need to skip to move to the next point in a particular direction within a flattened multi-dimensional array.
    const stride = this.getStride(direction);

    for (let i = 0; i <= degree; i++) {
      const point = this.vectorSpace.multiply(
        basis[i],
        points[span - degree + i]
      );
      newPoints.push(point);
    }

    return this.combinePoints(newPoints);
  }
}
```

## 2. Recursive Implementation

### Mathematical Background

The recursive approach follows the mathematical definition of tensor products:

For a 2D surface S(u,v): S(u,v) = ∑ᵢ (∑ⱼ Pᵢⱼ Nⱼ,ₖ(v)) Nᵢ,ₖ(u)

The recursion implements this as:

1. First evaluate in v direction: Qᵢ(v) = ∑ⱼ Pᵢⱼ Nⱼ,ₖ(v)
2. Then evaluate in u direction: S(u,v) = ∑ᵢ Qᵢ(v) Nᵢ,ₖ(u)

```typescript
class RecursiveTensorProduct<K extends Scalar, V extends Vector> {
  // Main evaluation entry point
  evaluate(parameters: number[]): V {
    // Start recursion with all control points and dimension 0
    return this.evaluateRecursive(
      this.controlNet.getControlPoints(),
      parameters,
      0
    );
  }

  private evaluateRecursive(
    points: V[],
    parameters: number[],
    currentDimension: number
  ): V {
    // Base case: when we've processed all dimensions
    if (currentDimension === parameters.length) {
      // Only one point should remain after all evaluations
      return points[0];
    }

    // Get current dimension's parameters
    const u = parameters[currentDimension];
    const degree = this.degrees[currentDimension];
    const knots = this.knots.getKnotSequence(currentDimension);

    // Evaluate along current dimension
    const newPoints = this.evaluateInDirection(
      points,
      u,
      degree,
      knots,
      currentDimension
    );

    // Recurse to next dimension with reduced point set
    return this.evaluateRecursive(newPoints, parameters, currentDimension + 1);
  }

  private evaluateInDirection(
    points: V[],
    u: number,
    degree: number,
    knots: number[],
    direction: number
  ): V[] {
    // Find the knot span containing u
    const span = this.findSpan(u, degree, knots);

    // Compute basis function values
    const basis = this.computeBasisValues(span, u, degree, knots);

    // Get the stride for this direction
    const stride = this.getDirectionalStride(direction);

    // Number of points to process in this direction
    const numPoints = this.getNumPointsInDirection(direction);

    // Result array for this direction's evaluation
    const result: V[] = [];

    // Process each group of points
    for (let i = 0; i < points.length / stride; i++) {
      let sum = this.vectorSpace.zero();

      // Apply basis functions to points in span
      for (let j = 0; j <= degree; j++) {
        const pointIndex =
          i * stride + (span - degree + j) * this.getStride(direction);
        const weighted = this.vectorSpace.multiply(
          basis[j],
          points[pointIndex]
        );
        sum = this.vectorSpace.add(sum, weighted);
      }

      result.push(sum);
    }

    return result;
  }

  // Helper method to get stride for a direction
  private getDirectionalStride(direction: number): number {
    const sizes = this.controlNet.getSizes();
    let stride = 1;
    for (let i = direction + 1; i < sizes.length; i++) {
      stride *= sizes[i];
    }
    return stride;
  }

  // Helper method to get number of points in a direction
  private getNumPointsInDirection(direction: number): number {
    return this.controlNet.getSizes()[direction];
  }
}
```

## 3. Basis Function Product Approach

```typescript
class BasisProductTensorProduct<K extends Scalar, V extends Vector> {
  evaluate(parameters: number[]): V {
    const dimensions = this.getDimension();
    const basisValues: number[][] = [];

    // Compute basis functions for each direction
    for (let dim = 0; dim < dimensions; dim++) {
      const u = parameters[dim];
      const degree = this.degrees[dim];
      const knots = this.knots.getKnotSequence(dim);
      const span = this.findSpan(u, degree, knots);

      basisValues[dim] = this.computeBasisValues(span, u, degree, knots);
    }

    // Compute tensor product using basis function products
    let result = this.vectorSpace.zero();
    this.iterateControlPoints((point: V, indices: number[]) => {
      let weight = 1.0;
      for (let dim = 0; dim < dimensions; dim++) {
        weight *= basisValues[dim][indices[dim]];
      }

      const weighted = this.vectorSpace.multiply(weight, point);
      result = this.vectorSpace.add(result, weighted);
    });

    return result;
  }

  private iterateControlPoints(
    callback: (point: V, indices: number[]) => void
  ): void {
    const dimensions = this.getDimension();
    const sizes = this.controlNet.getSizes();
    const indices = new Array(dimensions).fill(0);

    const iterate = (dim: number) => {
      if (dim === dimensions) {
        callback(this.controlNet.getPoint(...indices), indices);
        return;
      }

      for (let i = 0; i < sizes[dim]; i++) {
        indices[dim] = i;
        iterate(dim + 1);
      }
    };

    iterate(0);
  }
}
```

## 4. Optimized Block Processing

This method processes multiple dimensions simultaneously in blocks to improve cache utilization and reduce memory access overhead.

Instead of evaluating dimensions one at a time (like in recursive approach), block processing evaluates multiple dimensions simultaneously to improve performance.

// Traditional approach (one dimension at a time):
Dimension 1 → Dimension 2 → Dimension 3 → Dimension 4

// Block approach (with block size = 2):
[Dimensions 1,2] → [Dimensions 3,4]

```typescript
class BlockTensorProduct<K extends Scalar, V extends Vector> {
  // Block size determines how many dimensions we process at once
  private static readonly BLOCK_SIZE = 4;

  evaluate(parameters: number[]): V {
    const dimensions = this.getDimension();
    let points = this.controlNet.getControlPoints();

    // Process dimensions in blocks
    for (let dim = 0; dim < dimensions; dim += this.BLOCK_SIZE) {
      // Handle remaining dimensions if less than BLOCK_SIZE
      const blockSize = Math.min(this.BLOCK_SIZE, dimensions - dim);

      // Process current block of dimensions
      points = this.evaluateBlock(
        points,
        parameters.slice(dim, dim + blockSize),
        dim,
        blockSize
      );
    }

    return points[0];
  }

  private evaluateBlock(
    points: V[],
    blockParameters: number[],
    startDim: number,
    blockSize: number
  ): V[] {
    // Pre-compute all basis values and spans for the block
    const basisValues: number[][] = [];
    const spans: number[] = [];

    // Compute basis values for all dimensions in block
    for (let i = 0; i < blockSize; i++) {
      const dim = startDim + i;
      const u = blockParameters[i];
      const degree = this.degrees[dim];
      const knots = this.knots.getKnotSequence(dim);

      spans[i] = this.findSpan(u, degree, knots);
      basisValues[i] = this.computeBasisValues(spans[i], u, degree, knots);
    }

    // Process points using block-optimized algorithm
    return this.processPointsBlock(
      points,
      basisValues,
      spans,
      startDim,
      blockSize
    );
  }

  private processPointsBlock(
    points: V[],
    basisValues: number[][],
    spans: number[],
    startDim: number,
    blockSize: number
  ): V[] {
    // Calculate strides for each dimension in block
    const strides: number[] = [];
    const numPoints: number[] = [];
    for (let i = 0; i < blockSize; i++) {
      const dim = startDim + i;
      strides[i] = this.getDirectionalStride(dim);
      numPoints[i] = this.getNumPointsInDirection(dim);
    }

    // Initialize result array
    const result: V[] = [];

    // Create indices array for block processing
    const indices = new Array(blockSize).fill(0);
    const maxIndices = indices.map((_, i) => this.degrees[startDim + i] + 1);

    // Process all combinations of basis functions in block
    this.iterateBlockIndices(indices, maxIndices, (indices) => {
      // Calculate combined weight for current combination
      let weight = 1.0;
      for (let i = 0; i < blockSize; i++) {
        weight *= basisValues[i][indices[i]];
      }

      // Calculate point index for current combination
      let pointIndex = 0;
      for (let i = 0; i < blockSize; i++) {
        pointIndex +=
          (spans[i] - this.degrees[startDim + i] + indices[i]) * strides[i];
      }

      // Add weighted point to result
      const weighted = this.vectorSpace.multiply(weight, points[pointIndex]);

      if (result.length === 0) {
        result.push(weighted);
      } else {
        result[0] = this.vectorSpace.add(result[0], weighted);
      }
    });

    return result;
  }

  private iterateBlockIndices(
    indices: number[],
    maxIndices: number[],
    callback: (indices: number[]) => void,
    dimension: number = 0
  ) {
    if (dimension === indices.length) {
      callback(indices);
      return;
    }

    for (let i = 0; i < maxIndices[dimension]; i++) {
      indices[dimension] = i;
      this.iterateBlockIndices(indices, maxIndices, callback, dimension + 1);
    }
  }

  // Helper method to optimize memory access patterns
  private optimizeMemoryAccess(
    points: V[],
    blockSize: number,
    strides: number[]
  ): V[] {
    // Reorder points to improve cache locality
    const reorderedPoints: V[] = [];
    const blockStride = this.calculateBlockStride(strides, blockSize);

    for (let i = 0; i < points.length; i += blockStride) {
      for (let j = 0; j < blockStride; j++) {
        reorderedPoints.push(points[this.reorderIndex(i + j, blockStride)]);
      }
    }

    return reorderedPoints;
  }

  // Calculate optimal block stride
  private calculateBlockStride(strides: number[], blockSize: number): number {
    let stride = 1;
    for (let i = 0; i < blockSize; i++) {
      stride *= strides[i];
    }
    return stride;
  }

  // Example usage
  private example(): void {
    // Process a 4D tensor product with block size 2
    const parameters = [0.5, 0.3, 0.7, 0.2];
    let points = this.controlNet.getControlPoints();

    // First block: dimensions 0 and 1
    points = this.evaluateBlock(points, [0.5, 0.3], 0, 2);

    // Second block: dimensions 2 and 3
    points = this.evaluateBlock(points, [0.7, 0.2], 2, 2);

    // Result is in points[0]
  }
}
```

## 5. Cache-Aware Implementation

```typescript
class CacheAwareTensorProduct<K extends Scalar, V extends Vector> {
  private readonly basisCache: Map<string, number[]> = new Map();

  evaluate(parameters: number[]): V {
    const dimensions = this.getDimension();
    let points = this.controlNet.getControlPoints();

    // Reorder dimensions for cache efficiency
    const dimOrder = this.optimizeDimensionOrder();

    for (const dim of dimOrder) {
      const u = parameters[dim];
      points = this.evaluateWithCaching(points, u, dim);
    }

    return points[0];
  }

  private evaluateWithCaching(points: V[], u: number, dimension: number): V[] {
    const cacheKey = `${dimension}-${u.toFixed(6)}`;
    let basis = this.basisCache.get(cacheKey);

    if (!basis) {
      const degree = this.degrees[dimension];
      const knots = this.knots.getKnotSequence(dimension);
      const span = this.findSpan(u, degree, knots);

      basis = this.computeBasisValues(span, u, degree, knots);
      this.basisCache.set(cacheKey, basis);
    }

    return this.evaluateWithBasis(points, basis, dimension);
  }

  private optimizeDimensionOrder(): number[] {
    // Implement dimension reordering for cache efficiency
    // Consider control point layout and memory access patterns
    return [...Array(this.getDimension())].map((_, i) => i);
  }
}
```

## 6. Category Theory View

```typescript
// Functor representing our evaluation process
interface EvaluationFunctor<V> {
  map<A, B>(f: (a: A) => B, fa: (p: V) => A): (p: V) => B;
}

// Monoid for combining evaluation results
interface Monoid<V> {
  empty(): V;
  concat(a: V, b: V): V;
}

// Natural transformation for dimension reduction
interface DimensionReducer<V> {
  transform(points: V[], parameter: number): V[];
}
```

### B-Spline Evaluation as a Fold Operation

```typescript
class CategoryTheoreticBSpline<K extends Scalar, V extends Vector> {
  // Monoid instance for vector space operations
  private vectorMonoid: Monoid<V> = {
    empty: () => this.vectorSpace.zero(),
    concat: (a: V, b: V) => this.vectorSpace.add(a, b),
  };

  // Fold operation over dimensions
  evaluate(parameters: number[]): V {
    return this.fold(
      this.controlNet.getControlPoints(),
      parameters,
      this.vectorMonoid
    );
  }

  private fold(points: V[], parameters: number[], monoid: Monoid<V>): V {
    // Base case: single point
    if (parameters.length === 0) {
      return points[0];
    }

    // Current dimension reducer
    const reducer = this.createDimensionReducer(
      parameters[0],
      this.degrees[0],
      this.knots.getKnotSequence(0)
    );

    // Reduce current dimension
    const reducedPoints = reducer.transform(points, parameters[0]);

    // Recursive fold over remaining dimensions
    return this.fold(reducedPoints, parameters.slice(1), monoid);
  }

  // Create a dimension reducer (natural transformation)
  private createDimensionReducer(
    parameter: number,
    degree: number,
    knots: number[]
  ): DimensionReducer<V> {
    return {
      transform: (points: V[], u: number) => {
        const span = this.findSpan(u, degree, knots);
        const basis = this.computeBasisValues(span, u, degree, knots);

        // Map-Reduce operation
        return this.reduceWithBasis(points, basis, span, degree);
      },
    };
  }

  // Map-Reduce operation for a single dimension
  private reduceWithBasis(
    points: V[],
    basis: number[],
    span: number,
    degree: number
  ): V[] {
    return points.map((_, i) => {
      // Map operation
      const weightedPoints = basis.map((b, j) =>
        this.vectorSpace.multiply(b, points[span - degree + j])
      );

      // Reduce operation
      return weightedPoints.reduce(
        (acc, p) => this.vectorSpace.add(acc, p),
        this.vectorSpace.zero()
      );
    });
  }
}
```

### Categorical Composition of Operations

```typescript
// Kleisli composition for evaluation steps
class EvaluationComposition<V> {
  compose<A, B, C>(f: (a: A) => V[], g: (b: V[]) => C): (a: A) => C {
    return (a: A) => g(f(a));
  }
}

// Natural transformation between dimensions
interface DimensionalTransformation<V> {
  transform(source: V[], targetDimension: number): V[];
}

// Implementation using categorical concepts
class CategoryTheoreticEvaluation<V extends Vector> {
  private compose(
    transformations: DimensionalTransformation<V>[]
  ): (points: V[]) => V {
    return (points: V[]) =>
      transformations.reduce(
        (acc, transform, dim) => transform.transform(acc, dim),
        points
      )[0];
  }

  evaluate(parameters: number[]): V {
    const transformations = parameters.map((p, i) =>
      this.createTransformation(p, i)
    );

    return this.compose(transformations)(this.controlNet.getControlPoints());
  }
}
```

### Lens-like structure for parameter access

```typescript
// Lens-like structure for parameter access
interface ParameterLens<S, A> {
  get(s: S): A;
  set(s: S, a: A): S;
}

// Implementation for B-spline parameters
class BSplineParameterLens<V> {
  createLens(dimension: number): ParameterLens<number[], number> {
    return {
      get: (parameters: number[]) => parameters[dimension],
      set: (parameters: number[], value: number) => {
        const newParams = [...parameters];
        newParams[dimension] = value;
        return newParams;
      },
    };
  }

  // Use lens for parameter manipulation
  modifyParameter(
    parameters: number[],
    dimension: number,
    f: (n: number) => number
  ): number[] {
    const lens = this.createLens(dimension);
    return lens.set(parameters, f(lens.get(parameters)));
  }
}
```

### Functorial Nature of Evaluation

```typescript
class EvaluationFunctor<V extends Vector> {
  // Functor map operation
  map<A, B>(f: (a: A) => B, evaluation: (p: V[]) => A): (p: V[]) => B {
    return (points: V[]) => f(evaluation(points));
  }

  // Natural transformation between different vector spaces
  transform<W extends Vector>(points: V[], transformation: (v: V) => W): W[] {
    return points.map(transformation);
  }
}
```

This categorical view provides several benefits:

1. Clear separation of concerns

2. Composable operations

3. Mathematical rigor

4. Abstraction over implementation details

### The categorical approach can be made efficient through careful implementation

```typescript
class OptimizedCategoryEvaluation<K extends Scalar, V extends Vector> {
  // Cache for intermediate results
  private readonly cache: Map<string, V[]> = new Map();

  // Optimized fold operation with caching
  private foldWithCache(
    points: V[],
    parameters: number[],
    monoid: Monoid<V>
  ): V {
    const cacheKey = this.createCacheKey(parameters);

    // Check cache first
    if (this.cache.has(cacheKey)) {
      return this.cache.get(cacheKey)![0];
    }

    // Optimize memory allocation
    const workingBuffer = new Array<V>(points.length);
    let currentBuffer = points;

    // Process dimensions in optimal order
    const optimizedOrder = this.getOptimizedDimensionOrder(parameters);

    for (const dim of optimizedOrder) {
      const reducer = this.createOptimizedReducer(
        parameters[dim],
        this.degrees[dim],
        this.knots.getKnotSequence(dim)
      );

      // In-place reduction when possible
      currentBuffer = reducer.transformInPlace(
        currentBuffer,
        workingBuffer,
        this.getStride(dim)
      );
    }

    // Cache result for future use
    this.cache.set(cacheKey, [currentBuffer[0]]);

    return currentBuffer[0];
  }

  // Optimized reducer with SIMD operations when available
  private createOptimizedReducer(
    parameter: number,
    degree: number,
    knots: number[]
  ): OptimizedDimensionReducer<V> {
    return {
      transformInPlace: (points: V[], buffer: V[], stride: number) => {
        const span = this.findSpan(parameter, degree, knots);
        const basis = this.computeBasisValues(span, parameter, degree, knots);

        // Use SIMD operations when available
        if (this.vectorSpace.supportsSIMD()) {
          return this.reduceSIMD(points, basis, span, degree, stride, buffer);
        }

        return this.reduceStandard(points, basis, span, degree, stride, buffer);
      },
    };
  }

  // SIMD-enabled reduction
  private reduceSIMD(
    points: V[],
    basis: number[],
    span: number,
    degree: number,
    stride: number,
    buffer: V[]
  ): V[] {
    const vectorSize = this.vectorSpace.getSIMDSize();
    const alignedLength = Math.floor(points.length / vectorSize) * vectorSize;

    // Process SIMD-aligned chunks
    for (let i = 0; i < alignedLength; i += vectorSize) {
      let sum = this.vectorSpace.simdZero();

      for (let j = 0; j <= degree; j++) {
        const weight = basis[j];
        const pointIndex = span - degree + j;
        sum = this.vectorSpace.simdAdd(
          sum,
          this.vectorSpace.simdMultiply(
            weight,
            this.vectorSpace.simdLoad(points, i + pointIndex * stride)
          )
        );
      }

      this.vectorSpace.simdStore(buffer, i, sum);
    }

    // Handle remaining elements
    for (let i = alignedLength; i < points.length; i++) {
      buffer[i] = this.reducePoint(points, basis, span, degree, stride, i);
    }

    return buffer;
  }

  // Memory-efficient standard reduction
  private reduceStandard(
    points: V[],
    basis: number[],
    span: number,
    degree: number,
    stride: number,
    buffer: V[]
  ): V[] {
    // Process points in cache-friendly order
    const blockSize = 64 / this.vectorSpace.getElementSize(); // Typical cache line size

    for (let block = 0; block < points.length; block += blockSize) {
      const end = Math.min(block + blockSize, points.length);

      for (let i = block; i < end; i++) {
        buffer[i] = this.reducePoint(points, basis, span, degree, stride, i);
      }
    }

    return buffer;
  }

  // Optimized single point reduction
  private reducePoint(
    points: V[],
    basis: number[],
    span: number,
    degree: number,
    stride: number,
    index: number
  ): V {
    let sum = this.vectorSpace.zero();

    // Unroll small loops
    if (degree <= 4) {
      for (let j = 0; j <= degree; j++) {
        sum = this.vectorSpace.add(
          sum,
          this.vectorSpace.multiply(
            basis[j],
            points[index + (span - degree + j) * stride]
          )
        );
      }
    } else {
      // Use regular loop for higher degrees
      for (let j = 0; j <= degree; j++) {
        sum = this.vectorSpace.add(
          sum,
          this.vectorSpace.multiply(
            basis[j],
            points[index + (span - degree + j) * stride]
          )
        );
      }
    }

    return sum;
  }

  // Optimize dimension processing order
  private getOptimizedDimensionOrder(parameters: number[]): number[] {
    // Order dimensions based on memory access patterns and degrees
    return Array.from(parameters.keys()).sort((a, b) => {
      const costA = this.evaluationCost(a);
      const costB = this.evaluationCost(b);
      return costB - costA;
    });
  }

  private evaluationCost(dimension: number): number {
    return this.degrees[dimension] * this.getStride(dimension);
  }
}
```

### Key Optimization Techniques

1. Memory Management

```typescript
// Reuse buffers instead of creating new arrays
private readonly workingBuffer: V[];
private readonly resultBuffer: V[];
```

2. SIMD Operations

```typescript
if (this.vectorSpace.supportsSIMD()) {
  return this.reduceSIMD(points, basis, span, degree, stride, buffer);
}
```

3. Cache Optimization

```typescript
// Process in cache-friendly blocks
const blockSize = 64 / this.vectorSpace.getElementSize();
for (let block = 0; block < points.length; block += blockSize) {
  // Process block
}
```

4. Dimension Reordering

```typescript
private getOptimizedDimensionOrder(parameters: number[]): number[] {
    return Array.from(parameters.keys()).sort((a, b) => {
        const costA = this.evaluationCost(a);
        const costB = this.evaluationCost(b);
        return costB - costA;
    });
}
```

### A Categorical Approach to Pyramid Computations

```typescript
// Core abstractions for pyramid computations
interface PyramidLevel<V> {
  readonly level: number;
  readonly data: V[];
  readonly metadata: PyramidMetadata;
}

interface PyramidMetadata {
  readonly degree: number;
  readonly knots: number[];
  readonly dimension: number;
}

// Functor representing pyramid transformations
interface PyramidFunctor<V> {
  map<A, B>(f: (a: A) => B, fa: PyramidLevel<A>): PyramidLevel<B>;
  lift<A>(value: A): PyramidLevel<A>;
}

// Natural transformations between pyramid levels
interface PyramidTransformation<V> {
  transform(source: PyramidLevel<V>): PyramidLevel<V>;
}

// Generic pyramid computation framework
class PyramidComputation<V extends Vector> {
  // Unified interface for pyramid-like algorithms
  private compute<T extends PyramidTransformation<V>>(
    initial: PyramidLevel<V>,
    transformation: T
  ): PyramidLevel<V> {
    return transformation.transform(initial);
  }

  // Composition of transformations
  private compose<T extends PyramidTransformation<V>>(
    transformations: T[]
  ): PyramidTransformation<V> {
    return {
      transform: (source: PyramidLevel<V>) =>
        transformations.reduce((acc, t) => t.transform(acc), source),
    };
  }
}

// Implementation for B-spline evaluation
class BSplineEvaluationPyramid<V extends Vector>
  implements PyramidTransformation<V>
{
  constructor(private parameter: number) {}

  transform(source: PyramidLevel<V>): PyramidLevel<V> {
    const { degree, knots } = source.metadata;
    const span = this.findSpan(this.parameter, degree, knots);
    const basis = this.computeBasisValues(span, this.parameter, degree, knots);

    return {
      level: source.level + 1,
      data: this.computeLevel(source.data, basis, span, degree),
      metadata: {
        ...source.metadata,
        degree: degree - 1,
      },
    };
  }

  private computeLevel(
    points: V[],
    basis: number[],
    span: number,
    degree: number
  ): V[] {
    // Implementation of de Boor's algorithm
    return /* ... */;
  }
}

// Implementation for knot insertion
class KnotInsertionPyramid<V extends Vector>
  implements PyramidTransformation<V>
{
  constructor(private newKnot: number) {}

  transform(source: PyramidLevel<V>): PyramidLevel<V> {
    const { degree, knots } = source.metadata;
    const span = this.findSpan(this.newKnot, degree, knots);

    return {
      level: source.level + 1,
      data: this.computeLevel(source.data, span, degree),
      metadata: {
        ...source.metadata,
        knots: this.insertKnot(knots, this.newKnot),
      },
    };
  }

  private computeLevel(points: V[], span: number, degree: number): V[] {
    // Implementation of knot insertion algorithm
    return /* ... */;
  }
}

// Implementation for degree elevation
class DegreeElevationPyramid<V extends Vector>
  implements PyramidTransformation<V>
{
  transform(source: PyramidLevel<V>): PyramidLevel<V> {
    const { degree, knots } = source.metadata;

    return {
      level: source.level + 1,
      data: this.computeLevel(source.data, degree),
      metadata: {
        ...source.metadata,
        degree: degree + 1,
        knots: this.elevateKnots(knots),
      },
    };
  }

  private computeLevel(points: V[], degree: number): V[] {
    // Implementation of degree elevation algorithm
    return /* ... */;
  }
}

// Unified interface for all pyramid computations
class UnifiedBSplineComputation<V extends Vector> {
  evaluate(points: V[], parameter: number): V {
    const initial: PyramidLevel<V> = {
      level: 0,
      data: points,
      metadata: this.getInitialMetadata(),
    };

    const evaluation = new BSplineEvaluationPyramid(parameter);
    const result = this.compute(initial, evaluation);
    return result.data[0];
  }

  insertKnot(points: V[], newKnot: number): V[] {
    const initial: PyramidLevel<V> = {
      level: 0,
      data: points,
      metadata: this.getInitialMetadata(),
    };

    const insertion = new KnotInsertionPyramid(newKnot);
    const result = this.compute(initial, insertion);
    return result.data;
  }

  // Composition example
  evaluateWithRefinement(
    points: V[],
    parameter: number,
    refinementKnots: number[]
  ): V {
    const initial: PyramidLevel<V> = {
      level: 0,
      data: points,
      metadata: this.getInitialMetadata(),
    };

    // Create transformations
    const insertions = refinementKnots.map(
      (k) => new KnotInsertionPyramid<V>(k)
    );
    const evaluation = new BSplineEvaluationPyramid<V>(parameter);

    // Compose transformations
    const transformation = this.compose([...insertions, evaluation]);

    const result = this.compute(initial, transformation);
    return result.data[0];
  }

  // Optimized implementation
  private compute<T extends PyramidTransformation<V>>(
    initial: PyramidLevel<V>,
    transformation: T
  ): PyramidLevel<V> {
    // Add caching, memory management, etc.
    return transformation.transform(initial);
  }
}
```
