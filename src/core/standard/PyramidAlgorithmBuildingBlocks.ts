// 1. Core Concepts - The Building Blocks
interface PyramidLevel<V> {
    readonly data: V[];
    readonly metadata: LevelMetadata;
}

interface LevelMetadata {
    readonly degree: number;
    readonly knots: number[];
    readonly level: number;
}

// 2. Transformation Rules - How levels change
interface TransformationRule<V> {
    computeNewData(level: PyramidLevel<V>): V[];
    updateMetadata(metadata: LevelMetadata): LevelMetadata;
}

// 3. Level Generator - Creates new levels
class LevelGenerator<V> {
    createNext(
        current: PyramidLevel<V>, 
        rule: TransformationRule<V>
    ): PyramidLevel<V> {
        return {
            data: rule.computeNewData(current),
            metadata: rule.updateMetadata(current.metadata)
        };
    }
}

// 4. Specific Rules for Different Operations
class EvaluationRule<V> implements TransformationRule<V> {
    constructor(private parameter: number) {}

    computeNewData(level: PyramidLevel<V>): V[] {
        return this.computeEvaluationPoints(
            level.data, 
            this.parameter, 
            level.metadata
        );
    }

    updateMetadata(metadata: LevelMetadata): LevelMetadata {
        return {
            ...metadata,
            degree: metadata.degree - 1,
            level: metadata.level + 1
        };
    }
}

class KnotInsertionRule<V> implements TransformationRule<V> {
    constructor(private newKnot: number) {}

    computeNewData(level: PyramidLevel<V>): V[] {
        return this.computeRefinedPoints(
            level.data, 
            this.newKnot, 
            level.metadata
        );
    }

    updateMetadata(metadata: LevelMetadata): LevelMetadata {
        return {
            ...metadata,
            knots: [...metadata.knots, this.newKnot].sort(),
            level: metadata.level + 1
        };
    }
}

// 5. Pyramid Builder - Manages the computation process
class PyramidBuilder<V> {
    private generator = new LevelGenerator<V>();

    build(
        initial: PyramidLevel<V>,
        rule: TransformationRule<V>,
        predicate: (level: PyramidLevel<V>) => boolean
    ): PyramidLevel<V> {
        let current = initial;

        while (predicate(current)) {
            current = this.generator.createNext(current, rule);
        }

        return current;
    }
}

// 6. High-Level Operations
class BSplineOperations<V> {
    private builder = new PyramidBuilder<V>();

        const initial = this.createInitialLevel(points);
        const rule = new EvaluationRule<V>(parameter);
        
        const result = this.builder.build(
            initial,
            rule,
            level => level.metadata.degree > 0
        );

        return result.data[0];
    }

    insertKnot(points: V[], newKnot: number): V[] {
        const initial = this.createInitialLevel(points);
        const rule = new KnotInsertionRule<V>(newKnot);
        
        const result = this.builder.build(
            initial,
            rule,
            level => this.needsMoreRefinement(level)
        );

        return result.data;
    }
}

// 7. Composition of Rules
class CompositeRule<V> implements TransformationRule<V> {
    constructor(private rules: TransformationRule<V>[]) {}

    computeNewData(level: PyramidLevel<V>): V[] {
        return this.rules.reduce(
            (data, rule) => rule.computeNewData({
                data,
                metadata: level.metadata
            }),
            level.data
        );
    }

    updateMetadata(metadata: LevelMetadata): LevelMetadata {
        return this.rules.reduce(
            (meta, rule) => rule.updateMetadata(meta),
            metadata
        );
    }
}

// 8. Example Usage
class BSplineExample<V> {
    private operations = new BSplineOperations<V>();

    evaluateWithRefinement(
        points: V[],
        parameter: number,
        newKnots: number[]
    ): V {
        const initial = this.createInitialLevel(points);
        
        // Compose rules
        const rules = [
            ...newKnots.map(knot => new KnotInsertionRule<V>(knot)),
            new EvaluationRule<V>(parameter)
        ];
        
        const compositeRule = new CompositeRule<V>(rules);
        
        const result = new PyramidBuilder<V>().build(
            initial,
            compositeRule,
            level => level.metadata.degree > 0
        );

        return result.data[0];
    }
}

// 9. Optimization Support
interface OptimizedRule<V> extends TransformationRule<V> {
    canOptimize(level: PyramidLevel<V>): boolean;
    optimizedCompute(level: PyramidLevel<V>): V[];
}

class OptimizedPyramidBuilder<V> {
    build(
        initial: PyramidLevel<V>,
        rule: OptimizedRule<V>
    ): PyramidLevel<V> {
        if (rule.canOptimize(initial)) {
            return {
                data: rule.optimizedCompute(initial),
                metadata: rule.updateMetadata(initial.metadata)
            };
        }
        
        return new PyramidBuilder<V>().build(
            initial,
            rule,
            level => true
        );
    }
}

// 1. Core Blossom Concepts
interface Blossom<V> {
    // Blossom is symmetric and multi-affine
    evaluate(...parameters: number[]): V;
    degree: number;
}

// 2. Blossom-based Level Structure
interface BlossomLevel<V> {
    readonly blossom: Blossom<V>;
    readonly metadata: LevelMetadata;
}

interface LevelMetadata {
    readonly degree: number;
    readonly knots: number[];
    readonly level: number;
}

// 3. Blossom Transformations
interface BlossomTransformation<V> {
    transform(blossom: Blossom<V>): Blossom<V>;
    updateMetadata(metadata: LevelMetadata): LevelMetadata;
}

// 4. Implementation of Blossom
class BSplineBlossom<V> implements Blossom<V> {
    constructor(
        private controlPoints: V[],
        private knots: number[],
        readonly degree: number,
        private vectorSpace: VectorSpace<V>
    ) {}

    evaluate(...parameters: number[]): V {
        if (parameters.length !== this.degree) {
            throw new Error('Wrong number of parameters');
        }

        // Base case: degree 0
        if (this.degree === 0) {
            return this.controlPoints[0];
        }

        // Recursive case: de Boor's algorithm through blossoming
        const newBlossom = this.reduce(parameters[0]);
        return newBlossom.evaluate(...parameters.slice(1));
    }

    // Reduce degree by one parameter
    private reduce(parameter: number): Blossom<V> {
        const span = this.findSpan(parameter);
        const newPoints = this.computeNewPoints(span, parameter);
        
        return new BSplineBlossom(
            newPoints,
            this.knots,
            this.degree - 1,
            this.vectorSpace
        );
    }
}

// 5. Blossom-based Operations
class EvaluationByBlossom<V> implements BlossomTransformation<V> {
    constructor(private parameter: number) {}

    transform(blossom: Blossom<V>): Blossom<V> {
        return {
            evaluate: (...params) => 
                blossom.evaluate(this.parameter, ...params),
            degree: blossom.degree - 1
        };
    }

    updateMetadata(metadata: LevelMetadata): LevelMetadata {
        return {
            ...metadata,
            degree: metadata.degree - 1,
            level: metadata.level + 1
        };
    }
}

class KnotInsertionByBlossom<V> implements BlossomTransformation<V> {
    constructor(private newKnot: number) {}

    transform(blossom: Blossom<V>): Blossom<V> {
        return {
            evaluate: (...params) => {
                // Use blossoming property for knot insertion
                const leftParams = params.map(() => this.newKnot);
                const rightParams = params;
                
                return this.vectorSpace.combine(
                    blossom.evaluate(...leftParams),
                    blossom.evaluate(...rightParams),
                    this.computeWeight(params)
                );
            },
            degree: blossom.degree
        };
    }

    updateMetadata(metadata: LevelMetadata): LevelMetadata {
        return {
            ...metadata,
            knots: [...metadata.knots, this.newKnot].sort(),
            level: metadata.level + 1
        };
    }
}

// 6. Pyramid Builder with Blossoming
class BlossomPyramidBuilder<V> {
    build(
        initial: BlossomLevel<V>,
        transformation: BlossomTransformation<V>,
        predicate: (level: BlossomLevel<V>) => boolean
    ): BlossomLevel<V> {
        let current = initial;

        while (predicate(current)) {
            current = {
                blossom: transformation.transform(current.blossom),
                metadata: transformation.updateMetadata(current.metadata)
            };
        }

        return current;
    }
}

// 7. High-level Operations Using Blossoming
class BSplineOperations<V> {
    evaluate(points: V[], parameter: number): V {
        const initial = this.createInitialBlossom(points);
        const transformation = new EvaluationByBlossom<V>(parameter);
        
        // Use blossom to evaluate
        return this.builder
            .build(initial, transformation, level => level.metadata.degree > 0)
            .blossom
            .evaluate();
    }

    // Multiple evaluation using blossom properties
    evaluateMultiple(points: V[], parameters: number[]): V[] {
        const blossom = this.createInitialBlossom(points).blossom;
        
        return parameters.map(u => {
            // Reuse blossom for efficiency
            const params = Array(blossom.degree).fill(u);
            return blossom.evaluate(...params);
        });
    }
}

// 8. Composition with Blossoming
class CompositeBlossomTransformation<V> implements BlossomTransformation<V> {
    constructor(private transformations: BlossomTransformation<V>[]) {}

    transform(blossom: Blossom<V>): Blossom<V> {
        return this.transformations.reduce(
            (current, transform) => transform.transform(current),
            blossom
        );
    }

    updateMetadata(metadata: LevelMetadata): LevelMetadata {
        return this.transformations.reduce(
            (current, transform) => transform.updateMetadata(current),
            metadata
        );
    }
}

// 9. Example Usage with Blossoming
class BSplineExample<V> {
    evaluateWithBlossoming(
        points: V[],
        parameters: number[]
    ): V[] {
        const blossom = this.createInitialBlossom(points);
        
        // Efficient multiple evaluation using blossom properties
        return parameters.map(u => {
            const params = Array(blossom.degree).fill(u);
            return blossom.evaluate(...params);
        });
    }

    // Derivative computation using blossoming
    evaluateDerivative(
        points: V[],
        parameter: number,
        order: number
    ): V {
        const blossom = this.createInitialBlossom(points);
        const params = [
            ...Array(blossom.degree - order).fill(parameter),
            ...Array(order).fill(parameter + ε)
        ];
        
        return this.vectorSpace.scale(
            blossom.evaluate(...params),
            this.computeDerivativeScale(order)
        );
    }
}


// Simple example: Linear interpolation pyramid
interface Point {
    x: number;
    y: number;
}

// 1. Basic Level Structure
interface PyramidLevel {
    points: Point[];        // Data at this level
    parameter: number;      // Parameter value for interpolation
}

// 2. Simple Visualization
/*
Level 0:    P0 -------- P1 -------- P2 -------- P3    (Original points)
                  ↙           ↙           ↙
Level 1:        Q0 -------- Q1 -------- Q2           (Interpolated points)
                      ↙           ↙
Level 2:            R0 -------- R1                    (Further interpolated)
                          ↙
Level 3:                S0                            (Final point)
*/

// 3. Basic Implementation
class LinearInterpolationPyramid {
    // Linear interpolation between two points
    private interpolate(p1: Point, p2: Point, t: number): Point {
        return {
            x: p1.x * (1 - t) + p2.x * t,
            y: p1.y * (1 - t) + p2.y * t
        };
    }

    // Compute one level of the pyramid
    private computeNextLevel(level: PyramidLevel): PyramidLevel {
        const newPoints: Point[] = [];
        const points = level.points;
        const t = level.parameter;

        // Create new points by interpolating adjacent pairs
        for (let i = 0; i < points.length - 1; i++) {
            newPoints.push(
                this.interpolate(points[i], points[i + 1], t)
            );
        }

        return {
            points: newPoints,
            parameter: t
        };
    }

    // Main evaluation function
    evaluate(points: Point[], t: number): Point {
        let currentLevel: PyramidLevel = {
            points: points,
            parameter: t
        };

        // Build pyramid until we reach a single point
        while (currentLevel.points.length > 1) {
            currentLevel = this.computeNextLevel(currentLevel);
            
            // Optional: log for visualization
            console.log(
                "Level:", 
                currentLevel.points.map(p => `(${p.x}, ${p.y})`)
            );
        }

        return currentLevel.points[0];
    }
}

// 4. Usage Example
const pyramid = new LinearInterpolationPyramid();

const points: Point[] = [
    { x: 0, y: 0 },   // P0
    { x: 1, y: 2 },   // P1
    { x: 2, y: -1 },  // P2
    { x: 3, y: 1 }    // P3
];

const result = pyramid.evaluate(points, 0.5);
console.log("Final point:", result);

// Now with more dimensions !

// 1. Generic Vector representation
interface Vector {
    dimensions: number[];  // Array of values for each dimension
}

// 2. Enhanced Level Structure
interface PyramidLevel {
    points: Vector[];
    parameter: number;
}

// 3. Generic Implementation
class MultiDimensionalPyramid {
    // Generic interpolation between two vectors
    private interpolate(v1: Vector, v2: Vector, t: number): Vector {
        return {
            dimensions: v1.dimensions.map((value, index) => 
                value * (1 - t) + v2.dimensions[index] * t
            )
        };
    }

    // Compute next level (same logic, works for any dimension)
    private computeNextLevel(level: PyramidLevel): PyramidLevel {
        const newPoints: Vector[] = [];
        const points = level.points;
        const t = level.parameter;

        for (let i = 0; i < points.length - 1; i++) {
            newPoints.push(
                this.interpolate(points[i], points[i + 1], t)
            );
        }

        return { points: newPoints, parameter: t };
    }

    evaluate(points: Vector[], t: number): Vector {
        let currentLevel: PyramidLevel = {
            points: points,
            parameter: t
        };

        while (currentLevel.points.length > 1) {
            currentLevel = this.computeNextLevel(currentLevel);
            
            // Visualization helper
            this.visualizeLevel(currentLevel);
        }

        return currentLevel.points[0];
    }

    // Helper to visualize levels
    private visualizeLevel(level: PyramidLevel) {
        console.log("Level points:", level.points.map(p => 
            `(${p.dimensions.join(", ")})`
        ));
    }
}

// 4. Example Usage with different dimensions
function demonstrateMultiDimensional() {
    const pyramid = new MultiDimensionalPyramid();

    // 2D Example
    const points2D: Vector[] = [
        { dimensions: [0, 0] },    // (x, y)
        { dimensions: [1, 2] },
        { dimensions: [2, -1] },
        { dimensions: [3, 1] }
    ];

    console.log("2D Evaluation:");
    const result2D = pyramid.evaluate(points2D, 0.5);

    // 3D Example
    const points3D: Vector[] = [
        { dimensions: [0, 0, 0] }, // (x, y, z)
        { dimensions: [1, 2, 1] },
        { dimensions: [2, -1, 2] },
        { dimensions: [3, 1, 0] }
    ];

    console.log("\n3D Evaluation:");
    const result3D = pyramid.evaluate(points3D, 0.5);

    // 4D Example (e.g., space-time or color with alpha)
    const points4D: Vector[] = [
        { dimensions: [0, 0, 0, 1] }, // (x, y, z, w)
        { dimensions: [1, 2, 1, 0.5] },
        { dimensions: [2, -1, 2, 0.7] },
        { dimensions: [3, 1, 0, 1] }
    ];

    console.log("\n4D Evaluation:");
    const result4D = pyramid.evaluate(points4D, 0.5);
}

// 5. Specialized Visualization for different dimensions
class PyramidVisualizer {
    static visualize2D(points: Vector[]) {
        return points.map(p => 
            `(${p.dimensions[0]}, ${p.dimensions[1]})`
        ).join(" → ");
    }

    static visualize3D(points: Vector[]) {
        return points.map(p => 
            `(${p.dimensions[0]}, ${p.dimensions[1]}, ${p.dimensions[2]})`
        ).join(" → ");
    }

    static visualizeND(points: Vector[]) {
        return points.map(p => 
            `(${p.dimensions.join(", ")})`
        ).join(" → ");
    }
}

// 6. Example with specific applications
class ApplicationExamples {
    // Color interpolation (4D: R,G,B,A)
    static interpolateColors() {
        const pyramid = new MultiDimensionalPyramid();
        const colors: Vector[] = [
            { dimensions: [255, 0, 0, 1] },    // Red
            { dimensions: [0, 255, 0, 1] },    // Green
            { dimensions: [0, 0, 255, 1] }     // Blue
        ];

        return pyramid.evaluate(colors, 0.5);
    }

    // Animation keyframes (varying dimensions)
    static interpolateKeyframes() {
        const pyramid = new MultiDimensionalPyramid();
        const keyframes: Vector[] = [
            { dimensions: [0, 0, 0, 1, 0] },   // position_x, y, z, scale, rotation
            { dimensions: [10, 5, 2, 1.5, 45] },
            { dimensions: [20, 0, 0, 1, 90] }
        ];

        return pyramid.evaluate(keyframes, 0.5);


// Here's an example of pyramid computation for surfaces and volumes, starting with a simple bilinear/trilinear interpolation:

// 1. Core Structures
interface Vector {
    dimensions: number[];
}

interface GridLevel {
    points: Vector[][];      // 2D array for surface
    parameterU: number;      // Parameter in U direction
    parameterV: number;      // Parameter in V direction
}

interface VolumeLevel {
    points: Vector[][][];    // 3D array for volume
    parameterU: number;
    parameterV: number;
    parameterW: number;
}

// 2. Surface Pyramid Computation
class SurfacePyramid {
    // Interpolate in one direction (U or V)
    private interpolateRow(row: Vector[], t: number): Vector[] {
        const result: Vector[] = [];
        
        for (let i = 0; i < row.length - 1; i++) {
            result.push({
                dimensions: row[i].dimensions.map((value, dim) =>
                    value * (1 - t) + row[i + 1].dimensions[dim] * t
                )
            });
        }
        
        return result;
    }

    // Compute next level (alternating U and V directions)
    private computeNextLevel(
        level: GridLevel, 
        direction: 'U' | 'V'
    ): GridLevel {
        if (direction === 'U') {
            // Interpolate rows (U direction)
            const newPoints = level.points.map(row => 
                this.interpolateRow(row, level.parameterU)
            );

            return {
                points: newPoints,
                parameterU: level.parameterU,
                parameterV: level.parameterV
            };
        } else {
            // Interpolate columns (V direction)
            const newPoints: Vector[][] = [];
            for (let j = 0; j < level.points[0].length; j++) {
                const column = level.points.map(row => row[j]);
                newPoints.push(this.interpolateRow(column, level.parameterV));
            }

            return {
                points: newPoints,
                parameterU: level.parameterU,
                parameterV: level.parameterV
            };
        }
    }

    // Main evaluation function
    evaluate(
        points: Vector[][], 
        u: number, 
        v: number
    ): Vector {
        let currentLevel: GridLevel = {
            points: points,
            parameterU: u,
            parameterV: v
        };

        // Visualize initial grid
        console.log("Initial Grid:");
        this.visualizeGrid(currentLevel);

        // Alternate between U and V directions
        while (currentLevel.points.length > 1 || 
               currentLevel.points[0].length > 1) {
            
            if (currentLevel.points[0].length > 1) {
                currentLevel = this.computeNextLevel(currentLevel, 'U');
                console.log("After U interpolation:");
            } else {
                currentLevel = this.computeNextLevel(currentLevel, 'V');
                console.log("After V interpolation:");
            }

            this.visualizeGrid(currentLevel);
        }

        return currentLevel.points[0][0];
    }

    // Visualization helper
    private visualizeGrid(level: GridLevel) {
        level.points.forEach(row => {
            console.log(row.map(p => 
                `(${p.dimensions.join(",")})`
            ).join(" "));
        });
        console.log("");
    }
}

// 3. Volume Pyramid Computation
class VolumePyramid {
    // Interpolate in one direction (U, V, or W)
    private interpolate1D(points: Vector[], t: number): Vector[] {
        const result: Vector[] = [];
        
        for (let i = 0; i < points.length - 1; i++) {
            result.push({
                dimensions: points[i].dimensions.map((value, dim) =>
                    value * (1 - t) + points[i + 1].dimensions[dim] * t
                )
            });
        }
        
        return result;
    }

    // Compute next level for volume
    private computeNextLevel(
        level: VolumeLevel, 
        direction: 'U' | 'V' | 'W'
    ): VolumeLevel {
        switch (direction) {
            case 'U': {
                // Interpolate in U direction
                const newPoints = level.points.map(plane => 
                    plane.map(row => 
                        this.interpolate1D(row, level.parameterU)
                    )
                );
                return { ...level, points: newPoints };
            }
            case 'V': {
                // Interpolate in V direction
                const newPoints = level.points.map(plane => {
                    const newPlane: Vector[][] = [];
                    for (let j = 0; j < plane[0].length; j++) {
                        const column = plane.map(row => row[j]);
                        newPlane.push(
                            this.interpolate1D(column, level.parameterV)
                        );
                    }
                    return newPlane;
                });
                return { ...level, points: newPoints };
            }
            case 'W': {
                // Interpolate in W direction
                const newPoints: Vector[][][] = [];
                for (let i = 0; i < level.points[0].length; i++) {
                    const plane: Vector[][] = [];
                    for (let j = 0; j < level.points[0][0].length; j++) {
                        const line = level.points.map(p => p[i][j]);
                        plane.push(
                            this.interpolate1D(line, level.parameterW)
                        );
                    }
                    newPoints.push(plane);
                }
                return { ...level, points: newPoints };
            }
        }
    }

    // Example usage
    evaluate(
        points: Vector[][][], 
        u: number, 
        v: number, 
        w: number
    ): Vector {
        let currentLevel: VolumeLevel = {
            points,
            parameterU: u,
            parameterV: v,
            parameterW: w
        };

        // Reduce in each direction until single point
        while (currentLevel.points.length > 1 || 
               currentLevel.points[0].length > 1 || 
               currentLevel.points[0][0].length > 1) {
            
            if (currentLevel.points[0][0].length > 1) {
                currentLevel = this.computeNextLevel(currentLevel, 'U');
            } else if (currentLevel.points[0].length > 1) {
                currentLevel = this.computeNextLevel(currentLevel, 'V');
            } else {
                currentLevel = this.computeNextLevel(currentLevel, 'W');
            }

            this.visualizeVolume(currentLevel);
        }

        return currentLevel.points[0][0][0];
    }

    private visualizeVolume(level: VolumeLevel) {
        console.log("Volume Level:");
        level.points.forEach((plane, i) => {
            console.log(`Plane ${i}:`);
            plane.forEach(row => {
                console.log(row.map(p => 
                    `(${p.dimensions.join(",")})`
                ).join(" "));
            });
            console.log("");
        });
    }
}

// 4. Usage Examples
function demonstrateSurfaceInterpolation() {
    const surface = new SurfacePyramid();
    
    // Create a 3x3 control grid (2D surface)
    const controlPoints: Vector[][] = [
        [
            {dimensions: [0,0,0]}, 
            {dimensions: [1,0,1]}, 
            {dimensions: [2,0,0]}
        ],
        [
            {dimensions: [0,1,1]}, 
            {dimensions: [1,1,2]}, 
            {dimensions: [2,1,1]}
        ],
        [
            {dimensions: [0,2,0]}, 
            {dimensions: [1,2,1]}, 
            {dimensions: [2,2,0]}
        ]
    ];

    const result = surface.evaluate(controlPoints, 0.5, 0.5);
    console.log("Surface Result:", result);
}

function demonstrateVolumeInterpolation() {
    const volume = new VolumePyramid();
    
    // Create a 2x2x2 control grid (3D volume)
    const controlPoints: Vector[][][] = [
        [
            [
                {dimensions: [0,0,0]}, 
                {dimensions: [1,0,0]}
            ],
            [
                {dimensions: [0,1,0]}, 
                {dimensions: [1,1,0]}
            ]
        ],
        [
            [
                {dimensions: [0,0,1]}, 
                {dimensions: [1,0,1]}
            ],
            [
                {dimensions: [0,1,1]}, 
                {dimensions: [1,1,1]}
            ]
        ]
    ];

    const result = volume.evaluate(controlPoints, 0.5, 0.5, 0.5);
    console.log("Volume Result:", result);
}


// Here's a clearer abstraction focusing on the dimensional reduction concept:

// 1. Core Concept: Dimensional Structure
interface DimensionalStructure<T> {
    data: T;                     // Actual data (points, grid, volume)
    dimensions: number[];        // Size in each dimension
    visualize(): string;         // How to display this structure
}

// 2. Specific Dimensional Structures
class Point implements DimensionalStructure<number[]> {
    constructor(
        public data: number[],   // coordinates
        public dimensions: number[] = [1]
    ) {}

    visualize(): string {
        return `(${this.data.join(", ")})`;
    }
}

class Grid implements DimensionalStructure<Point[][]> {
    constructor(
        public data: Point[][],
        public dimensions: number[] = [data.length, data[0].length]
    ) {}

    visualize(): string {
        return this.data.map(row =>
            row.map(p => p.visualize()).join(" ")
        ).join("\n");
    }
}

class Volume implements DimensionalStructure<Point[][][]> {
    constructor(
        public data: Point[][][],
        public dimensions: number[] = [
            data.length, 
            data[0].length, 
            data[0][0].length
        ]
    ) {}

    visualize(): string {
        return this.data.map((plane, i) =>
            `Layer ${i}:\n` + new Grid(plane).visualize()
        ).join("\n\n");
    }
}

// 3. Dimensional Reduction Step
interface ReductionStep<T extends DimensionalStructure<any>> {
    reduce(structure: T, parameter: number): T;
    dimension: number;  // Which dimension we're reducing
}

// 4. Pyramid Process showing dimensional reduction
class DimensionalPyramid<T extends DimensionalStructure<any>> {
    constructor(private reductionSteps: ReductionStep<T>[]) {}

    evaluate(initial: T, parameters: number[]): T {
        let current = initial;
        console.log("\nInitial Structure:");
        console.log(current.visualize());

        // Apply each reduction step
        this.reductionSteps.forEach((step, index) => {
            console.log(`\nReducing dimension ${step.dimension}:`);
            console.log(`Before: ${current.dimensions.join(" × ")}`);
            
            current = step.reduce(current, parameters[index]);
            
            console.log(`After: ${current.dimensions.join(" × ")}`);
            console.log(current.visualize());
        });

        return current;
    }
}

// 5. Example: Surface Reduction
class SurfaceReduction {
    static uDirection: ReductionStep<Grid> = {
        dimension: 0,
        reduce(grid: Grid, u: number): Grid {
            const newData = grid.data.map(row => {
                const newRow: Point[] = [];
                for (let i = 0; i < row.length - 1; i++) {
                    newRow.push(new Point(
                        row[i].data.map((val, j) => 
                            val * (1-u) + row[i+1].data[j] * u
                        )
                    ));
                }
                return newRow;
            });
            return new Grid(newData);
        }
    };

    static vDirection: ReductionStep<Grid> = {
        dimension: 1,
        reduce(grid: Grid, v: number): Grid {
            const newData: Point[][] = [];
            for (let i = 0; i < grid.data.length - 1; i++) {
                const newRow: Point[] = [];
                for (let j = 0; j < grid.data[0].length; j++) {
                    newRow.push(new Point(
                        grid.data[i][j].data.map((val, k) =>
                            val * (1-v) + grid.data[i+1][j].data[k] * v
                        )
                    ));
                }
                newData.push(newRow);
            }
            return new Grid(newData);
        }
    };
}

// 6. Usage Example
function demonstrateSurfaceReduction() {
    // Create a 3×3 control grid
    const controlPoints: Point[][] = [
        [
            new Point([0, 0, 0]),
            new Point([1, 0, 1]),
            new Point([2, 0, 0])
        ],
        [
            new Point([0, 1, 1]),
            new Point([1, 1, 2]),
            new Point([2, 1, 1])
        ],
        [
            new Point([0, 2, 0]),
            new Point([1, 2, 1]),
            new Point([2, 2, 0])
        ]
    ];

    const initialGrid = new Grid(controlPoints);
    
    // Create pyramid with reduction steps
    const pyramid = new DimensionalPyramid<Grid>([
        SurfaceReduction.uDirection,
        SurfaceReduction.vDirection
    ]);

    // Evaluate at (u,v) = (0.5, 0.5)
    console.log("Surface Reduction Process:");
    const result = pyramid.evaluate(initialGrid, [0.5, 0.5]);
}

// 7. Example Output will look like:
/*
Initial Structure: (3 × 3 grid)
(0,0,0) (1,0,1) (2,0,0)
(0,1,1) (1,1,2) (2,1,1)
(0,2,0) (1,2,1) (2,2,0)

Reducing U dimension:
Before: 3 × 3
After: 3 × 2
(0.5,0,0.5) (1.5,0,0.5)
(0.5,1,1.5) (1.5,1,1.5)
(0.5,2,0.5) (1.5,2,0.5)

Reducing V dimension:
Before: 3 × 2
After: 2 × 2
(0.5,1,1) (1.5,1,1)
*/

// 8. More Complex Example: Volume Reduction
class VolumeReduction {
    static wDirection: ReductionStep<Volume> = {
        dimension: 2,
        reduce(volume: Volume, w: number): Volume {
            // Similar reduction logic for w direction
            // ...
            return new Volume(/* ... */);
        }
    };
}


