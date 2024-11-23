// Core abstractions and types
type Dimension = number;
type ParameterValue = number[];
type Index = number[];

// Base interfaces
interface Point {
    coordinates: number[];
    weight: number;
}

interface ParametricDomain {
    dimension: number;
    evaluate(parameters: ParameterValue): boolean;
    project(point: Point): ParameterValue;
}

// Abstract basis function interface
interface BasisFunction {
    evaluate(parameters: ParameterValue): number;
    support: ParametricDomain;
}

// Core space definition
abstract class ParametricSpace {
    abstract dimension: number;
    abstract domain: ParametricDomain;
    
    abstract evaluateBasis(
        parameters: ParameterValue,
        index: Index
    ): number;
    
    abstract getActiveBasis(
        parameters: ParameterValue
    ): Map<string, BasisFunction>;
}

// Base class for all spline types
abstract class SplineBase {
    protected space: ParametricSpace;
    protected controlPoints: Map<string, Point>;

    constructor(space: ParametricSpace) {
        this.space = space;
        this.controlPoints = new Map();
    }

    abstract evaluate(parameters: ParameterValue): Point;
}

// B-spline specific implementations
class BSplineBasis implements BasisFunction {
    private degree: number;
    private knots: number[];
    private index: number;

    constructor(degree: number, knots: number[], index: number) {
        this.degree = degree;
        this.knots = knots;
        this.index = index;
    }

    evaluate(parameters: ParameterValue): number {
        return this.deBoorCox(parameters[0], this.index);
    }

    private deBoorCox(t: number, i: number): number {
        // Implementation of de Boor-Cox algorithm
        return 0; // Placeholder
    }
}

// Tensor product space for standard B-splines
class TensorProductSpace extends ParametricSpace {
    private degrees: number[];
    private knotVectors: number[][];

    constructor(degrees: number[], knotVectors: number[][]) {
        super();
        this.degrees = degrees;
        this.knotVectors = knotVectors;
    }

    get dimension(): number {
        return this.degrees.length;
    }

        let result = 1;
        for (let i = 0; i < this.dimension; i++) {
            const basis = new BSplineBasis(
                this.degrees[i],
                this.knotVectors[i],
                index[i]
            );
            result *= basis.evaluate([parameters[i]]);
        }
        return result;
    }
}

// Extension for research on new spline types
namespace Research {
    // Triangular domain implementation
    class TriangularDomain implements ParametricDomain {
        dimension = 3;

        evaluate(parameters: ParameterValue): boolean {
            const [u, v, w] = parameters;
            return Math.abs(u + v + w - 1.0) < 1e-10;
        }

        project(point: Point): ParameterValue {
            // Project point onto triangular domain
            return [0, 0, 0]; // Placeholder
        }
    }

    // Base class for T-spline research
    abstract class TSplineBase extends SplineBase {
        protected tJunctions: Map<string, ParameterValue>;

        constructor(space: ParametricSpace) {
            super(space);
            this.tJunctions = new Map();
        }

        abstract computeLocalKnotStructure(
            parameters: ParameterValue
        ): Map<string, number[]>;
    }

    // Experimental triangular T-spline implementation
    class TriangularTSpline extends TSplineBase {
        private triangularSpace: TriangularBSplineSpace;

        constructor(degree: number) {
            const space = new TriangularBSplineSpace(degree);
            super(space);
            this.triangularSpace = space;
        }

        evaluate(parameters: ParameterValue): Point {
            // Implementation using triangular basis
            return { coordinates: [], weight: 0 }; // Placeholder
        }

        computeLocalKnotStructure(
            parameters: ParameterValue
        ): Map<string, number[]> {
            // Compute local knot structure considering T-junctions
            return new Map();
        }
    }

    // Research utilities
    class SplineAnalyzer {
        analyzeLinearIndependence(space: ParametricSpace): boolean {
            // Analyze linear independence of basis functions
            return true;
        }

        analyzeContinuity(spline: SplineBase): number {
            // Analyze continuity order
            return 0;
        }

        computeApproximationError(
            spline: SplineBase,
            targetFunction: (p: ParameterValue) => Point
        ): number {
            // Compute approximation error
            return 0;
        }
    }

    // Experimental features
    class AdaptiveRefinement {
        private errorTolerance: number;

        constructor(errorTolerance: number) {
            this.errorTolerance = errorTolerance;
        }

        refine(spline: TSplineBase, parameters: ParameterValue): void {
            // Implement adaptive refinement strategy
        }
    }
}

// Example of framework extension for new research
namespace ExperimentalSplines {
    // New basis function type
    class ExperimentalBasis implements BasisFunction {
        evaluate(parameters: ParameterValue): number {
            // New basis function definition
            return 0;
        }
    }

    // New spline space
    class ExperimentalSpace extends ParametricSpace {
        dimension = 0;
        domain!: ParametricDomain;

        evaluateBasis(parameters: ParameterValue, index: Index): number {
            // Experimental basis evaluation
            return 0;
        }
    }

    // Analysis tools for new spline types
    class ExperimentalAnalyzer {
        analyzeProperties(space: ExperimentalSpace): void {
            // Analyze mathematical properties
        }
    }
}

// Usage example for research
function researchExample() {
    // Standard B-spline setup
    const degrees = [3, 3];
    const knotVectors = [
        [0, 0, 0, 0, 1, 1, 1, 1],
        [0, 0, 0, 0, 1, 1, 1, 1]
    ];
    const standardSpace = new TensorProductSpace(degrees, knotVectors);

    // Experimental setup
    const experimentalSpace = new Research.TriangularTSpline(3);
    const analyzer = new Research.SplineAnalyzer();

    // Analysis
    const isLinearIndependent = analyzer.analyzeLinearIndependence(standardSpace);
    const continuityOrder = analyzer.analyzeContinuity(experimentalSpace);

    // Adaptive refinement experiment
    const refinement = new Research.AdaptiveRefinement(1e-6);
    refinement.refine(experimentalSpace, [0.5, 0.3, 0.2]);
}


///////////////////////////

/*

1. Abstraction Layers :

- Clear separation between topology and geometry

- Flexible knot structure handling

- Support for different domain types

2. Extension Points :

- Local refinement capabilities

- Custom basis functions

- Different mesh topologies

3. Analysis Tools :

- Linear independence checking

- Continuity analysis

- Error estimation

4. Research Support :

- Experimental implementations

- Analysis tools

- Refinement strategies


*/

// Core types
type Dimension = number;
type ParameterValue = number[];
type Index = number[];

// Fundamental interfaces
interface Point {
    coordinates: number[];
    weight: number;
}

// Domain abstractions
interface ParametricDomain {
    dimension: number;
    evaluate(parameters: ParameterValue): boolean;
    project(point: Point): ParameterValue;
}

// Knot structure abstractions
interface KnotStructure {
    getDegree(dimension: Dimension): number;
    getKnots(dimension: Dimension, index: Index): number[];
    isValid(): boolean;
}

// Local knot structure for adaptive refinement
interface LocalKnotStructure extends KnotStructure {
    getLocalKnots(parameters: ParameterValue): Map<string, number[]>;
    computeSupport(): ParametricDomain;
}

// Basis function abstraction
interface BasisFunction {
    evaluate(parameters: ParameterValue): number;
    support: ParametricDomain;
    knotStructure: KnotStructure;
}

// Space abstraction
abstract class ParametricSpace {
    abstract dimension: number;
    abstract domain: ParametricDomain;
    
    abstract evaluateBasis(
        parameters: ParameterValue,
        index: Index
    ): number;
    
    abstract getActiveBasis(
        parameters: ParameterValue
    ): Map<string, BasisFunction>;
}

// Mesh topology abstractions
interface MeshTopology {
    getNeighbors(index: Index): Index[];
    isValid(): boolean;
    refine(region: ParametricDomain): void;
}

// Base spline abstraction
abstract class SplineBase {
    protected space: ParametricSpace;
    protected controlPoints: Map<string, Point>;
    protected topology: MeshTopology;

    constructor(
        space: ParametricSpace,
        topology: MeshTopology
    ) {
        this.space = space;
        this.topology = topology;
        this.controlPoints = new Map();
    }

    abstract evaluate(parameters: ParameterValue): Point;
}

// Analysis framework
interface SplineAnalysis {
    checkLinearIndependence(): boolean;
    computeContinuity(): number;
    validatePartitionOfUnity(): boolean;
    computeApproximationError(
        targetFunction: (p: ParameterValue) => Point
    ): number;
}

// Refinement strategies
interface RefinementStrategy {
    refine(
        spline: SplineBase,
        region: ParametricDomain,
        errorTolerance: number
    ): void;
    
    validate(): boolean;
}

// Implementation examples
namespace Implementations {
    // Standard tensor-product B-spline implementation
    class TensorProductKnots implements KnotStructure {
        private degrees: number[];
        private knotVectors: number[][];

        constructor(degrees: number[], knotVectors: number[][]) {
            this.degrees = degrees;
            this.knotVectors = knotVectors;
        }

        getDegree(dimension: Dimension): number {
            return this.degrees[dimension];
        }

        getKnots(dimension: Dimension, index: Index): number[] {
            return this.knotVectors[dimension];
        }

        isValid(): boolean {
            // Validate knot vector properties
            return true;
        }
    }

    // Local refinement structure
    class LocalRefinementKnots implements LocalKnotStructure {
        private baseStructure: KnotStructure;
        private localModifications: Map<string, number[]>;

        constructor(
            baseStructure: KnotStructure,
            modifications: Map<string, number[]>
        ) {
            this.baseStructure = baseStructure;
            this.localModifications = modifications;
        }

        // Implementation of LocalKnotStructure interface
    }

    // Triangular domain implementation
    class TriangularDomain implements ParametricDomain {
        dimension = 3;

        evaluate(parameters: ParameterValue): boolean {
            const [u, v, w] = parameters;
            return Math.abs(u + v + w - 1.0) < 1e-10;
        }

        project(point: Point): ParameterValue {
            // Project point onto triangular domain
            return [0, 0, 0];
        }
    }
}

// Research extensions
namespace Research {
    // Experimental mesh topology
    class AdaptiveMeshTopology implements MeshTopology {
        private nodes: Map<string, Index>;
        private connections: Map<string, Set<string>>;

        constructor() {
            this.nodes = new Map();
            this.connections = new Map();
        }

        // Implementation of MeshTopology interface
    }

    // Analysis tools
    class SplineAnalyzer implements SplineAnalysis {
        private spline: SplineBase;

        constructor(spline: SplineBase) {
            this.spline = spline;
        }

        // Implementation of SplineAnalysis interface
    }

    // Refinement implementation
    class AdaptiveRefinement implements RefinementStrategy {
        private errorTolerance: number;
        private maxLevel: number;

        constructor(errorTolerance: number, maxLevel: number) {
            this.errorTolerance = errorTolerance;
            this.maxLevel = maxLevel;
        }

        // Implementation of RefinementStrategy interface
    }
}

// Example usage
function researchExample() {
    // Create base structures
    const knots = new Implementations.TensorProductKnots(
        [3, 3],
        [
            [0, 0, 0, 0, 1, 1, 1, 1],
            [0, 0, 0, 0, 1, 1, 1, 1]
        ]
    );

    const topology = new Research.AdaptiveMeshTopology();
    
    // Analysis setup
    const analyzer = new Research.SplineAnalyzer(/* spline instance */);
    
    // Refinement
    const refinement = new Research.AdaptiveRefinement(1e-6, 5);
}

/// Discussion over the parametric domain

interface Point {
    coordinates: number[];
    weight: number;
}

interface ParametricDomain {
    dimension: number;
    evaluate(parameters: ParameterValue): boolean;
    project(point: Point): ParameterValue;
    // New methods for better domain handling
    getClosestPoint(point: Point): ParameterValue;
    computeGradient(parameters: ParameterValue): number[][];
    computeBoundary(): BoundaryDescription;
}

interface BoundaryDescription {
    type: 'rectangular' | 'triangular' | 'simplex' | 'custom';
    constraints: ParametricConstraint[];
}

interface ParametricConstraint {
    evaluate(parameters: ParameterValue): number;
    gradient(parameters: ParameterValue): number[];
}

// Rectangular domain implementation
class RectangularDomain implements ParametricDomain {
    private ranges: [number, number][];

    constructor(ranges: [number, number][]) {
        this.ranges = ranges;
    }

    get dimension(): number {
        return this.ranges.length;
    }

    evaluate(parameters: ParameterValue): boolean {
        if (parameters.length !== this.dimension) return false;
        
        return parameters.every((p, i) => 
            p >= this.ranges[i][0] && p <= this.ranges[i][1]
        );
    }

    project(point: Point): ParameterValue {
        // Simple clamping for rectangular domains
        return this.ranges.map((range, i) => {
            const param = point.coordinates[i];
            return Math.max(range[0], Math.min(range[1], param));
        });
    }

    getClosestPoint(point: Point): ParameterValue {
        return this.project(point); // For rectangular, same as project
    }

    computeGradient(parameters: ParameterValue): number[][] {
        // Return identity matrix for rectangular domain
        return Array(this.dimension).fill(0).map((_, i) => 
            Array(this.dimension).fill(0).map((_, j) => i === j ? 1 : 0)
        );
    }

    computeBoundary(): BoundaryDescription {
        const constraints: ParametricConstraint[] = [];
        
        // Add min/max constraints for each dimension
        this.ranges.forEach(([min, max], dim) => {
            constraints.push({
                evaluate: (p: ParameterValue) => p[dim] - min,
                gradient: (p: ParameterValue) => 
                    Array(this.dimension).fill(0).map((_, i) => i === dim ? 1 : 0)
            });
            constraints.push({
                evaluate: (p: ParameterValue) => max - p[dim],
                gradient: (p: ParameterValue) => 
                    Array(this.dimension).fill(0).map((_, i) => i === dim ? -1 : 0)
            });
        });

        return {
            type: 'rectangular',
            constraints
        };
    }
}

// Triangular domain implementation
class TriangularDomain implements ParametricDomain {
    get dimension(): number {
        return 3; // Barycentric coordinates
    }

    evaluate(parameters: ParameterValue): boolean {
        const [u, v, w] = parameters;
        return (
            u >= 0 && v >= 0 && w >= 0 &&
            Math.abs(u + v + w - 1.0) < 1e-10
        );
    }

    project(point: Point): ParameterValue {
        // Project to barycentric coordinates and ensure sum = 1
        return this.projectToBarycentric(point);
    }

    private projectToBarycentric(point: Point): ParameterValue {
        // Convert to barycentric coordinates
        let [u, v, w] = this.computeInitialBarycentric(point);

        // Project to valid barycentric coordinates
        if (u < 0 || v < 0 || w < 0) {
            [u, v, w] = this.projectToTriangleBoundary(u, v, w);
        }

        // Ensure sum = 1
        const sum = u + v + w;
        return [u/sum, v/sum, w/sum];
    }

    private computeInitialBarycentric(point: Point): [number, number, number] {
        // Implement conversion to barycentric coordinates
        // This is a simplified version
        const [x, y] = point.coordinates;
        const u = x;
        const v = y;
        const w = 1 - u - v;
        return [u, v, w];
    }

        // Project point to nearest edge or vertex
        if (u < 0) {
            // Project to v-w edge
            const sum = v + w;
            return [0, v/sum, w/sum];
        }
        // Similar for v < 0 and w < 0
        return [u, v, w];
    }

    getClosestPoint(point: Point): ParameterValue {
        // Implement efficient closest point computation
        // Could use geometric algorithms for triangle
        return this.project(point);
    }

    computeGradient(parameters: ParameterValue): number[][] {
        // Compute gradient of barycentric coordinates
        const [u, v, w] = parameters;
        return [
            [1, 0, -1],
            [0, 1, -1],
            [-1, -1, 1]
        ];
    }

    computeBoundary(): BoundaryDescription {
        const constraints: ParametricConstraint[] = [
            // Sum to 1 constraint
            {
                evaluate: (p: ParameterValue) => p[0] + p[1] + p[2] - 1,
                gradient: (_: ParameterValue) => [1, 1, 1]
            },
            // Non-negativity constraints
            {
                evaluate: (p: ParameterValue) => p[0],
                gradient: (_: ParameterValue) => [1, 0, 0]
            },
            {
                evaluate: (p: ParameterValue) => p[1],
                gradient: (_: ParameterValue) => [0, 1, 0]
            },
            {
                evaluate: (p: ParameterValue) => p[2],
                gradient: (_: ParameterValue) => [0, 0, 1]
            }
        ];

        return {
            type: 'triangular',
            constraints
        };
    }
}

// Generic simplex domain
class SimplexDomain implements ParametricDomain {
    private dim: number;

    constructor(dimension: number) {
        this.dim = dimension;
    }

    get dimension(): number {
        return this.dim;
    }

    evaluate(parameters: ParameterValue): boolean {
        if (parameters.length !== this.dim) return false;
        
        const sum = parameters.reduce((a, b) => a + b, 0);
        return (
            parameters.every(p => p >= 0) &&
            Math.abs(sum - 1.0) < 1e-10
        );
    }

    project(point: Point): ParameterValue {
        // Project to n-dimensional simplex
        return this.projectToSimplex(point.coordinates);
    }

    private projectToSimplex(coords: number[]): ParameterValue {
        // Implement projection to n-dimensional simplex
        // This is a complex operation requiring optimization
        return this.simplexProjection(coords);
    }

    private simplexProjection(x: number[]): number[] {
        // Implement efficient simplex projection
        // Could use various algorithms (e.g., sorting-based)
        return [];  // Placeholder
    }

    // ... other methods
}

// Optimization utilities for projection
class ProjectionOptimizer {
    static projectToConstraints(
        point: Point,
        constraints: ParametricConstraint[]
    ): ParameterValue {
        // Implement constrained optimization
        // Could use various methods:
        // - Sequential quadratic programming
        // - Interior point methods
        // - Active set methods
        return [];  // Placeholder
    }
}

/// example project method

interface Point3D {
    x: number;
    y: number;
    z: number;
}

interface ParametricSurface {
    evaluate(u: number, v: number): Point3D;
    normal(u: number, v: number): Point3D;
    domain: ParametricDomain;
}

class ProjectionUtils {
    // Project point to parametric domain using Newton iteration
    static projectPointToSurface(
        surface: ParametricSurface,
        point: Point3D,
        tolerance: number = 1e-6,
        maxIterations: number = 20
    ): ParameterValue {
        // Initial guess (could be improved)
        
        for (let iter = 0; iter < maxIterations; iter++) {
            // Current point on surface
            const surfPoint = surface.evaluate(params[0], params[1]);
            
            // Compute distance vector
            const distance = {
                x: point.x - surfPoint.x,
                y: point.y - surfPoint.y,
                z: point.z - surfPoint.z
            };

            // Compute surface derivatives
            const [Su, Sv] = surface.derivatives(params[0], params[1]);
            
            // Form the system to solve
            // [Su·Su  Su·Sv] [Δu] = [Su·d]
            // [Sv·Su  Sv·Sv] [Δv] = [Sv·d]
            const a11 = dot(Su, Su);
            const a12 = dot(Su, Sv);
            const a21 = a12;
            const a22 = dot(Sv, Sv);
            
            const b1 = dot(Su, distance);
            const b2 = dot(Sv, distance);
            
            // Solve 2x2 system
            const det = a11 * a22 - a12 * a21;
            const du = (a22 * b1 - a12 * b2) / det;
            const dv = (-a21 * b1 + a11 * b2) / det;
            
            // Update parameters
            params[0] += du;
            params[1] += dv;
            
            // Project back to domain if needed
            params = surface.domain.project(params);
            
            // Check convergence
            if (Math.sqrt(du*du + dv*dv) < tolerance) {
                break;
            }
        }
        
        return params;
    }
}

// Example for a NURBS surface
class NURBSSurface implements ParametricSurface {
    domain: RectangularDomain;
    
    constructor() {
        this.domain = new RectangularDomain([
            [0, 1],
            [0, 1]
        ]);
    }

    evaluate(u: number, v: number): Point3D {
        // NURBS surface evaluation
        return { x: 0, y: 0, z: 0 }; // Placeholder
    }

    derivatives(u: number, v: number): [Point3D, Point3D] {
        // Compute partial derivatives ∂S/∂u and ∂S/∂v
        return [
            { x: 0, y: 0, z: 0 },
            { x: 0, y: 0, z: 0 }
        ]; // Placeholder
    }
}

// Example for cylindrical surface
class CylindricalSurface implements ParametricSurface {
    private radius: number;
    private height: number;
    domain: RectangularDomain;

    constructor(radius: number, height: number) {
        this.radius = radius;
        this.height = height;
        this.domain = new RectangularDomain([
            [0, 2 * Math.PI], // Angular parameter
            [0, height]       // Height parameter
        ]);
    }

    evaluate(u: number, v: number): Point3D {
        return {
            x: this.radius * Math.cos(u),
            y: this.radius * Math.sin(u),
            z: v
        };
    }

    derivatives(u: number, v: number): [Point3D, Point3D] {
        // Partial derivatives
        const Su = {
            x: -this.radius * Math.sin(u),
            y: this.radius * Math.cos(u),
            z: 0
        };
        
        const Sv = {
            x: 0,
            y: 0,
            z: 1
        };

        return [Su, Sv];
    }

    // Specialized projection for cylinder (more efficient than Newton)
    projectPoint(point: Point3D): ParameterValue {
        // Angular parameter
        const u = Math.atan2(point.y, point.x);
        
        // Height parameter (clamped to domain)
        
        return [u, v];
    }
}

// Example usage
function example() {
    // Create a cylindrical surface
    const cylinder = new CylindricalSurface(1.0, 2.0);
    
    // Point to project
    const point: Point3D = { x: 2.0, y: 1.0, z: 1.5 };
    
    // Project using general method
    const params = ProjectionUtils.projectPointToSurface(cylinder, point);
    
    // Or use specialized method
    const specializedParams = cylinder.projectPoint(point);
    
    // Evaluate projected point on surface
    const surfacePoint = cylinder.evaluate(params[0], params[1]);
}

// Helper function
function dot(a: Point3D, b: Point3D): number {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}


///// These are the core abstractions for spline-based geometric modeling. Let me break down these fundamental concepts and their relationships:

// 1. Point representation
interface Point {
    coordinates: number[];
    weight: number;  // For rational geometries (NURBS)
}

// 2. Parametric Domain - defines valid parameter space
interface ParametricDomain {
    dimension: number;
    evaluate(parameters: ParameterValue): boolean;
    project(point: Point): ParameterValue;
}

// 3. Basis Function - fundamental building block
interface BasisFunction {
    evaluate(parameters: ParameterValue): number;
    support: ParametricDomain;  // Where the basis is non-zero
}

// 4. Parametric Space - manages basis functions
abstract class ParametricSpace {
    abstract dimension: number;
    abstract domain: ParametricDomain;
    
    // Get value of specific basis function
    abstract evaluateBasis(
        parameters: ParameterValue,
        index: Index
    ): number;
    
    // Get all non-zero basis functions at parameters
    abstract getActiveBasis(
        parameters: ParameterValue
    ): Map<string, BasisFunction>;
}

// Example implementation for B-splines
class BSplineSpace extends ParametricSpace {
    private degrees: number[];
    private knots: number[][];

    constructor(degrees: number[], knots: number[][]) {
        super();
        this.degrees = degrees;
        this.knots = knots;
    }

    get dimension(): number {
        return this.degrees.length;
    }

    get domain(): ParametricDomain {
        return new RectangularDomain(
            this.knots.map(k => [k[0], k[k.length-1]])
        );
    }

    evaluateBasis(parameters: ParameterValue, index: Index): number {
        let result = 1.0;
        for (let i = 0; i < this.dimension; i++) {
            const basis = new BSplineBasis(
                this.degrees[i],
                this.knots[i],
                index[i]
            );
            result *= basis.evaluate([parameters[i]]);
        }
        return result;
    }

    getActiveBasis(parameters: ParameterValue): Map<string, BasisFunction> {
        const active = new Map<string, BasisFunction>();
        
        // Find active spans in each dimension
        const spans = parameters.map((p, i) => 
            this.findSpan(p, this.degrees[i], this.knots[i])
        );

        // Create basis functions for active region
        this.forEachActiveIndex(spans, this.degrees, (index) => {
            const basis = this.createBasisFunction(index);
            active.set(index.join(','), basis);
        });

        return active;
    }

    private createBasisFunction(index: Index): BasisFunction {
        return {
            evaluate: (params: ParameterValue) => 
                this.evaluateBasis(params, index),
            support: this.computeSupport(index)
        };
    }
}

// Usage example
function example() {
    // Create a B-spline space
    const degrees = [3, 3];  // Bicubic
    const knots = [
        [0, 0, 0, 0, 1, 1, 1, 1],
        [0, 0, 0, 0, 1, 1, 1, 1]
    ];
    const space = new BSplineSpace(degrees, knots);

    // Evaluate at a point
    const params = [0.5, 0.5];
    
    // Get active basis functions
    const activeBasis = space.getActiveBasis(params);
    
    // Evaluate geometry
    let point: Point = { coordinates: [0, 0, 0], weight: 0 };
    
    for (const [index, basis] of activeBasis) {
        const controlPoint = getControlPoint(index);
        const value = basis.evaluate(params);
        
        // Add contribution
        point.coordinates = point.coordinates.map(
            (c, i) => c + value * controlPoint.coordinates[i]
        );
        point.weight += value * controlPoint.weight;
    }
}

// Simplex Spline !!!!!!!!!!!!!!!!!!!!!!!!!

// Barycentric coordinates for triangular domain
interface BarycentricCoordinate {
    coordinates: [number, number, number]; // u + v + w = 1
}

class TriangularBSplineSpace extends ParametricSpace {
    private degree: number;
    private mesh: TriangularMesh;

    constructor(degree: number, mesh: TriangularMesh) {
        super();
        this.degree = degree;
        this.mesh = mesh;
    }

    get dimension(): number {
        return 3; // Barycentric coordinates
    }

    get domain(): ParametricDomain {
        return new TriangularDomain();
    }

    // Get active basis functions for given parameters
    getActiveBasis(parameters: ParameterValue): Map<string, BasisFunction> {
        const active = new Map<string, BasisFunction>();
        
        // Convert to barycentric coordinates
        const [u, v, w] = parameters;
        
        // Find containing triangle in mesh
        const triangle = this.mesh.findContainingTriangle(u, v, w);
        
        // Get neighboring triangles up to degree
        const neighborhood = this.mesh.getNeighborhood(triangle, this.degree);
        
        // Create basis functions for each control point in neighborhood
        for (const vertex of neighborhood.vertices) {
            const basis = this.createBasisFunction(vertex, triangle);
            active.set(vertex.index.toString(), basis);
        }
        
        return active;
    }

    evaluateBasis(
        parameters: ParameterValue,
        index: Index
    ): number {
        // Evaluate triangular B-spline basis function
        return this.evaluateTriangularBasis(
            parameters,
            index,
            this.degree
        );
    }

    private evaluateTriangularBasis(
        parameters: ParameterValue,
        index: Index,
        degree: number
    ): number {
        if (degree === 0) {
            // Base case: characteristic function of triangle
            return this.evaluateCharacteristicFunction(parameters, index);
        }

        // Recursive case using de Casteljau-like algorithm
        let result = 0;
        const weights = this.computeBarycentric(parameters, index);
        
        for (let i = 0; i <= degree; i++) {
            result += weights[i] * this.evaluateTriangularBasis(
                parameters,
                this.getSubIndex(index, i),
                degree - 1
            );
        }
        
        return result;
    }
}

// Helper classes for triangular mesh structure
class TriangularMesh {
    vertices: Vertex[];
    triangles: Triangle[];
    
    findContainingTriangle(u: number, v: number, w: number): Triangle {
        // Find triangle containing point using barycentric coordinates
        for (const triangle of this.triangles) {
            if (triangle.contains([u, v, w])) {
                return triangle;
            }
        }
        throw new Error("Point not in mesh");
    }

    getNeighborhood(
        center: Triangle,
        radius: number
    ): { vertices: Vertex[], triangles: Triangle[] } {
        const vertices = new Set<Vertex>();
        const triangles = new Set<Triangle>();
        
        // Start with center triangle
        let current = new Set([center]);
        
        // Expand neighborhood by radius steps
        for (let i = 0; i < radius; i++) {
            const next = new Set<Triangle>();
            
            // Add neighbors of current triangles
            for (const triangle of current) {
                for (const neighbor of triangle.neighbors) {
                    next.add(neighbor);
                    
                    // Collect vertices
                    for (const vertex of neighbor.vertices) {
                        vertices.add(vertex);
                    }
                }
                triangles.add(triangle);
            }
            
            current = next;
        }
        
        return {
            vertices: Array.from(vertices),
            triangles: Array.from(triangles)
        };
    }
}

class Triangle {
    vertices: Vertex[];
    neighbors: Triangle[];
    
    contains(barycentric: [number, number, number]): boolean {
        const [u, v, w] = barycentric;
        return (
            u >= 0 && v >= 0 && w >= 0 &&
            Math.abs(u + v + w - 1.0) < 1e-10
        );
    }
}

class Vertex {
    index: number;
    position: Point;
    adjacentTriangles: Triangle[];
}

// Example usage
function triangularExample() {
    // Create mesh structure
    const mesh = new TriangularMesh(/* ... */);
    
    // Create triangular B-spline space
    const space = new TriangularBSplineSpace(3, mesh); // Cubic
    
    // Point to evaluate
    const params: ParameterValue = [0.3, 0.3, 0.4]; // Barycentric
    
    // Get active basis functions
    const activeBasis = space.getActiveBasis(params);
    
    // Evaluate geometry
    let point: Point = { coordinates: [0, 0, 0], weight: 0 };
    
    for (const [index, basis] of activeBasis) {
        const controlPoint = getControlPoint(index);
        const value = basis.evaluate(params);
        
        // Add contribution
        point.coordinates = point.coordinates.map(
            (c, i) => c + value * controlPoint.coordinates[i]
        );
        point.weight += value * controlPoint.weight;
    }
}

/// Blossoming is indeed a fundamental concept that could be incorporated into the framework. It provides a more elegant and unified way to handle basis functions and their evaluation. Here's how we could extend the framework:

// Extended basis function interface with blossoming
interface BasisFunction {
    evaluate(parameters: ParameterValue): number;
    support: ParametricDomain;
    
    // Blossom evaluation
    blossom(parameters: ParameterValue[]): number;
    
    // Degree information (needed for blossoming)
    degree: number;
}

// Polar form representation
interface PolarForm {
    evaluate(...parameters: ParameterValue[]): number;
    degree: number;
    isSymmetric: boolean;
}

// B-spline basis with blossoming
class BSplineBasisWithBlossom implements BasisFunction {
    private knots: number[];
    degree: number;
    index: number;
    support: ParametricDomain;

    constructor(degree: number, knots: number[], index: number) {
        this.degree = degree;
        this.knots = knots;
        this.index = index;
        this.support = this.computeSupport();
    }

    // Standard evaluation (using blossom)
    evaluate(parameter: ParameterValue): number {
        // Evaluate blossom with repeated parameter
        return this.blossom(
            Array(this.degree + 1).fill(parameter)
        );
    }

    // Blossom evaluation
    blossom(parameters: ParameterValue[]): number {
        if (parameters.length !== this.degree + 1) {
            throw new Error("Invalid number of parameters for blossom");
        }

        return this.evaluateBlossom(
            parameters,
            this.index,
            this.knots
        );
    }

    private evaluateBlossom(
        parameters: ParameterValue[],
        index: number,
        knots: number[]
    ): number {
        if (parameters.length === 1) {
            // Base case: characteristic function
            return (
                knots[index] <= parameters[0][0] &&
                parameters[0][0] < knots[index + 1]
            ) ? 1 : 0;
        }

        // Recursive case using affine combination
        const alpha = this.computeBlossomWeight(
            parameters[0][0],
            knots[index],
            knots[index + parameters.length - 1]
        );

        return (
            alpha * this.evaluateBlossom(
                parameters.slice(1),
                index,
                knots
            ) +
            (1 - alpha) * this.evaluateBlossom(
                parameters.slice(1),
                index + 1,
                knots
            )
        );
    }
}

// Triangular B-spline basis with blossoming
class TriangularBasisWithBlossom implements BasisFunction {
    degree: number;
    support: ParametricDomain;

    blossom(parameters: ParameterValue[]): number {
        if (parameters.length !== this.degree + 1) {
            throw new Error("Invalid number of parameters for blossom");
        }

        return this.evaluateTriangularBlossom(parameters);
    }

    private evaluateTriangularBlossom(
        parameters: ParameterValue[]
    ): number {
        if (parameters.length === 1) {
            // Base case
            return this.evaluateCharacteristicFunction(parameters[0]);
        }

        // Recursive case using barycentric combination
        const weights = this.computeBarycentric(parameters[0]);
        let result = 0;

        for (let i = 0; i < 3; i++) {
            result += weights[i] * this.evaluateTriangularBlossom(
                parameters.slice(1)
            );
        }

        return result;
    }
}

// Parametric space with blossoming support
abstract class ParametricSpaceWithBlossom extends ParametricSpace {
    // Get polar form for specific index
    abstract getPolarForm(index: Index): PolarForm;

    // Evaluate blossom for specific index
    abstract evaluateBlossom(
        parameters: ParameterValue[],
        index: Index
    ): number;

    // Get derivatives using blossoming
    getDerivative(
        parameters: ParameterValue,
        index: Index,
        direction: number,
        order: number = 1
    ): number {
        // Use divided differences of blossom
        const h = 1e-6;
        const params = Array(this.degree + 1).fill(parameters);
        
        // Replace 'order' parameters with divided difference parameters
        for (let i = 0; i < order; i++) {
            params[i] = parameters.map((p, j) => 
                j === direction ? p + h : p
            );
        }

        return this.evaluateBlossom(params, index) / Math.pow(h, order);
    }
}

// Example usage with geometric operations
class GeometricOperations {
    // De Casteljau algorithm using blossoming
    static deCasteljau(
        controlPoints: Point[],
        basis: BasisFunction,
        parameter: ParameterValue
    ): Point {
        const n = controlPoints.length - 1;
        const params = Array(n + 1).fill(parameter);

        // Use blossoming for intermediate points
        const intermediate: Point[][] = Array(n + 1)
            .fill(null)
            .map(() => []);

        intermediate[0] = controlPoints;

        for (let r = 1; r <= n; r++) {
            for (let i = 0; i <= n - r; i++) {
                const weight = basis.blossom(params.slice(i, i + r + 1));
                
                intermediate[r][i] = interpolate(
                    intermediate[r-1][i],
                    intermediate[r-1][i+1],
                    weight
                );
            }
        }

        return intermediate[n][0];
    }

    // Degree elevation using blossoming
    static elevate(
        controlPoints: Point[],
        basis: BasisFunction
    ): Point[] {
        const elevated: Point[] = [];
        const n = controlPoints.length - 1;

        for (let i = 0; i <= n + 1; i++) {
            // Use blossoming to compute new control points
            const params = computeElevationParameters(i, n);
            const weight = basis.blossom(params);
            
            elevated[i] = computeElevatedPoint(
                controlPoints,
                weight,
                i
            );
        }

        return elevated;
    }
}

// Helper functions
function interpolate(p1: Point, p2: Point, t: number): Point {
    return {
        coordinates: p1.coordinates.map((c, i) => 
            c * (1 - t) + p2.coordinates[i] * t
        ),
        weight: p1.weight * (1 - t) + p2.weight * t
    };
}

// Blossoom and T-spline

interface TSplineBasis extends BasisFunction {
    // Local knot vectors at T-junctions
    localKnots: {
        u: number[],
        v: number[]
    };

    // Blossom using local knot structure
    blossom(parameters: ParameterValue[]): number {
        return this.evaluateLocalBlossom(
            parameters,
            this.localKnots
        );
    }
}

class TSplineSpace extends ParametricSpaceWithBlossom {
    // Extract local knot vectors at T-junction
    private getLocalKnots(index: Index): {u: number[], v: number[]} {
        // Navigate T-mesh to construct local knot vectors
        const tMeshNode = this.tMesh.getNode(index);
        return tMeshNode.extractLocalKnots();
    }

    // Evaluate basis using local structure
    evaluateBasis(parameters: ParameterValue, index: Index): number {
        const localKnots = this.getLocalKnots(index);
        return this.evaluateLocalBSplineBasis(parameters, localKnots);
    }
}

// 2. Anchors and T-junctions:

interface TMeshNode {
    index: Index;
    position: ParameterValue;
    // Local topology at T-junction
    neighbors: {
        north?: TMeshNode,
        south?: TMeshNode,
        east?: TMeshNode,
        west?: TMeshNode
    };
}

class TSplineBlossom {
    // Evaluate blossom at T-junction
    evaluateAtTJunction(
        parameters: ParameterValue[],
        node: TMeshNode
    ): number {
        // Get local knot structure
        const localKnots = this.extractLocalKnots(node);
        
        // Apply blossoming with local structure
        return this.evaluateLocalBlossom(
            parameters,
            localKnots,
            node.position
        );
    }

    // Extract local knot vectors considering T-junctions
    private extractLocalKnots(node: TMeshNode): LocalKnotStructure {
        const uKnots: number[] = [];
        const vKnots: number[] = [];

        // Navigate T-mesh to collect knots
        this.collectKnotsInDirection(node, 'u', uKnots);
        this.collectKnotsInDirection(node, 'v', vKnots);

        return { u: uKnots, v: vKnots };
    }
}

// 3. Analysis of T-junctions

class TSplineAnalysis {
    // Check linear independence using blossoming
    checkLinearIndependence(
        tMesh: TMesh,
        degree: number
    ): boolean {
        // Collect all basis functions
        const bases = this.getAllBases(tMesh);
        
        // Use blossoming to analyze linear independence
        return this.analyzeBlossomedBases(bases, degree);
    }

    // Analyze continuity at T-junction
    analyzeTJunctionContinuity(
        node: TMeshNode,
        degree: number
    ): number {
        // Use blossoming to analyze continuity
        const localBlossom = new TSplineBlossom(degree);
        return localBlossom.analyzeContinuity(node);
    }
}

// 4. Refinement Operations:

class TSplineRefinement {
    // Local refinement using blossoming
    refineAtTJunction(
        node: TMeshNode,
        degree: number
    ): TMeshNode[] {
        // Use blossoming to compute new control points
        const blossom = new TSplineBlossom(degree);
        
        // Compute new nodes
        return this.computeNewNodes(node, blossom);
    }

    private computeNewNodes(
        node: TMeshNode,
        blossom: TSplineBlossom
    ): TMeshNode[] {
        const newNodes: TMeshNode[] = [];
        
        // Use blossoming to ensure geometric continuity
        const parameters = this.computeRefinementParameters(node);
        
        for (const params of parameters) {
            const newPosition = blossom.evaluateAtTJunction(
                params,
                node
            );
            newNodes.push(this.createNewNode(newPosition));
        }
        
        return newNodes;
    }
}

// 5. Geometric Operations

class TSplineGeometry {
    // Evaluation using blossoming
    evaluate(
        parameters: ParameterValue,
        tMesh: TMesh
    ): Point {
        let result: Point = { coordinates: [0, 0, 0], weight: 0 };
        
        // Get active T-mesh nodes
        const activeNodes = tMesh.getActiveNodes(parameters);
        
        for (const node of activeNodes) {
            const blossom = new TSplineBlossom(this.degree);
            const value = blossom.evaluateAtTJunction(
                Array(this.degree + 1).fill(parameters),
                node
            );
            
            // Add contribution
            const controlPoint = node.controlPoint;
            result = this.addWeightedPoint(result, controlPoint, value);
        }
        
        return result;
    }
}

const tSpline = new TSplineSpace(3, [0, 0, 0], [1, 1, 1]);
export default tSpline;

