enum TopologyType {
    DISK,              // Simple planar domain
    SPHERE,            // No boundary
    TORUS,            // Genus 1
    MULTIPLE_TORI,    // Higher genus
    CUSTOM            // Arbitrary topology
}

interface RiemannSurface {
    genus: number;
    branchPoints: BranchPoint[];
    charts: Chart[];
    transitions: TransitionFunction[];
}

interface Chart {
    id: number;
    domain: ParametricDomain;
    localCoordinates: (u: number, v: number) => Complex;
}

interface TransitionFunction {
    sourceChart: Chart;
    targetChart: Chart;
    map: (z: Complex) => Complex;  // Holomorphic transition function
}

class TSplineRiemannApproximation {
    private tMeshes: TMesh[];
    private transitionMaps: Map<string, TransitionFunction>;

    constructor(
        private surface: RiemannSurface,
        private accuracy: number
    ) {
        this.tMeshes = [];
        this.transitionMaps = new Map();
    }

    approximate(): void {
        // 1. Decompose surface into charts
        const charts = this.decomposeIntoCharts();

        // 2. Create T-mesh for each chart
        for (const chart of charts) {
            const tMesh = this.createChartTMesh(chart);
            this.tMeshes.push(tMesh);
        }

        // 3. Handle transitions between charts
        this.setupTransitions();

        // 4. Refine mesh based on geometric features
        this.adaptiveMeshRefinement();
    }

    private decomposeIntoCharts(): Chart[] {
        const charts: Chart[] = [];

        switch (this.surface.genus) {
            case 0: // Sphere case
                return this.createSphereCharts();
            case 1: // Torus case
                return this.createTorusCharts();
            default: // Higher genus
                return this.createHigherGenusCharts();
        }
    }

    private createChartTMesh(chart: Chart): TMesh {
        // Create initial T-mesh structure
        const tMesh = new TMesh();

        // Handle different chart types
        switch (chart.domain.topology) {
            case TopologyType.DISK:
                this.setupDiskMesh(tMesh, chart);
                break;
            case TopologyType.SPHERE:
                this.setupSphereMesh(tMesh, chart);
                break;
            case TopologyType.TORUS:
                this.setupTorusMesh(tMesh, chart);
                break;
            // ... other cases
        }

        return tMesh;
    }

    private setupTransitions(): void {
        // Handle chart transitions
        for (let i = 0; i < this.tMeshes.length; i++) {
            for (let j = i + 1; j < this.tMeshes.length; j++) {
                if (this.chartsOverlap(i, j)) {
                    this.createTransitionMap(i, j);
                }
            }
        }
    }

    private adaptiveMeshRefinement(): void {
        let refined: boolean;
        do {
            refined = false;
            for (const tMesh of this.tMeshes) {
                // Refine based on:
                // 1. Curvature
                refined ||= this.refineByCurvature(tMesh);
                
                // 2. Branch points
                refined ||= this.refineNearBranchPoints(tMesh);
                
                // 3. Transition regions
                refined ||= this.refineTransitionRegions(tMesh);
            }
        } while (refined && !this.meetAccuracy());
    }

    private refineByCurvature(tMesh: TMesh): boolean {
        const refinements: TSplineRefinement[] = [];

        // Calculate Gaussian curvature at sample points
        for (const face of tMesh.faces) {
            const K = this.calculateGaussianCurvature(face);
            if (Math.abs(K) > this.accuracy) {
                refinements.push(this.createCurvatureBasedRefinement(face));
            }
        }

        return this.applyRefinements(tMesh, refinements);
    }

    private refineNearBranchPoints(tMesh: TMesh): boolean {
        const refinements: TSplineRefinement[] = [];

        for (const point of this.surface.branchPoints) {
            if (this.pointInChartDomain(point, tMesh)) {
                refinements.push(
                    this.createBranchPointRefinement(point, tMesh)
                );
            }
        }

        return this.applyRefinements(tMesh, refinements);
    }

    class ChartTransition {
        constructor(
            private source: TMesh,
            private target: TMesh,
            private transitionFunction: TransitionFunction
        ) {}

        // Ensure smooth transition between charts
        enforceCompatibility(): void {
            // 1. Match control points at boundaries
            this.matchBoundaryControlPoints();

            // 2. Ensure G1 continuity
            this.enforceG1Continuity();
        }

        private matchBoundaryControlPoints(): void {
            const boundaryVertices = this.source.getBoundaryVertices();
            for (const vertex of boundaryVertices) {
                const mappedPoint = this.transitionFunction.map(
                    vertex.position
                );
                // Create or adjust corresponding vertex in target mesh
                this.ensureMatchingVertex(vertex, mappedPoint);
            }
        }

        private enforceG1Continuity(): void {
            // Implement G1 continuity constraints at chart boundaries
            // This ensures smooth transitions between charts
        }
    }
}

class SphericalTMeshProjection {
    constructor(
        private tMesh: TMesh,
        private sphereRadius: number = 1.0
    ) {}

    project(): void {
        // 1. First move mesh towards convex configuration
        this.makeConvex();
        
        // 2. Project vertices onto sphere
        this.projectOntoSphere();
        
        // 3. Optimize vertex distribution
        this.optimizeSphericalDistribution();
    }

    private makeConvex(): void {
        // Use volumetric distance function gradient
        const distanceField = this.computeVolumetricDistance();
        let converged = false;
        
        while (!converged) {
            converged = true;
            
            for (const vertex of this.tMesh.vertices) {
                const gradient = distanceField.getGradient(vertex.position);
                const movement = this.computeConvexifyingMovement(vertex, gradient);
                
                if (movement.length() > this.tolerance) {
                    vertex.position = vertex.position.add(movement);
                    converged = false;
                }
            }
        }
    }

    private projectOntoSphere(): void {
        // Project each vertex onto sphere surface
        for (const vertex of this.tMesh.vertices) {
            const direction = vertex.position.normalize();
            vertex.position = direction.scale(this.sphereRadius);
        }
    }

    private optimizeSphericalDistribution(): void {
        // Minimize distortion while maintaining spherical shape
        const optimizer = new SphericalDistributionOptimizer(
            this.tMesh,
            this.sphereRadius
        );
        optimizer.optimize();
    }
}

class SphericalDistributionOptimizer {
    private energy: number = Infinity;
    
    constructor(
        private tMesh: TMesh,
        private radius: number,
        private maxIterations: number = 100,
        private tolerance: number = 1e-6
    ) {}

    optimize(): void {
        let iteration = 0;
        let prevEnergy = Infinity;
        
        while (iteration++ < this.maxIterations) {
            // 1. Compute current energy
            const currentEnergy = this.computeDistortionEnergy();
            
            // 2. Check convergence
            if (Math.abs(currentEnergy - prevEnergy) < this.tolerance) {
                break;
            }
            
            // 3. Optimize vertex positions
            this.optimizeVertexPositions();
            
            // 4. Project back to sphere
            this.projectToSphere();
            
            prevEnergy = currentEnergy;
        }
    }

    private computeDistortionEnergy(): number {
        let energy = 0;
        
        // Sum up various energy terms
        energy += this.computeAreaDistortion();
        energy += this.computeAngleDistortion();
        energy += this.computeEdgeLengthDistortion();
        
        return energy;
    }

    private optimizeVertexPositions(): void {
        for (const vertex of this.tMesh.vertices) {
            // Compute optimal position based on neighbors
            const newPosition = this.computeOptimalPosition(vertex);
            vertex.position = newPosition;
        }
    }

    private computeOptimalPosition(vertex: Vertex): Vector3 {
        // Get one-ring neighbors
        const neighbors = vertex.getNeighbors();
        
        // Compute Laplacian coordinates
        const laplacian = this.computeLaplacianCoordinates(
            vertex,
            neighbors
        );
        
        // Update position while considering spherical constraint
        return this.solveConstrainedPosition(vertex, laplacian);
    }

    private projectToSphere(): void {
        for (const vertex of this.tMesh.vertices) {
            const direction = vertex.position.normalize();
            vertex.position = direction.scale(this.radius);
        }
    }
}

class SphericalParameterization {
    static createForGenus0(mesh: TMesh): TMesh {
        // 1. Verify genus 0 topology
        if (!this.isGenus0(mesh)) {
            throw new Error("Mesh must be genus 0 for spherical parameterization");
        }

        // 2. Create initial projection
        const projection = new SphericalTMeshProjection(mesh);
        projection.project();

        // 3. Create chart structure
        return this.createSphericalCharts(mesh);
    }

    private static createSphericalCharts(mesh: TMesh): TMesh {
        const chartMesh = new TMesh();
        
        // Create two hemispheric charts
        const northChart = this.createHemisphereChart(mesh, "north");
        const southChart = this.createHemisphereChart(mesh, "south");
        
        // Handle transition region (equatorial band)
        this.setupChartTransition(northChart, southChart);
        
        return chartMesh;
    }

    private static createHemisphereChart(
        mesh: TMesh,
        hemisphere: "north" | "south"
    ): Chart {
        const chart = new Chart();
        
        // Select vertices for this hemisphere
        const vertices = this.selectHemisphereVertices(mesh, hemisphere);
        
        // Create chart mesh structure
        this.buildChartMesh(chart, vertices);
        
        // Setup parametric mapping
        this.setupHemisphericMapping(chart, hemisphere);
        
        return chart;
    }

    private static setupChartTransition(
        chart1: Chart,
        chart2: Chart
    ): void {
        // Create transition functions between charts
        const transition = new ChartTransition(chart1, chart2);
        
        // Ensure smooth transition near equator
        transition.enforceCompatibility();
    }
}

// Example usage
const genus0Surface = new RiemannSurface(/* ... */);
const tMesh = new TMesh(/* ... */);

// Project to sphere
const sphericalParam = SphericalParameterization.createForGenus0(tMesh);

// Verify quality
const qualityMetrics = new MeshQualityAnalyzer(sphericalParam);
console.log("Parameterization quality:", qualityMetrics.analyze());


interface DiscreteMobiusGroup {
    generators: MobiusTransformation[];
    fundamentalDomain: Domain;
}

class AutomorphicFunctionApproximation {
    constructor(
        private group: DiscreteMobiusGroup,
        private accuracy: number
    ) {}

    // Main approximation method
    approximate(): TSplineAutomorphicFunction {
        // 1. Set up fundamental domain
        const fundamentalMesh = this.createFundamentalDomainMesh();
        
        // 2. Enforce automorphic conditions
        this.enforceGroupInvariance(fundamentalMesh);
        
        // 3. Handle singularities (if any)
        this.handleSingularities(fundamentalMesh);
        
        return new TSplineAutomorphicFunction(fundamentalMesh, this.group);
    }

    private createFundamentalDomainMesh(): TMesh {
        // Create T-mesh for fundamental domain
        const mesh = new TMesh();
        
        // Set up initial mesh based on domain geometry
        if (this.group.fundamentalDomain.type === "HyperbolicPolygon") {
            this.setupHyperbolicMesh(mesh);
        } else {
            this.setupEuclideanMesh(mesh);
        }
        
        return mesh;
    }

    private enforceGroupInvariance(mesh: TMesh): void {
        // Ensure T-spline is invariant under group action
        for (const generator of this.group.generators) {
            this.enforceGeneratorInvariance(mesh, generator);
        }
    }

    private enforceGeneratorInvariance(
        mesh: TMesh,
        generator: MobiusTransformation
    ): void {
        // Set up constraints for control points
        const constraints = new GroupInvarianceConstraints(generator);
        
        // Apply constraints to mesh
        for (const vertex of mesh.vertices) {
            constraints.applyToVertex(vertex);
        }
    }
}

class GroupInvarianceConstraints {
    constructor(private transformation: MobiusTransformation) {}

    applyToVertex(vertex: TVertex): void {
        // Compute transformed position
        const transformedPos = this.transformation.apply(vertex.position);
        
        // Set up constraint equations
        this.addPositionConstraint(vertex, transformedPos);
        
        // Handle derivatives for smoothness
        this.addDerivativeConstraints(vertex, transformedPos);
    }

    private addPositionConstraint(vertex: TVertex, target: Complex): void {
        // Ensure vertex value matches transformed value
        vertex.addConstraint(new PositionConstraint(target));
    }

    private addDerivativeConstraints(vertex: TVertex, target: Complex): void {
        // Ensure derivatives transform correctly
        const jacobian = this.transformation.getJacobian(vertex.position);
        vertex.addConstraint(new DerivativeConstraint(jacobian));
    }
}

class TSplineAutomorphicFunction {
    constructor(
        private mesh: TMesh,
        private group: DiscreteMobiusGroup
    ) {}

    evaluate(z: Complex): Complex {
        // 1. Map z to fundamental domain if needed
        const w = this.mapToFundamentalDomain(z);
        
        // 2. Evaluate T-spline at mapped point
        return this.evaluateSpline(w);
    }

    private mapToFundamentalDomain(z: Complex): Complex {
        // Use group operations to map point to fundamental domain
        let w = z;
        let moved = true;
        
        while (moved && !this.group.fundamentalDomain.contains(w)) {
            moved = false;
            for (const generator of this.group.generators) {
                const newW = generator.apply(w);
                if (this.isCloserToFundamentalDomain(newW, w)) {
                    w = newW;
                    moved = true;
                    break;
                }
            }
        }
        
        return w;
    }

    private evaluateSpline(w: Complex): Complex {
        // Standard T-spline evaluation at point w
        return this.mesh.evaluate(w);
    }
}

// Example for modular functions
class ModularFunctionApproximation extends AutomorphicFunctionApproximation {
    constructor(accuracy: number) {
        // Set up modular group PSL(2,Z)
        const generators = [
            new MobiusTransformation([1, 1, 0, 1]),  // T: z → z + 1
            new MobiusTransformation([0, -1, 1, 0])  // S: z → -1/z
        ];
        
        const fundamentalDomain = new HyperbolicPolygon([
            // Standard fundamental domain for modular group
            new Complex(0.5, Math.sqrt(3)/2),
            new Complex(-0.5, Math.sqrt(3)/2),
            new Complex(Infinity)
        ]);
        
        super(new DiscreteMobiusGroup(generators, fundamentalDomain), accuracy);
    }

    protected handleSpecialPoints(): void {
        // Handle cusps and elliptic points
        this.handleCusps();
        this.handleEllipticPoints();
    }

    private handleCusps(): void {
        // Refine mesh near cusps (i∞, 0, 1)
        for (const cusp of this.getCusps()) {
            this.refineMeshNearCusp(cusp);
        }
    }

    private handleEllipticPoints(): void {
        // Handle special behavior near elliptic points
        // (i and ρ = e^{2πi/3})
        for (const point of this.getEllipticPoints()) {
            this.refineMeshNearEllipticPoint(point);
        }
    }
}

// Example usage
const modularFunction = new ModularFunctionApproximation(1e-6);
const jFunction = modularFunction.approximate();

// Evaluate at some points
console.log("j(i):", jFunction.evaluate(new Complex(0, 1)));
console.log("j(ρ):", jFunction.evaluate(
    new Complex(0.5, Math.sqrt(3)/2)
));


private setupFundamentalDomain(): void {
    // Create mesh covering fundamental domain
    // Handle boundary identifications
    // Deal with cusps and special points
}

private handleSpecialPoints(): void {
    // Refine mesh near:
    // - Cusps
    // - Elliptic points
    // - Other singularities
}

enum HyperbolicModel {
    UPPER_HALF_PLANE,  // H = {z ∈ C : Im(z) > 0}
    POINCARE_DISK      // D = {z ∈ C : |z| < 1}
}

class HyperbolicTMesh {
    constructor(
        private model: HyperbolicModel,
        private fundamentalDomain: HyperbolicPolygon
    ) {}

    // Hyperbolic metric tensor
    getMetric(z: Complex): Matrix2x2 {
        if (this.model === HyperbolicModel.UPPER_HALF_PLANE) {
            const y = z.imaginary;
            return new Matrix2x2(
                1/(y*y), 0,
                0, 1/(y*y)
            );
        } else {  // Poincaré disk
            const r2 = z.normSquared();
            return new Matrix2x2(
                4/((1-r2)*(1-r2)), 0,
                0, 4/((1-r2)*(1-r2))
            );
        }
    }

    // Hyperbolic distance
    distance(z1: Complex, z2: Complex): number {
        if (this.model === HyperbolicModel.UPPER_HALF_PLANE) {
            return this.upperHalfPlaneDistance(z1, z2);
        } else {
            return this.poincareDistance(z1, z2);
        }
    }

    private upperHalfPlaneDistance(z1: Complex, z2: Complex): number {
        // d(z1,z2) = arcosh(1 + |z1-z2|²/(2Im(z1)Im(z2)))
        const diff = z1.subtract(z2);
        return Math.acosh(1 + 
            diff.normSquared()/(2*z1.imaginary*z2.imaginary)
        );
    }

    private poincareDistance(z1: Complex, z2: Complex): number {
        // d(z1,z2) = 2arctanh(|z1-z2|/|1-z1z2*|)
        const num = z1.subtract(z2).abs();
        const den = Complex.ONE.subtract(
            z1.multiply(z2.conjugate())
        ).abs();
        return 2 * Math.atanh(num/den);
    }
}

class HyperbolicTMeshGenerator {
    constructor(
        private group: DiscreteMobiusGroup,
        private model: HyperbolicModel
    ) {}

    generateMesh(): TMesh {
        // 1. Create initial mesh in fundamental domain
        const mesh = this.createInitialMesh();

        // 2. Refine based on hyperbolic geometry
        this.refineHyperbolicMesh(mesh);

        // 3. Ensure group invariance
        this.enforceGroupInvariance(mesh);

        return mesh;
    }

    private createInitialMesh(): TMesh {
        const mesh = new TMesh();

        // Create vertices following hyperbolic geometry
        if (this.model === HyperbolicModel.UPPER_HALF_PLANE) {
            this.createUpperHalfPlaneMesh(mesh);
        } else {
            this.createPoincareDiskMesh(mesh);
        }

        return mesh;
    }

    private refineHyperbolicMesh(mesh: TMesh): void {
        const refinement = new HyperbolicMeshRefinement(mesh, this.model);
        
        // Refine based on:
        // 1. Hyperbolic curvature
        refinement.refineByCurvature();
        
        // 2. Distance to special points
        refinement.refineNearSpecialPoints();
        
        // 3. Group action boundaries
        refinement.refineNearGroupBoundaries();
    }
}

class HyperbolicMeshRefinement {
    constructor(
        private mesh: TMesh,
        private model: HyperbolicModel
    ) {}

    refineByCurvature(): void {
        for (const face of this.mesh.faces) {
            const K = this.hyperbolicCurvature(face);
            if (Math.abs(K) > this.threshold) {
                this.refineFace(face);
            }
        }
    }

    refineNearSpecialPoints(): void {
        // Handle different types of special points
        this.refineNearCusps();
        this.refineNearEllipticPoints();
        this.refineNearBranchPoints();
    }

    private refineNearCusps(): void {
        for (const cusp of this.cusps) {
            // Use horocyclic neighborhoods
            const neighborhood = this.getHorocyclicNeighborhood(cusp);
            this.refineInNeighborhood(neighborhood);
        }
    }

    private getHorocyclicNeighborhood(cusp: Complex): HyperbolicDomain {
        if (this.model === HyperbolicModel.UPPER_HALF_PLANE) {
            // For cusp at infinity, use horizontal strip
            return new HorizontalStrip(this.height);
        } else {
            // For finite cusps, use horocycle
            return new Horocycle(cusp, this.radius);
        }
    }
}

class FundamentalDomainMesh {
    constructor(
        private domain: HyperbolicPolygon,
        private model: HyperbolicModel
    ) {}

    createMesh(): TMesh {
        const mesh = new TMesh();

        // 1. Create vertices at corners
        this.addCornerVertices(mesh);

        // 2. Add boundary vertices
        this.addBoundaryVertices(mesh);

        // 3. Fill interior
        this.fillInterior(mesh);

        // 4. Create edges and faces
        this.createConnectivity(mesh);

        return mesh;
    }

    private addCornerVertices(mesh: TMesh): void {
        for (const corner of this.domain.vertices) {
            if (this.isFinitePoint(corner)) {
                mesh.addVertex(new TVertex(corner));
            } else {
                // Handle points at infinity
                this.handleInfiniteVertex(mesh, corner);
            }
        }
    }

    private fillInterior(mesh: TMesh): void {
        // Use hyperbolic-aware triangulation
        const triangulation = new HyperbolicDelaunay(
            this.domain,
            this.model
        );
        
        // Convert to T-mesh structure
        this.convertTriangulationToTMesh(triangulation, mesh);
    }
}

// Example usage for modular group
class ModularTMeshGenerator {
    static createModularMesh(): TMesh {
        // Standard fundamental domain for PSL(2,Z)
        const fundamentalDomain = new HyperbolicPolygon([
            new Complex(-0.5, Math.sqrt(3)/2),
            new Complex(0.5, Math.sqrt(3)/2),
            new Complex(Infinity)
        ]);

        const generator = new HyperbolicTMeshGenerator(
            ModularGroup.PSL2Z,
            HyperbolicModel.UPPER_HALF_PLANE
        );

        return generator.generateMesh();
    }
}

// Example usage
const modularMesh = ModularTMeshGenerator.createModularMesh();
const jInvariant = new ModularFunctionApproximation(modularMesh);

// Test evaluation at special points
console.log("j(i):", jInvariant.evaluate(new Complex(0, 1)));
console.log("j(ρ):", jInvariant.evaluate(
    new Complex(0.5, Math.sqrt(3)/2)
));

class HyperbolicEdge {
    constructor(
        private start: Complex,
        private end: Complex,
        private model: HyperbolicModel
    ) {}

    // Get geodesic parameters
    getGeodesicParameters(): GeodesicParameters {
        if (this.model === HyperbolicModel.POINCARE_DISK) {
            return this.getPoincareDiskGeodesic();
        } else {
            return this.getUpperHalfPlaneGeodesic();
        }
    }

    // Get point along geodesic at parameter t ∈ [0,1]
    getPointAtParameter(t: number): Complex {
        const params = this.getGeodesicParameters();
        
        if (this.model === HyperbolicModel.POINCARE_DISK) {
            return this.getPointOnCircularArc(t, params);
        } else {
            return this.getPointOnUpperHalfGeodesic(t, params);
        }
    }

    private getPoincareDiskGeodesic(): GeodesicParameters {
        // For Poincaré disk, compute center and radius of circular arc
        const z1 = this.start;
        const z2 = this.end;

        // If points are on a diameter, return straight line parameters
        if (this.isOnDiameter(z1, z2)) {
            return {
                type: 'LINE',
                direction: z2.subtract(z1)
            };
        }

        // Otherwise, compute circle parameters
        return this.computePerpendicularCircle(z1, z2);
    }

    private getUpperHalfPlaneGeodesic(): GeodesicParameters {
        const z1 = this.start;
        const z2 = this.end;

        // Check if points are on same vertical line
        if (Math.abs(z1.real - z2.real) < 1e-10) {
            return {
                type: 'VERTICAL',
                x: z1.real
            };
        }

        // Otherwise, compute semicircle parameters
        return this.computeSemicircle(z1, z2);
    }

    private computePerpendicularCircle(z1: Complex, z2: Complex): GeodesicParameters {
        // Compute center and radius of circle perpendicular to unit circle
        // passing through z1 and z2
        const a = z1.conjugate().multiply(z2).subtract(Complex.ONE);
        const b = z2.subtract(z1);
        const c = z1.conjugate().multiply(z2.multiply(z1)).subtract(z1);

        const center = c.divide(a);
        const radius = b.abs() / (2 * a.abs());

        return {
            type: 'CIRCLE',
            center,
            radius,
            startAngle: Math.atan2(z1.subtract(center).imaginary, 
                                 z1.subtract(center).real),
            endAngle: Math.atan2(z2.subtract(center).imaginary, 
                                z2.subtract(center).real)
        };
    }

    private computeSemicircle(z1: Complex, z2: Complex): GeodesicParameters {
        // Compute center and radius of semicircle in upper half-plane
        const x1 = z1.real;
        const x2 = z2.real;
        const y1 = z1.imaginary;
        const y2 = z2.imaginary;

        const center = new Complex(
            (x1 + x2)/2 + (y2 - y1)*(y1 + y2)/(2*(x2 - x1)),
            0
        );
        const radius = center.subtract(z1).abs();

        return {
            type: 'SEMICIRCLE',
            center,
            radius
        };
    }
}

class HyperbolicTMesh {
    private edges: HyperbolicEdge[] = [];
    private vertices: Complex[] = [];

    constructor(
        private model: HyperbolicModel
    ) {}

    addEdge(start: Complex, end: Complex): void {
        this.edges.push(new HyperbolicEdge(start, end, this.model));
    }

    // Subdivide edge using hyperbolic geodesic
    subdivideEdge(edge: HyperbolicEdge, t: number): Complex {
        return edge.getPointAtParameter(t);
    }

    // Get control points along geodesic
    getControlPoints(edge: HyperbolicEdge, numPoints: number): Complex[] {
        const points: Complex[] = [];
        for (let i = 0; i <= numPoints; i++) {
            const t = i / numPoints;
            points.push(edge.getPointAtParameter(t));
        }
        return points;
    }
}

class HyperbolicTSpline {
    constructor(
        private mesh: HyperbolicTMesh,
        private degree: number
    ) {}

    // Evaluate T-spline using hyperbolic geodesics
    evaluate(u: number, v: number): Complex {
        const controlPoints = this.getRelevantControlPoints(u, v);
        const weights = this.computeWeights(u, v);
        
        // Compute weighted sum using hyperbolic geometry
        return this.hyperbolicWeightedSum(controlPoints, weights);
    }

    private hyperbolicWeightedSum(
        points: Complex[],
        weights: number[]
    ): Complex {
        // Use hyperbolic barycenter computation
        let result = Complex.ZERO;
        let weightSum = 0;

        for (let i = 0; i < points.length; i++) {
            const geodesic = new HyperbolicEdge(
                result,
                points[i],
                this.mesh.model
            );
            result = geodesic.getPointAtParameter(weights[i]);
            weightSum += weights[i];
        }

        return result;
    }
}

// Example usage for modular forms
class ModularTSplineApproximation {
    constructor(
        private fundamentalDomain: HyperbolicPolygon,
        private model: HyperbolicModel = HyperbolicModel.UPPER_HALF_PLANE
    ) {}

    createMesh(): HyperbolicTMesh {
        const mesh = new HyperbolicTMesh(this.model);

        // Add edges following hyperbolic geodesics
        for (let i = 0; i < this.fundamentalDomain.vertices.length; i++) {
            const start = this.fundamentalDomain.vertices[i];
            const end = this.fundamentalDomain.vertices[
                (i + 1) % this.fundamentalDomain.vertices.length
            ];
            
            mesh.addEdge(start, end);
        }

        // Add interior edges and refinement
        this.addInteriorStructure(mesh);

        return mesh;
    }

    private addInteriorStructure(mesh: HyperbolicTMesh): void {
        // Add interior geodesic edges for T-mesh structure
        // Handle special points (cusps, elliptic points)
        // Ensure proper refinement near these points
    }
}

// Example visualization
class HyperbolicTMeshVisualizer {
    static drawMesh(
        mesh: HyperbolicTMesh,
        context: CanvasRenderingContext2D
    ): void {
        for (const edge of mesh.edges) {
            const params = edge.getGeodesicParameters();
            
            if (params.type === 'CIRCLE') {
                this.drawCircularArc(context, params);
            } else if (params.type === 'LINE') {
                this.drawLine(context, params);
            } else if (params.type === 'SEMICIRCLE') {
                this.drawSemicircle(context, params);
            }
        }
    }
}

class StarPointAnalysis {
    constructor(
        private centralVertex: TVertex,
        private valence: number,    // Number of edges meeting at star point
        private degree: number = 3  // Spline degree
    ) {}

    analyzeStarPoint(): StarPointProperties {
        // 1. Set up local parameterization
        const localParam = this.createLocalParameterization();
        
        // 2. Analyze geometric properties
        const properties = this.analyzeGeometry(localParam);
        
        // 3. Check continuity conditions
        const continuity = this.analyzeContinuity();

        return { localParam, properties, continuity };
    }

    private createLocalParameterization(): LocalParameterization {
        // Create polar-like structure around star point
        return new StarPointParameterization(
            this.centralVertex,
            this.valence
        );
    }
}

class StarPointParameterization {
    constructor(
        private center: TVertex,
        private valence: number
    ) {}

    // Map from polar-like coordinates to surface
    mapToSurface(r: number, θ: number): Point3D {
        // r: radial distance from star point
        // θ: angle parameter
        
        // Normalize angle to sector
        const sector = Math.floor(θ * this.valence / (2 * Math.PI));
        const localθ = θ - (2 * Math.PI * sector) / this.valence;
        
        return this.evaluateSector(r, localθ, sector);
    }

    // Create geodesic edges around star point
    createStarEdges(): HyperbolicEdge[] {
        const edges: HyperbolicEdge[] = [];
        
        for (let i = 0; i < this.valence; i++) {
            const angle = (2 * Math.PI * i) / this.valence;
            edges.push(this.createRadialEdge(angle));
        }
        
        return edges;
    }
}

class ExtendedTSplineSurface {
    constructor(
        private controlPoints: ControlPoint[],
        private knots: KnotStructure,
        private starPoints: StarPoint[]
    ) {}

    evaluate(u: number, v: number): Point3D {
        // Check if we're near a star point
        const nearestStar = this.findNearestStarPoint(u, v);
        
        if (nearestStar && this.isInStarRegion(u, v, nearestStar)) {
            return this.evaluateStarRegion(u, v, nearestStar);
        } else {
            return this.evaluateRegularRegion(u, v);
        }
    }

    private evaluateStarRegion(
        u: number, 
        v: number, 
        starPoint: StarPoint
    ): Point3D {
        // Convert to polar-like coordinates
        const { r, θ } = this.convertToStarCoordinates(u, v, starPoint);
        
        // Use special star point parameterization
        return starPoint.parameterization.mapToSurface(r, θ);
    }

    private isInStarRegion(
        u: number, 
        v: number, 
        starPoint: StarPoint
    ): boolean {
        // Check if point is within star point's influence region
        const dist = this.getParametricDistance(u, v, starPoint.center);
        return dist < starPoint.influenceRadius;
    }
}

class StarPoint {
    constructor(
        public center: TVertex,
        public valence: number,
        public parameterization: StarPointParameterization
    ) {
        this.influenceRadius = this.computeInfluenceRadius();
        this.transitionFunction = this.createTransitionFunction();
    }

    // Compute geometric properties
    computeProperties(): StarPointGeometry {
        return {
            tangentCone: this.computeTangentCone(),
            gaussianCurvature: this.computeGaussianCurvature(),
            continuity: this.analyzeLocalContinuity()
        };
    }

    // Create transition between star region and regular region
    private createTransitionFunction(): (r: number) => number {
        return (r: number) => {
            // Smooth transition function
            // 1 at center, 0 outside influence radius
            const t = r / this.influenceRadius;
            return t < 1 ? Math.pow(1 - t*t, 2) : 0;
        };
    }

    // Analyze local differential geometry
    analyzeDifferentialGeometry(): DifferentialProperties {
        const sectors = this.createSectorDecomposition();
        const properties: DifferentialProperties = {
            principalCurvatures: [],
            gaussianCurvature: 0,
            meanCurvature: 0
        };

        // Analyze each sector
        sectors.forEach(sector => {
            const sectorProps = this.analyzeSector(sector);
            this.updateGlobalProperties(properties, sectorProps);
        });

        return properties;
    }
}

class StarPointSector {
    constructor(
        private startAngle: number,
        private endAngle: number,
        private controlPoints: ControlPoint[]
    ) {}

    // Create geodesic boundaries
    createBoundaries(): HyperbolicEdge[] {
        return [
            this.createRadialEdge(),
            this.createCircularArcEdge()
        ];
    }

    // Evaluate points within sector
    evaluate(r: number, θ: number): Point3D {
        // Convert to local coordinates
        const { u, v } = this.sectorToLocal(r, θ);
        
        // Use modified basis functions near star point
        return this.evaluateWithModifiedBasis(u, v);
    }

    private evaluateWithModifiedBasis(u: number, v: number): Point3D {
        const basis = this.createModifiedBasisFunctions();
        return basis.evaluate(u, v, this.controlPoints);
    }
}

class ModifiedBasisFunctions {
    constructor(
        private valence: number,
        private degree: number
    ) {}

    evaluate(u: number, v: number, controlPoints: ControlPoint[]): Point3D {
        // Create special basis functions for star point
        const weights = this.computeModifiedWeights(u, v);
        
        // Blend control points using modified weights
        return this.blendWithWeights(controlPoints, weights);
    }

    private computeModifiedWeights(u: number, v: number): number[] {
        // Modify standard B-spline weights to handle star point
        const standardWeights = this.computeStandardWeights(u, v);
        return this.adjustWeightsForStarPoint(standardWeights);
    }
}

// Example usage
const starPoint = new StarPoint(
    centralVertex,
    5,  // valence
    new StarPointParameterization(centralVertex, 5)
);

const analysis = new StarPointAnalysis(
    centralVertex,
    5,  // valence
    3   // degree
);

const properties = analysis.analyzeStarPoint();
console.log("Star Point Analysis:", properties);

// Visualize star point region
const visualizer = new StarPointVisualizer(starPoint);
visualizer.drawParametricLines();
visualizer.drawGeometricProperties();


class ExtraordinaryPointAnalysis {
    constructor(
        private vertex: TVertex,
        private valence: number
    ) {}

    createHyperbolicVisualization(): HyperbolicTMesh {
        // Map the neighborhood of extraordinary point to hyperbolic plane
        const hyperbolicMesh = new HyperbolicTMesh(
            HyperbolicModel.POINCARE_DISK  // or UPPER_HALF_PLANE
        );

        // Create radial structure in hyperbolic plane
        this.createRadialStructure(hyperbolicMesh);

        return hyperbolicMesh;
    }

    private createRadialStructure(mesh: HyperbolicTMesh): void {
        // Create central point
        const center = new Complex(0, 0);  // Center of Poincaré disk

        // Create radial edges with equal hyperbolic angles
        for (let i = 0; i < this.valence; i++) {
            const angle = (2 * Math.PI * i) / this.valence;
            
            // Create vertices along geodesic rays
            const ray = this.createGeodesicRay(center, angle);
            this.subdivideRay(mesh, ray);
        }

        // Add connecting edges between rays
        this.createConnectingEdges(mesh);
    }

    private createGeodesicRay(
        start: Complex,
        angle: number
    ): HyperbolicEdge {
        // Create geodesic ray in chosen direction
        const endPoint = this.getIdealPoint(angle);
        return new HyperbolicEdge(start, endPoint, this.mesh.model);
    }

    private subdivideRay(
        mesh: HyperbolicTMesh,
        ray: HyperbolicEdge
    ): void {
        // Create vertices along ray at appropriate hyperbolic distances
        const points: Complex[] = [];
        
        // Use geometric progression for spacing in hyperbolic metric
        for (let i = 0; i < this.subdivisions; i++) {
            const t = this.getHyperbolicParameter(i);
            points.push(ray.getPointAtParameter(t));
        }

        // Add vertices and edges to mesh
        for (let i = 0; i < points.length; i++) {
            mesh.addVertex(points[i]);
            if (i > 0) {
                mesh.addEdge(points[i-1], points[i]);
            }
        }
    }
}

class HyperbolicMeshVisualizer {
    constructor(
        private mesh: HyperbolicTMesh,
        private canvas: HTMLCanvasElement
    ) {}

    draw(): void {
        const ctx = this.canvas.getContext('2d');
        if (!ctx) return;

        // Draw unit circle for Poincaré disk
        if (this.mesh.model === HyperbolicModel.POINCARE_DISK) {
            this.drawUnitCircle(ctx);
        }

        // Draw mesh edges as geodesics
        for (const edge of this.mesh.edges) {
            this.drawGeodesicEdge(ctx, edge);
        }

        // Draw vertices
        for (const vertex of this.mesh.vertices) {
            this.drawVertex(ctx, vertex);
        }

        // Optionally draw coordinate grid
        this.drawHyperbolicGrid(ctx);
    }

    private drawGeodesicEdge(
        ctx: CanvasRenderingContext2D,
        edge: HyperbolicEdge
    ): void {
        const params = edge.getGeodesicParameters();

        ctx.beginPath();
        if (params.type === 'CIRCLE') {
            // Draw circular arc
            const center = this.transformPoint(params.center!);
            const radius = params.radius! * this.scale;
            ctx.arc(
                center.x,
                center.y,
                radius,
                params.startAngle!,
                params.endAngle!
            );
        } else if (params.type === 'LINE') {
            // Draw straight line
            const start = this.transformPoint(edge.start);
            const end = this.transformPoint(edge.end);
            ctx.moveTo(start.x, start.y);
            ctx.lineTo(end.x, end.y);
        }
        ctx.stroke();
    }

    private drawHyperbolicGrid(ctx: CanvasRenderingContext2D): void {
        // Draw radial geodesics
        for (let angle = 0; angle < 2 * Math.PI; angle += Math.PI / 12) {
            const ray = new HyperbolicEdge(
                new Complex(0, 0),
                this.getIdealPoint(angle),
                this.mesh.model
            );
            this.drawGeodesicEdge(ctx, ray);
        }

        // Draw concentric circles
        for (let r = 0.2; r < 1; r += 0.2) {
            this.drawHyperbolicCircle(ctx, r);
        }
    }
}

class ValenceAnalyzer {
    constructor(private surface: TSplineSurface) {}

    analyzeExtraordinaryPoint(vertex: TVertex): ExtraordinaryPointAnalysis {
        const valence = this.getVertexValence(vertex);
        const analysis = new ExtraordinaryPointAnalysis(vertex, valence);

        // Create hyperbolic visualization
        const hyperbolicMesh = analysis.createHyperbolicVisualization();

        // Analyze local properties
        const properties = {
            geometricProperties: this.analyzeGeometry(vertex),
            continuityProperties: this.analyzeContinuity(vertex),
            parameterization: this.analyzeParameterization(vertex)
        };

        return {
            hyperbolicMesh,
            properties,
            vertex,
            valence
        };
    }

    private analyzeGeometry(vertex: TVertex): GeometricProperties {
        return {
            gaussianCurvature: this.computeGaussianCurvature(vertex),
            meanCurvature: this.computeMeanCurvature(vertex),
            principalDirections: this.computePrincipalDirections(vertex)
        };
    }

    private analyzeContinuity(vertex: TVertex): ContinuityProperties {
        return {
            tangentContinuity: this.checkTangentContinuity(vertex),
            curvatureContinuity: this.checkCurvatureContinuity(vertex)
        };
    }
}

// Example usage
const surface = new TSplineSurface(/* ... */);
const analyzer = new ValenceAnalyzer(surface);

// Analyze extraordinary point
const vertex = surface.findExtraordinaryPoint();
const analysis = analyzer.analyzeExtraordinaryPoint(vertex);

// Visualize in hyperbolic plane
const canvas = document.getElementById('canvas') as HTMLCanvasElement;
const visualizer = new HyperbolicMeshVisualizer(
    analysis.hyperbolicMesh,
    canvas
);
visualizer.draw();

console.log("Analysis results:", analysis.properties);


class HigherGenusAnalysis {
    constructor(
        private surface: TSplineSurface,
        private genus: number
    ) {
        if (genus <= 1) {
            throw new Error("This analysis is for genus > 1 surfaces");
        }
    }

    createHyperbolicRepresentation(): HyperbolicStructure {
        // By uniformization theorem, any Riemann surface of genus > 1
        // can be represented as H/Γ where H is hyperbolic plane
        // and Γ is a Fuchsian group

        // 1. Compute fundamental domain
        const domain = this.computeFundamentalDomain();

        // 2. Find generators of Fuchsian group
        const generators = this.computeFuchsianGenerators();

        // 3. Create T-mesh structure in hyperbolic plane
        const hyperbolicMesh = this.createHyperbolicTMesh(domain, generators);

        return {
            fundamentalDomain: domain,
            fuchsianGroup: generators,
            mesh: hyperbolicMesh
        };
    }

    private computeFundamentalDomain(): HyperbolicPolygon {
        // For genus g, create 4g-sided polygon in hyperbolic plane
        const vertices: Complex[] = [];
        const numSides = 4 * this.genus;

        for (let i = 0; i < numSides; i++) {
            const vertex = this.computeFundamentalVertex(i, numSides);
            vertices.push(vertex);
        }

        return new HyperbolicPolygon(vertices);
    }

    private computeFuchsianGenerators(): MobiusTransformation[] {
        // For genus g, need 2g generators
        const generators: MobiusTransformation[] = [];

        // Compute standard generators (a₁,b₁,...,aₘ,bₘ)
        // satisfying Π[aᵢ,bᵢ] = 1
        for (let i = 0; i < this.genus; i++) {
            const [a, b] = this.computeGeneratorPair(i);
            generators.push(a, b);
        }

        return generators;
    }

    private createHyperbolicTMesh(
        domain: HyperbolicPolygon,
        generators: MobiusTransformation[]
    ): HyperbolicTMesh {
        const mesh = new HyperbolicTMesh(HyperbolicModel.POINCARE_DISK);

        // Create base mesh in fundamental domain
        this.createBaseMesh(mesh, domain);

        // Apply group actions to create full mesh
        this.extendMeshByGroupAction(mesh, generators);

        return mesh;
    }
}

class HyperbolicSurfaceTMesh {
    constructor(
        private genus: number,
        private fuchsianGroup: FuchsianGroup
    ) {}

    createMesh(): TMesh {
        // Create T-mesh structure respecting hyperbolic geometry
        const mesh = new TMesh();

        // Set up fundamental domain
        const domain = this.createFundamentalDomain();

        // Create T-mesh structure in fundamental domain
        this.createTMeshStructure(domain, mesh);

        // Handle identifications by group action
        this.handleBoundaryIdentifications(mesh);

        return mesh;
    }

    private createTMeshStructure(domain: HyperbolicPolygon, mesh: TMesh): void {
        // Create vertices following hyperbolic geometry
        this.createVertices(domain, mesh);

        // Create edges as geodesic segments
        this.createGeodesicEdges(mesh);

        // Create faces respecting hyperbolic structure
        this.createHyperbolicFaces(mesh);
    }

    private handleBoundaryIdentifications(mesh: TMesh): void {
        // Implement side pairings according to Fuchsian group
        for (const generator of this.fuchsianGroup.generators) {
            this.identifyBoundaries(mesh, generator);
        }
    }
}

class HyperbolicGeometryAnalyzer {
    constructor(private surface: TSplineSurface) {}

    analyzeHyperbolicStructure(): HyperbolicProperties {
        // Compute fundamental geometric properties
        const gaussianCurvature = this.computeGaussianCurvature();
        const genus = this.computeGenus();
        
        // For genus > 1, analyze hyperbolic structure
        if (genus > 1) {
            return {
                genus,
                gaussianCurvature,
                hyperbolicMetric: this.computeHyperbolicMetric(),
                fuchsianGroup: this.computeFuchsianGroup(),
                fundamentalDomain: this.computeFundamentalDomain()
            };
        }

        throw new Error("Surface is not hyperbolic");
    }

    private computeHyperbolicMetric(): MetricTensor {
        // Compute hyperbolic metric induced on surface
        const metric = new MetricTensor();
        
        // Use uniformization theorem to find local hyperbolic structure
        this.computeUniformization();
        
        return metric;
    }

    private computeUniformization(): void {
        // Implement uniformization algorithm
        // This gives conformal map to H/Γ
    }
}

class UniformizationVisualizer {
    constructor(
        private surface: TSplineSurface,
        private canvas: HTMLCanvasElement
    ) {}

    visualizeHyperbolicStructure(): void {
        const analyzer = new HyperbolicGeometryAnalyzer(this.surface);
        const properties = analyzer.analyzeHyperbolicStructure();

        // Draw fundamental domain
        this.drawFundamentalDomain(properties.fundamentalDomain);

        // Draw T-mesh structure
        this.drawHyperbolicTMesh(properties);

        // Visualize group actions
        this.visualizeGroupActions(properties.fuchsianGroup);
    }

    private drawHyperbolicTMesh(properties: HyperbolicProperties): void {
        const ctx = this.canvas.getContext('2d');
        if (!ctx) return;

        // Draw geodesic edges
        for (const edge of properties.mesh.edges) {
            this.drawGeodesicEdge(ctx, edge);
        }

        // Draw vertices
        for (const vertex of properties.mesh.vertices) {
            this.drawVertex(ctx, vertex);
        }

        // Draw parameter lines
        this.drawParameterLines(ctx, properties);
    }
}

// Example usage
const surface = new TSplineSurface(/* ... */);
const genus = 2; // Example: genus 2 surface
const analyzer = new HigherGenusAnalysis(surface, genus);
const hyperbolicStructure = analyzer.createHyperbolicRepresentation();

// Visualize
const canvas = document.getElementById('canvas') as HTMLCanvasElement;
const visualizer = new UniformizationVisualizer(surface, canvas);
visualizer.visualizeHyperbolicStructure();

class SphereProjectionAnalyzer {
    constructor(
        private mesh: TMesh,
        private tolerance: number = 1e-6
    ) {}

    canProjectToSphere(): boolean {
        // 1. Check topological requirements
        if (!this.isGenus0()) {
            return false;
        }

        // 2. Check for topological obstructions
        return this.checkProjectability();
    }

    private isGenus0(): boolean {
        // Use Euler characteristic
        const V = this.mesh.vertices.length;
        const E = this.mesh.edges.length;
        const F = this.mesh.faces.length;
        
        // χ = V - E + F should be 2 for genus 0
        return V - E + F === 2;
    }

    private checkProjectability(): boolean {
        // Check for topological obstacles
        // Returns true if projection is possible
        return !this.hasTopologicalObstructions();
    }
}

class SphereProjection {
    constructor(
        private mesh: TMesh,
        private method: SphereProjectionMethod = SphereProjectionMethod.TUTTE
    ) {}

    project(): SphericalTMesh {
        switch (this.method) {
            case SphereProjectionMethod.TUTTE:
                return this.tutteEmbedding();
            case SphereProjectionMethod.PROGRESSIVE:
                return this.progressiveMapping();
            case SphereProjectionMethod.HARMONIC:
                return this.harmonicMapping();
        }
    }

    private tutteEmbedding(): SphericalTMesh {
        // 1. Initial projection to convex shape
        const initialMesh = this.createInitialConvexShape();
        
        // 2. Move towards sphere while maintaining no self-intersections
        return this.projectToSphere(initialMesh);
    }

    private progressiveMapping(): SphericalTMesh {
        const sphericalMesh = new SphericalTMesh();
        
        // 1. Start with a base triangle
        this.initializeBaseTriangle(sphericalMesh);
        
        // 2. Progressively add vertices
        for (const vertex of this.getOrderedVertices()) {
            this.addVertexToSphericalMesh(sphericalMesh, vertex);
        }
        
        return sphericalMesh;
    }

    private harmonicMapping(): SphericalTMesh {
        // Use harmonic map to sphere
        const sphericalMesh = new SphericalTMesh();
        
        // 1. Set up harmonic energy
        const energy = new HarmonicEnergy(this.mesh);
        
        // 2. Minimize energy while maintaining spherical constraint
        this.minimizeEnergyOnSphere(energy, sphericalMesh);
        
        return sphericalMesh;
    }

    private createInitialConvexShape(): TMesh {
        // Create initial convex embedding
        const convexMesh = new TMesh();
        
        // Use Tutte embedding for initial planar layout
        const tutteEmbedding = new TutteEmbedding(this.mesh);
        const planarMesh = tutteEmbedding.compute();
        
        // Lift to 3D convex shape
        return this.liftToPseudosphere(planarMesh);
    }

    private projectToSphere(mesh: TMesh): SphericalTMesh {
        const sphericalMesh = new SphericalTMesh();
        
        // Iteratively move vertices to sphere while:
        // 1. Maintaining no self-intersections
        // 2. Minimizing distortion
        let currentMesh = mesh;
        while (!this.isConverged(currentMesh)) {
            currentMesh = this.stepTowardsSphere(currentMesh);
            
            // Check and prevent self-intersections
            if (this.detectSelfIntersections(currentMesh)) {
                currentMesh = this.resolveIntersections(currentMesh);
            }
        }
        
        return sphericalMesh;
    }
}

class IntersectionChecker {
    constructor(private mesh: TMesh) {}

    checkSelfIntersections(): boolean {
        // Check for face-face intersections
        for (let i = 0; i < this.mesh.faces.length; i++) {
            for (let j = i + 1; j < this.mesh.faces.length; j++) {
                if (this.facesIntersect(
                    this.mesh.faces[i],
                    this.mesh.faces[j]
                )) {
                    return true;
                }
            }
        }
        return false;
    }

    private facesIntersect(face1: Face, face2: Face): boolean {
        // Implement triangle-triangle intersection test
        // Use geometric predicates for robustness
        return this.triangleIntersectionTest(face1, face2);
    }
}

class SphereProjectionOptimizer {
    constructor(
        private mesh: TMesh,
        private maxIterations: number = 1000
    ) {}

    optimize(): SphericalTMesh {
        let currentMesh = this.initializeProjection();
        
        for (let i = 0; i < this.maxIterations; i++) {
            // 1. Compute gradient of distortion energy
            const gradient = this.computeGradient(currentMesh);
            
            // 2. Project gradient to maintain spherical constraint
            const projectedGradient = this.projectToSphere(gradient);
            
            // 3. Take step while preventing self-intersections
            const nextMesh = this.takeStep(currentMesh, projectedGradient);
            
            // 4. Check convergence
            if (this.hasConverged(currentMesh, nextMesh)) {
                return this.createSphericalMesh(nextMesh);
            }
            
            currentMesh = nextMesh;
        }
        
        throw new Error("Failed to converge");
    }

    private computeGradient(mesh: TMesh): Vector3D[] {
        // Compute energy gradient for each vertex
        return mesh.vertices.map(vertex => 
            this.computeVertexGradient(vertex, mesh)
        );
    }

    private projectToSphere(gradient: Vector3D[]): Vector3D[] {
        // Project gradient to maintain spherical constraint
        return gradient.map(g => this.projectVectorToSphere(g));
    }

    private takeStep(
        mesh: TMesh,
        gradient: Vector3D[]
    ): TMesh {
        const newMesh = mesh.clone();
        
        // Move vertices while checking for intersections
        for (let i = 0; i < mesh.vertices.length; i++) {
            const newPosition = this.computeNewPosition(
                mesh.vertices[i],
                gradient[i]
            );
            
            if (!this.causesIntersection(newMesh, i, newPosition)) {
                newMesh.vertices[i].position = newPosition;
            }
        }
        
        return newMesh;
    }
}

// Example usage
const mesh = new TMesh(/* ... */);
const analyzer = new SphereProjectionAnalyzer(mesh);

if (analyzer.canProjectToSphere()) {
    const projector = new SphereProjection(mesh);
    const sphericalMesh = projector.project();
    
    // Optimize the projection
    const optimizer = new SphereProjectionOptimizer(sphericalMesh);
    const optimizedMesh = optimizer.optimize();
    
    console.log("Successfully projected to sphere!");
} else {
    console.log("Cannot project to sphere without self-intersections");
}

class CircularArcTMesh {
    constructor(
        private vertices: Vertex[],
        private circularEdges: CircularArc[],
        private patches: CircularPatch[]
    ) {}

    // Create circular arc grid
    createCircularGrid(): void {
        // Create grid of circular arcs that meet at right angles
        // (orthogonal circle patterns)
        this.createOrthogonalCirclePattern();
    }
}

class CircularArc {
    constructor(
        private center: Point2D,
        private radius: number,
        private startAngle: number,
        private endAngle: number
    ) {}

    // Check if arc preserves angles at intersection
    isConformal(other: CircularArc): boolean {
        // Two circles intersect at equal angles
        return this.computeIntersectionAngle(other) === Math.PI/2;
    }

    // Get intersection points with another arc
    intersectWith(other: CircularArc): Point2D[] {
        const d = this.center.distanceTo(other.center);
        
        // Using circle intersection formula
        if (d > this.radius + other.radius) return []; // No intersection
        
        // Compute intersection points
        return this.computeCircleIntersections(other);
    }
}

class CircularPatch {
    constructor(
        private boundaryArcs: CircularArc[],
        private innerStructure: CircularArcGrid
    ) {}

    // Create conformal structure inside patch
    createConformalStructure(): void {
        // Use circle patterns that preserve angles
        this.createOrthogonalCirclePattern();
    }

    // Check angle preservation
    checkAnglePreservation(): AngleMetrics {
        const metrics = new AngleMetrics();
        
        // Check all intersections except star points
        for (const intersection of this.getIntersections()) {
            if (!this.isStarPoint(intersection)) {
                metrics.addMeasurement(
                    this.measureAngleAt(intersection)
                );
            }
        }
        
        return metrics;
    }
}

class CircularArcGrid {
    constructor(
        private domain: Domain,
        private gridSpacing: number
    ) {}

    generateGrid(): CircularArc[] {
        const arcs: CircularArc[] = [];
        
        // Create two families of circular arcs
        this.createHorizontalCircleFamily(arcs);
        this.createVerticalCircleFamily(arcs);
        
        return arcs;
    }

    private createHorizontalCircleFamily(arcs: CircularArc[]): void {
        // Create family of circles that will become "horizontal" lines
        for (let y = this.domain.minY; y <= this.domain.maxY; y += this.gridSpacing) {
            const circle = this.createHorizontalCircle(y);
            arcs.push(circle);
        }
    }

    private createVerticalCircleFamily(arcs: CircularArc[]): void {
        // Create orthogonal family of circles
        for (let x = this.domain.minX; x <= this.domain.maxX; x += this.gridSpacing) {
            const circle = this.createVerticalCircle(x);
            arcs.push(circle);
        }
    }
}

class StarPointHandler {
    constructor(
        private starPoint: Point2D,
        private valence: number
    ) {}

    createCircularStructure(): CircularArc[] {
        const arcs: CircularArc[] = [];
        
        // Create circular arcs meeting at star point
        const angleStep = 2 * Math.PI / this.valence;
        
        for (let i = 0; i < this.valence; i++) {
            const angle = i * angleStep;
            arcs.push(this.createRadialArc(angle));
        }
        
        return arcs;
    }

    private createRadialArc(angle: number): CircularArc {
        // Create circular arc emanating from star point
        // at given angle
        return new CircularArc(
            this.computeArcCenter(angle),
            this.computeArcRadius(angle),
            angle,
            this.computeEndAngle(angle)
        );
    }
}

class ConformalTMeshOptimizer {
    constructor(
        private mesh: CircularArcTMesh,
        private tolerance: number = 1e-6
    ) {}

    optimize(): void {
        let improved = true;
        while (improved) {
            improved = false;
            
            // Optimize circle positions and radii
            improved = this.optimizeCircleParameters() || improved;
            
            // Maintain angle preservation where required
            improved = this.adjustForAnglePreservation() || improved;
            
            // Handle star points specially
            this.handleStarPoints();
        }
    }

    private optimizeCircleParameters(): boolean {
        let improved = false;
        
        // Optimize each circular arc while maintaining
        // orthogonality constraints
        for (const arc of this.mesh.getCircularArcs()) {
            if (this.optimizeArc(arc)) {
                improved = true;
            }
        }
        
        return improved;
    }

    private adjustForAnglePreservation(): boolean {
        // Adjust circles to maintain orthogonality
        // except at star points
        return this.enforceOrthogonality();
    }
}

// Example usage
const domain = new Domain(/* bounds */);
const gridSpacing = 0.1;

// Create circular arc grid
const grid = new CircularArcGrid(domain, gridSpacing);
const arcs = grid.generateGrid();

// Create T-mesh with circular arcs
const mesh = new CircularArcTMesh(/* vertices */, arcs, /* patches */);

// Handle star points
const starPoints = mesh.findStarPoints();
for (const starPoint of starPoints) {
    const handler = new StarPointHandler(
        starPoint.position,
        starPoint.valence
    );
    handler.createCircularStructure();
}

// Optimize for conformality
const optimizer = new ConformalTMeshOptimizer(mesh);
optimizer.optimize();

// Check angle preservation
for (const patch of mesh.getPatches()) {
    const angleMetrics = patch.checkAnglePreservation();
    console.log("Angle preservation metrics:", angleMetrics);
}


class MobiusTransformation {
    constructor(
        private matrix: [[Complex, Complex], [Complex, Complex]]
    ) {}

    // (az + b)/(cz + d) form
    apply(z: Complex): Complex {
        const [[a, b], [c, d]] = this.matrix;
        return a.multiply(z).add(b).divide(
            c.multiply(z).add(d)
        );
    }

    // Transform a circular arc
    transformArc(arc: CircularArc): CircularArc {
        // Möbius transformation maps circles to circles
        const startPoint = this.apply(arc.getStartPoint());
        const endPoint = this.apply(arc.getEndPoint());
        const midPoint = this.apply(arc.getMidPoint());
        
        return CircularArc.fromThreePoints(
            startPoint, 
            midPoint, 
            endPoint
        );
    }
}

class CircularArcTMeshTransformer {
    constructor(
        private mesh: CircularArcTMesh,
        private mobius: MobiusTransformation
    ) {}

    transform(): CircularArcTMesh {
        const transformedMesh = new CircularArcTMesh();

        // Transform regular vertices
        this.transformRegularVertices(transformedMesh);

        // Keep star points fixed or handle specially
        this.handleStarPoints(transformedMesh);

        // Transform circular arcs
        this.transformArcs(transformedMesh);

        return transformedMesh;
    }

    private transformRegularVertices(transformedMesh: CircularArcTMesh): void {
        for (const vertex of this.mesh.vertices) {
            if (!this.isStarPoint(vertex)) {
                const transformedPos = this.mobius.apply(vertex.position);
                transformedMesh.addVertex(new Vertex(transformedPos));
            } else {
                // Keep star point as is or handle specially
                transformedMesh.addVertex(vertex.clone());
            }
        }
    }

    private transformArcs(transformedMesh: CircularArcTMesh): void {
        for (const arc of this.mesh.circularArcs) {
            if (!this.isConnectedToStarPoint(arc)) {
                // Regular arc transformation
                const transformedArc = this.mobius.transformArc(arc);
                transformedMesh.addArc(transformedArc);
            } else {
                // Special handling for arcs connected to star points
                this.handleStarPointArc(arc, transformedMesh);
            }
        }
    }

    private handleStarPointArc(
        arc: CircularArc, 
        transformedMesh: CircularArcTMesh
    ): void {
        // For arcs connected to star points:
        // 1. Keep star point end fixed
        // 2. Transform only the regular end
        // 3. Create new arc maintaining angle conditions
        const starPoint = this.getStarPointEnd(arc);
        const regularEnd = this.getRegularEnd(arc);
        const transformedEnd = this.mobius.apply(regularEnd);

        const newArc = this.createTransitionArc(
            starPoint,
            transformedEnd
        );
        transformedMesh.addArc(newArc);
    }

    private createTransitionArc(
        starPoint: Point2D,
        transformedEnd: Point2D
    ): CircularArc {
        // Create arc that:
        // 1. Maintains original angle at star point
        // 2. Connects smoothly to transformed region
        return new CircularArc(
            this.computeTransitionCenter(starPoint, transformedEnd),
            this.computeTransitionRadius(starPoint, transformedEnd),
            this.computeStartAngle(starPoint),
            this.computeEndAngle(transformedEnd)
        );
    }
}

class ConformallyCriticalTMesh {
    constructor(
        private baseMesh: CircularArcTMesh,
        private starPoints: StarPoint[]
    ) {}

    applyMobiusTransformation(mobius: MobiusTransformation): void {
        const transformer = new CircularArcTMeshTransformer(
            this.baseMesh,
            mobius
        );

        // Transform mesh while preserving conformal structure
        this.baseMesh = transformer.transform();

        // Verify angle preservation
        this.verifyConformality();
    }

    private verifyConformality(): void {
        const angleChecker = new AnglePreservationChecker(this.baseMesh);
        
        for (const intersection of this.baseMesh.getIntersections()) {
            if (!this.isStarPoint(intersection)) {
                angleChecker.checkAngleAt(intersection);
            }
        }
    }
}

// Example usage
const mesh = new CircularArcTMesh(/* ... */);

// Create a Möbius transformation
const mobius = new MobiusTransformation([
    [new Complex(1, 0), new Complex(0, 1)],
    [new Complex(0, -1), new Complex(1, 0)]
]);

// Transform the mesh
const transformer = new CircularArcTMeshTransformer(mesh, mobius);
const transformedMesh = transformer.transform();

// Verify properties
const conformality = new ConformallyCriticalTMesh(
    transformedMesh,
    mesh.getStarPoints()
);
conformality.verifyConformality();

// Visualize
const visualizer = new CircularArcMeshVisualizer(transformedMesh);
visualizer.draw();


class HyperbolicTSplineProjection {
    constructor(
        private surface: TSplineSurface,
        private genus: number,
        private model: HyperbolicModel = HyperbolicModel.POINCARE_DISK
    ) {
        if (genus < 2) {
            throw new Error("Surface must have genus ≥ 2 for hyperbolic projection");
        }
    }

    project(): HyperbolicTMesh {
        // 1. Compute fundamental domain
        const fundamentalDomain = this.computeFundamentalDomain();

        // 2. Find Fuchsian group generators
        const generators = this.computeFuchsianGroup();

        // 3. Project T-spline mesh to hyperbolic plane
        return this.projectToHyperbolicPlane(fundamentalDomain, generators);
    }

    private computeFundamentalDomain(): HyperbolicPolygon {
        // For genus g, create 4g-sided polygon
        const vertices: Complex[] = [];
        const numSides = 4 * this.genus;

        // Compute vertices in hyperbolic plane
        for (let i = 0; i < numSides; i++) {
            vertices.push(this.computeFundamentalVertex(i, numSides));
        }

        return new HyperbolicPolygon(vertices);
    }

    private projectToHyperbolicPlane(
        domain: HyperbolicPolygon,
        generators: MobiusTransformation[]
    ): HyperbolicTMesh {
        const mesh = new HyperbolicTMesh(this.model);

        // Project vertices while handling star points specially
        this.projectVertices(mesh);

        // Create circular arc edges
        this.createHyperbolicEdges(mesh);

        // Handle identifications by Fuchsian group
        this.handleFuchsianIdentifications(mesh, generators);

        return mesh;
    }
}

class HyperbolicStarPointHandler {
    constructor(
        private starPoint: StarPoint,
        private model: HyperbolicModel
    ) {}

    createLocalStructure(): HyperbolicCircularArcs[] {
        const arcs: HyperbolicCircularArcs[] = [];
        const valence = this.starPoint.valence;

        // Create geodesic rays in hyperbolic plane
        for (let i = 0; i < valence; i++) {
            const angle = (2 * Math.PI * i) / valence;
            arcs.push(this.createHyperbolicRay(angle));
        }

        return arcs;
    }

    private createHyperbolicRay(angle: number): HyperbolicCircularArc {
        if (this.model === HyperbolicModel.POINCARE_DISK) {
            return this.createPoincareDiskRay(angle);
        } else {
            return this.createUpperHalfPlaneRay(angle);
        }
    }

    private createPoincareDiskRay(angle: number): HyperbolicCircularArc {
        // Create circular arc perpendicular to boundary circle
        const center = this.computeCircleCenter(angle);
        const radius = this.computeCircleRadius(angle);
        
        return new HyperbolicCircularArc(
            center,
            radius,
            this.starPoint.position,
            this.computeIdealEndpoint(angle)
        );
    }
}

class HyperbolicMobiusTransformation {
    constructor(
        private matrix: [[Complex, Complex], [Complex, Complex]],
        private model: HyperbolicModel
    ) {}

    // Preserve hyperbolic distance and angles
    apply(z: Complex): Complex {
        const transformed = this.applyMobius(z);
        return this.normalizeToModel(transformed);
    }

    // Transform hyperbolic geodesic
    transformGeodesic(arc: HyperbolicCircularArc): HyperbolicCircularArc {
        const start = this.apply(arc.getStartPoint());
        const end = this.apply(arc.getEndPoint());
        const mid = this.apply(arc.getMidPoint());

        return HyperbolicCircularArc.throughThreePoints(
            start, mid, end, this.model
        );
    }

    private normalizeToModel(z: Complex): Complex {
        if (this.model === HyperbolicModel.POINCARE_DISK) {
            // Ensure point stays in unit disk
            if (z.abs() >= 1) {
                return z.scale(0.99 / z.abs());
            }
        } else {
            // Ensure point stays in upper half-plane
            if (z.imaginary <= 0) {
                return new Complex(z.real, 0.01);
            }
        }
        return z;
    }
}

class HyperbolicTMeshProjector {
    constructor(
        private mesh: TSplineSurface,
        private genus: number,
        private model: HyperbolicModel
    ) {}

    project(): HyperbolicTMesh {
        // 1. Initialize hyperbolic structure
        const hyperbolicMesh = new HyperbolicTMesh(this.model);

        // 2. Handle regular vertices
        this.projectRegularVertices(hyperbolicMesh);

        // 3. Handle star points specially
        this.handleStarPoints(hyperbolicMesh);

        // 4. Create hyperbolic geodesic edges
        this.createGeodesicEdges(hyperbolicMesh);

        return hyperbolicMesh;
    }

    private handleStarPoints(mesh: HyperbolicTMesh): void {
        for (const starPoint of this.mesh.getStarPoints()) {
            const handler = new HyperbolicStarPointHandler(
                starPoint,
                this.model
            );
            
            // Create local hyperbolic structure around star point
            const localStructure = handler.createLocalStructure();
            mesh.addStarPointStructure(starPoint, localStructure);
        }
    }

    private createGeodesicEdges(mesh: HyperbolicTMesh): void {
        for (const edge of this.mesh.edges) {
            if (this.isRegularEdge(edge)) {
                // Create hyperbolic geodesic
                mesh.addGeodesic(
                    this.createHyperbolicGeodesic(edge)
                );
            } else {
                // Handle edges connected to star points
                this.handleStarPointEdge(edge, mesh);
            }
        }
    }
}

// Example usage
const surface = new TSplineSurface(/* ... */);
const genus = 2;

const projector = new HyperbolicTSplineProjection(
    surface,
    genus,
    HyperbolicModel.POINCARE_DISK
);

// Project to hyperbolic plane
const hyperbolicMesh = projector.project();

// Apply Möbius transformation
const mobius = new HyperbolicMobiusTransformation(
    [[new Complex(1, 0.5), new Complex(0, 1)],
     [new Complex(0, -1), new Complex(1, -0.5)]],
    HyperbolicModel.POINCARE_DISK
);

// Transform while preserving hyperbolic structure
const transformer = new HyperbolicMeshTransformer(
    hyperbolicMesh,
    mobius
);
const transformedMesh = transformer.transform();

// Visualize
const visualizer = new HyperbolicMeshVisualizer(transformedMesh);
visualizer.draw();


class FuchsianIdentificationHandler {
    constructor(
        private mesh: HyperbolicTMesh,
        private generators: MobiusTransformation[]
    ) {}

    handleIdentifications(): void {
        // 1. Identify paired edges of fundamental domain
        this.identifyPairedEdges();

        // 2. Create orbit data structure
        this.createVertexOrbits();

        // 3. Handle transition functions between patches
        this.createTransitionMaps();
    }

    private identifyPairedEdges(): void {
        // For genus g surface, we have 4g edges paired as (a₁,b₁,...,aₘ,bₘ)
        for (let i = 0; i < this.generators.length; i += 2) {
            const generatorA = this.generators[i];
            const generatorB = this.generators[i + 1];

            // Handle edge pairs
            this.identifyEdgePair(
                this.mesh.getBoundaryEdge(2*i),
                this.mesh.getBoundaryEdge(2*i + 1),
                generatorA
            );

            this.identifyEdgePair(
                this.mesh.getBoundaryEdge(2*i + 2),
                this.mesh.getBoundaryEdge(2*i + 3),
                generatorB
            );
        }
    }

    private createVertexOrbits(): void {
        const orbitMap = new Map<Vertex, VertexOrbit>();

        // For each vertex, compute its orbit under the group action
        for (const vertex of this.mesh.vertices) {
            if (!orbitMap.has(vertex)) {
                const orbit = this.computeVertexOrbit(vertex);
                for (const v of orbit.vertices) {
                    orbitMap.set(v, orbit);
                }
            }
        }

        this.mesh.setVertexOrbits(orbitMap);
    }

    private computeVertexOrbit(vertex: Vertex): VertexOrbit {
        const orbit = new VertexOrbit();
        const visited = new Set<Vertex>();
        const queue: Vertex[] = [vertex];

        while (queue.length > 0) {
            const current = queue.shift()!;
            if (visited.has(current)) continue;
            
            visited.add(current);
            orbit.addVertex(current);

            // Apply all generators and their inverses
            for (const generator of this.generators) {
                const transformed = this.applyGenerator(current, generator);
                if (!visited.has(transformed)) {
                    queue.push(transformed);
                }
            }
        }

        return orbit;
    }

    private createTransitionMaps(): void {
        // Create transition functions between adjacent fundamental domains
        const transitions = new Map<EdgePair, MobiusTransformation>();

        for (const edge of this.mesh.getBoundaryEdges()) {
            const pairedEdge = this.findPairedEdge(edge);
            if (pairedEdge) {
                const transition = this.computeTransition(edge, pairedEdge);
                transitions.set({edge, pairedEdge}, transition);
            }
        }

        this.mesh.setTransitionMaps(transitions);
    }

    private computeTransition(
        edge1: HyperbolicEdge,
        edge2: HyperbolicEdge
    ): MobiusTransformation {
        // Find Möbius transformation taking edge1 to edge2
        const points1 = edge1.getEndpoints();
        const points2 = edge2.getEndpoints();

        return MobiusTransformation.fromPointPairs(
            points1,
            points2,
            this.mesh.model
        );
    }
}

class VertexOrbit {
    private vertices: Set<Vertex> = new Set();
    private transformations: Map<Vertex, MobiusTransformation> = new Map();

    addVertex(vertex: Vertex, transformation?: MobiusTransformation): void {
        this.vertices.add(vertex);
        if (transformation) {
            this.transformations.set(vertex, transformation);
        }
    }

    getTransformationTo(from: Vertex, to: Vertex): MobiusTransformation {
        // Get transformation taking one vertex to another in orbit
        const fromTrans = this.transformations.get(from);
        const toTrans = this.transformations.get(to);
        
        if (!fromTrans || !toTrans) {
            throw new Error("Vertices not in orbit");
        }

        return toTrans.compose(fromTrans.inverse());
    }
}

class HyperbolicPatchTransition {
    constructor(
        private source: HyperbolicPatch,
        private target: HyperbolicPatch,
        private transformation: MobiusTransformation
    ) {}

    mapPoint(point: Complex): Complex {
        return this.transformation.apply(point);
    }

    mapTangentVector(
        point: Complex,
        vector: Complex
    ): Complex {
        // Transform tangent vector according to Möbius transformation
        const derivative = this.transformation.derivative(point);
        return derivative.multiply(vector);
    }
}

class HyperbolicTMesh {
    private vertexOrbits: Map<Vertex, VertexOrbit> = new Map();
    private transitions: Map<EdgePair, MobiusTransformation> = new Map();

    setVertexOrbits(orbits: Map<Vertex, VertexOrbit>): void {
        this.vertexOrbits = orbits;
    }

    setTransitionMaps(transitions: Map<EdgePair, MobiusTransformation>): void {
        this.transitions = transitions;
    }

    // Get all equivalent points under group action
    getEquivalentPoints(point: Complex): Complex[] {
        const equivalentPoints: Complex[] = [];

        for (const generator of this.getGenerators()) {
            equivalentPoints.push(generator.apply(point));
        }

        return equivalentPoints;
    }

    // Check if two points are equivalent under group action
    areEquivalentPoints(p1: Complex, p2: Complex): boolean {
        for (const generator of this.getGenerators()) {
            if (generator.apply(p1).equals(p2, 1e-10)) {
                return true;
            }
        }
        return false;
    }
}

// Example usage
const genus = 2;
const mesh = new HyperbolicTMesh(HyperbolicModel.POINCARE_DISK);
const generators = computeFuchsianGenerators(genus);

const handler = new FuchsianIdentificationHandler(mesh, generators);
handler.handleIdentifications();

// Verify identifications
for (const edge of mesh.getBoundaryEdges()) {
    const pairedEdge = mesh.getPairedEdge(edge);
    if (pairedEdge) {
        console.log(
            "Edge pair identified by transformation:",
            mesh.getTransformation(edge, pairedEdge)
        );
    }
}

class RiemannSurfaceAnalogy {
    // Key concepts from Riemann surface theory applicable to T-splines:
    
    // 1. Local Charts and Transitions
    class LocalParametrization {
        constructor(
            private patch: TSplinePatch,
            private transitionMaps: Map<PatchBoundary, TransitionFunction>
        ) {}

        // Analogous to complex charts on Riemann surface
        getLocalParameters(point: Point3D): UV {
            return this.patch.getParameters(point);
        }

        // Similar to transition functions between charts
        transitionToPatch(
            uv: UV, 
            targetPatch: TSplinePatch
        ): UV {
            const transition = this.transitionMaps.get(
                this.getBoundaryWith(targetPatch)
            );
            return transition ? transition(uv) : null;
        }
    }

    // 2. Branching Behavior (Star Points)
    class BranchPointAnalysis {
        constructor(
            private starPoint: StarPoint,
            private valence: number
        ) {}

        // Analogous to branch points in Riemann surfaces
        analyzeLocalStructure(): BranchStructure {
            return {
                branchOrder: this.valence - 4, // Regular points have valence 4
                localMonodromy: this.computeMonodromy(),
                ramificationIndex: this.computeRamificationIndex()
            };
        }

        private computeMonodromy(): MobiusTransformation {
            // Compute how parameters transform when going around star point
            // Similar to monodromy around branch points
            return this.calculateLocalParameterTransformation();
        }
    }

    // 3. Global Topology
    class GlobalStructureAnalysis {
        constructor(
            private surface: TSplineSurface,
            private genus: number
        ) {}

        analyzeGlobalStructure(): GlobalStructure {
            return {
                eulerCharacteristic: this.computeEulerCharacteristic(),
                fundamentalGroup: this.computeFundamentalGroup(),
                canonicalForm: this.computeCanonicalForm()
            };
        }

        private computeEulerCharacteristic(): number {
            const V = this.surface.vertices.length;
            const E = this.surface.edges.length;
            const F = this.surface.faces.length;
            return V - E + F;
        }
    }
}

class ConformalTSplineTheory {
    // Applications of conformal theory to T-splines

    class ConformalStructure {
        constructor(
            private surface: TSplineSurface,
            private starPoints: StarPoint[]
        ) {}

        analyzeConformalStructure(): ConformalAnalysis {
            return {
                // Similar to complex structure on Riemann surface
                localParameters: this.computeConformalParameters(),
                holomorphicForms: this.findHolomorphicForms(),
                periodMatrix: this.computePeriodMatrix()
            };
        }

        private computeConformalParameters(): LocalParametrization[] {
            // Create conformal parameterization away from star points
            const params: LocalParametrization[] = [];
            
            for (const patch of this.surface.patches) {
                if (!this.containsStarPoint(patch)) {
                    params.push(this.createConformalParam(patch));
                }
            }
            
            return params;
        }
    }

    class HolomorphicStructure {
        // Study holomorphic properties of T-spline maps

        analyzeHolomorphicBehavior(): HolomorphicProperties {
            return {
                branchPoints: this.findBranchPoints(),
                zeros: this.findZeros(),
                poles: this.findPoles(),
                residues: this.computeResidues()
            };
        }
    }
}

class UniformizationApplication {
    // Apply uniformization theorem concepts to T-splines

    constructor(
        private surface: TSplineSurface,
        private genus: number
    ) {}

    createUniformStructure(): UniformStructure {
        if (this.genus === 0) {
            return this.createSphericalStructure();
        } else if (this.genus === 1) {
            return this.createToroidalStructure();
        } else {
            return this.createHyperbolicStructure();
        }
    }

    private createHyperbolicStructure(): HyperbolicStructure {
        // For genus > 1:
        // Use Fuchsian group representation
        const fuchsianGroup = this.computeFuchsianGroup();
        return new HyperbolicStructure(
            this.surface,
            fuchsianGroup
        );
    }
}

class TeichmullerTheoryApplication {
    // Apply Teichmüller theory to T-spline deformations

    analyzeDeformationSpace(): TeichmullerSpace {
        return {
            moduliSpace: this.computeModuliSpace(),
            beltramiDifferentials: this.computeBeltramiDiffs(),
            quadraticDifferentials: this.findQuadraticDiffs()
        };
    }

    optimizeParametrization(): OptimalParametrization {
        // Use Teichmüller theory to find optimal parameterization
        return this.minimizeTeichmullerDistortion();
    }
}

class LocalStructureAnalogy {
    // Teaches us about:
    // - Local parameter behavior
    // - Transition functions between patches
    // - Regularity conditions
}


class BranchPointTheory {
    // Insights about:
    // - Local structure around star points
    // - Ramification and branching order
    // - Monodromy transformations
}


class GlobalTopologyInsights {
    // Understanding of:
    // - Topological constraints
    // - Global existence conditions
    // - Period relations
}


class UniformizationPrinciples {
    // Provides:
    // - Canonical forms for different genera
    // - Universal covering spaces
    // - Group actions and quotients
}

// 1. Specific applications of conformal theory to T-splines?

// 2. Star point analysis using branching theory?

// 3. Global structure theorems for T-splines?

// 4. Deformation theory applications?

class ConeApexStructure {
    constructor(
        private apex: Point3D,
        private baseRadius: number,
        private height: number,
        private numSectors: number = 4  // Can be adjusted for different refinement levels
    ) {}

    createTSplineStructure(): TSplineMesh {
        const mesh = new TSplineMesh();
        
        // 1. Create apex as special control point
        const apexVertex = this.createApexVertex();
        mesh.addVertex(apexVertex);

        // 2. Create base control points
        const baseVertices = this.createBaseVertices();
        baseVertices.forEach(v => mesh.addVertex(v));

        // 3. Create intermediate rings for better shape control
        const intermediateRings = this.createIntermediateRings();
        intermediateRings.forEach(ring => 
            ring.forEach(v => mesh.addVertex(v))
        );

        // 4. Create T-spline connectivity
        this.createTSplineTopology(mesh);

        return mesh;
    }

    private createApexVertex(): TSplineVertex {
        // Apex needs special weights and knot intervals
        return new TSplineVertex({
            position: this.apex,
            weight: 1.0,  // Special weight for apex
            // Use multiple coincident knots at apex
            knotIntervals: new KnotIntervals({
                u: [0, 0, 0, 0],  // Multiplicity 4 for cubic
                v: [0, 0, 0, 0]
            })
        });
    }

    private createBaseVertices(): TSplineVertex[] {
        const vertices: TSplineVertex[] = [];
        const angleStep = (2 * Math.PI) / this.numSectors;

        for (let i = 0; i < this.numSectors; i++) {
            const angle = i * angleStep;
            const position = new Point3D(
                this.baseRadius * Math.cos(angle),
                this.baseRadius * Math.sin(angle),
                0
            );

            vertices.push(new TSplineVertex({
                position: position,
                weight: 1.0,
                knotIntervals: this.computeBaseKnotIntervals(i)
            }));
        }

        return vertices;
    }

    private createIntermediateRings(): TSplineVertex[][] {
        const rings: TSplineVertex[][] = [];
        const numRings = 2;  // Can be adjusted for refinement

        for (let ring = 1; ring <= numRings; ring++) {
            const ringRadius = this.baseRadius * (1 - ring/(numRings + 1));
            const ringHeight = this.height * ring/(numRings + 1);
            rings.push(this.createRingVertices(ringRadius, ringHeight));
        }

        return rings;
    }

    private createRingVertices(
        radius: number, 
        height: number
    ): TSplineVertex[] {
        const vertices: TSplineVertex[] = [];
        const angleStep = (2 * Math.PI) / this.numSectors;

        for (let i = 0; i < this.numSectors; i++) {
            const angle = i * angleStep;
            const position = new Point3D(
                radius * Math.cos(angle),
                radius * Math.sin(angle),
                height
            );

            vertices.push(new TSplineVertex({
                position: position,
                weight: 1.0,
                knotIntervals: this.computeRingKnotIntervals(i)
            }));
        }

        return vertices;
    }

    private createTSplineTopology(mesh: TSplineMesh): void {
        // Create special connectivity around apex
        this.createApexConnectivity(mesh);

        // Create regular connectivity for rest of cone
        this.createRegularConnectivity(mesh);
    }

    private createApexConnectivity(mesh: TSplineMesh): void {
        const apex = mesh.getApexVertex();
        const firstRing = mesh.getFirstRingVertices();

        // Create special T-mesh structure around apex
        for (let i = 0; i < firstRing.length; i++) {
            const next = firstRing[(i + 1) % firstRing.length];
            
            // Create face with degenerate edge at apex
            mesh.createFace([
                apex,
                firstRing[i],
                next
            ]);
        }
    }

    private computeKnotVectorsAroundApex(): KnotVectors {
        // Special knot vector configuration for apex
        return new KnotVectors({
            u: [0, 0, 0, 0, 1],  // Multiple knots at apex
            v: [0, 0, 0, 0, 1]
        });
    }
}

class ConeParameterization {
    constructor(
        private mesh: TSplineMesh,
        private apex: Point3D
    ) {}

    computeParameters(point: Point3D): UV {
        // Special handling near apex
        if (this.isNearApex(point)) {
            return this.computeApexParameters(point);
        }

        // Regular parameterization elsewhere
        return this.computeRegularParameters(point);
    }

    private isNearApex(point: Point3D): boolean {
        const distance = point.distanceTo(this.apex);
        return distance < this.apexThreshold;
    }

    private computeApexParameters(point: Point3D): UV {
        // Project point onto cone surface and compute parameters
        const projectedPoint = this.projectOntoConeSurface(point);
        
        // Use polar-like coordinates near apex
        const angle = this.computeAngle(projectedPoint);
        const height = this.computeHeight(projectedPoint);

        return new UV(
            angle / (2 * Math.PI),
            height / this.mesh.height
        );
    }
}

class ConeEvaluation {
    constructor(private mesh: TSplineMesh) {}

    evaluate(u: number, v: number): Point3D {
        if (this.isAtApex(u, v)) {
            return this.mesh.getApexPosition();
        }

        // Regular T-spline evaluation with special handling near apex
        return this.evaluateWithApexBlending(u, v);
    }

    private evaluateWithApexBlending(u: number, v: number): Point3D {
        const regularEval = this.regularTSplineEval(u, v);
        const apexInfluence = this.computeApexInfluence(u, v);

        if (apexInfluence > 0) {
            // Blend between regular evaluation and apex
            return this.blendWithApex(regularEval, apexInfluence);
        }

        return regularEval;
    }
}

// Example usage
const coneStructure = new ConeApexStructure(
    new Point3D(0, 0, 10),  // apex
    5,                      // base radius
    10,                     // height
    8                       // number of sectors
);

const coneMesh = coneStructure.createTSplineStructure();

// Create evaluator
const evaluator = new ConeEvaluation(coneMesh);

// Test evaluation
const point = evaluator.evaluate(0.5, 0.5);
console.log("Evaluated point:", point);

class SphericalTSplineSurface {
    constructor(
        private radius: number = 1.0,
        private refinementLevel: number = 2  // Controls mesh density
    ) {}

    createSphericalMesh(): TSplineMesh {
        // Start with cube topology and project to sphere
        const mesh = new TSplineMesh();
        
        // 1. Create initial cube-based control net
        this.createInitialControlNet(mesh);
        
        // 2. Project and adjust control points
        this.projectToSphere(mesh);
        
        // 3. Assign appropriate weights
        this.computeRationalWeights(mesh);
        
        // 4. Refine if needed
        if (this.refinementLevel > 1) {
            this.refineSphericalMesh(mesh);
        }

        return mesh;
    }

    private createInitialControlNet(mesh: TSplineMesh): void {
        // Create cube-based topology with 8 patches
        const controlPoints: TSplineVertex[] = [
            // Top square
            new TSplineVertex({ x: -1, y: -1, z: 1 }),
            new TSplineVertex({ x: 1, y: -1, z: 1 }),
            new TSplineVertex({ x: 1, y: 1, z: 1 }),
            new TSplineVertex({ x: -1, y: 1, z: 1 }),
            
            // Middle ring
            new TSplineVertex({ x: -1, y: -1, z: 0 }),
            new TSplineVertex({ x: 1, y: -1, z: 0 }),
            new TSplineVertex({ x: 1, y: 1, z: 0 }),
            new TSplineVertex({ x: -1, y: 1, z: 0 }),
            
            // Bottom square
            new TSplineVertex({ x: -1, y: -1, z: -1 }),
            new TSplineVertex({ x: 1, y: -1, z: -1 }),
            new TSplineVertex({ x: 1, y: 1, z: -1 }),
            new TSplineVertex({ x: -1, y: 1, z: -1 })
        ];

        // Add vertices to mesh
        controlPoints.forEach(v => mesh.addVertex(v));
        
        // Create connectivity
        this.createPatchConnectivity(mesh);
    }

    /*
    private projectToSphere(mesh: TSplineMesh): void {
        for (const vertex of mesh.vertices) {
            // Project to sphere and adjust position
            const pos = vertex.position;
            const length = Math.sqrt(
                pos.x * pos.x + 
                pos.y * pos.y + 
                pos.z * pos.z
            );
            
            // Project to sphere surface
            vertex.position = new Point3D(
                this.radius * pos.x / length,
                this.radius * pos.y / length,
                */

    private projectToSphere(mesh: TSplineMesh): void {
        // First normalize all vertices to unit sphere
        for (const vertex of mesh.vertices) {
            // Get the position vector
            const position = vertex.position;
            
            // Skip if at origin to avoid division by zero
            if (position.length() < 1e-10) {
                continue;
            }
    
            // Project onto unit sphere by normalizing
            vertex.position = position.normalize();
        }
    
        // Handle star points specially
        for (const starPoint of mesh.getStarPoints()) {
            this.handleStarPointProjection(starPoint, mesh);
        }
    
        // Adjust edges to be great circle arcs
        this.projectEdgesToGreatCircles(mesh);
    
        // Verify and adjust for spherical consistency
        this.enforceSphericalConstraints(mesh);
    }
    
    private handleStarPointProjection(
        starPoint: StarPoint, 
        mesh: TSplineMesh
    ): void {
        // Get all edges connected to star point
        const connectedEdges = mesh.getEdgesFromVertex(starPoint);
        
        // Compute the average direction for the star point
        const avgDirection = this.computeAverageDirection(connectedEdges);
        
        // Project star point to sphere along average direction
        starPoint.position = avgDirection.normalize();
        
        // Adjust neighboring vertices to maintain mesh structure
        this.adjustStarPointNeighbors(starPoint, connectedEdges);
    }
    
    private computeAverageDirection(edges: Edge[]): Vector3D {
        let avgDirection = new Vector3D(0, 0, 0);
        
        for (const edge of edges) {
            const direction = edge.end.position.subtract(edge.start.position);
            avgDirection = avgDirection.add(direction);
        }
        
        return avgDirection.length() > 1e-10 ? 
            avgDirection.normalize() : 
            new Vector3D(0, 0, 1);
    }
    
    private projectEdgesToGreatCircles(mesh: TSplineMesh): void {
        for (const edge of mesh.edges) {
            // Get endpoints on sphere
            const start = edge.start.position;
            const end = edge.end.position;
            
            // Create great circle interpolation
            const interpolator = new GreatCircleInterpolator(start, end);
            
            // Adjust interior edge points to lie on great circle
            if (edge.hasInteriorPoints()) {
                for (const point of edge.interiorPoints) {
                    const t = this.computeParametricPosition(
                        point.position, 
                        start, 
                        end
                    );
                    point.position = interpolator.interpolate(t);
                }
            }
        }
    }
    
    private computeParametricPosition(
        point: Point3D,
        start: Point3D,
        end: Point3D
    ): number {
        // Project point onto line between start and end
        const v = end.subtract(start);
        const w = point.subtract(start);
        
        // Compute parametric position
        return w.dot(v) / v.dot(v);
    }
    
    private enforceSphericalConstraints(mesh: TSplineMesh): void {
        // Ensure all points lie exactly on unit sphere
        for (const vertex of mesh.vertices) {
            if (Math.abs(vertex.position.length() - 1.0) > 1e-10) {
                vertex.position = vertex.position.normalize();
            }
        }
        
        // Check and adjust for spherical patch consistency
        for (const patch of mesh.patches) {
            this.adjustSphericalPatch(patch);
        }
    }
    
    private adjustSphericalPatch(patch: TSplinePatch): void {
        // Ensure patch boundaries follow great circles
        for (const boundary of patch.boundaries) {
            this.adjustBoundaryToGreatCircle(boundary);
        }
        
        // Adjust interior control points to maintain shape
        this.adjustPatchInteriorPoints(patch);
    }
    
    private adjustBoundaryToGreatCircle(boundary: PatchBoundary): void {
        const start = boundary.startPoint.position;
        const end = boundary.endPoint.position;
        
        // Create great circle interpolator
        const greatCircle = new GreatCircleInterpolator(start, end);
        
        // Adjust boundary points to lie on great circle
        for (const point of boundary.points) {
            const t = this.computeParametricPosition(
                point.position,
                start,
                end
            );
            point.position = greatCircle.interpolate(t);
        }
    }
    
    private adjustPatchInteriorPoints(patch: TSplinePatch): void {
        // Use spherical barycentric coordinates to adjust interior points
        for (const point of patch.interiorPoints) {
            const coords = this.computeSphericalBarycentricCoords(
                point,
                patch.getBoundaryPoints()
            );
            point.position = this.interpolateSphericalPosition(
                coords,
                patch.getBoundaryPoints()
            );
        }
    }
    
    class GreatCircleInterpolator {
        constructor(
            private start: Point3D,
            private end: Point3D
        ) {
            // Normalize points to ensure they're on unit sphere
            this.start = start.normalize();
            this.end = end.normalize();
        }
    
        interpolate(t: number): Point3D {
            // Use spherical linear interpolation (SLERP)
            const omega = Math.acos(
                Math.max(-1, Math.min(1, this.start.dot(this.end)))
            );
            
            if (Math.abs(omega) < 1e-10) {
                return this.start;
            }
    
            const sinOmega = Math.sin(omega);
            const scale1 = Math.sin((1 - t) * omega) / sinOmega;
            const scale2 = Math.sin(t * omega) / sinOmega;
    
            return this.start.scale(scale1)
                .add(this.end.scale(scale2))
                .normalize();
        }
    }

    // Project vertices to unit sphere by normalizing
    vertex.position = position.normalize();
    
    private handleStarPointProjection(starPoint: StarPoint): void {
        // Special handling for star points to maintain mesh structure
        // Uses average direction of connected edges
    }
    
    private projectEdgesToGreatCircles(mesh: TSplineMesh): void {
        // Ensure edges follow great circles on sphere
        // Uses spherical linear interpolation (SLERP)
    }
    
    private enforceSphericalConstraints(mesh: TSplineMesh): void {
        // Ensure all points lie exactly on sphere
        // Maintain mesh structure and shape
    }

    class SphericalInterpolation {
        constructor(
            private start: Point3D,
            private end: Point3D,
            private radius: number = 1.0
        ) {
            // Normalize points to sphere surface
            this.start = start.normalize().scale(radius);
            this.end = end.normalize().scale(radius);
        }
    
        slerp(t: number): Point3D {
            // Compute angle between vectors
            const cosOmega = Math.max(-1, Math.min(1, 
                this.start.dot(this.end) / 
                (this.start.length() * this.end.length())
            ));
            const omega = Math.acos(cosOmega);
    
            // Handle special cases
            if (Math.abs(omega) < 1e-10) {
                // Points are very close - linear interpolation
                return this.start.scale(1 - t).add(this.end.scale(t));
            }
            if (Math.abs(Math.abs(cosOmega) - 1) < 1e-10) {
                // Points are antipodal - choose arbitrary great circle path
                return this.handleAntipodal(t);
            }
    
            // Compute SLERP
            const sinOmega = Math.sin(omega);
            const scale1 = Math.sin((1 - t) * omega) / sinOmega;
            const scale2 = Math.sin(t * omega) / sinOmega;
    
            return this.start.scale(scale1)
                .add(this.end.scale(scale2));
        }
    
        // Extended SLERP with acceleration control
        slerpWithAcceleration(
            t: number, 
            acceleration: number = 1.0
        ): Point3D {
            // Modify t with acceleration parameter
            const modifiedT = this.accelerationFunction(t, acceleration);
            return this.slerp(modifiedT);
        }
    
        // Quaternion SLERP for better numerical stability
        quaternionSlerp(t: number): Point3D {
            const startQuat = this.vectorToQuaternion(this.start);
            const endQuat = this.vectorToQuaternion(this.end);
            
            const interpolatedQuat = this.quaternionInterpolation(
                startQuat,
                endQuat,
                t
            );
            
            return this.quaternionToVector(interpolatedQuat);
        }
    
        // Double SLERP for smooth interpolation between multiple points
        doubleSlerp(
            middle: Point3D,
            t: number
        ): Point3D {
            if (t <= 0.5) {
                // First half: interpolate between start and middle
                return this.slerp(2 * t);
            } else {
                // Second half: interpolate between middle and end
                const secondHalf = new SphericalInterpolation(
                    middle,
                    this.end
                );
                return secondHalf.slerp(2 * (t - 0.5));
            }
        }
    
        private handleAntipodal(t: number): Point3D {
            // Handle antipodal points by choosing arbitrary great circle
            const arbitrary = this.findPerpendicularVector(this.start);
            const rotationAxis = this.start.cross(arbitrary).normalize();
            
            // Rotate around axis by 180 degrees * t
            return this.rotateAroundAxis(
                this.start,
                rotationAxis,
                Math.PI * t
            );
        }
    
        private findPerpendicularVector(v: Point3D): Point3D {
            // Find a vector perpendicular to v
            if (Math.abs(v.x) < Math.abs(v.y)) {
    
    







