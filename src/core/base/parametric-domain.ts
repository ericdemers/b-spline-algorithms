

interface ParametricDomain {
    dimension: number;          // Topological dimension
    
    // Core topology properties
    isPeriodic: boolean[];     // Per dimension periodicity
    hasBoundary: boolean[];    // Per dimension boundary existence
    
    // Fundamental operations
    isValid(parameter: number[]): boolean;
    getDistance(p1: number[], p2: number[]): number;
    
    // Topology-related operations
    getEquivalentPoints(parameter: number[]): number[][];
    mapToFundamentalDomain(parameter: number[]): number[];

    // Singularity handling
    getSingularities(): Singularity[];
    isRegularPoint(parameter: number[]): boolean;
    getNearestSingularity(parameter: number[]): Singularity | null;
    
    // Local structure analysis
    getLocalStructure(parameter: number[]): LocalStructure;

    // Add decompose method
    decompose(): Cell[];
}

// Example: Standard parameter space [0,1]^n
class StandardParametricDomain implements ParametricDomain {
    constructor(
        public dimension: number,
        public isPeriodic: boolean[] = new Array(dimension).fill(false),
        public hasBoundary: boolean[] = new Array(dimension).fill(true)
    ) {}

    isValid(parameter: number[]): boolean {
        return parameter.length === this.dimension &&
            parameter.every((p, i) => 
                this.isPeriodic[i] ? true : (p >= 0 && p <= 1));
    }

    getDistance(p1: number[], p2: number[]): number {
        return Math.sqrt(p1.reduce((sum, val, i) => {
            let diff = p2[i] - val;
            if (this.isPeriodic[i]) {
                // Handle periodic case - take shortest path
                diff = Math.min(Math.abs(diff), 1 - Math.abs(diff));
            }
            return sum + diff * diff;
        }, 0));
    }

    decompose(): Cell[] {
        // For standard domain, just return a single cell
        return [{
            bounds: {
                min: [0, 0],
                max: [1, 1]
            },
            type: CellType.REGULAR
        }];
    }

}

// Example: Periodic parameter space (e.g., for periodic B-splines)
class PeriodicParametricSpace implements ParametricSpace {
    constructor(
        public dimension: number,
        public periods: number[] // Length per dimension
    ) {
        this.isPeriodic = new Array(dimension).fill(true);
        this.hasBoundary = new Array(dimension).fill(false);
    }

    mapToFundamentalDomain(parameter: number[]): number[] {
        return parameter.map((p, i) => {
            const period = this.periods[i];
            return ((p % period) + period) % period;
        });
    }

    getEquivalentPoints(parameter: number[]): number[][] {
        // Returns equivalent points under periodicity
        // Useful for basis function evaluation
    }
}




// For topology-aware parameter spaces (e.g., T-splines with extraordinary points)
interface TopologicalParametricDomain extends ParametricDomain {
    // Local structure information
    getLocalChart(parameter: number[]): LocalChart;
    getTransitionMaps(from: LocalChart, to: LocalChart): TransitionMap;
}

interface LocalChart {
    center: number[];
    neighborhood: Map<Direction, number[]>;
    isRegular: boolean;
    valence?: number;  // For extraordinary points
}

interface TransitionMap {
    forward: (p: number[]) => number[];
    inverse: (p: number[]) => number[];
    isSmooth: boolean;
}


export {
    ParametricDomain as ParametricDomain,
    TopologicalParametricDomain as TopologicalParametricDomainSpace,
    LocalChart,
    TransitionMap,
    StandardParametric
}

// Types of singularities we might encounter
enum SingularityType {
    EXTRAORDINARY_POINT,    // Vertex with valence != 4
    T_JUNCTION,            // T-spline T-junction
    POLAR,                 // Polar singularity (e.g., sphere pole)
    DEGENERATE            // Collapsed edge/point (like cone apex)
}

interface Singularity {
    type: SingularityType;
    position: number[];           // Location in parameter space
    valence?: number;            // For extraordinary points
    tangentCone?: TangentCone;   // Local geometric structure
    chart: LocalChart;           // Local parameterization
}

interface LocalChart {
    // Maps from local to global coordinates
    toLocal(parameter: number[]): number[];
    toGlobal(localParam: number[]): number[];
    
    // Describes valid region around singularity
    isInValidRegion(localParam: number[]): boolean;
    
    // For basis function evaluation
    getLocalBasisFunctions(): BasisFunction[];
}


// Implementation example
class SingularParametricDomain implements ParametricDomain {
    private singularities: Map<string, Singularity>;

    constructor(
        public dimension: number,
        public bounds: { min: number[], max: number[] },
        singularPoints: Singularity[]
    ) {
        this.singularities = new Map(
            singularPoints.map(s => [this.hashPosition(s.position), s])
        );
    }

    isRegularPoint(parameter: number[]): boolean {
        // Check if point is away from all singularities
        return !this.getNearestSingularity(parameter);
    }

    getLocalStructure(parameter: number[]): LocalStructure {
        const singularity = this.getNearestSingularity(parameter);
        
        if (!singularity) {
            return this.getRegularStructure(parameter);
        }

        switch (singularity.type) {
            case SingularityType.EXTRAORDINARY_POINT:
                return this.getExtraordinaryStructure(parameter, singularity);
            case SingularityType.T_JUNCTION:
                return this.getTJunctionStructure(parameter, singularity);
            case SingularityType.POLAR:
                return this.getPolarStructure(parameter, singularity);
            default:
                throw new Error('Unknown singularity type');
        }
    }

    private getExtraordinaryStructure(
        parameter: number[], 
        singularity: Singularity
    ): LocalStructure {
        const localParam = singularity.chart.toLocal(parameter);
        
        return {
            type: 'extraordinary',
            valence: singularity.valence!,
            sectors: this.computeSectors(localParam, singularity),
            basisFunctions: this.computeExtraordinaryBasis(localParam, singularity),
            jacobian: this.computeLocalJacobian(localParam, singularity)
        };
    }

    private computeSectors(
        localParam: number[], 
        singularity: Singularity
    ): Sector[] {
        // Divide region around extraordinary point into sectors
        const valence = singularity.valence!;
        const sectorAngle = 2 * Math.PI / valence;
        
        return Array.from({ length: valence }, (_, i) => ({
            index: i,
            angle: i * sectorAngle,
            basis: this.computeSectorBasis(localParam, i, valence)
        }));
    }

    private computeLocalJacobian(
        localParam: number[], 
        singularity: Singularity
    ): Matrix {
        // Compute transformation derivatives for geometric mapping
        // This is crucial for maintaining surface continuity
        // ...
    }
}

// Usage example
const domainWithSingularities = new SingularParametricDomain(
    2, // 2D domain
    { min: [0, 0], max: [1, 1] },
    [
        {
            type: SingularityType.EXTRAORDINARY_POINT,
            position: [0.5, 0.5],
            valence: 3,
            chart: new ExtraordinaryChart(3),
            tangentCone: new TangentCone(3)
        },
        {
            type: SingularityType.T_JUNCTION,
            position: [0.25, 0.5],
            chart: new TJunctionChart()
        }
    ]
);

interface TSplineDomain extends ParametricDomain {
    // T-mesh structure
    tMesh: TMesh;
    
    // Core mapping functionality
    mapToParameterSpace(index: Index2D): number[];
    mapFromParameterSpace(parameter: number[]): Index2D;
    
    // T-junction handling
    getTJunctions(): TJunction[];
    getLocalKnotVectors(parameter: number[]): LocalKnotStructure;
}

interface TMesh {
    // Topological structure
    vertices: TVertex[];
    edges: TEdge[];
    faces: TFace[];
    
    // T-mesh specific operations
    findCell(parameter: number[]): TCell;
    getLocalRefinement(cell: TCell): RefinementInfo;
    
    // Extraordinary point handling
    getExtraordinaryVertices(): TVertex[];
    getSectorDecomposition(vertex: TVertex): Sector[];
}

class TSplineParametricDomain implements ParametricDomain {
    constructor(
        private tMesh: TMesh,
        public dimension: number,
        public topology: TopologyType
    ) {}

    decompose(): Cell[] {
        const cells: Cell[] = [];
        const vertices = this.tMesh.getVertices();
        
        // Get all extraordinary points and T-junctions
        const extraordinaryPoints = vertices.filter(v => v.valence !== 4);
        const tJunctions = vertices.filter(v => v.isTJunction());
        
        // Create cells based on T-mesh faces
        for (const face of this.tMesh.getFaces()) {
            const cell = this.createCellFromFace(face);
            
            // Check if cell contains special points
            const containedExtraordinary = extraordinaryPoints.filter(
                p => this.cellContainsPoint(cell, p.position)
            );
            
            const containedTJunctions = tJunctions.filter(
                p => this.cellContainsPoint(cell, p.position)
            );
            
            if (containedExtraordinary.length > 0) {
                cell.type = CellType.EXTRAORDINARY;
                cell.features = {
                    extraordinaryPoints: containedExtraordinary
                };
            } else if (containedTJunctions.length > 0) {
                cell.type = CellType.T_JUNCTION;
                cell.features = {
                    tJunctions: containedTJunctions
                };
            } else {
                cell.type = CellType.REGULAR;
            }
            
            cells.push(cell);
        }
        
        return cells;
    }

    private createCellFromFace(face: Face): Cell {
        const vertices = face.getVertices();
        const bounds = this.computeBounds(vertices);
        
        return {
            bounds: bounds,
            type: CellType.REGULAR, // Initial type, may be updated
            features: {}
        };
    }

    private computeBounds(vertices: Vertex[]): { min: number[], max: number[] } {
        const coords = vertices.map(v => v.position);
        
        const min = coords.reduce((acc, curr) => [
            Math.min(acc[0], curr[0]),
            Math.min(acc[1], curr[1])
        ], [Infinity, Infinity]);
        
        const max = coords.reduce((acc, curr) => [
            Math.max(acc[0], curr[0]),
            Math.max(acc[1], curr[1])
        ], [-Infinity, -Infinity]);
        
        return { min, max };
    }

    private cellContainsPoint(cell: Cell, point: number[]): boolean {
        return point[0] >= cell.bounds.min[0] && 
               point[0] <= cell.bounds.max[0] &&
               point[1] >= cell.bounds.min[1] && 
               point[1] <= cell.bounds.max[1];
    }
}

class TSplineParameterSpace implements TSplineDomain {
    constructor(
        private tMesh: TMesh,
        private degree: [number, number] // [p,q] degrees
    ) {}

    mapToParameterSpace(index: Index2D): number[] {
        const cell = this.tMesh.findCell(index);
        return cell.computeParameterCoordinates(index);
    }

    getLocalKnotVectors(parameter: number[]): LocalKnotStructure {
        const cell = this.tMesh.findCell(parameter);
        
        if (cell.hasExtraordinaryPoint()) {
            return this.getExtraordinaryKnotStructure(cell);
        }

        if (cell.hasTJunction()) {
            return this.getTJunctionKnotStructure(cell);
        }

        return this.getRegularKnotStructure(cell);
    }

    private getExtraordinaryKnotStructure(
        cell: TCell
    ): ExtraordinaryKnotStructure {
        const vertex = cell.getExtraordinaryVertex();
        const sectors = this.tMesh.getSectorDecomposition(vertex);
        
        return new ExtraordinaryKnotStructure(
            vertex,
            sectors,
            this.degree
        );
    }
}

// Handle local parameter space near T-junctions
class TJunctionParameterSpace {
    constructor(
        private junction: TJunction,
        private transitionFunction: BlendingFunction
    ) {}

    evaluateParameter(u: number, v: number): number {
        const distance = this.junction.getDistance([u, v]);
        const blend = this.transitionFunction.evaluate(distance);
        
        return this.blendKnotIntervals(
            this.junction.getMainInterval(),
            this.junction.getSubInterval(),
            blend
        );
    }
}

// Handle extraordinary points in T-splines
class ExtraordinaryParameterSpace {
    constructor(
        private vertex: TVertex,
        private valence: number,
        private sectors: Sector[]
    ) {}

    mapToLocal(parameter: number[]): LocalParameters {
        const sector = this.findSector(parameter);
        return sector.mapToLocal(parameter);
    }

    evaluateBasis(
        parameter: number[], 
        basisIndex: number
    ): number {
        const localParams = this.mapToLocal(parameter);
        
        if (this.isAtExtraordinaryPoint(parameter)) {
            return this.evaluateExtraordinaryBasis(
                localParams,
                basisIndex
            );
        }

        // Blend between extraordinary and regular basis
        const distance = this.getDistanceToExtraordinary(parameter);
        const blend = this.getBlendingFactor(distance);
        
        return (1 - blend) * this.evaluateRegularBasis(parameter, basisIndex) +
               blend * this.evaluateExtraordinaryBasis(localParams, basisIndex);
    }
}

// Example usage
class TSplineSurface {
    constructor(
        private domain: TSplineDomain,
        private controlPoints: ControlPoint[],
        private weights: number[]
    ) {}

    evaluate(u: number, v: number): Point3D {
        const parameter = [u, v];
        const cell = this.domain.tMesh.findCell(parameter);
        
        if (cell.isRegular()) {
            return this.evaluateRegular(parameter);
        }
        
        if (cell.hasExtraordinaryPoint()) {
            return this.evaluateExtraordinary(parameter);
        }
        
        if (cell.hasTJunction()) {
            return this.evaluateTJunction(parameter);
        }
    }

    private evaluateExtraordinary(parameter: number[]): Point3D {
        const extraSpace = new ExtraordinaryParameterSpace(
            cell.getExtraordinaryVertex(),
            cell.getValence(),
            this.domain.tMesh.getSectorDecomposition(cell)
        );

        return extraSpace.evaluate(
            parameter,
            this.controlPoints,
            this.weights
        );
    }
}

interface TVertex {
    id: number;
    position: number[];    // [u, v] in parameter space
    valence: number;       // Number of incident edges
    
    // Topology relationships
    edges: TEdge[];       // Connected edges
    faces: TFace[];       // Adjacent faces
    
    // Helper methods
    isTJunction(): boolean {
        // T-junction has 3 incident edges
        return this.edges.length === 3;
    }

    isExtraordinary(): boolean {
        // Extraordinary point has valence != 4
        return this.valence !== 4;
    }
}

interface TEdge {
    id: number;
    start: TVertex;       // Start vertex
    end: TVertex;         // End vertex
    faces: TFace[];       // Adjacent faces (1 or 2)
    
    // Helper methods
    isOnBoundary(): boolean {
        // Boundary edge has only one adjacent face
        return this.faces.length === 1;
    }

    getOtherVertex(vertex: TVertex): TVertex {
        return vertex === this.start ? this.end : this.start;
    }
}

interface TFace {
    id: number;
    vertices: TVertex[];  // Vertices in CCW order
    edges: TEdge[];       // Bounding edges
    neighbors: TFace[];   // Adjacent faces
    
    // Helper methods
    getVertexLoop(): TVertex[] {
        return this.vertices;
    }

    getEdgeLoop(): TEdge[] {
        return this.edges;
    }

    isNeighbor(face: TFace): boolean {
        return this.neighbors.includes(face);
    }
}

// Example implementation of mesh traversal
class TMeshTraversal {
    // Get all faces connected to a vertex
    getFacesAtVertex(vertex: TVertex): TFace[] {
        return vertex.faces;
    }

    // Get vertices connected to a vertex through edges
    getVertexNeighbors(vertex: TVertex): TVertex[] {
        return vertex.edges.map(edge => edge.getOtherVertex(vertex));
    }

    // Get faces sharing an edge
    getFacesAtEdge(edge: TEdge): TFace[] {
        return edge.faces;
    }

    // Get edges around a face
    getEdgesAroundFace(face: TFace): TEdge[] {
        return face.edges;
    }

    // Get vertex ring (1-neighborhood) around a vertex
    getVertexRing(vertex: TVertex): TVertex[] {
        const ring: TVertex[] = [];
        const visited = new Set<number>();
        
        // Start with first edge
        let currentEdge = vertex.edges[0];
        let startVertex = vertex;
        
        do {
            // Get other vertex of current edge
            const nextVertex = currentEdge.getOtherVertex(startVertex);
            
            if (!visited.has(nextVertex.id)) {
                ring.push(nextVertex);
                visited.add(nextVertex.id);
            }
            
            // Find next edge in CCW order
            currentEdge = this.getNextEdgeInRing(currentEdge, startVertex);
            startVertex = nextVertex;
            
        } while (currentEdge !== vertex.edges[0]);
        
        return ring;
    }

    // Get face ring around a vertex
    getFaceRing(vertex: TVertex): TFace[] {
        return vertex.faces;
    }
}

// Example usage
class TMeshOperations {
    // Split an edge
    splitEdge(edge: TEdge, parameter: number): TVertex {
        const newPosition = [
            edge.start.position[0] * (1 - parameter) + edge.end.position[0] * parameter,
            edge.start.position[1] * (1 - parameter) + edge.end.position[1] * parameter
        ];
        
        const newVertex = new TVertex(newPosition);
        const newEdges = [
            new TEdge(edge.start, newVertex),
            new TEdge(newVertex, edge.end)
        ];
        
        // Update face connectivity
        for (const face of edge.faces) {
            this.updateFaceAfterEdgeSplit(face, edge, newVertex, newEdges);
        }
        
        return newVertex;
    }

    // Insert a T-junction
    insertTJunction(edge: TEdge, parameter: number): TVertex {
        const tJunction = this.splitEdge(edge, parameter);
        tJunction.valence = 3;
        return tJunction;
    }
}


interface BSplineMesh {
    vertices: BSplineVertex[];
    edges: BSplineEdge[];
    faces: BSplineFace[];
    
    // Grid structure for regular topology
    numU: number;  // Number of vertices in U direction
    numV: number;  // Number of vertices in V direction
}

class BSplineVertex {
    constructor(
        public id: number,
        public position: number[],    // [u, v] in parameter space
        public indices: [number, number]  // [i, j] grid indices
    ) {}

    // Always 4 for regular B-spline (except boundary)
    getValence(): number {
        return this.edges.length;
    }

    // Topological relationships
    edges: BSplineEdge[] = [];
    faces: BSplineFace[] = [];
}

class BSplineEdge {
    constructor(
        public id: number,
        public start: BSplineVertex,
        public end: BSplineVertex
    ) {}

    faces: BSplineFace[] = [];  // 1 or 2 faces

    isUDirection(): boolean {
        return this.start.indices[1] === this.end.indices[1];
    }

    isVDirection(): boolean {
        return this.start.indices[0] === this.end.indices[0];
    }
}

class BSplineFace {
    constructor(
        public id: number,
        public indices: [number, number]  // Bottom-left corner indices
    ) {}

    vertices: BSplineVertex[] = [];  // Always 4 vertices
    edges: BSplineEdge[] = [];      // Always 4 edges
}

class RegularBSplineMesh implements BSplineMesh {
    vertices: BSplineVertex[] = [];
    edges: BSplineEdge[] = [];
    faces: BSplineFace[] = [];
    
    constructor(
        public numU: number,
        public numV: number,
        private uKnots: number[],
        private vKnots: number[]
    ) {
        this.constructMesh();
    }

    private constructMesh(): void {
        // 1. Create vertices
        this.createVertices();
        
        // 2. Create edges
        this.createEdges();
        
        // 3. Create faces
        this.createFaces();
        
        // 4. Set up connectivity
        this.establishConnectivity();
    }

    private createVertices(): void {
        let id = 0;
        
        // Create vertices in grid pattern
        for (let j = 0; j < this.numV; j++) {
            for (let i = 0; i < this.numU; i++) {
                const u = this.uKnots[i + 3]; // Assuming degree 3
                const v = this.vKnots[j + 3];
                
                this.vertices.push(new BSplineVertex(
                    id++,
                    [u, v],
                    [i, j]
                ));
            }
        }
    }

    private createEdges(): void {
        let id = 0;
        
        // Create horizontal edges
        for (let j = 0; j < this.numV; j++) {
            for (let i = 0; i < this.numU - 1; i++) {
                const start = this.getVertex(i, j);
                const end = this.getVertex(i + 1, j);
                
                const edge = new BSplineEdge(id++, start, end);
                this.edges.push(edge);
                
                start.edges.push(edge);
                end.edges.push(edge);
            }
        }
        
        // Create vertical edges
        for (let j = 0; j < this.numV - 1; j++) {
            for (let i = 0; i < this.numU; i++) {
                const start = this.getVertex(i, j);
                const end = this.getVertex(i, j + 1);
                
                const edge = new BSplineEdge(id++, start, end);
                this.edges.push(edge);
                
                start.edges.push(edge);
                end.edges.push(edge);
            }
        }
    }

    private createFaces(): void {
        let id = 0;
        
        // Create quad faces
        for (let j = 0; j < this.numV - 1; j++) {
            for (let i = 0; i < this.numU - 1; i++) {
                const face = new BSplineFace(id++, [i, j]);
                
                // Add vertices (counter-clockwise)
                face.vertices = [
                    this.getVertex(i, j),
                    this.getVertex(i + 1, j),
                    this.getVertex(i + 1, j + 1),
                    this.getVertex(i, j + 1)
                ];
                
                // Add edges
                face.edges = [
                    this.getHorizontalEdge(i, j),
                    this.getVerticalEdge(i + 1, j),
                    this.getHorizontalEdge(i, j + 1),
                    this.getVerticalEdge(i, j)
                ];
                
                this.faces.push(face);
                
                // Update vertex-face connectivity
                face.vertices.forEach(v => v.faces.push(face));
                
                // Update edge-face connectivity
                face.edges.forEach(e => e.faces.push(face));
            }
        }
    }

    private getVertex(i: number, j: number): BSplineVertex {
        return this.vertices[j * this.numU + i];
    }

    private getHorizontalEdge(i: number, j: number): BSplineEdge {
        return this.edges[j * (this.numU - 1) + i];
    }

    private getVerticalEdge(i: number, j: number): BSplineEdge {
        const horizontalEdgeCount = this.numU * (this.numV - 1);
        return this.edges[horizontalEdgeCount + j * this.numU + i];
    }

    // Utility methods for mesh traversal
    getParametricCoordinates(i: number, j: number): number[] {
        return this.getVertex(i, j).position;
    }

    getFaceParametricBounds(face: BSplineFace): {
        min: number[];
        max: number[];
    } {
        const [i, j] = face.indices;
        return {
            min: this.getParametricCoordinates(i, j),
            max: this.getParametricCoordinates(i + 1, j + 1)
        };
    }
}

// Example usage
const mesh = new RegularBSplineMesh(
    5,  // numU
    4,  // numV
    [0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1],  // uKnots
    [0, 0, 0, 0, 0.33, 0.67, 1, 1, 1, 1]        // vKnots
);


class BSplineSurfaceRenderer {
    private gl: WebGLRenderingContext;
    private program: WebGLProgram;
    private buffers: {
        position: WebGLBuffer;
        normal: WebGLBuffer;
        uv: WebGLBuffer;
        index: WebGLBuffer;
    };

    constructor(canvas: HTMLCanvasElement) {
        this.gl = canvas.getContext('webgl')!;
        this.program = this.createShaderProgram();
        this.buffers = this.createBuffers();
    }

    updateMesh(mesh: BSplineMesh): void {
        const gl = this.gl;
        
        // Update vertex positions
        gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.position);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(mesh.vertices), gl.STATIC_DRAW);
        
        // Update normals
        gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.normal);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(mesh.normals), gl.STATIC_DRAW);
        
        // Update indices
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.buffers.index);
        gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(mesh.indices), gl.STATIC_DRAW);
    }

    render(): void {
        const gl = this.gl;
        gl.useProgram(this.program);
        
        // Set up attributes and uniforms
        this.setupAttributes();
        this.setupUniforms();
        
        // Draw triangles
        gl.drawElements(gl.TRIANGLES, this.indexCount, gl.UNSIGNED_SHORT, 0);
    }
}

interface BSplineCurveMesh {
    vertices: CurveVertex[];
    edges: CurveEdge[];
    isClosed: boolean;
}

class CurveVertex {
    constructor(
        public id: number,
        public position: number,    // u in parameter space
        public index: number        // i index in control point sequence
    ) {}

    // For open curve: end points have 1 edge, interior have 2
    // For closed curve: all vertices have 2 edges
    edges: CurveEdge[] = [];

    isEndPoint(): boolean {
        return this.edges.length === 1;
    }
}

class CurveEdge {
    constructor(
        public id: number,
        public start: CurveVertex,
        public end: CurveVertex
    ) {}

    getParametricInterval(): [number, number] {
        return [this.start.position, this.end.position];
    }
}

class BSplineCurveMeshBuilder {
    constructor(
        private numPoints: number,
        private knots: number[],
        private isClosed: boolean
    ) {}

    build(): BSplineCurveMesh {
        const mesh: BSplineCurveMesh = {
            vertices: [],
            edges: [],
            isClosed: this.isClosed
        };

        this.createVertices(mesh);
        this.createEdges(mesh);

        return mesh;
    }

    private createVertices(mesh: BSplineCurveMesh): void {
        // For closed curve, we need to handle periodic knot vector
        const n = this.isClosed ? this.numPoints + 1 : this.numPoints;
        
        for (let i = 0; i < n; i++) {
            const vertex = new CurveVertex(
                i,
                this.knots[i + 3], // Assuming cubic (degree 3)
                i % this.numPoints  // Wrap around for closed curve
            );
            mesh.vertices.push(vertex);
        }
    }

    private createEdges(mesh: BSplineCurveMesh): void {
        const n = mesh.vertices.length;
        
        for (let i = 0; i < (this.isClosed ? n : n - 1); i++) {
            const start = mesh.vertices[i];
            const end = mesh.vertices[(i + 1) % n];
            
            const edge = new CurveEdge(i, start, end);
            mesh.edges.push(edge);
            
            start.edges.push(edge);
            end.edges.push(edge);
        }
    }
}

class RegularBSplineCurve {
    private mesh: BSplineCurveMesh;

    constructor(
        numPoints: number,
        knots: number[],
        isClosed: boolean = false
    ) {
        const builder = new BSplineCurveMeshBuilder(
            numPoints,
            knots,
            isClosed
        );
        this.mesh = builder.build();
    }

    // Mesh traversal methods
    getVertex(index: number): CurveVertex {
        return this.mesh.vertices[index];
    }

    getEdge(index: number): CurveEdge {
        return this.mesh.edges[index];
    }

    // Topology queries
    findEdgeContainingParameter(u: number): CurveEdge | null {
        for (const edge of this.mesh.edges) {
            const [start, end] = edge.getParametricInterval();
            if (u >= start && u <= end) {
                return edge;
            }
        }
        return null;
    }

    // Iterator for traversing the curve
    *traverseVertices(): Generator<CurveVertex> {
        for (const vertex of this.mesh.vertices) {
            yield vertex;
        }
    }

    *traverseEdges(): Generator<CurveEdge> {
        for (const edge of this.mesh.edges) {
            yield edge;
        }
    }
}

// Example usage for open curve
const openCurve = new RegularBSplineCurve(
    5,  // numPoints
    [0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1], // knots
    false // open
);

// Example usage for closed curve
const closedCurve = new RegularBSplineCurve(
    5,  // numPoints
    [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], // periodic knots
    true  // closed
);

// Example of mesh operations
class CurveMeshOperations {
    // Insert a new vertex by splitting an edge
    insertVertex(
        curve: RegularBSplineCurve,
        parameter: number
    ): CurveVertex | null {
        const edge = curve.findEdgeContainingParameter(parameter);
        if (!edge) return null;

        const t = this.computeRelativeParameter(
            parameter,
            edge.getParametricInterval()
        );

        return this.splitEdge(edge, t);
    }

    private computeRelativeParameter(
        u: number,
        [start, end]: [number, number]
    ): number {
        return (u - start) / (end - start);
    }

    private splitEdge(edge: CurveEdge, t: number): CurveVertex {
        const newPosition = edge.start.position * (1 - t) + 
                          edge.end.position * t;
        
        const newVertex = new CurveVertex(
            -1, // temporary id
            newPosition,
            -1  // temporary index
        );

        // Create new edges
        const newEdge1 = new CurveEdge(-1, edge.start, newVertex);
        const newEdge2 = new CurveEdge(-1, newVertex, edge.end);

        // Update connectivity
        newVertex.edges = [newEdge1, newEdge2];
        
        // Update original vertices' edge lists
        this.replaceEdge(edge.start.edges, edge, newEdge1);
        this.replaceEdge(edge.end.edges, edge, newEdge2);

        return newVertex;
    }

    private replaceEdge(
        edges: CurveEdge[],
        oldEdge: CurveEdge,
        newEdge: CurveEdge
    ): void {
        const index = edges.indexOf(oldEdge);
        if (index !== -1) {
            edges[index] = newEdge;
        }
    }
}

interface BSplineCurveMesh {
    vertices: CurveVertex[];
    edges: CurveEdge[];
    isClosed: boolean;
    // Add reference to peripheral vertices
    peripheralVertices: CurveVertex[];
}

class CurveVertex {
    constructor(
        public id: number,
        public position: number,    // u in parameter space
        public index: number,       // i index in control point sequence
        public isPeripheral: boolean = false  // Flag for peripheral vertices
    ) {}

    edges: CurveEdge[] = [];

    isEndPoint(): boolean {
        return !this.isPeripheral && this.edges.length === 1;
    }
}

class BSplineCurveMeshBuilder {
    private degree: number;

    constructor(
        private numPoints: number,
        private knots: number[],
        private isClosed: boolean,
        degree: number = 3  // Default to cubic
    ) {
        this.degree = degree;
    }

    build(): BSplineCurveMesh {
        const mesh: BSplineCurveMesh = {
            vertices: [],
            edges: [],
            isClosed: this.isClosed,
            peripheralVertices: []
        };

        if (this.isClosed) {
            this.createPeriodicMesh(mesh);
        } else {
            this.createOpenMeshWithPeripheral(mesh);
        }

        return mesh;
    }

    private createOpenMeshWithPeripheral(mesh: BSplineCurveMesh): void {
        // Create peripheral vertices at start
        for (let i = 0; i < this.degree; i++) {
            const vertex = new CurveVertex(
                -i - 1,              // Negative IDs for peripheral
                this.knots[i],       // Use first knots
                -i - 1,             // Negative indices for peripheral
                true                // Mark as peripheral
            );
            mesh.peripheralVertices.push(vertex);
        }

        // Create regular vertices
        for (let i = 0; i < this.numPoints; i++) {
            const vertex = new CurveVertex(
                i,
                this.knots[i + this.degree],  // Offset by degree
                i
            );
            mesh.vertices.push(vertex);
        }

        // Create peripheral vertices at end
        for (let i = 0; i < this.degree; i++) {
            const vertex = new CurveVertex(
                this.numPoints + i,
                this.knots[this.numPoints + this.degree + i],
                this.numPoints + i,
                true
            );
            mesh.peripheralVertices.push(vertex);
        }

        // Create edges including peripheral vertices
        this.createEdgesWithPeripheral(mesh);
    }

    private createEdgesWithPeripheral(mesh: BSplineCurveMesh): void {
        let id = 0;
        const allVertices = [
            ...mesh.peripheralVertices.slice(0, this.degree),  // Start peripheral
            ...mesh.vertices,                                  // Regular vertices
            ...mesh.peripheralVertices.slice(this.degree)      // End peripheral
        ];

        // Create edges connecting all vertices in sequence
        for (let i = 0; i < allVertices.length - 1; i++) {
            const start = allVertices[i];
            const end = allVertices[i + 1];
            
            const edge = new CurveEdge(id++, start, end);
            mesh.edges.push(edge);
            
            start.edges.push(edge);
            end.edges.push(edge);
        }
    }

    private createPeriodicMesh(mesh: BSplineCurveMesh): void {
        // For periodic case, we wrap around without peripheral vertices
        for (let i = 0; i < this.numPoints; i++) {
            const vertex = new CurveVertex(
                i,
                this.knots[i + this.degree],
                i
            );
            mesh.vertices.push(vertex);
        }

        // Create edges with wrap-around
        let id = 0;
        for (let i = 0; i < this.numPoints; i++) {
            const start = mesh.vertices[i];
            const end = mesh.vertices[(i + 1) % this.numPoints];
            
            const edge = new CurveEdge(id++, start, end);
            mesh.edges.push(edge);
            
            start.edges.push(edge);
            end.edges.push(edge);
        }
    }
}

class RegularBSplineCurve {
    private mesh: BSplineCurveMesh;

    constructor(
        numPoints: number,
        knots: number[],
        isClosed: boolean = false,
        degree: number = 3
    ) {
        const builder = new BSplineCurveMeshBuilder(
            numPoints,
            knots,
            isClosed,
            degree
        );
        this.mesh = builder.build();
    }

    // Get all vertices including peripheral
    getAllVertices(): CurveVertex[] {
        if (this.mesh.isClosed) {
            return this.mesh.vertices;
        }
        return [
            ...this.mesh.peripheralVertices.slice(0, this.degree),
            ...this.mesh.vertices,
            ...this.mesh.peripheralVertices.slice(this.degree)
        ];
    }

    // Get only the regular vertices
    getRegularVertices(): CurveVertex[] {
        return this.mesh.vertices;
    }

    // Get peripheral vertices
    getPeripheralVertices(): CurveVertex[] {
        return this.mesh.peripheralVertices;
    }
}

// Example usage
const openCurve = new RegularBSplineCurve(
    5,  // numPoints
    [0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1], // knots
    false, // open
    3      // cubic
);

const closedCurve = new RegularBSplineCurve(
    5,  // numPoints
    [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], // periodic knots
    true,  // closed
    3      // cubic
);


class BSplineBasisCalculator {
    constructor(
        private degree: number,
        private knots: number[],
        private isClosed: boolean
    ) {}

    // Calculate basis function N_i,p(u) using Cox-de Boor recursion
    calculateBasis(i: number, p: number, u: number): number {
        // Base case for degree 0
        if (p === 0) {
            return (u >= this.knots[i] && u < this.knots[i + 1]) ? 1.0 : 0.0;
        }

        // Cox-de Boor recursion formula
        let value = 0.0;

        // First term
        const d1 = this.knots[i + p] - this.knots[i];
        if (d1 !== 0) {
            value += ((u - this.knots[i]) / d1) * 
                    this.calculateBasis(i, p - 1, u);
        }

        // Second term
        const d2 = this.knots[i + p + 1] - this.knots[i + 1];
        if (d2 !== 0) {
            value += ((this.knots[i + p + 1] - u) / d2) * 
                    this.calculateBasis(i + 1, p - 1, u);
        }

        return value;
    }

    // Calculate all non-zero basis functions at u
    calculateAllNonzeroBasis(u: number): Map<number, number> {
        const result = new Map<number, number>();
        
        // Find knot span
        const span = this.findKnotSpan(u);
        
        // Calculate basis functions that are non-zero at u
        for (let i = span - this.degree; i <= span; i++) {
            const value = this.calculateBasis(i, this.degree, u);
            if (value > 0) {
                result.set(i, value);
            }
        }

        return result;
    }
}

class BSplineCurveEvaluator {
    private basisCalculator: BSplineBasisCalculator;

    constructor(private curve: RegularBSplineCurve) {
        this.basisCalculator = new BSplineBasisCalculator(
            curve.degree,
            curve.knots,
            curve.isClosed
        );
    }

    // Evaluate curve point at parameter u
    evaluate(u: number): number[] {
        // Get all vertices including peripheral
        const vertices = this.curve.getAllVertices();
        
        // Calculate non-zero basis functions at u
        const basisValues = this.basisCalculator.calculateAllNonzeroBasis(u);
        
        // Initialize result
        const result = [0.0, 0.0];  // Assuming 2D points
        
        // Sum up contributions from all non-zero basis functions
        for (const [i, basis] of basisValues) {
            // Include peripheral vertices in calculation
            const vertex = vertices[i + this.curve.degree];  // Offset by degree
            const point = vertex.controlPoint;
            
            result[0] += basis * point[0];
            result[1] += basis * point[1];
        }
        
        return result;
    }

    // Demonstrate effect of peripheral vertices
    demonstrateBoundaryBehavior(): void {
        // Example showing basis functions near boundary
        const u = 0.1;  // Near start of curve
        
        console.log("Basis functions near start:");
        
        // Without peripheral vertices (incorrect)
        const regularBasis = this.basisCalculator
            .calculateAllNonzeroBasis(u);
        
        console.log("Without peripheral vertices:");
        for (const [i, value] of regularBasis) {
            console.log(`N_{${i},${this.curve.degree}}(${u}) = ${value}`);
        }
        
        // With peripheral vertices (correct)
        const withPeripheral = this.basisCalculator
            .calculateAllNonzeroBasis(u);
            
        console.log("\nWith peripheral vertices:");
        for (const [i, value] of withPeripheral) {
            const isPeripheral = i < 0;
            console.log(
                `N_{${i},${this.curve.degree}}(${u}) = ${value}` +
                (isPeripheral ? " (peripheral)" : "")
            );
        }
    }
}

// Example showing the importance of peripheral vertices
class BoundaryBehaviorDemo {
    static demonstrate(): void {
        // Create open curve with and without peripheral vertices
        const knotsWithPeripheral = [0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4];
        const curveWithPeripheral = new RegularBSplineCurve(
            5,              // numPoints
            knotsWithPeripheral,
            false,          // open
            3              // cubic
        );

        const evaluator = new BSplineCurveEvaluator(curveWithPeripheral);

        // Demonstrate basis functions near boundary
        console.log("Evaluating near start of curve (u = 0.1):");
        const p1 = evaluator.evaluate(0.1);
        console.log("Point with peripheral vertices:", p1);

        // Show basis functions
        evaluator.demonstrateBoundaryBehavior();

        // Key effects of peripheral vertices:
        // 1. Proper partition of unity at boundaries
        // 2. Correct end derivatives
        // 3. Proper interpolation of end points
    }
}

// Example usage
BoundaryBehaviorDemo.demonstrate();
