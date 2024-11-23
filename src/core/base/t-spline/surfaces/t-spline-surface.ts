interface SurfaceSampler {
    // Generate points covering the domain
    generateSamplePoints(): ParameterPoint[];
    // Adapt sampling density based on local features
    getLocalDensity(parameter: number[]): number;
}

class TSplineSampler implements SurfaceSampler {
    constructor(
        private domain: ParametricDomain,
        private baseResolution: number,
        private adaptiveThreshold: number
    ) {}

    generateSamplePoints(): ParameterPoint[] {
        const points: ParameterPoint[] = [];
        const cells = this.decomposeDomain();

        for (const cell of cells) {
            if (cell.hasExtraordinaryPoint()) {
                points.push(...this.sampleExtraordinaryCell(cell));
            } else if (cell.hasTJunction()) {
                points.push(...this.sampleTJunctionCell(cell));
            } else {
                points.push(...this.sampleRegularCell(cell));
            }
        }

        return points;
    }

    private decomposeDomain(): Cell[] {
        // Decompose domain into cells based on T-mesh structure
        // Each cell is either regular, contains T-junction, or extraordinary point
        return this.domain.decompose();
    }

    private sampleExtraordinaryCell(cell: Cell): ParameterPoint[] {
        const center = cell.getExtraordinaryPoint();
        const valence = cell.getValence();
        
        // Use polar-like coordinates around extraordinary point
        return this.generatePolarSamples(
            center,
            cell.getRadius(),
            valence,
            (r, theta) => this.adaptSamplingDensity(r, theta, cell)
        );
    }

    private sampleTJunctionCell(cell: Cell): ParameterPoint[] {
        // Handle transition near T-junction
        const samples: ParameterPoint[] = [];
        const transitions = cell.getTransitionRegions();

        for (const region of transitions) {
            const localDensity = this.getLocalDensity(region.center);
            samples.push(...this.generateTransitionSamples(
                region,
                localDensity
            ));
        }

        return samples;
    }

    private generatePolarSamples(
        center: number[],
        radius: number,
        valence: number,
        densityFn: (r: number, theta: number) => number
    ): ParameterPoint[] {
        const samples: ParameterPoint[] = [];
        const sectors = valence;
        
        // Radial sampling
        for (let r = 0; r <= radius; r += this.getRadialStep(r)) {
            // Angular sampling
            for (let s = 0; s < sectors; s++) {
                const theta = (2 * Math.PI * s) / sectors;
                const density = densityFn(r, theta);
                const angleStep = this.getAngularStep(r, density);
                
                for (let t = 0; t < 2 * Math.PI / sectors; t += angleStep) {
                    const angle = theta + t;
                    samples.push({
                        parameters: this.polarToParameter(r, angle, center),
                        isExtraordinary: r === 0,
                        sectorIndex: s
                    });
                }
            }
        }

        return samples;
    }

    private adaptSamplingDensity(
        r: number, 
        theta: number, 
        cell: Cell
    ): number {
        const base = this.baseResolution;
        const curvature = cell.getLocalCurvature(r, theta);
        const featureSize = cell.getLocalFeatureSize(r, theta);
        
        return base * Math.min(
            1.0,
            Math.max(
                curvature * this.adaptiveThreshold,
                1.0 / featureSize
            )
        );
    }

    private generateTransitionSamples(
        region: TransitionRegion,
        density: number
    ): ParameterPoint[] {
        const samples: ParameterPoint[] = [];
        const [uMin, uMax] = region.uRange;
        const [vMin, vMax] = region.vRange;
        
        const stepU = (uMax - uMin) / Math.max(2, density);
        const stepV = (vMax - vMin) / Math.max(2, density);

        for (let u = uMin; u <= uMax; u += stepU) {
            for (let v = vMin; v <= vMax; v += stepV) {
                samples.push({
                    parameters: [u, v],
                    isTransition: true,
                    transitionBlend: region.getBlendFactor([u, v])
                });
            }
        }

        return samples;
    }
}

// Usage example
const sampler = new TSplineSampler(
    domain,
    baseResolution: 32,
    adaptiveThreshold: 0.1
);

const visualizationPoints = sampler.generateSamplePoints();

interface SurfaceTriangulator {
    // Generate triangulation covering the domain
    generateMesh(): TriangularMesh;
    // Get triangulation topology
    getTopology(): MeshTopology;
}

interface TriangularMesh {
    vertices: number[];      // Flattened array of vertex positions
    normals: number[];      // Flattened array of vertex normals
    uvs: number[];          // Texture coordinates/parameter space coords
    indices: number[];      // Triangle indices
    features: MeshFeatures; // Special features for rendering
}

interface MeshFeatures {
    extraordinaryPoints: number[];  // Indices of extraordinary vertices
    tJunctions: number[];          // Indices of T-junction vertices
    boundaries: number[];          // Indices of boundary vertices
}

class TSplineTriangulator implements SurfaceTriangulator {
    constructor(
        private domain: ParametricDomain,
        private baseResolution: number,
        private adaptiveThreshold: number
    ) {}

    generateMesh(): TriangularMesh {
        const cells = this.decomposeDomain();
        const meshBuilder = new MeshBuilder();

        for (const cell of cells) {
            if (cell.hasExtraordinaryPoint()) {
                this.triangulateExtraordinaryCell(cell, meshBuilder);
            } else if (cell.hasTJunction()) {
                this.triangulateTJunctionCell(cell, meshBuilder);
            } else {
                this.triangulateRegularCell(cell, meshBuilder);
            }
        }

        return meshBuilder.build();
    }

    private triangulateExtraordinaryCell(
        cell: Cell, 
        builder: MeshBuilder
    ): void {
        const center = cell.getExtraordinaryPoint();
        const valence = cell.getValence();
        
        // Create center vertex
        const centerIndex = builder.addVertex({
            position: center,
            isExtraordinary: true,
            valence: valence
        });

        // Create vertices in sectors around extraordinary point
        const sectorVertices = this.generateSectorVertices(
            cell,
            builder
        );

        // Create triangles connecting sectors
        this.generateSectorTriangles(
            centerIndex,
            sectorVertices,
            builder
        );
    }

    private generateSectorVertices(
        cell: Cell,
        builder: MeshBuilder
    ): number[][] {
        const sectors: number[][] = [];
        const valence = cell.getValence();
        
        for (let s = 0; s < valence; s++) {
            const sectorVerts = this.generatePolarVertices(
                cell,
                s,
                builder
            );
            sectors.push(sectorVerts);
        }

        return sectors;
    }

    private generatePolarVertices(
        cell: Cell,
        sectorIndex: number,
        builder: MeshBuilder
    ): number[] {
        const vertices: number[] = [];
        const radialSteps = this.getRadialResolution(cell);
        const angularSteps = this.getAngularResolution(cell);
        
        for (let r = 0; r <= radialSteps; r++) {
            const radius = (r / radialSteps) * cell.getRadius();
            
            for (let a = 0; a <= angularSteps; a++) {
                const angle = this.getSectorAngle(sectorIndex, a, angularSteps);
                const params = this.polarToParameter(radius, angle);
                
                vertices.push(builder.addVertex({
                    position: params,
                    polar: [radius, angle],
                    sectorIndex: sectorIndex
                }));
            }
        }

        return vertices;
    }

    private triangulateTJunctionCell(
        cell: Cell, 
        builder: MeshBuilder
    ): void {
        const transitions = cell.getTransitionRegions();
        
        for (const region of transitions) {
            // Create conforming triangulation at T-junction
            this.generateTransitionTriangles(region, builder);
        }
    }

    private generateTransitionTriangles(
        region: TransitionRegion,
        builder: MeshBuilder
    ): void {
        // Create vertices with proper density adaptation
        const vertices = this.generateTransitionVertices(region, builder);
        
        // Create triangles ensuring conforming mesh
        this.createConformingTriangles(vertices, region, builder);
    }

    ////////

    private decomposeDomain(): Cell[] {
        const cells: Cell[] = [];
        const tMesh = this.domain.getTMesh();
        
        // Get critical points (extraordinary points and T-junctions)
        const extraordinaryPoints = tMesh.getExtraordinaryPoints();
        const tJunctions = tMesh.getTJunctions();
        
        // Create initial cell decomposition
        const initialCells = this.createInitialDecomposition(
            extraordinaryPoints,
            tJunctions
        );

        // Refine cells based on features
        for (const cell of initialCells) {
            cells.push(...this.refineCellIfNeeded(cell));
        }

        return cells;
    }

    private createInitialDecomposition(
        extraordinaryPoints: ExtraordinaryPoint[],
        tJunctions: TJunction[]
    ): Cell[] {
        // Create Quad-tree like structure
        const quadTree = new QuadTree(this.domain.getBounds());
        
        // Insert critical points
        for (const point of extraordinaryPoints) {
            quadTree.insert({
                position: point.position,
                type: 'extraordinary',
                data: point
            });
        }

        for (const junction of tJunctions) {
            quadTree.insert({
                position: junction.position,
                type: 't-junction',
                data: junction
            });
        }

        // Initial subdivision based on feature proximity
        return this.subdivideQuadTree(quadTree);
    }

    private subdivideQuadTree(quadTree: QuadTree): Cell[] {
        const cells: Cell[] = [];
        const leaves = quadTree.getLeafNodes();

        for (const leaf of leaves) {
            const features = this.analyzeCellFeatures(leaf);
            
            if (this.needsSubdivision(leaf, features)) {
                // Subdivide further if features are too close
                const subCells = this.subdivideCell(leaf);
                cells.push(...subCells);
            } else {
                // Create cell with identified features
                cells.push({
                    bounds: leaf.bounds,
                    type: this.determineCellType(features),
                    features: features
                });
            }
        }

        return cells;
    }

    private analyzeCellFeatures(
        leaf: QuadTreeNode
    ): CellFeatures {
        const features: CellFeatures = {};
        const points = leaf.getPoints();

        // Group points by type
        features.extraordinaryPoints = points
            .filter(p => p.type === 'extraordinary')
            .map(p => p.data);

        features.tJunctions = points
            .filter(p => p.type === 't-junction')
            .map(p => p.data);

        // Check for boundary intersection
        if (this.domain.hasBoundary()) {
            features.boundaries = this.findBoundaryIntersections(leaf.bounds);
        }

        return features;
    }

    private determineCellType(features: CellFeatures): CellType {
        const hasExtraordinary = features.extraordinaryPoints?.length > 0;
        const hasTJunction = features.tJunctions?.length > 0;
        
        if (hasExtraordinary && hasTJunction) {
            return CellType.MIXED;
        }
        if (hasExtraordinary) {
            return CellType.EXTRAORDINARY;
        }
        if (hasTJunction) {
            return CellType.T_JUNCTION;
        }
        return CellType.REGULAR;
    }

    private needsSubdivision(
        leaf: QuadTreeNode, 
        features: CellFeatures
    ): boolean {
        // Check feature density
        const featureCount = (features.extraordinaryPoints?.length ?? 0) +
                           (features.tJunctions?.length ?? 0);
                           
        if (featureCount > 1) {
            return true;
        }

        // Check cell size vs feature influence radius
        if (features.extraordinaryPoints?.length === 1) {
            const point = features.extraordinaryPoints[0];
            return this.cellSizeTooLargeForFeature(leaf.bounds, point);
        }

        if (features.tJunctions?.length === 1) {
            const junction = features.tJunctions[0];
            return this.cellSizeTooLargeForFeature(leaf.bounds, junction);
        }

        return false;
    }

    private refineCellIfNeeded(cell: Cell): Cell[] {
        if (cell.type === CellType.REGULAR) {
            return [cell];
        }

        // Refine based on feature type
        switch (cell.type) {
            case CellType.EXTRAORDINARY:
                return this.refineExtraordinaryCell(cell);
            case CellType.T_JUNCTION:
                return this.refineTJunctionCell(cell);
            case CellType.MIXED:
                return this.refineMixedCell(cell);
            default:
                return [cell];
        }
    }

    private refineExtraordinaryCell(cell: Cell): Cell[] {
        const point = cell.features.extraordinaryPoints![0];
        const valence = point.valence;
        
        // Create polar-like subdivision around extraordinary point
        return this.createPolarSubdivision(
            point.position,
            valence,
            cell.bounds
        );
    }

    private createPolarSubdivision(
        center: number[],
        valence: number,
        bounds: { min: number[]; max: number[] }
    ): Cell[] {
        const cells: Cell[] = [];
        const radius = this.computeInfluenceRadius(bounds);
        
        // Create sectors around extraordinary point
        for (let i = 0; i < valence; i++) {
            
            cells.push({
                bounds: this.computeSectorBounds(
                    center, 
                    radius, 
                    angle, 
                    nextAngle
                ),
                type: CellType.EXTRAORDINARY,
                features: {
                    extraordinaryPoints: [{ 
                        position: center,
                        valence: valence,
                        sectorIndex: i 
                    }]
                }
            });
        }

        return cells;
    }
}

// WebGL rendering setup
class TSplineSurfaceRenderer {
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

    updateMesh(mesh: TriangularMesh): void {
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

// Usage example
const triangulator = new TSplineTriangulator(
    domain,
    baseResolution: 32,
    adaptiveThreshold: 0.1
);

const mesh = triangulator.generateMesh();
const renderer = new TSplineSurfaceRenderer(canvas);
renderer.updateMesh(mesh);
renderer.render();

