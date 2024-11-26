interface KnotVertex {
    connectedEdges: Set<number>;  // Indices into edges array
}

interface KnotEdge {
    startVertex: number;  // Index into vertices array
    endVertex: number;    // Index into vertices array
    length: number;       // Zero length for multiplicity
}

class GraphKnotStructure {
    protected readonly vertices: KnotVertex[];
    protected readonly edges: KnotEdge[];
    
    constructor(vertices: KnotVertex[], edges: KnotEdge[]) {
        this.vertices = vertices;
        this.edges = edges;
        this.validateStructure();
    }

    private validateStructure(): void {
        // Ensure the graph is properly connected
        // Check that all edge references are valid
        // Verify that zero-length edges are properly used for multiplicities
    }

    getParameterAtVertex(vertexIndex: number): number {
        let parameter = 0;
        // Sum lengths of all edges before this vertex
        for (let i = 0; i < vertexIndex; i++) {
            parameter += this.edges[i].length;
        }
        return parameter;
    }


    // Get multiplicity at a vertex by counting zero-length edges
    getMultiplicityAtVertex(vertexIndex: number): number {
        return 0
    }

    insertKnot(parameter: number): void {
        // Find the edge where the parameter falls
        // Update edges...
    }
    

    // Convert to traditional knot vector
    toKnotVector(): number[] {
        const knots: number[] = [];
        for (let i = 0; i < this.vertices.length; i++) {
            const parameter = this.getParameterAtVertex(i);
            const multiplicity = this.getMultiplicityAtVertex(i);
            
            // Add knot value according to multiplicity
            for (let j = 0; j < multiplicity; j++) {
                knots.push(parameter);
            }
        }
        return knots;
    }
}

class BezierKnotStructure extends GraphKnotStructure {
    static createForDegree(degree: number): BezierKnotStructure {
        // For cubic (degree 3), we need 8 vertices and 7 edges
        const vertexCount = 2 * (degree + 1);
        const edgeCount = vertexCount - 1;
        
        const vertices: KnotVertex[] = [];
        const edges: KnotEdge[] = [];

        // Create all vertices first
        for (let i = 0; i < vertexCount; i++) {
            vertices.push({ connectedEdges: new Set() });
        }

        // Create edges - following the chain structure
        for (let i = 0; i < edgeCount; i++) {
            const edgeIndex = edges.length;
            const isMiddleEdge = i === degree; // The one edge with non-zero length
            
            edges.push({
                startVertex: i,
                endVertex: i + 1,
                length: isMiddleEdge ? 1 : 0
            });
            
            // Connect edges to vertices
            vertices[i].connectedEdges.add(edgeIndex);
            vertices[i + 1].connectedEdges.add(edgeIndex);
        }

        return new BezierKnotStructure(vertices, edges);
    }
}

// For a cubic Bezier curve (degree 3):
/*
const vertices = [
    { connectedEdges: new Set([0]) },       // Vertex 0
    { connectedEdges: new Set([0, 1]) },    // Vertex 1
    { connectedEdges: new Set([1, 2]) },    // Vertex 2
    { connectedEdges: new Set([2, 3]) },    // Vertex 3
    { connectedEdges: new Set([3, 4]) },    // Vertex 4
    { connectedEdges: new Set([4, 5]) },    // Vertex 5
    { connectedEdges: new Set([5, 6]) },    // Vertex 6
    { connectedEdges: new Set([6]) }        // Vertex 7
];


const edges = [
    // First three edges have zero length (start multiplicity)
    { startVertex: 0, endVertex: 1, length: 0 },  // Edge 0
    { startVertex: 1, endVertex: 2, length: 0 },  // Edge 1
    { startVertex: 2, endVertex: 3, length: 0 },  // Edge 2
    
    // Middle edge with unit length
    { startVertex: 3, endVertex: 4, length: 1 },  // Edge 3
    
    // Last three edges have zero length (end multiplicity)
    { startVertex: 4, endVertex: 5, length: 0 },  // Edge 4
    { startVertex: 5, endVertex: 6, length: 0 },  // Edge 5
    { startVertex: 6, endVertex: 7, length: 0 }   // Edge 6
];
*/


class ClosedBSplineKnotStructure extends GraphKnotStructure {
    static createClosedCurve(degree: number, numControlPoints: number): ClosedBSplineKnotStructure {
        // For a cubic closed curve with n control points
        // We need n vertices and n edges
        // Each vertex connects to the next, and the last connects back to first
        
        const vertices: KnotVertex[] = [];
        const edges: KnotEdge[] = [];

        // Create vertices
        for (let i = 0; i < numControlPoints; i++) {
            vertices.push({ connectedEdges: new Set() });
        }

        // Create edges - all with equal length for uniform parameterization
        for (let i = 0; i < numControlPoints; i++) {
            const edgeIndex = edges.length;
            edges.push({
                startVertex: i,
                endVertex: (i + 1) % numControlPoints,  // Wrap around to 0 at the end
                length: 1  // Unit length for uniform parameterization
            });
            
            vertices[i].connectedEdges.add(edgeIndex);
            vertices[(i + 1) % numControlPoints].connectedEdges.add(edgeIndex);
        }

        return new ClosedBSplineKnotStructure(vertices, edges);
    }

    // Override to handle periodic parameter space
    getParameterAtVertex(vertexIndex: number): number {
        let parameter = 0;
        for (let i = 0; i < vertexIndex; i++) {
            parameter += this.edges[i].length;
        }
        return parameter % this.getTotalLength();
    }

    private getTotalLength(): number {
        return this.edges.reduce((sum, edge) => sum + edge.length, 0);
    }

    // Convert to traditional periodic knot vector
    toKnotVector(): number[] {
        const knots: number[] = [];
        const totalLength = this.getTotalLength();
        
        // For a cubic curve, we need to repeat the first 3 control points
        // at the end to maintain C2 continuity
        for (let i = -3; i <= this.vertices.length; i++) {
            knots.push((i * totalLength) / this.vertices.length);
        }
        
        return knots;
    }

    // Get the n vertices around a given vertex index
    getNeighboringVertices(centerIndex: number, count: number): number[] {
        const neighbors: number[] = [];
        for (let i = -Math.floor(count/2); i <= Math.floor(count/2); i++) {
            let index = (centerIndex + i) % this.vertices.length;
            if (index < 0) index += this.vertices.length;
            neighbors.push(index);
        }
        return neighbors;
    }

    // Get parameter value considering periodicity
    getPeriodicParameter(parameter: number): number {
        const totalLength = this.getTotalLength();
        return parameter - Math.floor(parameter / totalLength) * totalLength;
    }
}

/*
// Example usage for a cubic closed B-spline with 5 control points:
const closedCurve = ClosedBSplineKnotStructure.createClosedCurve(3, 5);

// This would create:
const vertices = [
    { connectedEdges: new Set([0, 4]) },    // Vertex 0
    { connectedEdges: new Set([0, 1]) },    // Vertex 1
    { connectedEdges: new Set([1, 2]) },    // Vertex 2
    { connectedEdges: new Set([2, 3]) },    // Vertex 3
    { connectedEdges: new Set([3, 4]) }     // Vertex 4
];

const edges = [
    { startVertex: 0, endVertex: 1, length: 1 },  // Edge 0
    { startVertex: 1, endVertex: 2, length: 1 },  // Edge 1
    { startVertex: 2, endVertex: 3, length: 1 },  // Edge 2
    { startVertex: 3, endVertex: 4, length: 1 },  // Edge 3
    { startVertex: 4, endVertex: 0, length: 1 }   // Edge 4 (closing edge)
];
*/



