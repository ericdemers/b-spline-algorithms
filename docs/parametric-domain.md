# Parametric Domain Space

- Definition: A parametric domain refers to the open connected sets in a closed, compact, orientable manifold that can be represented parametrically by an atlas.

## The concept represents:

- A domain for parameter values
- Has a topological structure
- Contains singular points
- Manages local charts and transitions

## It is not a vector space because

The parameteric domain space is a manifold-like structure, but with some special properties due to the spline nature.

- It can be periodic
- It can have singularities
- It might not have global coordinates
- The topology might be non-trivial

1. Explore how to handle domain mapping for T-splines?

2. Discuss how the domain interacts with basis functions?

3. Look at how to handle domain subdivision for adaptive refinement?

4. Consider special cases like singular points in the domain?

- Local Charts:

  - Provide smooth transition between regular and singular regions

  - Define valid neighborhoods around singularities

  - Handle coordinate transformations

- Basis Functions:

  - Need special treatment near singularities

  - Must maintain desired continuity

  - May require blending between regular and singular regions

- Geometric Properties:

  - Tangent cone structure at extraordinary points

  - Controlled degradation of smoothness

  - Preservation of important geometric properties

## For T-splines with general topology, generating points for visualization is more complex than regular B-splines due to extraordinary points and T-junctions. Here's an approach:

Key aspects of this approach:

1. Domain Decomposition:

- Split domain into cells based on T-mesh structure

- Handle each cell type differently

- Maintain proper transitions between cells

2. Extraordinary Points:

- Use polar-like sampling around extraordinary points

- Adapt sampling density based on valence

- Ensure proper sector transitions

3. T-junctions:

- Handle transition regions carefully

- Adapt sampling density near T-junctions

- Maintain smooth transitions

### For WebGL rendering, we need to create a triangulation of the surface. For T-splines with general topology, we need to handle different cases carefully:

1. Mesh Generation:

- Adaptive triangulation based on local features

- Special handling of extraordinary points

- Conforming mesh at T-junctions

2. Rendering Features:

- WebGL-ready mesh format

- Normal computation for smooth shading

- Parameter space coordinates for texturing

3. Quality Control:

- Adaptive resolution based on curvature

- Proper transitions between different regions

- Crack-free triangulation

4. Other aspect to consider

- More sophisticated triangulation strategies?

- Specific WebGL shader implementations?

- How to handle texture mapping?

- Adaptive refinement based on view parameters?

## Vertex Edge Face relationships

1. Vertex-Edge:

- Each vertex knows its incident edges

- Each edge knows its endpoint vertices

- Vertices can be regular (valence 4), extraordinary (valence ≠ 4), or T-junctions (valence 3)

2. Edge-Face:

- Each edge is bounded by 1 face (boundary) or 2 faces (interior)

- Each face has a loop of edges forming its boundary [1]

- Edges connect adjacent faces

3. Vertex-Face:

- Each vertex knows its adjacent faces

- Each face has an ordered list of vertices

- Faces share vertices along common edges

4. Face-Face:

- Faces are connected through shared edges

- Each face knows its neighboring faces

- T-junctions create special face relationships
