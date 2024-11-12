// import { Vector, VectorSpace } from "./vector-space";



// // control-points.ts
// export class ControlPoints<K extends number, V extends Vector> {
//     private readonly points: V[];
//     private readonly isPeriodic: boolean;
//     private readonly vectorSpace: VectorSpace<K, V>;

//     constructor(
//         points: V[], 
//         vectorSpace: VectorSpace<K, V>,
//         isPeriodic: boolean = false
//     ) {
//         this.points = [...points];
//         this.vectorSpace = vectorSpace;
//         this.isPeriodic = isPeriodic;
//     }

//     public getPoint(index: number): V {
//         if (this.isPeriodic) {
//             return this.getPeriodicPoint(index);
//         }
//         return this.points[index];
//     }

//     public insert(
//         span: number, 
//         alpha: number
//     ): ControlPoints<K, V> {
//         if (this.isPeriodic) {
//             return this.insertPeriodic(span, alpha);
//         }
//         return this.insertSingle(span, alpha);
//     }

//     private getPeriodicPoint(index: number): V {
//         const baseIndex = index % this.points.length;
//         return this.points[baseIndex];
//     }

//     // ... other control point methods
// }

// // blossom.ts
// export class BSplineBlossom<K extends number, V extends Vector> {
//     private readonly degree: number;
//     private readonly knots: Knots;
//     private readonly controlPoints: ControlPoints<K, V>;
//     private readonly vectorSpace: VectorSpace<K, V>;

//     constructor(
//         degree: number,
//         knots: number[],
//         controlPoints: V[],
//         vectorSpace: VectorSpace<K, V>,
//         isPeriodic: boolean = false
//     ) {
//         this.degree = degree;
//         this.controlPoints = new ControlPoints(
//             controlPoints,
//             vectorSpace,
//             isPeriodic
//         );
//         this.knots = new Knots(knots, isPeriodic);
//         this.vectorSpace = vectorSpace;
//     }

//     public blossomValue(parameters: number[]): V {
//         if (parameters.length !== this.degree) {
//             throw new Error('Number of parameters must equal degree');
//         }

//         return this.computeBlossom(
//             0,
//             this.controlPoints.length - 1,
//             parameters
//         );
//     }

//     private computeBlossom(
//         left: number,
//         right: number,
//         parameters: number[]
//     ): V {
//         // Base case: degree 0
//         if (parameters.length === 0) {
//             return this.controlPoints.getPoint(left);
//         }

//         // Recursive case
//         const u = parameters[0];
//         const remainingParams = parameters.slice(1);
        
//         const span = this.knots.getSpan(u);
//         const alpha = this.computeAlpha(u, span);

//         const left_result = thi        const right_result = this.computeBlossom(
//             span,
//             right,
//             remainingParams
//         );

//         return this.vectorSpace.add(
//             this.vectorSpace.scale((1 - alpha) as K, left_result),
//             this.vectorSpace.scale(alpha as K, right_result)
//         );
//     }

//     public evaluate(u: number): V {
//         const parameters = new Array(this.degree).fill(u);
//         return this.blossomValue(parameters);
//     }

//     public derivative(u: number, order: number): V {
//         if (order > this.degree) {
//             throw new Error('Order exceeds degree');
//         }

//         const parameters: number[] = [];
        
//         for (let i = 0; i < this.degree - order; i++) {
//             parameters.push(u);
//         }
        
//         for (let i = 0; i < order; i++) {
//             parameters.push(u + i * 1e-6);
//         }

//         return this.blossomValue(parameters);
//     }

//     public insertKnot(u: number): BSplineBlossom<K, V> {
//         const newKnots = this.knots.insert(u);
//         const span = this.knots.getSpan(u);
//         const alpha = this.computeAlpha(u, span);
//         const newPoints = this.controlPoints.insert(span, alpha);

//         return new BSplineBlossom(
//             this.degree,
//             newKnots.values,
//             newPoints.points,
//             this.vectorSpace,
//             this.knots.isPeriodic
//         );
//     }

//     public elevateDegree(): BSplineBlossom<K, V> {
//         // ... degree elevation implementation
//     }

//     private computeAlpha(u: number, span: number): number {
//         return (u - this.knots.getValue(span)) / 
//                (this.knots.getValue(span + 1) - this.knots.getValue(span));
//     }
// }

// // Example usage:
// const bspline = new BSplineBlossom(
//     3, // degree
//     [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], // knots
//     [[0,0], [1,1], [2,-1], [3,0], [0,0], [1,1], [2,-1]], // points
//     vectorSpace,
//     true // isPeriodic
// );




//  /**
//  * Finds the knot span index using binary search
//  */
//  function findSpanIndex(
//     t: number,
//     degree: number,
//     knots: ReadonlyArray<number>
// ): number {
//     const n = knots.length - degree - 2 // number of control points
//     if (t >= knots[n + 1]) return n;
//     if (t <= knots[degree]) return degree;
//     let low = degree;
//     let high = n + 1;
//     let mid = Math.floor((low + high) / 2);
//     while (t < knots[mid] || t >= knots[mid + 1]) {
//         if (t < knots[mid]) {
//             high = mid;
//         } else {
//             low = mid;
//         }
//         mid = Math.floor((low + high) / 2);
//     }
//     return mid;
// }
