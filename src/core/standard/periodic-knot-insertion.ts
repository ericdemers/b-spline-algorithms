// //Boehm's key insights for periodic cases:
// /**
//  * 1. Periodic Knot Vector Structure:
//    Original: [..., u₋₁, u₀, u₁, ..., uₙ, uₙ₊₁, ...]
//    Must satisfy: uᵢ₊ₖ = uᵢ + T (where T is the period)

//     2. Periodic Control Point Structure:
//    Pᵢ₊ₖ = Pᵢ (where k is number of control points)

//     3. Modified Insertion Algorithm:
//  */


// class PeriodicBoehm<K extends number, V> {
//     insertKnot(u: number): void {
//         // Map u to fundamental domain
//         const uMapped = this.mapToFundamentalDomain(u);
        
//         // Find affected control points
//         const span = this.findSpan(uMapped);
//         const affected = this.getAffectedPoints(span);
        
//         // Apply Boehm's formula with periodic indices
//         for (let i = span - this.degree + 1; i <= span; i++) {
//             const alpha = (u - this.knots[i]) / 
//                          (this.knots[i + this.degree] - this.knots[i]);
            
//             const periodicIdx = this.getPeriodicIndex(i);
//             const prevPeriodicIdx = this.getPeriodicIndex(i - 1);
            
//             this.newPoints[periodicIdx] = 
//                 this.vectorSpace.add(
//                     this.vectorSpace.scale(alpha as K, this.points[periodicIdx]),
//                     this.vectorSpace.scale(
//                         (1 - alpha) as K, 
//                         this.points[prevPeriodicIdx]
//                     )
//                 );
//         }
//     }

//     private getPeriodicIndex(i: number): number {
//         const n = this.controlPoints.length;
//         return ((i % n) + n) % n;  // Proper modulo for negative numbers
//     }
// }

// /*
// // Oslo algorithm for periodic cases
// class PeriodicOslo<K extends number, V> {
//     insertKnot(u: number): void {
//         // Map u to fundamental domain
//         const uMapped = this.mapToFundamentalDomain(u);

//         // Find affected control points
//         const span = this.findSpan(uMapped);
//         const affected = this.getAffectedPoints(span);

//         // Apply Oslo's formula with periodic indices
//         for (let i = span - this.degree + 1; i <= span; i++) {
//             const alpha = (u - this.knots[i]) /
//                          (this.knots[i + this.degree] - this.knots[i]);

//             const periodicIdx = this.getPeriodicIndex(i);
//             const prevPeriodicIdx = this.getPeriodicIndex(i - 1);

//             this.newPoints[periodicIdx] =
//                 this.vectorSpace.add(
//                     this.vectorSpace.scale(alpha as K, this.points[periodicIdx]),
//                     this.vectorSpace.scale(
//                         (1 - alpha) as K,
//                         this.points[prevPeriodicIdx]
//                     )
//                 );
//         }
//     }

//     private getPeriodicIndex(i: number): number {
//         const n = this.controlPoints.length;
//         return ((i % n) + n) % n;  // Proper modulo for negative numbers
//     }
// }
//     */

// /*
// The Oslo Algorithm's adaptation for periodic cases:

// 1. Discrete B-spline Functions:
//    aᵢ,ⱼ(u) represents the coefficient of the j-th new basis 
//    function in terms of the i-th original basis function

// 2. Periodic Modification:
//    aᵢ₊ₖ,ⱼ₊ₖ(u) = aᵢ,ⱼ(u - T)
// */

// class PeriodicOsloAlgorithm<K extends number, V> {
//     computeOsloCoefficients(
//         oldKnots: number[], 
//         newKnots: number[], 
//         degree: number
//     ): number[][] {
//         const coefficients: number[][] = [];
        
//         // Initialize with identity matrix
//         for (let i = 0; i < oldKnots.length - degree - 1; i++) {
//             coefficients[i] = new Array(newKnots.length - degree - 1).fill(0);
//             coefficients[i][i] = 1;
//         }
        
//         // Compute coefficients with periodicity
//         for (let r = 1; r <= degree; r++) {
//             for (let i = 0; i < oldKnots.length - degree - 1; i++) {
//                 for (let j = 0; j < newKnots.length - degree - 1; j++) {
//                     const periodicI = this.getPeriodicIndex(i);
//                     const periodicJ = this.getPeriodicIndex(j);
                    
//                     coefficients[periodicI][periodicJ] = 
//                         this.computePeriodicCoefficient(
//                             oldKnots, 
//                             newKnots, 
//                             i, 
//                             j, 
//                             r
//                         );
//                 }
//             }
//         }
        
//         return coefficients;
//     }
// }

// /*
// 1. Knot Vector Constraints:
//    uᵢ₊ₖ = uᵢ + T

// 2. Control Point Constraints:
//    Pᵢ₊ₖ = Pᵢ

// 3. Basis Function Constraints:
//    Nᵢ₊ₖ,ₚ(u + T) = Nᵢ,ₚ(u)
// */

// class PeriodicConstraints<K extends number, V> {
//     validatePeriodicConstraints(): boolean {
//         return (
//             this.checkParametricPeriodicity() &&
//             this.checkGeometricContinuity() &&
//             this.checkKnotPeriodicity()
//         );
//     }

//     enforcePeriodicConstraints(): void {
//         // Enforce parametric periodicity
//         this.enforceParametricPeriodicity();
        
//         // Enforce geometric continuity
//         this.enforceGeometricContinuity();
        
//         // Enforce knot periodicity
//         this.enforceKnotPeriodicity();
//     }

//     private enforceParametricPeriodicity(): void {
//         const period = this.getPeriod();
        
//         // Ensure curve is periodic
//         for (let i = 0; i < this.degree; i++) {
//             const t = i / (this.degree - 1);
//             const p1 = this.evaluate(t);
//             const p2 = this.evaluate(t + period);
            
//             if (!this.vectorSpace.equals(p1, p2)) {
//                 this.adjustControlPoints(t, p1, p2);
//             }
//         }
//     }
// }







// class PeriodicKnotInsertion<K extends number, V> {
//     insertKnot(u: number): void {
//         // Map u to fundamental domain
//         const uMapped = this.mapToFundamentalDomain(u);

//         // Find affected control points
//         const span = this.findSpan(uMapped);
//         const affected = this.getAffectedPoints(span);

//         // Apply Oslo's formula with periodic indices
//         for (let i = span - this.degree + 1; i <= span; i++) {
//             const alpha = (u - this.knots[i]) /
//                          (this.knots[i + this.degree] - this.knots[i]);

//             const periodicIdx = this.getPeriodicIndex(i);
//             const prevPeriodicIdx = this.getPeriodicIndex(i - 1);

//             this.newPoints[periodicIdx] =
//                 this.vectorSpace.add(
//                     this.vectorSpace.scale(alpha as K, this.points[periodicIdx]),
//                     this.vectorSpace.scale(
//                         (1 - alpha) as K,
//                         this.points[prevPeriodicIdx]
//                     )
//                 );
//         }
//     }

//     private getPeriodicIndex(i: number): number {
//         const n = this.controlPoints.length;
//         return ((i % n) + n) % n;  // Proper modulo for negative numbers
//     }
// }