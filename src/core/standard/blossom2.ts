/* Mathematical Foundation of Blossoms (also known as polar form)
 Properties of blossoms:
 1. Symmetry: b(u₁,u₂,u₃) = b(u₂,u₁,u₃) = b(u₃,u₁,u₂)
 2. Multi-affine: b(αu₁+(1-α)v₁,u₂,u₃) = αb(u₁,u₂,u₃) + (1-α)b(v₁,u₂,u₃)
 3. Diagonal property: b(u,u,u) = p(u) where p is the B-spline curve
 */

import { Vector, VectorSpace } from "./vector-space";

 class BSplineCurveBlossom<K extends number, V extends Vector> {
     private degree: number;
     private knots: number[];
     private controlPoints: V[];
     private vectorSpace: VectorSpace<K, V>;
 
     constructor(
         degree: number,
         knots: number[],
         controlPoints: V[],
         vectorSpace: VectorSpace<K, V>
     ) {
         this.degree = degree;
         this.knots = knots;
         this.controlPoints = controlPoints;
         this.vectorSpace = vectorSpace;
     }
 
     // Compute blossom value b(u₁,u₂,...,uₚ)
     public blossomValue(parameters: number[]): V {
         if (parameters.length !== this.degree) {
             throw new Error('Number of parameters must equal degree');
         }
 
         return this.computeBlossom(
             0,
             this.controlPoints.length - 1,
             parameters
         );
     }
 
     private computeBlossom(
         left: number,
         right: number,
         parameters: number[]
     ): V {
         // Base case: degree 0
         if (parameters.length === 0) {
             return this.controlPoints[left];
         }
 
         // Recursive case
         const u = parameters[0];
         const remainingParams = parameters.slice(1);
         
         // Find knot span
         let span = right;
         while (span > left && this.knots[span] > u) {
             span--;
         }
 
         // Compute alpha
         const alpha = (u - this.knots[span]) / 
                      (this.knots[span + 1] - this.knots[span]);
 
         // Recursive calls
         const left_result = this.computeBlossom(
             left,
             span,
             remainingParams
         );
         const right_result = this.computeBlossom(
             span,
             right,
             remainingParams
         );
 
         // Linear interpolation
         return this.vectorSpace.add(
             this.vectorSpace.scale((1 - alpha) as K, left_result),
             this.vectorSpace.scale(alpha as K, right_result)
         );
     }

     public evaluate(u: number): V {
        // Create parameter array with u repeated degree times
        const parameters = new Array(this.degree).fill(u);
        return this.blossomValue(parameters);
    }

    public derivative(u: number, order: number): V {
        if (order > this.degree) {
            throw new Error('Order exceeds degree');
        }

        // Parameters for derivative computation
        const parameters: number[] = [];
        
        // Add u (degree-order) times
        for (let i = 0; i < this.degree - order; i++) {
            parameters.push(u);
        }
        
        // Add parameters for derivative
        for (let i = 0; i < order; i++) {
            parameters.push(u + i * 1e-6); // Small offset
        }

        return this.blossomValue(parameters);
    }

    public elevateDegree(): {
        newPoints: V[],
        newKnots: number[]
    } {
        const newPoints: V[] = [];
        
        // Compute new control points using blossoms
        for (let i = 0; i < this.controlPoints.length + 1; i++) {
            // Parameters for elevated degree
            const parameters: number[] = [];
            
            // Fill parameters for degree elevation
            for (let j = 0; j <= this.degree; j++) {
                const t = j / this.degree;
                const u = this.parameterAtIndex(i, t);
                parameters.push(u);
            }

            newPoints.push(this.blossomValue(parameters));
        }

        // Create new knot vector with increased multiplicity
        const newKnots = this.elevateKnots();

        return { newPoints, newKnots };
    }

    private parameterAtIndex(i: number, t: number): number {
        const span = this.findKnotSpan(i);
        return (1 - t) * this.knots[span] + t * this.knots[span + 1];
    }

    private findKnotSpan(i: number): number {
        let span = this.knots.length - 2;
        while (span >= 0 && this.knots[span + 1] <= i) {
            span--;
        }
        return span;
    }

    private elevateKnots(): number[] {
        const newKnots: number[] = [];
        const m = this.knots.length;
        const n = this.controlPoints.length;
        const p = this.degree;

        // 1. Handle start of knot vector
        // Increase multiplicity at start by 1
        for (let i = 0; i <= p; i++) {
            newKnots.push(this.knots[0]);
        }

        // 2. Handle internal knots
        // Keep track of multiplicities
        let i = p + 1;
        while (i < m - p - 1) {
            // Find multiplicity of current knot
            let mult = 1;
            const currentKnot = this.knots[i];
            
            while (i + 1 < m && 
                   Math.abs(this.knots[i + 1] - currentKnot) < 1e-10) {
                mult++;
                i++;
            }

            // Add knot with increased multiplicity
            for (let j = 0; j < mult + 1; j++) {
                newKnots.push(currentKnot);
            }

            i++;
        }

        // 3. Handle end of knot vector
        // Increase multiplicity at end by 1
        for (let i = 0; i <= p + 1; i++) {
            newKnots.push(this.knots[m - 1]);
        }

        return newKnots;
    }

    // For periodic case
    private elevatePeriodicKnots(): number[] {
        const newKnots: number[] = [];
        const period = this.getPeriod();
        
        // 1. Handle base domain knots
        const baseKnots = this.elevateKnotsInDomain(0, period);
        newKnots.push(...baseKnots);

        // 2. Add periodic copies
        const numPeriods = Math.ceil(
            (this.degree + 2) / this.controlPoints.length
        );
        
        for (let i = 1; i <= numPeriods; i++) {
            const periodicKnots = baseKnots.map(k => k + i * period);
            newKnots.push(...periodicKnots);
        }

        return newKnots;
    }

    private elevateKnotsInDomain(
        start: number,
        end: number
    ): number[] {
        const domainKnots: number[] = [];
        const p = this.degree;

        // Find knots in domain
        const domainIndices: number[] = [];
        for (let i = 0; i < this.knots.length; i++) {
            if (this.knots[i] >= start && this.knots[i] <= end) {
                domainIndices.push(i);
            }
        }

        // Process knots in domain
        for (let i = 0; i < domainIndices.length; i++) {
            const idx = domainIndices[i];
            const currentKnot = this.knots[idx];

            // Find multiplicity
            let mult = 1;
            while (i + 1 < domainIndices.length && 
                   Math.abs(this.knots[domainIndices[i + 1]] - currentKnot) < 1e-10) {
                mult++;
                i++;
            }

            // Add knot with increased multiplicity
            for (let j = 0; j < mult + 1; j++) {
                domainKnots.push(currentKnot);
            }
        }

        return domainKnots;
    }

    // Helper method to get period
    private getPeriod(): number {
        const distinctKnots = Array.from(new Set(this.knots));
        return distinctKnots[distinctKnots.length - 1] - distinctKnots[0];
    }
 }
 


 /**
 * Finds the knot span index using binary search
 */
function findSpanIndex(
    t: number,
    degree: number,
    knots: ReadonlyArray<number>
): number {
    const n = knots.length - degree - 2 // number of control points
    if (t >= knots[n + 1]) return n;
    if (t <= knots[degree]) return degree;
    let low = degree;
    let high = n + 1;
    let mid = Math.floor((low + high) / 2);
    while (t < knots[mid] || t >= knots[mid + 1]) {
        if (t < knots[mid]) {
            high = mid;
        } else {
            low = mid;
        }
        mid = Math.floor((low + high) / 2);
    }
    return mid;
}

 
 