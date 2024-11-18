import { PeriodicControlPoints } from "./periodicControlPoints";
import { PeriodicKnots } from "./periodicKnots2";

export class PeriodicBspline {


    constructor(
        private readonly knots: PeriodicKnots, 
        private readonly controlPoints: PeriodicControlPoints, 
        private readonly degree: number) {
        this.validateBSpline();
    }
    validateBSpline() {
        if(this.knots.getPatternLength() !== this.controlPoints.getNumControlPoints()) {
            throw new Error('The number of control points must be equal to the number of knots');
        }

        if(this.degree < 1) {
            throw new Error('Degree must be greater than or equal to 1');
        }
    }

    //Unwind b-spline curve to uses standard open curve knot insertion algorithms
    unwind(insertKnotPosition: number) {
        const relevantKnotInfo = this.knots.getRelevantKnotsAndIndices(
            insertKnotPosition,
            this.degree
        );
        
        const relevantControlPoints = this.controlPoints.getRelevantControlPoints(
            this.knots,
            insertKnotPosition,
            this.degree
        );

        // Create temporary non-periodic B-spline for insertion using
        // relevantKnotInfo.knots and relevantControlPoints
    }

    

    elevateDegree() {
        const newKnots = this.knots.elevateDegree();
        const newControlPoints = new PeriodicControlPoints(this.controlPoints.getControlPoints());
        return new PeriodicBspline(newKnots, newControlPoints, this.degree + 1);
    }
}

function boehmKnotInsertion(
    knots: number[],
    controlPoints: Point[],
    degree: number,
    insertKnot: number
): { newKnots: number[], newControlPoints: Point[] } {
    const n = controlPoints.length - 1;
    const r = findKnotSpan(knots, insertKnot);
    
    const newControlPoints = [...controlPoints];
    const newKnots = [...knots.slice(0, r + 1), insertKnot, ...knots.slice(r + 1)];
    
    for (let i = r - degree + 1; i <= r; i++) {
        const alpha = (insertKnot - knots[i]) / (knots[i + degree] - knots[i]);
        newControlPoints[i] = interpolate(controlPoints[i - 1], controlPoints[i], alpha);
    }
    
    return { newKnots, newControlPoints };
}

function unrollPeriodicKnotVector(
    periodicKnots: number[],  // The base periodic knot pattern
    degree: number,
    insertPosition: number    // The parameter value where we want to insert
): number[] {
    // Get the pattern length (assuming periodicKnots contains one complete pattern)
    const patternLength = periodicKnots.length;
    
    // Find which periodic segment contains our insertion point
    const segment = Math.floor(insertPosition / patternLength);
    
    // We need (degree + 1) knots before the insertion point and
    // (degree + 1) knots after for the basis functions
    const knotsNeeded = 2 * (degree + 1);
    
    // Calculate how many pattern repetitions we need
    const startSegment = segment - 1;  // One segment before
    const endSegment = segment + 1;    // One segment after
    
    const unrolledKnots: number[] = [];
    
    // Generate the necessary knots
    for (let seg = startSegment; seg <= endSegment; seg++) {
        for (let i = 0; i < patternLength; i++) {
            const shiftedKnot = periodicKnots[i] + (seg * patternLength);
            unrolledKnots.push(shiftedKnot);
        }
    }
    
    // Find the relevant range of knots needed for the insertion
    const insertSpan = findKnotSpan(unrolledKnots, insertPosition);
    const startIndex = Math.max(0, insertSpan - degree);
    const endIndex = Math.min(unrolledKnots.length - 1, insertSpan + degree + 1);
    
    // Extract only the needed portion
    return unrolledKnots.slice(startIndex, endIndex + 1);
}

// Helper function to find the knot span
function findKnotSpan(knots: number[], u: number): number {
    // Find the first knot that is greater than u
    for (let i = 0; i < knots.length - 1; i++) {
        if (knots[i] <= u && u < knots[i + 1]) {
            return i;
        }
    }
    return knots.length - 2;
}

