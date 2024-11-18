/**
 * Represents periodic control points for closed B-spline curves.
 * Periodic control points repeat 
 * Note : The number of periodic knot in the pattern is equal
 * the number of periodic control point
 */

import { PeriodicKnots } from "./periodicKnots2";
import { Vector2D } from "./vector-space";

export class PeriodicControlPoints {


    constructor(private readonly points: Vector2D[]) {
        this.validateControlPoints(points);
    }

    validateControlPoints(points: Vector2D[]) {
        if (points.length < 2) {
            throw new Error('At least 2 control points are required');
        }
    }

    getControlPoint(index: number): Vector2D {
        const baseIndex = ((index % this.points.length) + this.points.length) % this.points.length;
        return this.points[baseIndex];
    }

    // Return a deep copy of the points
    getControlPoints(): Vector2D[] {
        return this.points.map(point => [...point]);
    }

    getNumControlPoints() {
        return this.points.length;
    }

   
    getRelevantControlPoints(periodicKnots: PeriodicKnots, insertPosition: number, degree: number): Point[] {
        // Get the relevant knot information from the PeriodicKnot object
        const {knots, startIndex, endIndex} = 
            periodicKnots.getRelevantKnotsAndIndices(insertPosition, degree);
        
        // Use this information to determine which control points are needed
        // The number of control points needed is (endIndex - startIndex + 1)
        return this.extractControlPointsForRange(startIndex, endIndex);
    }
   
    

}