interface IKnotVector {
    getKnots(): number[];
    findSpan(u: number): number;
    getDegree(): number;
}

interface IControlPoints {
    getPoints(): Point[];
    getPoint(index: number): Point;
    count(): number;
}

interface IBSplineCurve {
    evaluate(u: number): Point;
    derivative(u: number, order: number): Point;
    insertKnot(u: number): void;
}



class BSpline {
    private knots: IKnotVector;
    private controlPoints: IControlPoints;
    private degree: number;

    constructor(knots: IKnotVector, controlPoints: IControlPoints, degree: number) {
        // Validate relationship constraints
        this.validateConfiguration(knots, controlPoints, degree);
        
        this.knots = knots;
        this.controlPoints = controlPoints;
        this.degree = degree;
    }

    protected validateConfiguration(knots: IKnotVector, controlPoints: IControlPoints, degree: number) {
        // Relationship constraints:
        // m = n + p + 1
        // where m = number of knots
        //       n = number of control points
        //       p = degree
        
        const expectedKnots = controlPoints.count() + degree + 1;
        if (knots.getKnots().length !== expectedKnots) {
            throw new Error('Invalid knot/control point configuration');
        }
    }

    evaluate(u: number): Point {
        // Evaluate the B-spline at parameter u
        // Uses both knots and control points
    }
}

// Specific implementations
class PeriodicBSpline extends BSpline {
    constructor(periodicKnots: PeriodicKnotVector, periodicControlPoints: PeriodicControlPoints, degree: number) {
        super(periodicKnots, periodicControlPoints, degree);
    }

    insertKnot(u: number): void {
        // 1. Get relevant portion of curve
        const relevantKnots = this.knots.getRelevantKnotsAndIndices(u, this.degree);
        const relevantControlPoints = this.controlPoints.getRelevantControlPoints(u, this.degree);

        // 2. Perform insertion
        // 3. Update periodic structure
    }
}

class PeriodicKnotVector implements IKnotVector {
    private knotPattern: number[];
    
    getRelevantKnotsAndIndices(u: number, degree: number) {
        // Implementation from previous discussion
    }
}

class PeriodicControlPoints implements IControlPoints {
    private basePoints: Point[];
    
    getRelevantControlPoints(u: number, degree: number) {
        // Get relevant control points based on parameter value
    }
}
