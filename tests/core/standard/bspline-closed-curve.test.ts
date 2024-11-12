import { Vector2DSpace } from '../../../src/core/standard/vector-space';
import { ClosedBSplineCurve } from '../../../src/core/standard/bspline-closed-curve';

describe('ClosedBSplineCurve', () => {
    let vectorSpace: Vector2DSpace;
    let controlPoints: [number, number][];
    let curve: ClosedBSplineCurve<number, [number, number]>;
    
    beforeEach(() => {
        vectorSpace = new Vector2DSpace();
        // Square-like control points
        controlPoints = [
            [1, 1],
            [1, -1],
            [-1, -1],
            [-1, 1]
        ];
        curve = new ClosedBSplineCurve(vectorSpace, controlPoints, 3); // cubic
    });

    describe('Constructor', () => {
        it('should create a valid closed B-spline curve', () => {
            expect(curve).toBeDefined();
            expect(curve.getDegree()).toBe(3);
            expect(curve.getControlPoints()).toHaveLength(4);
        });

        it('should throw error for degree less than 1', () => {
            expect(() => {
                new ClosedBSplineCurve(vectorSpace, controlPoints, 0);
            }).toThrow('Degree must be positive');
        });

        it('should throw error when not enough control points', () => {
            expect(() => {
                new ClosedBSplineCurve(vectorSpace, [[0, 0], [1, 1]], 3);
            }).toThrow('Number of control points must be greater than degree');
        });
    });

    /*
    describe('Evaluation', () => {
        it('should evaluate at t = 0', () => {
            const point = curve.evaluate(0);
            // First point should be influenced by last and first control points
            expect(point[0]).toBeCloseTo(0.6667, 4);
            expect(point[1]).toBeCloseTo(1.0000, 4);
        });

        it('should evaluate at t = 0.5', () => {
            const point = curve.evaluate(0.5);
            // Middle point should be influenced by middle control points
            expect(point[0]).toBeCloseTo(-0.6667, 4);
            expect(point[1]).toBeCloseTo(-1.0000, 4);
        });

        it('should evaluate at t = 1', () => {
            const point = curve.evaluate(1);
            // Should return to starting point
            expect(point[0]).toBeCloseTo(0.6667, 4);
            expect(point[1]).toBeCloseTo(1.0000, 4);
        });

        it('should handle parameter wrapping', () => {
            const point1 = curve.evaluate(0.25);
            const point2 = curve.evaluate(1.25);
            // Same parameter modulo 1 should give same point
            expect(point1[0]).toBeCloseTo(point2[0], 4);
            expect(point1[1]).toBeCloseTo(point2[1], 4);
        });
    });

    describe('Derivative', () => {
        it('should compute first derivative', () => {
            const derivative = curve.getDerivative();
            expect(derivative.getDegree()).toBe(2);
            
            // Test derivative at t = 0
            const derPoint = derivative.evaluate(0);
            expect(derPoint[0]).toBeDefined();
            expect(derPoint[1]).toBeDefined();
        });

        it('should throw error for derivative of linear curve', () => {
            const linearCurve = new ClosedBSplineCurve(vectorSpace, controlPoints, 1);
            expect(() => {
                linearCurve.getDerivative();
            }).toThrow('Degree must be at least 2 for derivative of closed curve');
        });
    });

    describe('Knot Insertion', () => {
        it('should insert knot at t = 0.5', () => {
            const refinedCurve = curve.insertKnot(0.5);
            
            // Check curve shape preservation
            const testParams = [0, 0.25, 0.5, 0.75, 1];
            testParams.forEach(t => {
                const originalPoint = curve.evaluate(t);
                const newPoint = refinedCurve.evaluate(t);
                expect(newPoint[0]).toBeCloseTo(originalPoint[0], 4);
                expect(newPoint[1]).toBeCloseTo(originalPoint[1], 4);
            });
        });

        it('should maintain closure after knot insertion', () => {
            const refinedCurve = curve.insertKnot(0.5);
            const startPoint = refinedCurve.evaluate(0);
            const endPoint = refinedCurve.evaluate(1);
            
            expect(endPoint[0]).toBeCloseTo(startPoint[0], 4);
            expect(endPoint[1]).toBeCloseTo(startPoint[1], 4);
        });

        it('should handle multiple knot insertions', () => {
            const refinedCurve = curve
                .insertKnot(0.3)
                .insertKnot(0.6);
            
            // Check curve shape preservation
            const testParams = [0, 0.3, 0.6, 1];
            testParams.forEach(t => {
                const originalPoint = curve.evaluate(t);
                const newPoint = refinedCurve.evaluate(t);
                expect(newPoint[0]).toBeCloseTo(originalPoint[0], 4);
                expect(newPoint[1]).toBeCloseTo(originalPoint[1], 4);
            });
        });
    });

    describe('Geometric Properties', () => {
        it('should be periodic', () => {
            const point1 = curve.evaluate(0);
            const point2 = curve.evaluate(1);
            const point3 = curve.evaluate(2);
            
            expect(point1[0]).toBeCloseTo(point2[0], 4);
            expect(point1[1]).toBeCloseTo(point2[1], 4);
            expect(point1[0]).toBeCloseTo(point3[0], 4);
            expect(point1[1]).toBeCloseTo(point3[1], 4);
        });

        it('should preserve shape under parameter transformation', () => {
            const testParams = [0.1, 0.4, 0.7];
            testParams.forEach(t => {
                const point1 = curve.evaluate(t);
                const point2 = curve.evaluate(t + 1);
                const point3 = curve.evaluate(t - 1);
                
                expect(point1[0]).toBeCloseTo(point2[0], 4);
                expect(point1[1]).toBeCloseTo(point2[1], 4);
                expect(point1[0]).toBeCloseTo(point3[0], 4);
                expect(point1[1]).toBeCloseTo(point3[1], 4);
            });
        });

        it('should maintain continuity at joints', () => {
            const h = 0.001; // Small step for derivative approximation
            const t = 1.0; // Test at curve closure
            
            // Approximate derivative from left
            const leftPoint = curve.evaluate(t - h);
            const midPoint = curve.evaluate(t);
            const leftDerivative = vectorSpace.scale(
                1/h, 
                vectorSpace.subtract(midPoint, leftPoint)
            );
            
            // Approximate derivative from right
            const rightPoint = curve.evaluate(t + h);
            const rightDerivative = vectorSpace.scale(
                1/h, 
                vectorSpace.subtract(rightPoint, midPoint)
            );
            
            // Compare derivatives
            expect(rightDerivative[0]).toBeCloseTo(leftDerivative[0], 2);
            expect(rightDerivative[1]).toBeCloseTo(leftDerivative[1], 2);
        });
    }); */
});
