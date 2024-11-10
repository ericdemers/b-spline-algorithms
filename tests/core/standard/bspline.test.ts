//import { describe, it, expect, beforeEach } from 'vitest';
import { BSpline } from '../../../src/core/standard/bspline';
import { Complex, ComplexVector2DSpace, Vector2D, Vector2DSpace } from '../../../src/core/standard/vector-space';
import { ControlPolygonCurve } from '../../../src/core/standard/control-net';
import { Knots } from '../../../src/core/standard/knot-structure';

describe('BSpline', () => {
    let vectorSpace: Vector2DSpace;
    let controlPoints: Vector2D[];
    let controlNet: ControlPolygonCurve<Vector2D>;
    let knots: Knots;
    let bspline: BSpline<number, Vector2D>;

    beforeEach(() => {
        // Setup common test fixtures
        vectorSpace = new Vector2DSpace();
        controlPoints = [
            [0, 0],
            [1, 1],
            [2, 0],
            [3, 1]
        ];
        controlNet = new ControlPolygonCurve(controlPoints);
        knots = new Knots([0, 0, 0, 0, 1, 1, 1, 1]); // Cubic B-spline
        bspline = new BSpline(vectorSpace, controlNet, knots, [3]);
    });

    describe('Constructor', () => {
        it('should create a valid B-spline instance', () => {
            expect(bspline).toBeDefined();
        });

        it('should throw error for invalid degree', () => {
            expect(() => {
                new BSpline(vectorSpace, controlNet, knots, [4]); // Degree too high
            }).toThrow();
        });

        it('should throw error for mismatched dimensions', () => {
            expect(() => {
                new BSpline(vectorSpace, controlNet, knots, [3, 2]); // Wrong dimension
            }).toThrow();
        });
    });

    describe('evaluate', () => {
        it('should evaluate at start point (t=0)', () => {
            const point = bspline.evaluate([0]);
            expect(point).toEqual([0, 0]);
        });

        it('should evaluate at end point (t=1)', () => {
            const point = bspline.evaluate([1]);
            expect(point).toEqual([3, 1]);
        });

        it('should evaluate at middle point (t=0.5)', () => {
            const point = bspline.evaluate([0.5]);
            // Expected value should be calculated based on known B-spline properties
            expect(point[0]).toBeCloseTo(1.5, 5);
            expect(point[1]).toBeCloseTo(0.5, 5);
        });

        it('should throw error for parameter out of range', () => {
            expect(() => {
                bspline.evaluate([1.5]);
            }).toThrow();
        });

        it('should throw error for wrong number of parameters', () => {
            expect(() => {
                bspline.evaluate([0.5, 0.5]);
            }).toThrow();
        });
    });

    describe('evaluateRecursive', () => {
        it('should handle single parameter evaluation', () => {
            const point = bspline['evaluateRecursive']([0.5], 0);
            expect(point[0]).toBeCloseTo(1.5, 5);
            expect(point[1]).toBeCloseTo(0.5, 5);
        });

        it('should maintain partition of unity property', () => {
            // Test that basis functions sum to 1
            const t = 0.5;
            const span = bspline['findSpanIndex'](t, 3, knots.getKnotSequence(0));
            const basisFunctions = bspline['computeBasisFunctions'](t, span, 3, knots.getKnotSequence(0));
            const sum = basisFunctions.reduce((a, b) => a + b, 0);
            expect(sum).toBeCloseTo(1, 10);
        });
    });

    describe('Complex vector space', () => {
        let complexControlPoints: [Complex, Complex][];
        let complexBspline: BSpline<Complex, [Complex, Complex]>;

        beforeEach(() => {
            const complexVectorSpace = new ComplexVector2DSpace();
            complexControlPoints = [
                [[0, 0], [0, 0]],
                [[1, 1], [0, 1]],
                [[2, -1], [1, 0]],
                [[3, 0], [1, -1]]
            ];
            const complexControlNet = new ControlPolygonCurve(complexControlPoints);
            complexBspline = new BSpline(
                complexVectorSpace,
                complexControlNet,
                knots,
                [3]
            );
        });

        it('should evaluate complex B-spline at t=0.5', () => {
            const point = complexBspline.evaluate([0.5]);
            // Expected values should be calculated based on complex arithmetic
            expect(point[0][0]).toBeCloseTo(1.5, 5); // Real part
            expect(point[0][1]).toBeCloseTo(0, 5); // Imaginary part
            expect(point[1][0]).toBeCloseTo(0.5, 5);
            expect(point[1][1]).toBeCloseTo(0.25, 5);
        });
    });

    describe('Basis Functions', () => {
        it('should compute correct basis functions', () => {
            const t = 0.5;
            const span = bspline['findSpanIndex'](t, 3, knots.getKnotSequence(0));
            const basisFunctions = bspline['computeBasisFunctions'](t, span, 3, knots.getKnotSequence(0));
            
            // Test properties of basis functions
            expect(basisFunctions.length).toBe(4); // Degree + 1
            basisFunctions.forEach(value => {
                expect(value).toBeGreaterThanOrEqual(0); // Non-negative
                expect(value).toBeLessThanOrEqual(1); // Less than or equal to 1
            });
        });

        it('should satisfy partition of unity', () => {
            const testPoints = [0, 0.25, 0.5, 0.75, 1];
            testPoints.forEach(t => {
                const span = bspline['findSpanIndex'](t, 3, knots.getKnotSequence(0));
                const basisFunctions = bspline['computeBasisFunctions'](t, span, 3, knots.getKnotSequence(0));
                const sum = basisFunctions.reduce((a, b) => a + b, 0);
                expect(sum).toBeCloseTo(1, 10);
            });
        });
    });

    describe('Edge Cases', () => {
        it('should handle repeated knots', () => {
            // For a cubic B-spline (degree 3) with 4 control points
            // We need n + p + 1 knots where:
            // n = number of control points - 1 = 3
            // p = degree = 3
            // So we need 3 + 3 + 1 = 7 + 1 = 8 knots total
            
            const repeatedKnots = new Knots([0, 0, 0, 0, 1, 1, 1, 1]); // 8 knots total
            const bsplineRepeated = new BSpline(vectorSpace, controlNet, repeatedKnots, [3]);
            const point = bsplineRepeated.evaluate([0]);
            expect(point).toEqual(controlPoints[0]);
        });

        it('should handle minimal number of control points', () => {
            const minimalControlPoints: Vector2D[] = [[0, 0], [1, 1], [2, 0], [3, 1]];
            const minimalControlNet = new ControlPolygonCurve(minimalControlPoints);
            const minimalBspline = new BSpline(vectorSpace, minimalControlNet, knots, [3]);
            expect(() => {
                minimalBspline.evaluate([0.5]);
            }).not.toThrow();
        });
    });

    describe('Performance', () => {
        it('should evaluate quickly for many points', () => {
            const startTime = performance.now();
            const numPoints = 1000;
            
            for (let i = 0; i < numPoints; i++) {
                bspline.evaluate([i / (numPoints - 1)]);
            }
            
            const endTime = performance.now();
            const timePerEval = (endTime - startTime) / numPoints;
            
            expect(timePerEval).toBeLessThan(1); // Less than 1ms per evaluation
        });
    });
});
