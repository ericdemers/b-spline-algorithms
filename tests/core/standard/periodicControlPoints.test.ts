import { PeriodicControlPoints } from '../../../src/core/standard/periodicControlPoints';
import { Vector2D } from '../../../src/core/standard/vector-space';

describe('PeriodicControlPoints', () => {
    let periodicPoints: PeriodicControlPoints;
    
    describe('basic functionality', () => {
        beforeEach(() => {
            // Create a simple set of control points forming a triangle
            const points = [ 
                [0, 0], // P0
                [1, 1], // P1
                [2, 0] // P2  
            ] as Vector2D[];
            periodicPoints = new PeriodicControlPoints(points);
        });

        it('should return correct control point for positive indices within range', () => {
            expect(periodicPoints.getControlPoint(0)).toEqual([0, 0]);
            expect(periodicPoints.getControlPoint(1)).toEqual([1, 1]);
            expect(periodicPoints.getControlPoint(2)).toEqual([2, 0]);
        });

        it('should handle periodic wrapping for indices beyond array length', () => {
            expect(periodicPoints.getControlPoint(3)).toEqual([0, 0]);
            expect(periodicPoints.getControlPoint(4)).toEqual([1, 1]);
            expect(periodicPoints.getControlPoint(5)).toEqual([2, 0]);
        });

        it('should handle negative indices with periodic wrapping', () => {
            expect(periodicPoints.getControlPoint(-1)).toEqual([2, 0]);
            expect(periodicPoints.getControlPoint(-2)).toEqual([1, 1]);
            expect(periodicPoints.getControlPoint(-3)).toEqual([0, 0]);
        });

        it('should handle large positive indices', () => {
            expect(periodicPoints.getControlPoint(99)).toEqual(periodicPoints.getControlPoint(99 % 3));
            expect(periodicPoints.getControlPoint(100)).toEqual(periodicPoints.getControlPoint(100 % 3));
            expect(periodicPoints.getControlPoint(101)).toEqual(periodicPoints.getControlPoint(101 % 3));
        });

        it('should handle large negative indices', () => {
            expect(periodicPoints.getControlPoint(-99)).toEqual(periodicPoints.getControlPoint((-99 % 3 + 3) % 3));
            expect(periodicPoints.getControlPoint(-100)).toEqual(periodicPoints.getControlPoint((-100 % 3 + 3) % 3));
            expect(periodicPoints.getControlPoint(-101)).toEqual(periodicPoints.getControlPoint((-101 % 3 + 3) % 3));
        });
    });


    describe('complex patterns', () => {
        beforeEach(() => {
            // Create a more complex set of control points
            const points = [
                [0, 0],  // P0
                [1, 2],   // P1
                [2, -1],   // P2
                [3, 1],   // P3
                [4, 0]   // P4
            ] as Vector2D[];
            periodicPoints = new PeriodicControlPoints(points);
        });

        it('should maintain periodicity over multiple cycles', () => {
            const basePoint = periodicPoints.getControlPoint(2); // P2
            expect(periodicPoints.getControlPoint(2 + 5)).toEqual(basePoint);  // One cycle
            expect(periodicPoints.getControlPoint(2 + 10)).toEqual(basePoint); // Two cycles
            expect(periodicPoints.getControlPoint(2 - 5)).toEqual(basePoint);  // One cycle backward
        });

        it('should preserve point coordinates exactly', () => {
            expect(periodicPoints.getControlPoint(3)).toEqual([3, 1]);
            expect(periodicPoints.getControlPoint(8)).toEqual([3, 1]);
            expect(periodicPoints.getControlPoint(-2)).toEqual([3, 1]);
        });
    });
});
