import { extendPeriodicArray, Knots } from '../../../src/core/standard/knots';

describe('Knots', () => {
    // Test data
    const standardKnots = [0, 0, 0, 1, 2, 3, 4, 4, 4];
    const periodicKnots = [0, 1, 2, 3, 4, 5];
    const degree = 2;

    describe('Constructor', () => {
        test('should create non-periodic knot vector', () => {
            const knots = new Knots(standardKnots, degree);
            expect(knots.getDegree()).toBe(degree);
        });

        test('should throw error for invalid non-decreasing knots', () => {
            const invalidKnots = [0, 2, 1, 3, 4];
            expect(() => new Knots(invalidKnots, degree))
                .toThrow();
        });

        test('should throw error for periodic knots without numControlPoints', () => {
            expect(() => new Knots(periodicKnots, degree, true))
                .toThrow('Number of control points must be specified for periodic knot vector');
        });
    });

    describe('getValue', () => {
        
        test('should return correct knot value for non-periodic case', () => {
            const knots = new Knots(standardKnots, degree);
            expect(knots.getValue(3)).toBe(1);
        });
    
        test('should return correct periodic knot value', () => {
            // For degree = 2 and numControlPoints = 3, we need at least 6 knots
            const periodicKnots = [0, 1, 2, 3, 4, 5];  // 6 knots
            const numControlPoints = 3;
            const knots = new Knots(periodicKnots, degree, true, numControlPoints);
            
            // Test basic value retrieval within the original range
            expect(knots.getValue(0)).toBe(0);
            expect(knots.getValue(1)).toBe(1);
            expect(knots.getValue(2)).toBe(2);
    
            // Calculate the period (difference between the knot a numControlPoints position and first knot)
            const period = periodicKnots[numControlPoints] - periodicKnots[0];
    
            // Test periodic behavior by checking values beyond the original range
            // When we go beyond the knot vector length, it should wrap around
            // and add the period
            expect(knots.getValue(3)).toBe(0 + period);  // wraps to index 0 + period
            expect(knots.getValue(4)).toBe(1 + period);  // wraps to index 1 + period
            expect(knots.getValue(5)).toBe(2 + period);  // wraps to index 2 + period
    
            // Test negative indices (should wrap backwards and subtract period)
            expect(knots.getValue(-1)).toBe(2 - period);
            expect(knots.getValue(-2)).toBe(1 - period);
        });
        
    });
    
    
    // Add more comprehensive periodic tests
    describe('Periodic Knots Behavior', () => {
        const degree = 2;
        const numControlPoints = 3;
        const periodicKnots = [0, 1, 2, 3, 4, 5];
    
        test('should handle periodic value retrieval correctly', () => {
            const knots = new Knots(periodicKnots, degree, true, numControlPoints);
            const period = periodicKnots[numControlPoints] - periodicKnots[0]; // 3
    
            // Test values within one complete period
            for (let i = 0; i < periodicKnots.length; i++) {
                expect(knots.getValue(i)).toBe(periodicKnots[i]);
            }
    
            // Test values in the next period
            for (let i = 0; i < numControlPoints; i++) {
                const wrappedIndex = i % numControlPoints;
                expect(knots.getValue(i + numControlPoints))
                    .toBe(periodicKnots[wrappedIndex] + period);
            }
    
            // Test values in the previous period
            for (let i = 0; i < numControlPoints; i++) {
                const wrappedIndex = i % numControlPoints;
                expect(knots.getValue(i - numControlPoints))
                    .toBe(periodicKnots[wrappedIndex] - period);
            }
        });

    
        test('should maintain uniform spacing', () => {
            const knots = new Knots(periodicKnots, degree, true, numControlPoints);
            const spacing = 1; // Known spacing for our test data
    
            // Check spacing within original knot vector
            for (let i = 1; i < periodicKnots.length; i++) {
                expect(knots.getValue(i) - knots.getValue(i - 1)).toBe(spacing);
            }
    
            // Check spacing across period boundary
            expect(knots.getValue(periodicKnots.length) - knots.getValue(periodicKnots.length - 1))
                .toBe(spacing);
        });
    
        test('should handle multiple period shifts', () => {
            const knots = new Knots(periodicKnots, degree, true, numControlPoints);
            const period = periodicKnots[numControlPoints] - periodicKnots[0];
    
            // Test double period shift forward
            expect(knots.getValue(numControlPoints * 2)).toBe(0 + period * 2);
            expect(knots.getValue(numControlPoints * 2 + 1)).toBe(1 + period * 2);
    
            // Test double period shift backward
            expect(knots.getValue(-numControlPoints * 2)).toBe(0 - period * 2);
            expect(knots.getValue(-numControlPoints * 2 + 1)).toBe(1 - period * 2);
        });
        
    });
    
    
    
    describe('getSpan', () => {
        test('should find correct span in non-periodic case', () => {
            const knots = new Knots(standardKnots, degree);
            expect(knots.getSpan(1.5)).toBe(3);
        });

        

        test('should handle boundary cases in span finding', () => {
            const knots = new Knots(standardKnots, degree);
            expect(knots.getSpan(0)).toBe(0);
            expect(knots.getSpan(4)).toBe(standardKnots.length - 2);
        });
        
    });

    
    describe('insert', () => {
        test('should insert knot in non-periodic case', () => {
            const knots = new Knots(standardKnots, degree);
            const newKnots = knots.insert(1.5);
            expect(newKnots['values']).toContain(1.5);
            // Check if order is maintained
            expect(newKnots['values']).toEqual([...newKnots['values']].sort((a, b) => a - b));
        });

        test('should maintain periodicity when inserting knot', () => {
            const knots = new Knots(periodicKnots, degree, true, 3);
            const newKnots = knots.insert(1.5);
            // Check if inserted value exists
            expect(newKnots['values']).toContain(1.5);
            // Check if periodic copy exists
            expect(newKnots['values']).toContain(1.5 + newKnots['period']);
        });
    });

    

    describe('elevate', () => {
        test('should elevate degree for non-periodic case', () => {
            const knots = new Knots(standardKnots, degree);
            const elevated = knots.elevate();
            expect(elevated.getDegree()).toBe(degree + 1);
            // Check if multiplicities increased
            expect(elevated['values'].length).toBeGreaterThan(knots['values'].length);
        });

        
        
        test('should maintain periodicity when elevating degree', () => {
            const knots = new Knots(periodicKnots, degree, true, 3);
            const elevated = knots.elevate();
            expect(elevated.getDegree()).toBe(degree + 1);
            // Check if still periodic
            expect(elevated['_isPeriodic']).toBe(true);
        });

        

        test('should obtain this specific result when elevating degree for the uniform knot periodic case', () => {
            // Initialize a periodic knot vector with:
            // - degree 2 (quadratic)
            // - 3 periodic control points
            // - 5 wrapped control points
            // - 8 knots
            // - The domain is [2, 5)
            // - periodic flag set to true
            const periodicKnots = [0, 1, 2, 3, 4, 5, 6, 7];
            const degree = 2;
            const knots = new Knots(periodicKnots, degree, true, 3);
            expect(knots.getDegree()).toBe(2)
            expect(knots.getPeriod()).toBe(3)
            expect(knots.getDomain()).toEqual({"maxValue": 5, "minValue": 2})
            expect(knots.getNumberOfWrappedControlPoints()).toBe(5)
            expect(knots.getNumberOfControlPoints()).toBe(3)

            // Elevate degree from 2 to 3 (cubic)
            const elevated1 = knots.elevate();
            
            // After first elevation:
            // - degree is now 3
            // - number of periodic control points increases to 6
            // - number of wrapped control points increases to 9
            // - number of knots increases to 13
            // - each interior knot has multiplicity 2 for C¹ continuity
            // - domain remains [2,5)
            expect(elevated1['values']).toEqual([0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6]);
            expect(elevated1.getDegree()).toBe(3)
            expect(elevated1.getPeriod()).toBe(3)
            expect(elevated1.getDomain()).toEqual({"maxValue": 5, "minValue": 2})
            expect(elevated1.getNumberOfWrappedControlPoints()).toBe(9)
            expect(elevated1.getNumberOfControlPoints()).toBe(6)

            // Elevate again from degree 3 to 4
            const elevated2 = elevated1.elevate()
            
            // After second elevation:
            // - degree is now 4
            // - number of periodic control points increases to 9
            // - number of wrapped control points increases to 13
            // - number of knots increases to 18
            // - interior knots have multiplicity 3 for C¹ continuity
            // - domain remains [2,5)
            expect(elevated2['values']).toEqual([0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6]);
            expect(elevated2.getDegree()).toBe(4)
            expect(elevated2.getPeriod()).toBe(3)
            expect(elevated2.getDomain()).toEqual({"maxValue": 5, "minValue": 2})
            expect(elevated2.getNumberOfWrappedControlPoints()).toBe(13)
            expect(elevated2.getNumberOfControlPoints()).toBe(9)
        
        })

        test('should obtain this specific results when elevating degree for the non uniform knot periodic case', () => {
            // Initialize a periodic knot vector with:
            // - degree 2 (quadratic)
            // - 4 periodic control points
            // - 6 wrapped control points
            // - 9 knots
            // - The domain is [2, 5)
            // - periodic flag set to true
            const periodicKnots = [-0.2, 1, 2, 2.8, 2.8, 4, 5, 5.8, 5.8];
            const degree = 2;
            const knots = new Knots(periodicKnots, degree, true, 4);
            expect(knots.getDegree()).toBe(2)
            expect(knots.getPeriod()).toBe(3)
            expect(knots.getDomain()).toEqual({"maxValue": 5, "minValue": 2})
            expect(knots.getNumberOfWrappedControlPoints()).toBe(6)
            expect(knots.getNumberOfControlPoints()).toBe(4)

            // Elevate degree from 2 to 3 (cubic)
            const elevated1 = knots.elevate();
            
            // After first elevation:
            // - degree is now 3
            // - number of periodic control points increases to 7
            // - number of wrapped control points increases to 10
            // - number of knots increases to 14
            // - domain remains [2,5)

            const expected = [-0.2, 1, 1, 2, 2, 2.8, 2.8, 2.8, 4, 4, 5, 5, 5.8, 5.8];
            const result = elevated1['values']
            expect(result.length).toBe(expected.length);
            result.forEach((value, index) => {
                expect(value).toBeCloseTo(expected[index], 0.0001);
            });
            expect(elevated1.getDegree()).toBe(3)
            expect(elevated1.getPeriod()).toBe(3)
            expect(elevated1.getDomain()).toEqual({"maxValue": 5, "minValue": 2})
            expect(elevated1.getNumberOfWrappedControlPoints()).toBe(10)
            expect(elevated1.getNumberOfControlPoints()).toBe(7)

            // Elevate again from degree 3 to 4
            const elevated2 = elevated1.elevate()
            
            // After second elevation:
            // - degree is now 4
            // - number of periodic control points increases to 10
            // - number of wrapped control points increases to 14
            // - number of knots increases to 19
            // - domain remains [2,5)

            const expected2 = [-0.2, 1, 1, 1, 2, 2, 2, 2.8, 2.8, 2.8, 2.8, 4, 4, 4, 5, 5, 5, 5.8, 5.8];
            const result2 = elevated2['values']
            expect(result2.length).toBe(expected2.length);
            result2.forEach((value, index) => {
                expect(value).toBeCloseTo(expected2[index], 0.0001);
            });
            expect(elevated2.getDegree()).toBe(4)
            expect(elevated2.getPeriod()).toBe(3)
            expect(elevated2.getDomain()).toEqual({"maxValue": 5, "minValue": 2})
            expect(elevated2.getNumberOfWrappedControlPoints()).toBe(14)
            expect(elevated2.getNumberOfControlPoints()).toBe(10)
        
        })
            
        
    });

    describe('Edge cases', () => {
        test('should handle empty knot vector', () => {
            expect(() => new Knots([], degree))
                .toThrow();
        });

        test('should handle single knot', () => {
            expect(() => new Knots([1], degree))
                .toThrow();
        });
    });
    
    
});



describe('Periodic Knots', () => {
    const degree = 2;
    const numControlPoints = 3;
    const minKnots = degree + numControlPoints + 1; // 6 knots minimum
    const periodicKnots = [0, 1, 2, 3, 4, 5];

    test('should create valid periodic knot vector', () => {
        const knots = new Knots(periodicKnots, degree, true, numControlPoints);
        expect(knots.getDegree()).toBe(degree);
    });

    test('should throw error if not enough knots provided', () => {
        const insufficientKnots = [0, 1, 2, 3, 4]; // Only 5 knots
        expect(() => new Knots(insufficientKnots, degree, true, numControlPoints))
            .toThrow(`Periodic knot vector must have at least ${minKnots} knots`);
    });

    test('should correctly calculate period', () => {
        const knots = new Knots(periodicKnots, degree, true, numControlPoints);
        const period = periodicKnots[numControlPoints] - periodicKnots[0];
        expect(knots['period']).toBe(period);
    });

    test('should handle periodic value retrieval across multiple periods', () => {
        const knots = new Knots(periodicKnots, degree, true, numControlPoints);
        const period = periodicKnots[numControlPoints] - periodicKnots[0];
        
        // Test positive periods
        expect(knots.getValue(numControlPoints)).toBe(periodicKnots[0] + period);
        expect(knots.getValue(numControlPoints + 1)).toBe(periodicKnots[1] + period);
        
        // Test negative periods
        expect(knots.getValue(-1)).toBe(periodicKnots[numControlPoints - 1] - period);
        expect(knots.getValue(-2)).toBe(periodicKnots[numControlPoints - 2] - period);
    });

    test('should maintain uniform spacing in periodic knots', () => {
        const knots = new Knots(periodicKnots, degree, true, numControlPoints);
        const spacing = periodicKnots[1] - periodicKnots[0];
        
        for (let i = 1; i < periodicKnots.length; i++) {
            expect(periodicKnots[i] - periodicKnots[i-1]).toBe(spacing);
        }
    });
   
});

describe('extendPeriodicArray', () => {
    test('should handle empty array', () => {
        expect(extendPeriodicArray([], 3, 2, 2)).toEqual([]);
    });

    test('should handle single value array', () => {
        expect(extendPeriodicArray([2], 3, 2, 2))
            .toEqual([-4, -1, 2, 5, 8]);
    });

    test('should handle uniform spacing', () => {
        const result = extendPeriodicArray([2, 3, 4], 3, 2, 2);
        expect(result).toEqual([0, 1, 2, 3, 4, 5, 6]);
    });

    test('should handle non-uniform spacing', () => {
        const result = extendPeriodicArray([2, 2.8, 4], 3, 2, 2);
        const expected = [-0.2, 1, 2, 2.8, 4, 5, 5.8];
        
        // Compare arrays with floating point tolerance
        expect(result.length).toBe(expected.length);
        result.forEach((value, index) => {
            expect(value).toBeCloseTo(expected[index], 10);
        });
    });

    test('should handle repeated values', () => {
        const result = extendPeriodicArray([2, 2.8, 2.8, 4], 3, 2, 2);
        const expected = [-0.2, 1, 2, 2.8, 2.8, 4, 5, 5.8];
        expect(result.length).toBe(expected.length);
        result.forEach((value, index) => {
            expect(value).toBeCloseTo(expected[index], 10);
        });
    });

    test('should handle multiple periods', () => {
        const result = extendPeriodicArray([2, 3, 4], 3, 6, 6);
        expect(result).toEqual([
            -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
        ]);
    });

    test('should preserve periodicity', () => {
        const result = extendPeriodicArray([2, 3, 4], 3, 3, 3);
        
        // Check periodicity before
        expect(result[0] + 3).toBeCloseTo(result[3]);
        expect(result[1] + 3).toBeCloseTo(result[4]);
        expect(result[2] + 3).toBeCloseTo(result[5]);

        // Check periodicity after
        expect(result[3] + 3).toBeCloseTo(result[6]);
        expect(result[4] + 3).toBeCloseTo(result[7]);
        expect(result[5] + 3).toBeCloseTo(result[8]);
    });

    

    test('should handle zero period', () => {
        const result = extendPeriodicArray([1, 2, 3], 0, 2, 2);
        expect(result).toEqual([2, 3, 1, 2, 3, 1, 2]);
    });

    test('should handle negative period', () => {
        const result = extendPeriodicArray([1, 2, 3], -3, 2, 2);
        expect(result).toEqual([5, 6, 1, 2, 3, -2, -1]);
    });

    test('should handle zero extensions', () => {
        const original = [1, 2, 3];
        const result = extendPeriodicArray(original, 3, 0, 0);
        expect(result).toEqual(original);
    });

    test('should not modify original array', () => {
        const original = [1, 2, 3];
        const originalCopy = [...original];
        extendPeriodicArray(original, 3, 2, 2);
        expect(original).toEqual(originalCopy);
    });
    
});

 


