import { PeriodicKnots } from '../../../src/core/standard/periodicKnots';

describe('PeriodicKnots', () => {
    // Constructor tests
    describe('constructor', () => {
        describe('constructor', () => {
        it('should create a valid periodic knot vector when pattern range is less than period', () => {
            const pattern = [0, 0, 1, 2, 2];
            const period = 3;
            expect(() => new PeriodicKnots(pattern, period)).not.toThrow();
        });

        it('should throw error when pattern range equals period', () => {
            const pattern = [0, 1, 2, 3];
            const period = 3;
            expect(() => new PeriodicKnots(pattern, period))
                .toThrow('Pattern range (3) exceeds specified period (3)');
        });

        it('should throw error when pattern range exceeds period', () => {
            const pattern = [0, 1, 2, 4];
            const period = 3;
            expect(() => new PeriodicKnots(pattern, period))
                .toThrow('Pattern range (4) exceeds specified period (3)');
        });

        it('should throw error for non-decreasing pattern', () => {
            const pattern = [0, 2, 1];
            const period = 2;
            expect(() => new PeriodicKnots(pattern, period))
                .toThrow('Knot vector values must be in non-decreasing order');
        });

        it('should throw error when pattern range exceeds period', () => {
            const pattern = [0, 1, 2, 3];
            const period = 2;
            expect(() => new PeriodicKnots(pattern, period))
                .toThrow('Pattern range (3) exceeds specified period (2)');
        });
      });
    });

    
    // Static factory method tests
    describe('fromDistinctKnots', () => {
        it('should create correct periodic knots from distinct knots', () => {
            const values = {
                knots: [0, 1, 2],
                multiplicities: [2, 1, 2]
            };
            const period = 3;
            const periodicKnots = PeriodicKnots.fromDistinctKnots(values, period);
            expect(periodicKnots.getDistinctKnots()).toEqual([0, 1, 2]);
            expect(periodicKnots.getMultiplicities()).toEqual([2, 1, 2]);
        });
    });

    

    // Query methods tests
    describe('query methods', () => {
        let periodicKnots: PeriodicKnots;

        beforeEach(() => {
            periodicKnots = new PeriodicKnots([0, 0, 1, 2, 2], 3);
        });

        it('should get correct domain', () => {
            const domain = periodicKnots.getDomain();
            expect(domain).toEqual({ minValue: 0, maxValue: 3 });
        });

        it('should get distinct knots', () => {
            const distinctKnots = periodicKnots.getDistinctKnots();
            expect(distinctKnots).toEqual([0, 1, 2]);
        });

        it('should get correct multiplicities', () => {
            const multiplicities = periodicKnots.getMultiplicities();
            expect(multiplicities).toEqual([2, 1, 2]);
        });

        it('should get correct knot value for any index', () => {
            // Test within pattern
            expect(periodicKnots.getKnotValue(0)).toBe(0);
            expect(periodicKnots.getKnotValue(2)).toBe(1);
            
            // Test periodic behavior
            expect(periodicKnots.getKnotValue(5)).toBe(3);
            expect(periodicKnots.getKnotValue(10)).toBe(6);
            
            // Test negative indices
            expect(periodicKnots.getKnotValue(-5)).toBe(-3);
        });
    });

    
    
    // Modification methods tests
    describe('modification methods', () => {
        let periodicKnots: PeriodicKnots;

        beforeEach(() => {
            periodicKnots = new PeriodicKnots([0, 0, 1, 2, 2], 3);
        });

        it('should insert knot correctly', () => {
            const newKnots = periodicKnots.insertKnot(1.5);
            expect(newKnots.getDistinctKnots()).toContain(1.5);
        });

        it('should elevate degree correctly', () => {
            const elevated = periodicKnots.elevateDegree();
            const originalMults = periodicKnots.getMultiplicities();
            const newMults = elevated.getMultiplicities();
            
            // Check that each multiplicity increased by 1
            expect(newMults.length).toBe(originalMults.length);
            newMults.forEach((mult, i) => {
                expect(mult).toBe(originalMults[i] + 1);
            });
        });
    });

    
    // Utility methods tests
    describe('utility methods', () => {
        let periodicKnots: PeriodicKnots;

        beforeEach(() => {
            periodicKnots = new PeriodicKnots([0, 0, 1, 2, 2], 3);
        });

        it('should find correct span index', () => {
            expect(periodicKnots.findSpanIndex(0.5)).toBe(1);
            expect(periodicKnots.findSpanIndex(1.5)).toBe(2);
            expect(periodicKnots.findSpanIndex(2.5)).toBe(4);
            expect(periodicKnots.findSpanIndex(3.5)).toBe(6);
            expect(periodicKnots.findSpanIndex(-0.5)).toBe(-1);
        });

        it('should shift values inside domain correctly', () => {
            expect(periodicKnots.shiftInsideDomain(2.5)).toBe(2.5);
            expect(periodicKnots.shiftInsideDomain(-0.5)).toBe(2.5);
            expect(periodicKnots.shiftInsideDomain(4.5)).toBe(1.5);
        });
    });

    

    describe('unrollForKnotInsertion', () => {
        it('should correctly unroll knots for insertion', () => {
            const periodicKnots = new PeriodicKnots([0, 0, 1, 2, 2], 3);
            const degree = 3;
            const u = 1.5;
            
            const result = periodicKnots.unrollForKnotInsertion(degree, u);
            
            // Check that the unrolled knots span enough periods
            expect(result.knots).toBeDefined();
            expect(result.knotsToBeInserted).toBeDefined();
            
            // Verify the knots to be inserted are periodic
            const knotsToBeInserted = result.knotsToBeInserted;
            for (let i = 1; i < knotsToBeInserted.length; i++) {
                expect(knotsToBeInserted[i] - knotsToBeInserted[i-1])
                    .toBeCloseTo(periodicKnots.getDomain().maxValue - 
                                periodicKnots.getDomain().minValue);
            }
        });
    });
    


    
});
