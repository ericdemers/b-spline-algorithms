// 1. Mathematical Correctness
/*
interface ITest {
    evaluate(position: number): number;
    evaluate(position: number[]): number;
}

class BaseTest implements ITest {
    
    evaluate(position: number): number;
    evaluate(position: number[]): number;
    evaluate(position: number | number[]): number {
        if (Array.isArray(position)) {
            return position[0];
        }
        return position;
    }
}

class Test extends BaseTest {
    evaluate(position: number): number;
    evaluate(position: number[]): number;
    evaluate(position: number | number[]): number {
        if (Array.isArray(position)) {
            return super.evaluate(position);
        }
        return super.evaluate([position]);
    }
}
*/


function adaptParameter() {
    return function (target: any, propertyKey: string, descriptor: PropertyDescriptor) {
        const originalMethod = descriptor.value;
        descriptor.value = function (position: number | number[]) {
            const input = Array.isArray(position) ? position : [position];
            return originalMethod.call(this, input);
        };
    };
}
    

/*
function adaptParameter() {
    return function (
        target: Object,
        propertyKey: string | symbol,
        descriptor: TypedPropertyDescriptor<any>
    ): TypedPropertyDescriptor<any> | void {
        const originalMethod = descriptor.value;
        
        descriptor.value = function(position: number | number[]) {
            const input = Array.isArray(position) ? position : [position];
            return originalMethod.call(this, input);
        };

        return descriptor;
    };
}
    */

interface ITest {
    evaluate(position: number[]): number;
}

class BaseTest implements ITest {
    evaluate(position: number[]): number {
        return position[0];
    }
}

class Test extends BaseTest {
    @adaptParameter()
    evaluate(position: number | number[]): number {
        return super.evaluate(position as number[]);
    }
}

const test = new Test()
console.log(test.evaluate(1));


interface ITest2 {
    getKnotSequence(direction: number): number[];
}

class BaseTest2 implements ITest2 {
    getKnotSequence(direction: number = 0): number[] {
        return [1, 2, 3];
    }
}

class Test2 extends BaseTest2 {
    getKnotSequence(direction: number = 0): number[] {
        return super.getKnotSequence(direction);
    }
}

const test2 = new Test2()
console.log(test2.getKnotSequence());



// Start with abstract base that captures the essence of B-splines
export class BSpline<K extends Scalar, V extends Vector> {
    constructor(
        protected readonly vectorSpace: VectorSpace<K, V>,
        protected readonly controlNet: ControlNet<V>,
        protected readonly knots: KnotStructure,
        protected readonly degrees: ReadonlyArray<number>
    ) {
        // Core B-spline properties independent of dimension
    }
}

// Then specialize for curves with simplified interface
export class BSplineCurve<K extends Scalar, V extends Vector> extends BSpline<K, V> {
    constructor(
        vectorSpace: VectorSpace<K, V>,
        controlPoints: ReadonlyArray<V>,
        degree: number,
        knots?: ReadonlyArray<number>
    ) {
        // Convert to general form
        super(
            vectorSpace,
            new ControlPointArray(controlPoints),
            new UnivariateKnots(knots ?? generateUniformKnots(controlPoints.length, degree)),
            [degree]
        );
    }

    // Simplified curve-specific methods
    evaluate(t: number): V {
        return super.evaluate([t]);
    }

    derivative(order: number = 1): BSplineCurve<K, V> {
        // Curve-specific implementation using general machinery
    }
}


// 2. Algorithm Implementation

// Generic blossom implementation
class BlossomContext<K extends Scalar, V extends Vector> 
    readonly vectorSpace: VectorSpace<K, V>;
    readonly dimension: number;
    
    evaluate(parameters: number[]): V;
}

// Specialized for curves
class CurveBlossomContext<K extends Scalar, V extends Vector> 
    implements BlossomContext<K, V> {
    
    readonly dimension = 1;

    constructor(
        readonly vectorSpace: VectorSpace<K, V>,
        private readonly controlPoints: ReadonlyArray<V>,
        private readonly knots: ReadonlyArray<number>,
        private readonly degree: number
    ) {}

    evaluate(parameters: number[]): V {
        // Simplified implementation for curves
    }
}

// Can be extended to surfaces while reusing core logic
class SurfaceBlossomContext<K extends Scalar, V extends Vector>
    implements BlossomContext<K, V> {
    
    readonly dimension = 2;
    
    constructor(
        readonly vectorSpace: VectorSpace<K, V>,
        private readonly controlNet: ControlGrid<V>,
        private readonly knots: [ReadonlyArray<number>, ReadonlyArray<number>],
        private readonly degrees: [number, number]
    ) {}

    evaluate(parameters: number[]): V {
        // Reuse curve logic in each direction
    }
}

// 3 Type System Benefits

// Generic interfaces that work for all dimensions
class EvaluationAlgorithm<K extends Scalar, V extends Vector> 
    evaluate(
        bspline: BSpline<K, V>,
        parameters: ReadonlyArray<number>
    ): V;
}

// Specialized implementations with better ergonomics
class CurveEvaluation<K extends Scalar, V extends Vector> 
    implements EvaluationAlgorithm<K, V> {
    
    evaluate(curve: BSplineCurve<K, V>, t: number): V;
    evaluate(bspline: BSpline<K, V>, parameters: ReadonlyArray<number>): V;
    evaluate(input: BSpline<K, V> | BSplineCurve<K, V>, param: number | number[]): V {
        // Type-safe implementation that works with both interfaces
    }
}

// 4. Factory Methods

class BSplineFactory {
    // Generic factories
    static createFromNet<K extends Scalar, V extends Vector>(
        controlNet: ControlNet<V>,
        degrees: number[],
        vectorSpace: VectorSpace<K, V>
    ): BSpline<K, V>;

    // Specialized factories with simpler interfaces
    static createCurve2D(
        controlPoints: ReadonlyArray<[number, number]>,
        degree: number
    ): BSplineCurve<number, [number, number]> {
        return new BSplineCurve(
            vector2DSpace,
            controlPoints,
            degree
        );
    }
}


// Generic usage
const bspline = new BSpline(
    vector3DSpace,
    controlNet,
    knotStructure,
    [3, 2]  // degrees for each direction
);

// Specialized curve usage
const curve = new BSplineCurve(
    vector2DSpace,
    [[0,0], [1,1], [2,0]],
    2  // just one degree
);

// Factory usage
const curve2D = BSplineFactory.createCurve2D(
    [[0,0], [1,1], [2,0]],
    2
);





