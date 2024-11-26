import { Vector } from "./vector-space";
import { adaptParameter } from '../utils/decorators';

/**
 * Control net structure for B-spline definition
 * Represents the control points that define the B-spline shape
 */
export interface ControlNet<V extends Vector> {
    /** Get control point at given parametric indices */
    getPoint(indices: number[]): V;
    
    /** Get the dimension of the control net */
    getDimension(): number;
    
    /** Get size in each parametric direction */
    getSize(direction: number): number;
}


export class ControlPolygon<V extends Vector> implements ControlNet<V> {
    
    constructor(private readonly controlPoints: ReadonlyArray<V>,
                private readonly weights?: ReadonlyArray<V> ) {}

    @adaptParameter()
    getPoint(position: number | number[]): V {
        const index = position as number[]
        return this.controlPoints[index[0]]
    }

    getDimension() {
        return 1
    }

    getSize(direction: number = 0) {
        return this.controlPoints.length
    }
}

