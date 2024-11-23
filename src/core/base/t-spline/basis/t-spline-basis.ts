// src/core/base/t-spline/basis/t-spline-basis.ts

import { Point, Parameter, Index, Weight } from '../../../base/types';
import { BasisEvaluator, Domain } from '../../../base/parametric-space';
import { StarPoint } from '../star-point/star-point';

/**
 * Represents a T-spline basis function anchored at a T-junction
 */
interface TJunctionAnchor {
    position: Point;
    localKnots: number[][];
    transitionFunction?: (parameter: Parameter) => number;
}

/**
 * T-spline basis function configuration
 */
interface TSplineBasisConfig {
    index: Index;
    degree: number[];
    localKnots: number[][];
    starPoint?: StarPoint;
    tJunctions?: TJunctionAnchor[];
}

/**
 * Implementation of T-spline basis functions
 */
export class TSplineBasis implements BasisEvaluator {
    private index: Index;
    private degree: number[];
    private localKnots: number[][];
    private starPoint?: StarPoint;
    private tJunctions: TJunctionAnchor[];
    private support: TSplineDomain;

    constructor(config: TSplineBasisConfig) {
        this.index = config.index;
        this.degree = config.degree;
        this.localKnots = config.localKnots;
        this.starPoint = config.starPoint;
        this.tJunctions = config.tJunctions || [];
        this.support = this.computeSupport();
    }

    /**
     * Evaluate basis function at parameter
     */
    evaluate(parameter: Parameter): number {
        if (!this.support.contains(parameter)) {
            return 0;
        }

        if (this.isNearStarPoint(parameter)) {
            return this.evaluateNearStarPoint(parameter);
        }

        if (this.isNearTJunction(parameter)) {
            return this.evaluateNearTJunction(parameter);
        }

        return this.evaluateRegular(parameter);
    }

    /**
     * Evaluate basis function derivative
     */
    evaluateDerivative(parameter: Parameter, derivativeOrder: number[]): number {
        if (!this.support.contains(parameter)) {
            return 0;
        }

        if (this.isNearStarPoint(parameter)) {
            return this.evaluateStarPointDerivative(parameter, derivativeOrder);
        }

        if (this.isNearTJunction(parameter)) {
            return this.evaluateTJunctionDerivative(parameter, derivativeOrder);
        }

        return this.evaluateRegularDerivative(parameter, derivativeOrder);
    }

    /**
     * Get support domain of basis function
     */
    getSupport(): Domain {
        return this.support;
    }

    /**
     * Check if parameter is near a star point
     */
    private isNearStarPoint(parameter: Parameter): boolean {
        if (!this.starPoint) {
            return false;
        }

        const distance = this.computeDistanceToStarPoint(parameter);
        return distance < this.starPoint.getInfluenceRadius();
    }

    /**
     * Check if parameter is near a T-junction
     */
    private isNearTJunction(parameter: Parameter): boolean {
        return this.tJunctions.some(tj => 
            this.computeDistanceToTJunction(parameter, tj) < this.getTJunctionThreshold());
    }

    /**
     * Evaluate basis function in regular region
     */
    private evaluateRegular(parameter: Parameter): number {
        return this.degree.reduce((acc, deg, dim) => 
            acc * this.evaluateUnivariateBasis(
                parameter[dim],
                this.index[dim],
                deg,
                this.localKnots[dim]
            ), 1);
    }

    /**
     * Evaluate basis function near star point
     */
    private evaluateNearStarPoint(parameter: Parameter): number {
        if (!this.starPoint) {
            return 0;
        }

        const [r, theta] = this.starPoint.parameterization.mapToLocal(parameter);
        const regularValue = this.evaluateRegular(parameter);
        const starPointValue = this.starPoint.evaluateBasis(r, theta);
        
        const blend = this.computeBlendingFunction(r / this.starPoint.getInfluenceRadius());
        return regularValue * (1 - blend) + starPointValue * blend;
    }

    /**
     * Evaluate basis function near T-junction
     */
    private evaluateNearTJunction(parameter: Parameter): number {
        const regularValue = this.evaluateRegular(parameter);
        
        let tJunctionInfluence = 0;
        for (const tj of this.tJunctions) {
            const distance = this.computeDistanceToTJunction(parameter, tj);
            if (distance < this.getTJunctionThreshold()) {
                const blend = tj.transitionFunction?.(parameter) ?? 
                    this.computeBlendingFunction(distance / this.getTJunctionThreshold());
                tJunctionInfluence += blend;
            }
        }

        return regularValue * (1 - tJunctionInfluence);
    }

    /**
     * Evaluate univariate basis function using Cox-de Boor recursion
     */
    private evaluateUnivariateBasis(
        t: number,
        i: number,
        p: number,
        knots: number[]
    ): number {
        if (p === 0) {
            return (knots[i] <= t && t < knots[i + 1]) ? 1 : 0;
        }

        let value = 0;

        const d1 = knots[i + p] - knots[i];
        if (d1 > 0) {
            value += ((t - knots[i]) / d1) * 
                this.evaluateUnivariateBasis(t, i, p - 1, knots);
        }

        const d2 = knots[i + p + 1] - knots[i + 1];
        if (d2 > 0) {
            value += ((knots[i + p + 1] - t) / d2) * 
                this.evaluateUnivariateBasis(t, i + 1, p - 1, knots);
        }

        return value;
    }

    /**
     * Compute support domain
     */
    private computeSupport(): TSplineDomain {
        const min: Parameter = [];
        const max: Parameter = [];

        for (let i = 0; i < this.degree.length; i++) {
            min[i] = this.localKnots[i][0];
            max[i] = this.localKnots[i][this.localKnots[i].length - 1];
        }

        return new TSplineDomain(min, max);
    }

    /**
     * Compute blending function for smooth transitions
     */
    private computeBlendingFunction(t: number): number {
        // Cubic Hermite blending
        return t * t * (3 - 2 * t);
    }

    /**
     * Compute distance to star point
     */
    private computeDistanceToStarPoint(parameter: Parameter): number {
        if (!this.starPoint) {
            return Infinity;
        }

        const [r] = this.starPoint.parameterization.mapToLocal(parameter);
        return r;
    }

    /**
     * Compute distance to T-junction
     */
    private computeDistanceToTJunction(
        parameter: Parameter, 
        tJunction: TJunctionAnchor
    ): number {
        return Math.sqrt(parameter.reduce((sum, p, i) => 
            sum + Math.pow(p - tJunction.position[i], 2), 0));
    }

    /**
     * Get T-junction influence threshold
     */
    private getTJunctionThreshold(): number {
        return 0.1; // Configurable
    }

    /**
     * Evaluate derivative near star point
     */
    private evaluateStarPointDerivative(
        parameter: Parameter,
        derivativeOrder: number[]
    ): number {
        // Implement star point derivative computation
        return 0; // Placeholder
    }

    /**
     * Evaluate derivative near T-junction
     */
    private evaluateTJunctionDerivative(
        parameter: Parameter,
        derivativeOrder: number[]
    ): number {
        // Implement T-junction derivative computation
        return 0; // Placeholder
    }

    /**
     * Evaluate regular derivative
     */
    private evaluateRegularDerivative(
        parameter: Parameter,
        derivativeOrder: number[]
    ): number {
        // Implement regular derivative computation
        return 0; // Placeholder
    }
}

/**
 * T-spline domain implementation
 */
class TSplineDomain implements Domain {
    constructor(
        private min: Parameter,
        private max: Parameter
    ) {}

    get dimension(): number {
        return this.min.length;
    }

    contains(parameter: Parameter): boolean {
        return parameter.every((p, i) => 
            p >= this.min[i] && p <= this.max[i]);
    }

    project(parameter: Parameter): Parameter {
        return parameter.map((p, i) => 
            Math.max(this.min[i], Math.min(p, this.max[i])));
    }

    getBounds(): [Parameter, Parameter] {
        return [[...this.min], [...this.max]];
    }
}
