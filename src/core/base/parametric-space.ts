// src/core/base/parametric-space.ts

import { Point, Parameter, Index } from './types';

/**
 * Represents a domain in parametric space
 */
export interface Domain {
    /** Number of parametric dimensions */
    dimension: number;

    /** Check if parameter is within domain bounds */
    contains(parameter: Parameter): boolean;

    /** Project parameter onto valid domain */
    project(parameter: Parameter): Parameter;

    /** Get domain bounds */
    getBounds(): [Parameter, Parameter];
}

/**
 * Interface for basis function evaluation
 */
export interface BasisEvaluator {
    /** Evaluate basis function at parameter */
    evaluate(parameter: Parameter): number;

    /** Evaluate basis function derivative */
    evaluateDerivative(parameter: Parameter, derivativeOrder: number[]): number;

    /** Get support domain of basis function */
    getSupport(): Domain;
}

/**
 * Abstract class representing a parametric space
 */
export abstract class ParametricSpace {
    protected domain: Domain;
    protected degree: number[];
    
    constructor(domain: Domain, degree: number[]) {
        if (domain.dimension !== degree.length) {
            throw new Error('Domain dimension must match degree array length');
        }
        this.domain = domain;
        this.degree = degree;
    }

    /**
     * Evaluate all non-zero basis functions at parameter
     */
    abstract getActiveBasis(parameter: Parameter): Map<string, BasisEvaluator>;

    /**
     * Evaluate specific basis function at parameter
     */
    abstract evaluateBasis(parameter: Parameter, index: Index): number;

    /**
     * Evaluate basis function derivative
     */
    abstract evaluateBasisDerivative(
        parameter: Parameter, 
        index: Index, 
        derivativeOrder: number[]
    ): number;

    /**
     * Get indices of all basis functions affecting parameter
     */
    abstract getActiveIndices(parameter: Parameter): Index[];

    /**
     * Check if refinement is possible at parameter
     */
    abstract canRefine(parameter: Parameter): boolean;

    /**
     * Perform local refinement at parameter
     */
    abstract refine(parameter: Parameter): void;

    /**
     * Get parametric domain
     */
    getDomain(): Domain {
        return this.domain;
    }

    /**
     * Get degree in each parametric direction
     */
    getDegree(): number[] {
        return [...this.degree];
    }

    /**
     * Validate parameter
     */
    protected validateParameter(parameter: Parameter): Parameter {
        if (parameter.length !== this.domain.dimension) {
            throw new Error('Invalid parameter dimension');
        }
        return this.domain.contains(parameter) ? 
            parameter : 
            this.domain.project(parameter);
    }
}

/**
 * Implementation of a tensor product parametric space
 */
export class TensorProductSpace extends ParametricSpace {
    private knots: number[][];

    constructor(degree: number[], knots: number[][]) {
        const domain = new TensorProductDomain(knots);
        super(domain, degree);
        this.knots = knots;
    }

    getActiveBasis(parameter: Parameter): Map<string, BasisEvaluator> {
        const validated = this.validateParameter(parameter);
        const result = new Map<string, BasisEvaluator>();
        
        const activeIndices = this.getActiveIndices(validated);
        for (const index of activeIndices) {
            result.set(
                index.join(','),
                new TensorProductBasis(this.degree, this.knots, index)
            );
        }
        
        return result;
    }

    evaluateBasis(parameter: Parameter, index: Index): number {
        const validated = this.validateParameter(parameter);
        return new TensorProductBasis(
            this.degree, 
            this.knots, 
            index
        ).evaluate(validated);
    }

    evaluateBasisDerivative(
        parameter: Parameter, 
        index: Index, 
        derivativeOrder: number[]
    ): number {
        const validated = this.validateParameter(parameter);
        return new TensorProductBasis(
            this.degree, 
            this.knots, 
            index
        ).evaluateDerivative(validated, derivativeOrder);
    }

    getActiveIndices(parameter: Parameter): Index[] {
        const validated = this.validateParameter(parameter);
        // Implementation depends on knot structure
        // Returns array of multi-indices for active basis functions
        return [];  // Placeholder
    }

    canRefine(parameter: Parameter): boolean {
        const validated = this.validateParameter(parameter);
        // Check if refinement is possible at parameter
        return true;  // Placeholder
    }

    refine(parameter: Parameter): void {
        const validated = this.validateParameter(parameter);
        // Implement knot insertion
    }
}

/**
 * Implementation of tensor product domain
 */
class TensorProductDomain implements Domain {
    private bounds: [Parameter, Parameter];

    constructor(knots: number[][]) {
        this.bounds = [
            knots.map(k => k[0]),
            knots.map(k => k[k.length - 1])
        ];
    }

    get dimension(): number {
        return this.bounds[0].length;
    }

    contains(parameter: Parameter): boolean {
        return parameter.every((p, i) => 
            p >= this.bounds[0][i] && p <= this.bounds[1][i]);
    }

    project(parameter: Parameter): Parameter {
        return parameter.map((p, i) => 
            Math.max(this.bounds[0][i], 
                Math.min(p, this.bounds[1][i])));
    }

    getBounds(): [Parameter, Parameter] {
        return [[...this.bounds[0]], [...this.bounds[1]]];
    }
}

/**
 * Implementation of tensor product basis
 */
class TensorProductBasis implements BasisEvaluator {
    constructor(
        private degree: number[],
        private knots: number[][],
        private index: Index
    ) {}

    evaluate(parameter: Parameter): number {
        return this.degree.reduce((acc, deg, dim) => 
            acc * this.evaluateUnivariate(
                parameter[dim],
                this.index[dim],
                deg,
                this.knots[dim]
            ), 1);
    }

    evaluateDerivative(parameter: Parameter, derivativeOrder: number[]): number {
        // Implement derivative computation
        return 0;  // Placeholder
    }

    getSupport(): Domain {
        // Compute support domain from knots and index
        return new TensorProductDomain(this.knots);  // Placeholder
    }

    private evaluateUnivariate(
        t: number,
        i: number,
        p: number,
        knots: number[]
    ): number {
        // Implement Cox-de Boor recursion
        return 0;  // Placeholder
    }
}
