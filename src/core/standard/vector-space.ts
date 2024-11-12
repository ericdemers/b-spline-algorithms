/**
 * Core types and interfaces for B-spline implementation
 * Based on the mathematical theory from "The NURBS Book" by Les Piegl and Wayne Tiller
 */

// ------------ Type Definitions ------------

/** Real numbers (ℝ) */
export type Real = number;

/** Complex numbers (ℂ) represented as [real, imaginary] */
export type Complex = [number, number];

/** Scalar types supported in calculations */
export type Scalar = Real | Complex;

/** Generic vector type for n-dimensional space */
export type Vector = number[] | Complex[];

/** Specific vector type */
export type Vector1D = [number]
export type Vector2D = [number, number]
export type Vector3D = [number, number, number]
export type Vector4D = [number, number, number, number]
export type ComplexVector1D = [number, number]
export type ComplexVector2D = [[number, number], [number, number]]

// ------------ Mathematical Structures ------------

/**
 * Vector Space interface following mathematical axioms
 * V is a vector space over field K if it satisfies the vector space axioms
 */
export interface VectorSpace<K extends Scalar, V extends Vector> {
    /** Additive identity element (zero vector) */
    zero(): V;
    
    /** Vector addition (commutative group operation) */
    add(a: V, b: V): V;
    
    /** Scalar multiplication */
    scale(scalar: K, v: V): V;
    
    /** Vector subtraction (derived operation) */
    subtract(a: V, b: V): V;
    
    /** Dimension of the vector space */
    dimension(): number;
}

/**
 * Mathematical structures for B-spline computation
 */
interface Field<K> {
    zero(): K;
    one(): K;
    add(a: K, b: K): K;
    multiply(a: K, b: K): K;
    inverse(a: K): K;
}

/**
 * Real field implementation
 */
class RealField implements Field<number> {
    zero(): number { return 0; }
    one(): number { return 1; }
    add(a: number, b: number): number { return a + b; }
    multiply(a: number, b: number): number { return a * b; }
    inverse(a: number): number { return 1 / a; }
}

/**
 * Complex field implementation
 */
class ComplexField implements Field<Complex> {
    zero(): Complex { return [0, 0]; }
    one(): Complex { return [1, 0]; }
    add(a: Complex, b: Complex): Complex {
        return [a[0] + b[0], a[1] + b[1]];
    }
    multiply(a: Complex, b: Complex): Complex {
        return [
            a[0] * b[0] - a[1] * b[1],
            a[0] * b[1] + a[1] * b[0]
        ];
    }
    inverse(a: Complex): Complex {
        const denom = a[0] * a[0] + a[1] * a[1];
        return [a[0] / denom, -a[1] / denom];
    }
}


export class Vector1DSpace implements VectorSpace<Real, Vector1D> {
    zero(): Vector1D {
        return [0];
    }

    add(a: Vector1D, b: Vector1D): Vector1D {
        return [a[0] + b[0]];
    }

    scale(scalar: Real, v: Vector1D): Vector1D {
        return [scalar * v[0]];
    }

    subtract(a: Vector1D, b: Vector1D): Vector1D {
        return [a[0] - b[0]];
    }

    dimension() {
        return 1
    }
}

export class Vector2DSpace implements VectorSpace<Real, Vector2D> {
    zero(): Vector2D {
        return [0, 0];
    }

    add(a: Vector2D, b: Vector2D): Vector2D {
        return [a[0] + b[0], a[1] + b[1]];
    }

    scale(scalar: Real, v: Vector2D): Vector2D {
        return [scalar * v[0], scalar * v[1]];
    }

    subtract(a: Vector2D, b: Vector2D): Vector2D {
        return [a[0] - b[0], a[1] - b[1]];
    }

    dimension() {
        return 2
    }
}

export class Vector3DSpace implements VectorSpace<Real, Vector3D> {
    zero(): Vector3D {
        return [0, 0, 0];
    }

    add(a: Vector3D, b: Vector3D): Vector3D {
        return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
    }

    scale(scalar: Real, v: Vector3D): Vector3D {
        return [scalar * v[0], scalar * v[1], scalar * v[2]];
    }

    subtract(a: Vector3D, b: Vector3D): Vector3D {
        return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    }

    dimension() {
        return 3
    }
}

export class ComplexVector2DSpace implements VectorSpace<Complex, ComplexVector2D> {
    zero(): ComplexVector2D {
        return [[0, 0], [0, 0]];
    }

    add(a: ComplexVector2D, b: ComplexVector2D): ComplexVector2D {
        return [
            [a[0][0] + b[0][0], a[0][1] + b[0][1]],
            [a[1][0] + b[1][0], a[1][1] + b[1][1]]
        ];
    }

    scale(scalar: Complex, v: ComplexVector2D): ComplexVector2D {
        return [
            [scalar[0] * v[0][0] - scalar[1]*v[0][1], scalar[0] * v[0][1] + scalar[1] * v[0][0]],
            [scalar[0] * v[1][0] - scalar[1]*v[1][1], scalar[0] * v[1][1] + scalar[1] * v[1][0]]
        ];
    }

    subtract(a: ComplexVector2D, b: ComplexVector2D): ComplexVector2D {
        return [
            [a[0][0] - b[0][0], a[0][1] - b[0][1]],
            [a[1][0] - b[1][0], a[1][1] - b[1][1]]
        ];
    }

    dimension() {
        return 2
    }
}