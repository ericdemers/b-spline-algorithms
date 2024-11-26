/**
 * Core types and interfaces for vector space operations
 * Implements mathematical vector space axioms for real and complex numbers
 */

// ------------ Type Definitions ------------

/** Real numbers (ℝ) */
export type Real = number;

/** Complex numbers (ℂ) represented as [real, imaginary] */
export type Complex = [number, number];

/** Scalar types supported in calculations */
export type Scalar = Real | Complex;

/** Generic vector type for n-dimensional space */
export type RealVector = number | number[];
export type ComplexVector = Complex | Complex[];
export type Vector = RealVector | ComplexVector;

/** Specific vector type */
export type Vector1D = number
export type Vector2D = [number, number]
export type Vector3D = [number, number, number]
export type Vector4D = [number, number, number, number]
export type ComplexVector1D = Complex
export type ComplexVector2D = [Complex, Complex]

// ------------ Type Guards ------------

/**
 * Checks if the vector is one-dimensional
 * @param v Vector to check
 * @returns True if vector is 1D (number or Complex)
 */
export function isVector1D(v: Vector): v is number | Complex {
    return typeof v === 'number' || 
            (Array.isArray(v) && v.length === 2 && !Array.isArray(v[0]));
}

/**
 * Checks if the vector contains real numbers
 * @param v Vector to check
 * @returns True if vector contains real numbers
 */
export function isRealVector(v: Vector): v is RealVector {
    return typeof v === 'number' || (Array.isArray(v) && typeof v[0] === 'number');
}

/**
 * Checks if the vector contains complex numbers
 * @param v Vector to check
 * @returns True if vector contains complex numbers
 */
export function isComplexVector(v: Vector): v is ComplexVector {
    return Array.isArray(v) && (v.length === 2 || Array.isArray(v[0]));
}

// ------------ Complex Number Operations ------------

/**
 * Complex number operations helper class
 */
export class ComplexOps {
    /**
     * Adds two complex numbers
     */
    static add(a: Complex, b: Complex): Complex {
        return [a[0] + b[0], a[1] + b[1]];
    }

    /**
     * Multiplies two complex numbers
     */
    static multiply(a: Complex, b: Complex): Complex {
        return [
            a[0] * b[0] - a[1] * b[1],
            a[0] * b[1] + a[1] * b[0]
        ];
    }

    /**
     * Subtracts two complex numbers
     */
    static subtract(a: Complex, b: Complex): Complex {
        return [a[0] - b[0], a[1] - b[1]];
    }

    /**
     * Returns the complex conjugate
     */
    static conjugate(a: Complex): Complex {
        return [a[0], -a[1]];
    }
}

// ------------ Vector Space Interface ------------

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

// ------------ Vector Space Implementations ------------

/**
 * Implementation of a real vector space
 */
export class RealVectorSpace implements VectorSpace<Real, RealVector> {
    private readonly dim: number;

    constructor(dimension: number) {
        if (dimension < 1) {
            throw new Error('Dimension must be positive');
        }
        this.dim = dimension;
    }

    zero(): RealVector {
        return this.dim === 1 ? 0 : new Array(this.dim).fill(0);
    }

    add(a: RealVector, b: RealVector): RealVector {
        if (typeof a === 'number' && typeof b === 'number') {
            return a + b;
        }
        return (a as number[]).map((val, i) => val + (b as number[])[i]);
    }

    scale(scalar: Real, v: RealVector): RealVector {
        if (typeof v === 'number') {
            return scalar * v;
        }
        return (v as number[]).map(val => scalar * val);
    }

    subtract(a: RealVector, b: RealVector): RealVector {
        if (typeof a === 'number' && typeof b === 'number') {
            return a - b;
        }
        return (a as number[]).map((val, i) => val - (b as number[])[i]);
    }

    dimension(): number {
        return this.dim;
    }
}

export class ComplexVectorSpace implements VectorSpace<Complex, ComplexVector> {
    private readonly dim: number;

    constructor(dimension: number) {
        if (dimension < 1) {
            throw new Error('Dimension must be positive');
        }
        this.dim = dimension;
    }

    zero(): ComplexVector {
        return this.dim === 1 ? [0, 0] : Array(this.dim).fill([0, 0]);
    }

    add(a: ComplexVector, b: ComplexVector): ComplexVector {
        if (isVector1D(a) && isVector1D(b)) {
            return ComplexOps.add(a as Complex, b as Complex);
        }
        return (a as Complex[]).map((val, i) => 
            ComplexOps.add(val, (b as Complex[])[i])
        );
    }

    // Overloaded scale method
    scale(scalar: Complex, vector: Complex[]): Complex[];
    scale(scalar: number, vector: Complex[]): Complex[];
    // Implementation of the scale method
    scale(scalar: Complex | number, vector: Complex[]): Complex[] {
        if (typeof scalar === 'number') {
        // Directly scale the real part of the complex vector
        return vector.map(([real, imag]) => [
            real * scalar, // Scale real part
            imag * scalar, // Scale imaginary part
        ]);
        } else {
        // Handle the case where scalar is a complex number
        return vector.map(([real, imag]) => [
            real * scalar[0] - imag * scalar[1], // Real part
            real * scalar[1] + imag * scalar[0], // Imaginary part
        ]);
        }
    }



    subtract(a: ComplexVector, b: ComplexVector): ComplexVector {
        if (isVector1D(a) && isVector1D(b)) {
            return ComplexOps.subtract(a as Complex, b as Complex);
        }
        return (a as Complex[]).map((val, i) => 
            ComplexOps.subtract(val, (b as Complex[])[i])
        );
    }

    dimension() {
        return this.dim
    }
}