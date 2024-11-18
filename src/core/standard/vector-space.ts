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
    private dim: number;

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
    private dim: number;

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

    scale(scalar: Complex, v: ComplexVector): ComplexVector {
        if (isVector1D(v)) {
            return ComplexOps.multiply(scalar, v as Complex);
        }
        return (v as Complex[]).map(component => 
            ComplexOps.multiply(scalar, component)
        );
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




// ------------ Helper Functions ------------

/**
 * Creates a vector of specified dimension with optional initial value
 */
export function createVector(dimension: number, value: number = 0): Vector {
    if (dimension < 1) {
        throw new Error('Dimension must be positive');
    }
    return dimension === 1 ? value : new Array(dimension).fill(value);
}

/**
 * Gets a component of a vector at specified index
 */
export function getVectorComponent(v: Vector, index: number): number | Complex {
    if (isVector1D(v)) {
        if (index !== 0) {
            throw new Error('Index out of bounds for 1D vector');
        }
        return v;
    }
    if (index < 0 || index >= (v as any[]).length) {
        throw new Error('Index out of bounds');
    }
    return v[index];
}

export function setVectorComponent<T extends Vector>(v: T, index: number, value: number | Complex): T {
    if (isVector1D(v)) {
        if (index !== 0) {
            throw new Error('Index out of bounds for 1D vector');
        }
        return value as T;
    }

    if (index < 0 || index >= (v as any[]).length) {
        throw new Error('Index out of bounds');
    }
    
    if (isRealVector(v)) {
        const newV = Array.from(v as number[]);
        newV[index] = value as number;
        return newV as T;
    }
    
    if (isComplexVector(v)) {
        const newV = Array.from(v as Complex[]);
        newV[index] = value as Complex;
        return newV as T;
    }
    
    throw new Error('Unsupported vector type');
}





// /**
//  * Mathematical structures for B-spline computation
//  */
// interface Field<K> {
//     zero(): K;
//     one(): K;
//     add(a: K, b: K): K;
//     multiply(a: K, b: K): K;
//     inverse(a: K): K;
// }


// /**
//  * Real field implementation
//  */
// class RealField implements Field<number> {
//     zero(): number { return 0; }
//     one(): number { return 1; }
//     add(a: number, b: number): number { return a + b; }
//     multiply(a: number, b: number): number { return a * b; }
//     inverse(a: number): number { return 1 / a; }
// }

// /**
//  * Complex field implementation
//  */
// class ComplexField implements Field<Complex> {
//     zero(): Complex { return [0, 0]; }
//     one(): Complex { return [1, 0]; }
//     add(a: Complex, b: Complex): Complex {
//         return [a[0] + b[0], a[1] + b[1]];
//     }
//     multiply(a: Complex, b: Complex): Complex {
//         return [
//             a[0] * b[0] - a[1] * b[1],
//             a[0] * b[1] + a[1] * b[0]
//         ];
//     }
//     inverse(a: Complex): Complex {
//         const denom = a[0] * a[0] + a[1] * a[1];
//         return [a[0] / denom, -a[1] / denom];
//     }
// }

// export class Vector1DSpace implements VectorSpace<Real, Vector1D> {
//     zero(): Vector1D {
//         return 0;
//     }

//     add(a: Vector1D, b: Vector1D): Vector1D {
//         return a + b;
//     }

//     scale(scalar: Real, v: Vector1D): Vector1D {
//         return scalar * v;
//     }

//     subtract(a: Vector1D, b: Vector1D): Vector1D {
//         return a - b;
//     }

//     dimension() {
//         return 1
//     }
// }

// export class Vector2DSpace implements VectorSpace<Real, Vector2D> {
//     zero(): Vector2D {
//         return [0, 0];
//     }

//     add(a: Vector2D, b: Vector2D): Vector2D {
//         return [a[0] + b[0], a[1] + b[1]];
//     }

//     scale(scalar: Real, v: Vector2D): Vector2D {
//         return [scalar * v[0], scalar * v[1]];
//     }

//     subtract(a: Vector2D, b: Vector2D): Vector2D {
//         return [a[0] - b[0], a[1] - b[1]];
//     }

//     dimension() {
//         return 2
//     }
// }

// export class Vector3DSpace implements VectorSpace<Real, Vector3D> {
//     zero(): Vector3D {
//         return [0, 0, 0];
//     }

//     add(a: Vector3D, b: Vector3D): Vector3D {
//         return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
//     }

//     scale(scalar: Real, v: Vector3D): Vector3D {
//         return [scalar * v[0], scalar * v[1], scalar * v[2]];
//     }

//     subtract(a: Vector3D, b: Vector3D): Vector3D {
//         return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
//     }

//     dimension() {
//         return 3
//     }
// }

// export class ComplexVector2DSpace implements VectorSpace<Complex, ComplexVector2D> {
//     zero(): ComplexVector2D {
//         return [[0, 0], [0, 0]];
//     }

//     add(a: ComplexVector2D, b: ComplexVector2D): ComplexVector2D {
//         return [
//             [a[0][0] + b[0][0], a[0][1] + b[0][1]],
//             [a[1][0] + b[1][0], a[1][1] + b[1][1]]
//         ];
//     }

//     scale(scalar: Complex, v: ComplexVector2D): ComplexVector2D {
//         return [
//             [scalar[0] * v[0][0] - scalar[1]*v[0][1], scalar[0] * v[0][1] + scalar[1] * v[0][0]],
//             [scalar[0] * v[1][0] - scalar[1]*v[1][1], scalar[0] * v[1][1] + scalar[1] * v[1][0]]
//         ];
//     }

//     subtract(a: ComplexVector2D, b: ComplexVector2D): ComplexVector2D {
//         return [
//             [a[0][0] - b[0][0], a[0][1] - b[0][1]],
//             [a[1][0] - b[1][0], a[1][1] - b[1][1]]
//         ];
//     }

//     dimension() {
//         return 2
//     }
// }