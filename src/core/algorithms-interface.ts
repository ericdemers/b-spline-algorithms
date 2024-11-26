import { Vector } from "./vector-space";

export interface BSplineData<V extends Vector> {
  controlPoints: ReadonlyArray<V>;
  knots: ReadonlyArray<number>;
  degree: number;
}

// Algorithm interfaces
export interface EvaluationAlgorithm<V extends Vector> {
  evaluate(bspline: BSplineData<V>, parameter: number): V;
}

export interface KnotInsertionAlgorithm<V extends Vector> {
  insertKnot(bspline: BSplineData<V>, parameter: number): BSplineData<V>;
}

export interface MultiplicationAlgorithm<V extends Vector> {
  multiply(f: BSplineData<V>, g: BSplineData<V>): BSplineData<V>;
}

export interface BasisFunctionAlgorithm {
  basisFunction(
    knots: ReadonlyArray<number>,
    index: number,
    degree: number,
    u: number
  ): number;
}
export interface BasisFunctionDerivativeAlgorithm {
  basisFunctionDerivative(
    knots: ReadonlyArray<number>,
    index: number,
    degree: number,
    u: number
  ): number;
}
export interface KnotRemovalAlgorithm<V extends Vector> {
  removeKnot(bspline: BSplineData<V>, parameter: number): BSplineData<V>;
}
export interface KnotRefinementAlgorithm<V extends Vector> {
  refineKnots(
    bspline: BSplineData<V>,
    parameters: ReadonlyArray<number>
  ): BSplineData<V>;
}
export interface DegreeElevationAlgorithm<V extends Vector> {
  elevateDegree(bspline: BSplineData<V>, degree: number): BSplineData<V>;
}
export interface DegreeReductionAlgorithm<V extends Vector> {
  reduceDegree(bspline: BSplineData<V>, degree: number): BSplineData<V>;
}
