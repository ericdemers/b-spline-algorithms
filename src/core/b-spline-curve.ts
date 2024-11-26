import {
  BasisFunctionAlgorithm,
  BasisFunctionDerivativeAlgorithm,
  DegreeElevationAlgorithm,
  DegreeReductionAlgorithm,
  EvaluationAlgorithm,
  KnotInsertionAlgorithm,
  KnotRefinementAlgorithm,
  KnotRemovalAlgorithm,
  MultiplicationAlgorithm,
} from "./algorithms-interface";
import { Vector, VectorSpace } from "./vector-space";

/**
 * Implementation of a B-spline curve
 * Represents a parametric curve defined by control points, knot vector, and degree
 */
export class BSplineCurve<K extends number, V extends Vector> {
  //private readonly vectorSpace: VectorSpace<K, V>;
  //private readonly controlPoints: V[];
  //private readonly knots: number[];
  //private readonly degree: number;

  /**
   * Constructs a new B-spline curve
   *
   * @param vectorSpace - Vector space for control points
   * @param controlPoints - Array of control points
   * @param knots - Knot vector
   * @param degree - Degree of the B-spline
   */
  constructor(
    vectorSpace: VectorSpace<K, V>,
    controlPoints: V[],
    knots: number[],
    degree: number,
    options?: {
      evaluation?: EvaluationAlgorithm<V>;
      basisFunction?: BasisFunctionAlgorithm;
      basisFunctionDerivative?: BasisFunctionDerivativeAlgorithm;
      knotInsertion?: KnotInsertionAlgorithm<V>;
      knotRemoval?: KnotRemovalAlgorithm<V>;
      knotRefinement?: KnotRefinementAlgorithm<V>;
      degreeElevation?: DegreeElevationAlgorithm<V>;
      degreeReduction?: DegreeReductionAlgorithm<V>;
      multiplication?: MultiplicationAlgorithm<V>;
    }
  ) {}
}
