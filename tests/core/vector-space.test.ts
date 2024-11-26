import { RealVectorSpace, ComplexVectorSpace, ComplexVector2D, Complex, Real } from '../../src/core/vector-space';

describe('RealVectorSpace', () => {
  describe('Vector Space Axioms', () => {
    const V = new RealVectorSpace(3);
    
    test('additive identity', () => {
      const v = [1, 2, 3];
      const zero = V.zero();
      expect(V.add(v, zero)).toEqual(v);
      expect(V.add(zero, v)).toEqual(v);
    });

    test('multiplicative identity', () => {
      const v = [1, 2, 3];
      expect(V.scale(1, v)).toEqual(v);
    });
  });
});

describe('ComplexVectorSpace for 2D Complex Vectors', () => {
  const vectorSpace = new ComplexVectorSpace(2); // 2D complex vector space

  test('addition of complex vectors', () => {
    const v1 = [[1, 2], [3, 4]] as ComplexVector2D; // Representing (1 + 2i) and (3 + 4i)
    const v2 = [[5, 6], [7, 8]] as ComplexVector2D; // Representing (5 + 6i) and (7 + 8i)
    
    const result = vectorSpace.add(v1 , v2);
    expect(result).toEqual([[6, 8], [10, 12]]); // (1 + 2i) + (5 + 6i) = 6 + 8i and (3 + 4i) + (7 + 8i) = 10 + 12i
  });

  test('scaling a complex vector', () => {
    const v = [[1, 2], [3, 4]] as ComplexVector2D; // Representing (1 + 2i) and (3 + 4i)
    const scalar = [2, 0] as Complex;
    
    const scaled = vectorSpace.scale(scalar, v);
    expect(scaled).toEqual([[2, 4], [6, 8]]); // 2 * (1 + 2i) = 2 + 4i and 2 * (3 + 4i) = 6 + 8i
  });

  test('additive identity', () => {
    const v = [[1, 2], [3, 4]] as ComplexVector2D; // Representing (1 + 2i) and (3 + 4i)
    const zero = vectorSpace.zero(); // Assuming zero() returns [[0, 0], [0, 0]]
    
    expect(vectorSpace.add(v, zero)).toEqual(v);
    expect(vectorSpace.add(zero, v)).toEqual(v);
  });

  test('multiplicative identity', () => {
    const v = [[1, 2], [3, 4]] as ComplexVector2D; // Representing (1 + 2i) and (3 + 4i)
    expect(vectorSpace.scale([1, 0] as Complex, v)).toEqual(v);
    expect(vectorSpace.scale(1 as Real, v)).toEqual(v);
  });

  test('scaling a complex vector by a real number', () => {
    const vector = [[1, 2], [3, 4]] as ComplexVector2D; // Representing (1 + 2i) and (3 + 4i)
    const scalar = 2; // Real scalar
    const scaled = vectorSpace.scale(scalar, vector);
    expect(scaled).toEqual([[2, 4], [6, 8]]); // 2 * (1 + 2i) = 2 + 4i and 2 * (3 + 4i) = 6 + 8i
  });

  test('scaling a complex vector by a complex number', () => {
    const vector = [[1, 2], [3, 4]] as ComplexVector2D; // Representing (1 + 2i) and (3 + 4i)
    const scalar = [2, 3] as Complex; // Complex scalar (2 + 3i)
    const scaled = vectorSpace.scale(scalar, vector);
    expect(scaled).toEqual([[-4, 7], [-6, 17]]); 
  });

});

