import { ControlPolygon } from '../../src/core/control-net';
import { Vector2D } from '../../src/core/vector-space';

describe('ControlPolygon', () => {
  let controlPolygon: ControlPolygon<Vector2D>;

  beforeEach(() => {
    // Initialize a ControlPolygon with some control points before each test
    controlPolygon = new ControlPolygon<Vector2D>([
      [0, 0],      // Control point 1
      [1, 1],      // Control point 2
      [2, 0],      // Control point 3
    ]);
  });

  test('getPoint returns the correct control point', () => {
    const point = controlPolygon.getPoint(1); // Get the control point at index 1
    expect(point).toEqual([1, 1]); // Should return the second control point
  });

  test('getDimension returns the correct dimension', () => {
    const dimension = controlPolygon.getDimension();
    expect(dimension).toBe(1); // Should return 1
  });

  test('getSize returns the correct size', () => {
    const size = controlPolygon.getSize(); // Size in the first direction
    expect(size).toBe(3); // Should return the number of control points
  });
});