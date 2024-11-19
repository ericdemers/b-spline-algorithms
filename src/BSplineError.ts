export class BSplineError extends Error {
    constructor(message: string, public readonly code: string) {
        super(message);
        this.name = 'BSplineError';
    }
}

export function validateInputs(knots: number[], degree: number): void {
    if (knots.length < degree + 1) {
        throw new BSplineError(
            'Not enough knots for specified degree',
            'INVALID_KNOT_COUNT'
        );
    }
}
