// knot-structure.ts

import { BaseKnotStructure, Knots, KnotStructure } from "./knot-structure2";
import { PeriodicKnots } from "./periodicKnots2";

/**
 * Enum to specify the type of knot sequence in each dimension
 */
export enum KnotType {
    REGULAR = 'regular',
    PERIODIC = 'periodic'
}

/**
 * Configuration for a knot sequence dimension
 */
export interface KnotDimensionConfig {
    type: KnotType;
    knots: ReadonlyArray<number>;
    period?: number;  // Required for periodic dimensions
}

/**
 * Mixed knot structure supporting both periodic and regular dimensions
 */
export class MixedKnotStructure extends BaseKnotStructure {
    private readonly dimensionStructures: KnotStructure[];

    /**
     * Creates a mixed knot structure
     * @param configs - Configuration for each dimension
     */
    constructor(private readonly configs: ReadonlyArray<KnotDimensionConfig>) {
        super();
        this.dimensionStructures = this.initializeDimensions(configs);
    }

    getDimension(): number {
        return this.configs.length;
    }

    getKnotSequence(direction: number): ReadonlyArray<number> {
        this.validateDirection(direction);
        return this.dimensionStructures[direction].getKnotSequence(0);
    }

    /**
     * Gets the type of knot sequence for a given dimension
     */
    getKnotType(direction: number): KnotType {
        this.validateDirection(direction);
        return this.configs[direction].type;
    }

    withInsertedKnot(dimension: number, u: number): KnotStructure {
        this.validateDirection(dimension);
        
        const newConfigs = [...this.configs];
        const newDimension = this.dimensionStructures[dimension].withInsertedKnot(0, u);
        
        // Update the configuration for the modified dimension
        if (this.configs[dimension].type === KnotType.PERIODIC) {
            newConfigs[dimension] = {
                type: KnotType.PERIODIC,
                knots: (newDimension as PeriodicKnots).getPattern(),
                period: (newDimension as PeriodicKnots).getPeriod()
            };
        } else {
            newConfigs[dimension] = {
                type: KnotType.REGULAR,
                knots: newDimension.getKnotSequence(0)
            };
        }

        return new MixedKnotStructure(newConfigs);
    }

    withRemovedKnots(dimension: number): KnotStructure {
        this.validateDirection(dimension);
        
        const newConfigs = [...this.configs];
        const newDimension = this.dimensionStructures[dimension].withRemovedKnots(0);
        
        // Update the configuration for the modified dimension
        if (this.configs[dimension].type === KnotType.PERIODIC) {
            newConfigs[dimension] = {
                type: KnotType.PERIODIC,
                knots: (newDimension as PeriodicKnots).getKnotPattern(),
                period: (newDimension as PeriodicKnots).getPeriod()
            };
        } else {
            newConfigs[dimension] = {
                type: KnotType.REGULAR,
                knots: newDimension.getKnotSequence(0)
            };
        }

        return new MixedKnotStructure(newConfigs);
    }

    /**
     * Gets domain information for each dimension
     */
    getDomains(): Array<{ min: number; max: number; isPeriodic: boolean }> {
        return this.dimensionStructures.map((structure, index) => {
            const isPeriodic = this.configs[index].type === KnotType.PERIODIC;
            if (isPeriodic) {
                const periodicStruct = structure as PeriodicKnots;
                return {
                    ...periodicStruct.getDomain(),
                    isPeriodic: true
                };
            } else {
                const knots = structure.getKnotSequence(0);
                return {
                    min: knots[0],
                    max: knots[knots.length - 1],
                    isPeriodic: false
                };
            }
        });
    }

    private initializeDimensions(configs: ReadonlyArray<KnotDimensionConfig>): KnotStructure[] {
        return configs.map(config => {
            if (config.type === KnotType.PERIODIC) {
                if (config.period === undefined) {
                    throw new Error('Period must be specified for periodic dimensions');
                }
                return new PeriodicKnots(config.knots, config.period);
            } else {
                return new Knots(config.knots);
            }
        });
    }
}

/**
 * Builder class for creating mixed knot structures
 */
export class MixedKnotStructureBuilder {
    private configs: KnotDimensionConfig[] = [];

    /**
     * Adds a regular (non-periodic) dimension
     */
    addRegularDimension(knots: number[]): this {
        this.configs.push({
            type: KnotType.REGULAR,
            knots: knots
        });
        return this;
    }

    /**
     * Adds a periodic dimension
     */
    addPeriodicDimension(knots: number[], period: number): this {
        this.configs.push({
            type: KnotType.PERIODIC,
            knots: knots,
            period: period
        });
        return this;
    }

    /**
     * Builds the mixed knot structure
     */
    build(): MixedKnotStructure {
        if (this.configs.length === 0) {
            throw new Error('At least one dimension must be specified');
        }
        return new MixedKnotStructure(this.configs);
    }
}

// Create a mixed knot structure using the builder
const mixedKnots = new MixedKnotStructureBuilder()
    // Regular dimension (e.g., x direction)
    .addRegularDimension([0, 0, 0, 1, 2, 3, 3, 3])
    // Periodic dimension (e.g., y direction)
    .addPeriodicDimension([0, 0, 1, 2, 2], 2)
    // Another regular dimension (e.g., z direction)
    .addRegularDimension([0, 0, 1, 1])
    .build();

// Get knot sequence for each dimension
const xKnots = mixedKnots.getKnotSequence(0);  // Regular
const yKnots = mixedKnots.getKnotSequence(1);  // Periodic
const zKnots = mixedKnots.getKnotSequence(2);  // Regular

// Check knot type for a dimension
const isYPeriodic = mixedKnots.getKnotType(1) === KnotType.PERIODIC;

// Get domain information
const domains = mixedKnots.getDomains();
// domains[0] = { min: 0, max: 3, isPeriodic: false }
// domains[1] = { min: 0, max: 2, isPeriodic: true }
// domains[2] = { min: 0, max: 1, isPeriodic: false }

// Insert a knot in the regular dimension
const newMixedKnots = mixedKnots.withInsertedKnot(0, 1.5);

