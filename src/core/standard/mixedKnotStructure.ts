import { BaseKnotStructure, Domain, KnotStructure, KnotValue } from './knotStructure';
import { PeriodicKnots } from './periodicKnots';
import { Knots } from './knotStructure';

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
    readonly type: KnotType;
    readonly knots: ReadonlyArray<number>;
    readonly period?: number;  // Required for periodic dimensions
    readonly degree?: number;  // Optional: can be used for validation
}

/**
 * Domain information with periodicity
 */
export interface ExtendedDomain extends Domain {
    readonly isPeriodic: boolean;
    readonly period?: number;
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
        this.validateConfigs(configs);
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
            const periodicDim = newDimension as PeriodicKnots;
            newConfigs[dimension] = {
                ...this.configs[dimension],
                knots: periodicDim.getPattern(),
            };
        } else {
            newConfigs[dimension] = {
                ...this.configs[dimension],
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
            const periodicDim = newDimension as PeriodicKnots;
            newConfigs[dimension] = {
                ...this.configs[dimension],
                knots: periodicDim.getPattern(),
            };
        } else {
            newConfigs[dimension] = {
                ...this.configs[dimension],
                knots: newDimension.getKnotSequence(0)
            };
        }

        return new MixedKnotStructure(newConfigs);
    }

    getDomain(direction: number): ExtendedDomain {
        this.validateDirection(direction);
        const config = this.configs[direction];
        const baseDomain = this.dimensionStructures[direction].getDomain(0);

        return {
            ...baseDomain,
            isPeriodic: config.type === KnotType.PERIODIC,
            period: config.period
        };
    }

    getDistinctKnots(direction: number): ReadonlyArray<KnotValue> {
        this.validateDirection(direction);
        return this.dimensionStructures[direction].getDistinctKnots(0);
    }

    /**
     * Gets domain information for all dimensions
     */
    getAllDomains(): ReadonlyArray<ExtendedDomain> {
        return Array.from({ length: this.getDimension() }, (_, i) => this.getDomain(i));
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

    private validateConfigs(configs: ReadonlyArray<KnotDimensionConfig>): void {
        if (!configs || configs.length === 0) {
            throw new Error('At least one dimension configuration must be provided');
        }

        configs.forEach((config, i) => {
            if (!config.knots || config.knots.length < 2) {
                throw new Error(`Dimension ${i}: Knot vector must have at least 2 values`);
            }

            if (config.type === KnotType.PERIODIC && config.period === undefined) {
                throw new Error(`Dimension ${i}: Period must be specified for periodic dimensions`);
            }

            if (config.type === KnotType.PERIODIC && config.period! <= 0) {
                throw new Error(`Dimension ${i}: Period must be positive`);
            }

            if (config.degree !== undefined && config.degree < 0) {
                throw new Error(`Dimension ${i}: Degree must be non-negative`);
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
     * @param knots - Knot vector
     * @param degree - Optional degree for validation
     */
    addRegularDimension(knots: ReadonlyArray<number>, degree?: number): this {
        this.configs.push({
            type: KnotType.REGULAR,
            knots: knots,
            degree
        });
        return this;
    }

    /**
     * Adds a periodic dimension
     * @param knots - Knot pattern
     * @param period - Period length
     * @param degree - Optional degree for validation
     */
    addPeriodicDimension(
        knots: ReadonlyArray<number>, 
        period: number,
        degree?: number
    ): this {
        this.configs.push({
            type: KnotType.PERIODIC,
            knots: knots,
            period: period,
            degree
        });
        return this;
    }

    /**
     * Builds the mixed knot structure
     * @throws Error if no dimensions are specified
     */
    build(): MixedKnotStructure {
        if (this.configs.length === 0) {
            throw new Error('At least one dimension must be specified');
        }
        return new MixedKnotStructure(this.configs);
    }

    /**
     * Clears all configurations
     */
    clear(): this {
        this.configs = [];
        return this;
    }
}

// Example usage:
/*
const mixedKnots = new MixedKnotStructureBuilder()
    // Regular dimension (e.g., x direction)
    .addRegularDimension([0, 0, 0, 1, 2, 3, 3, 3], 2)
    // Periodic dimension (e.g., y direction)
    .addPeriodicDimension([0, 0, 1, 2, 2], 2, 2)
    // Another regular dimension (e.g., z direction)
    .addRegularDimension([0, 0, 1, 1], 1)
    .build();
*/
