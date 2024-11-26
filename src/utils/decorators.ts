export function adaptParameter() {
    return function (
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) {
      const originalMethod = descriptor.value;
      descriptor.value = function (position: number | number[]) {
        const input = Array.isArray(position) ? position : [position];
        return originalMethod.call(this, input);
      };
    };
  }