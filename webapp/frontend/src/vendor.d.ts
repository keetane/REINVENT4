declare module "process" {
  const processShim: {
    env: Record<string, string | undefined>;
    nextTick: (callback: (...args: unknown[]) => void, ...args: unknown[]) => void;
    [key: string]: unknown;
  };
  export default processShim;
}

