import { defineConfig } from "vitest/config";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  define: {
    global: "globalThis"
  },
  resolve: {
    alias: {
      events: "events/events.js"
    },
    dedupe: ["react", "react-dom", "@emotion/react", "@emotion/styled"]
  },
  build: {
    commonjsOptions: {
      // ketcher-core 3.18 ships an ESM bundle containing require("raphael").
      // Convert CommonJS calls inside otherwise-ESM dependencies as well.
      transformMixedEsModules: true
    }
  },
  server: {
    host: "127.0.0.1",
    port: 5173,
    proxy: {
      "/api": "http://127.0.0.1:8000"
    }
  },
  test: {
    environment: "jsdom",
    setupFiles: "./src/test/setup.ts",
    include: ["src/**/*.test.{ts,tsx}"],
    css: true
  }
});
