import ReactDOM from "react-dom/client";
import processShim from "process";

import "./i18n";
import { App } from "./App";
import "./styles.css";

Object.assign(globalThis, { process: processShim });

ReactDOM.createRoot(document.getElementById("root")!).render(<App />);
