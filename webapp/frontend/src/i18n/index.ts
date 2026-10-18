import i18n from "i18next";
import { initReactI18next } from "react-i18next";

import { resources } from "./resources";

const saved = window.localStorage.getItem("molecular-generator.language");

void i18n.use(initReactI18next).init({
  resources,
  lng: saved === "en" ? "en" : "ja",
  fallbackLng: "en",
  interpolation: { escapeValue: false },
  returnNull: false
});

i18n.on("languageChanged", (language) => {
  window.localStorage.setItem("molecular-generator.language", language === "en" ? "en" : "ja");
  document.documentElement.lang = language === "en" ? "en" : "ja";
});

document.documentElement.lang = i18n.language === "en" ? "en" : "ja";

export default i18n;
