import { useEffect, useRef, useState, type InputHTMLAttributes } from "react";

interface Props extends Omit<InputHTMLAttributes<HTMLInputElement>, "type" | "value" | "onChange"> {
  value: number | null | undefined;
  onValueChange: (value: number) => void;
}

/** Keeps the user's draft (including an empty value or trailing decimal) intact until blur. */
export function NumericInput({ value, onValueChange, onBlur, onFocus, ...props }: Props) {
  const [draft, setDraft] = useState(value == null ? "" : String(value));
  const focused = useRef(false);

  useEffect(() => {
    if (!focused.current) setDraft(value == null ? "" : String(value));
  }, [value]);

  return (
    <input
      {...props}
      type="number"
      value={draft}
      onFocus={(event) => {
        focused.current = true;
        onFocus?.(event);
      }}
      onChange={(event) => {
        const raw = event.target.value;
        setDraft(raw);
        if (raw.trim() === "") return;
        const parsed = Number(raw);
        if (Number.isFinite(parsed)) onValueChange(parsed);
      }}
      onBlur={(event) => {
        focused.current = false;
        const raw = event.target.value.trim();
        if (raw === "") {
          setDraft(value == null ? "" : String(value));
        } else {
          const parsed = Number(raw);
          if (Number.isFinite(parsed)) {
            onValueChange(parsed);
            setDraft(String(parsed));
          } else {
            setDraft(value == null ? "" : String(value));
          }
        }
        onBlur?.(event);
      }}
    />
  );
}
