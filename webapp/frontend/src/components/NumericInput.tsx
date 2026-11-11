import { useEffect, useRef, useState, type InputHTMLAttributes } from "react";

interface Props extends Omit<InputHTMLAttributes<HTMLInputElement>, "type" | "value" | "onChange"> {
  value: number | null | undefined;
  onValueChange: (value: number) => void;
}

/** Keeps the user's draft (including an empty value or trailing decimal) intact until blur. */
export function NumericInput({ value, onValueChange, onBlur, onFocus, min, max, ...props }: Props) {
  const minimum = min == null ? null : Number(min);
  const maximum = max == null ? null : Number(max);
  const initialValue = value != null && maximum !== null && Number.isFinite(maximum)
    ? Math.min(value, maximum)
    : value;
  const [draft, setDraft] = useState(initialValue == null ? "" : String(initialValue));
  const focused = useRef(false);

  const clamp = (next: number) => {
    const withMinimum = minimum !== null && Number.isFinite(minimum) ? Math.max(minimum, next) : next;
    return maximum !== null && Number.isFinite(maximum) ? Math.min(maximum, withMinimum) : withMinimum;
  };

  useEffect(() => {
    if (!focused.current) setDraft(value == null ? "" : String(value));
  }, [value]);

  useEffect(() => {
    if (value != null && maximum !== null && Number.isFinite(maximum) && value > maximum) {
      setDraft(String(maximum));
      onValueChange(maximum);
    }
  }, [maximum, onValueChange, value]);

  return (
    <input
      {...props}
      type="number"
      min={min}
      max={max}
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
        if (Number.isFinite(parsed)) {
          const capped = maximum !== null && Number.isFinite(maximum) ? Math.min(maximum, parsed) : parsed;
          if (capped !== parsed) setDraft(String(capped));
          onValueChange(capped);
        }
      }}
      onBlur={(event) => {
        focused.current = false;
        const raw = event.target.value.trim();
        if (raw === "") {
          setDraft(value == null ? "" : String(value));
        } else {
          const parsed = Number(raw);
          if (Number.isFinite(parsed)) {
            const normalized = clamp(parsed);
            onValueChange(normalized);
            setDraft(String(normalized));
          } else {
            setDraft(value == null ? "" : String(value));
          }
        }
        onBlur?.(event);
      }}
    />
  );
}
