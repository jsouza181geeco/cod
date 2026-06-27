from pathlib import Path
import json
import re
import numpy as np
import pandas as pd

from models import CandidateInput


def parse_xye(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = Path(path)
    if not path.exists():
        raise ValueError(f"File not found: {path}")

    rows = []
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                rows.append([float(p) for p in parts[:3]])
            except ValueError:
                continue

    if len(rows) < 10:
        raise ValueError(f"Too few valid data points ({len(rows)}) in {path}")

    data = np.array(rows)
    data = data[np.argsort(data[:, 0])]
    tth = data[:, 0]
    Iobs = data[:, 1]

    if data.shape[1] >= 3:
        sigma = data[:, 2]
    else:
        sigma = np.sqrt(np.maximum(Iobs, 1.0))

    bad = sigma <= 0
    if bad.any():
        sigma[bad] = np.sqrt(np.maximum(Iobs[bad], 1.0))

    return tth, Iobs, sigma


def parse_asc(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parse Bruker/Siemens/PANalytical .ASC diffractogram (2-column: 2theta, counts).

    No sigma column in .asc format — computed as sqrt(max(Iobs, 1)).
    Skips non-numeric header lines. Handles European comma decimal separator.
    """
    path = Path(path)
    if not path.exists():
        raise ValueError(f"File not found: {path}")

    rows = []
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            line = line.replace(',', '.')
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                rows.append([float(parts[0]), float(parts[1])])
            except ValueError:
                continue

    if len(rows) < 10:
        raise ValueError(f"Too few valid data points ({len(rows)}) in {path}")

    data = np.array(rows)
    data = data[np.argsort(data[:, 0])]
    tth = data[:, 0]
    Iobs = data[:, 1]
    sigma = np.sqrt(np.maximum(Iobs, 1.0))
    return tth, Iobs, sigma


def parse_diffractogram(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Dispatch to parse_xye or parse_asc by file extension (case-insensitive)."""
    p = Path(path)
    if p.suffix.lower() == '.asc':
        return parse_asc(p)
    return parse_xye(p)


def load_candidates_csv(path: str | Path) -> list[CandidateInput]:
    path = Path(path)
    df = pd.read_csv(path)

    if df.empty:
        raise ValueError(f"Empty CSV: {path}")
    missing = {'cod_id', 'reflections'} - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns {missing} in {path}")

    def _parse_reflections(s: str) -> list[dict]:
        try:
            return json.loads(s)
        except json.JSONDecodeError:
            # fix spurious spaces inside numbers (data quality issue in some CSV rows)
            fixed = re.sub(r'(\d) (\d)', r'\1\2', s)
            return json.loads(fixed)

    candidates = []
    for _, row in df.iterrows():
        reflections = _parse_reflections(str(row['reflections']))
        peak_matches = int(row['peak_matches']) if 'peak_matches' in df.columns else 0
        candidates.append(CandidateInput(
            cod_id=int(row['cod_id']),
            reflections=reflections,
            peak_matches=peak_matches,
        ))
    return candidates


if __name__ == '__main__':
    import sys
    import matplotlib.pyplot as plt

    path = sys.argv[1] if len(sys.argv) > 1 else 'Cu_synthetic.xye'
    tth, Iobs, sigma = parse_diffractogram(path)

    print(f"Arquivo : {path}")
    print(f"Pontos  : {len(tth)}")
    print(f"2θ range: {tth[0]:.2f}° — {tth[-1]:.2f}°")
    print(f"I max   : {Iobs.max():.1f}  em 2θ = {tth[Iobs.argmax()]:.3f}°")
    print(f"σ range : {sigma.min():.2f} — {sigma.max():.2f}")

    plt.figure(figsize=(12, 4))
    plt.plot(tth, Iobs, 'k-', lw=0.8, label='Iobs')
    plt.fill_between(tth, Iobs - sigma, Iobs + sigma, alpha=0.2, color='gray')
    plt.xlabel('2θ (graus)')
    plt.ylabel('Intensidade')
    plt.title(f'Difratograma experimental — {path}')
    plt.legend()
    plt.tight_layout()
    plt.show()

    csv_path = sys.argv[2] if len(sys.argv) > 2 else 'data-1782394014136.csv'
    candidates = load_candidates_csv(csv_path)

    print(f"\n{'cod_id':>10}  {'peak_matches':>12}  {'n_reflections':>13}  {'2θ min':>8}  {'2θ max':>8}")
    print('-' * 58)
    for c in candidates:
        tth_vals = [r['two_theta'] for r in c.reflections]
        print(f"{c.cod_id:>10}  {c.peak_matches:>12}  {len(c.reflections):>13}  "
              f"{min(tth_vals):>8.2f}  {max(tth_vals):>8.2f}")
