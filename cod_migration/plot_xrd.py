#!/usr/bin/env python3
"""
Plot theoretical XRD pattern for a COD entry.

Usage:
    python plot_xrd.py 1006056
    python plot_xrd.py 1006056 --fwhm 0.1
    python plot_xrd.py 1006056 --save xrd.png
"""
import argparse
import json
import math
import os
import sys
from pathlib import Path

from dotenv import load_dotenv
load_dotenv(Path(__file__).parent / '.env')

import psycopg2
import psycopg2.extras
import matplotlib.pyplot as plt
import numpy as np


def get_connection():
    return psycopg2.connect(
        host=os.environ.get('PG_HOST', 'localhost'),
        port=int(os.environ.get('PG_PORT', 5432)),
        dbname=os.environ.get('PG_DB', 'cod'),
        user=os.environ.get('PG_USER', 'cod_admin'),
        password=os.environ.get('PG_PASSWORD', ''),
    )


def fetch_pattern(conn, cod_id: int) -> tuple[list, dict]:
    with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
        cur.execute("""
            SELECT rp.reflections, rp.wavelength, rp.rad_symbol,
                   rp.has_intensities, rp.two_theta_min, rp.two_theta_max,
                   d.formula, d.mineral, d."sgNumber", d.sg
            FROM xrd_analysis.reference_patterns rp
            JOIN data d ON d.file = rp.cod_id
            WHERE rp.cod_id = %s
            ORDER BY rp.calculated_at DESC
            LIMIT 1
        """, (cod_id,))
        row = cur.fetchone()
    if not row:
        print(f'COD {cod_id} not in reference_patterns. Run xrd_schema_setup.py first.')
        sys.exit(1)
    reflections = row['reflections']
    if isinstance(reflections, str):
        reflections = json.loads(reflections)
    return reflections, dict(row)


def gaussian_profile(two_theta_grid, peaks_tt, peaks_I, fwhm):
    """Convolve stick pattern with Gaussian peaks."""
    sigma = fwhm / (2 * math.sqrt(2 * math.log(2)))
    y = np.zeros_like(two_theta_grid)
    for tt, I in zip(peaks_tt, peaks_I):
        y += I * np.exp(-0.5 * ((two_theta_grid - tt) / sigma) ** 2)
    # Renormalise to 100
    if y.max() > 0:
        y = 100.0 * y / y.max()
    return y


def plot(cod_id: int, fwhm: float, save: str | None):
    conn = get_connection()
    reflections, meta = fetch_pattern(conn, cod_id)
    conn.close()

    peaks_tt = [r['two_theta'] for r in reflections]
    peaks_I  = [r['intensity_rel'] if r['intensity_rel'] is not None else 1.0
                for r in reflections]
    hkl      = [f"{r['h']}{r['k']}{r['l']}" for r in reflections]

    tt_min = meta['two_theta_min']
    tt_max = meta['two_theta_max']
    grid   = np.linspace(tt_min, tt_max, 4000)
    profile = gaussian_profile(grid, peaks_tt, peaks_I, fwhm)

    has_I    = meta['has_intensities']
    rad      = meta['rad_symbol'] or f"{meta['wavelength']:.5f} Å"
    formula  = meta['formula'] or ''
    mineral  = meta['mineral'] or ''
    sg       = meta['sg'] or f"SG {meta['sgNumber']}"
    title    = f"COD {cod_id}  {formula}  {mineral}  [{sg}]"
    subtitle = f"{rad}   {'intensidades calculadas' if has_I else 'somente posições (sem atomic_sites)'}"

    fig, ax = plt.subplots(figsize=(12, 5))

    # Profile
    ax.plot(grid, profile, color='steelblue', lw=1.2, label='perfil')

    # Sticks below the profile
    for tt, I in zip(peaks_tt, peaks_I):
        ax.vlines(tt, -12, -2, color='dimgray', lw=0.8, alpha=0.7)

    ax.set_xlim(tt_min, tt_max)
    ax.set_ylim(-15, 110)
    ax.set_xlabel('2θ (°)', fontsize=12)
    ax.set_ylabel('Intensidade relativa', fontsize=12)
    ax.set_title(f'{title}\n{subtitle}', fontsize=11)
    ax.axhline(0, color='black', lw=0.5)
    ax.grid(axis='x', ls='--', alpha=0.3)

    # Label top-10 peaks by intensity
    top = sorted(zip(peaks_I, peaks_tt, hkl), reverse=True)[:10]
    for I, tt, h in top:
        ax.text(tt, I + 3, h, ha='center', va='bottom', fontsize=6.5, color='navy')

    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150)
        print(f'Saved: {save}')
    else:
        plt.show()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('cod_id', type=int)
    p.add_argument('--fwhm', type=float, default=0.15,
                   help='Peak width FWHM in degrees 2θ (default 0.15)')
    p.add_argument('--save', metavar='FILE',
                   help='Save to file instead of showing (png/pdf/svg)')
    args = p.parse_args()
    plot(args.cod_id, args.fwhm, args.save)


if __name__ == '__main__':
    main()
