import argparse
import sys

from data_loader import load_candidates_csv
from db_client import DBClient
from pipeline import run_pipeline, run_multiphase, prefilter_candidates, candidates_from_db


def main():
    p = argparse.ArgumentParser(description='Rietveld Phase Identification + QPA')
    p.add_argument('xye',              help='Arquivo .xye experimental')
    p.add_argument('candidates_csv',   nargs='?', default=None,
                   help='CSV com cod_id, peak_matches, reflections '
                        '(omitir se usar --from-db)')
    p.add_argument('--from-db',  action='store_true',
                   help='Modo DB-only: descobre candidatos do COD via MV peak_fingerprints '
                        '(detecta picos -> Hanawalt MV -> fetch_reflections). '
                        'Nao precisa de candidates_csv. Exige DB ativo (incompativel com --no-db).')
    p.add_argument('--no-db',    action='store_true',
                   help='Nao conectar ao PostgreSQL (sem metadados)')
    p.add_argument('--top',      type=int, default=5,
                   help='Numero de candidatos a exibir (default: 5)')
    p.add_argument('--phases',   type=int, default=None, metavar='K',
                   help='Modo multi-fase: ajusta as K melhores fases e quantifica (QPA). Exige DB.')
    p.add_argument('--prefilter', type=int, default=None, metavar='N',
                   help='Em modo CSV: pre-filtra top-N por matches de pico (Hanawalt) antes do '
                        'Rietveld. Em modo --from-db: controla top_n da busca na MV (default 50).')
    p.add_argument('--wavelength', type=float, default=1.54056, metavar='LAMBDA',
                   help='Comprimento de onda da amostra em Angstrom (default: 1.54056 = CuKa1). '
                        'Usar CoKa=1.7902 ou CrKa=2.2909 para amostras com Fe.')
    p.add_argument('--plot',     action='store_true',
                   help='Plot (melhor candidato; ou pizza QPA + decomposicao em --phases)')
    p.add_argument('--env',      default=None,
                   help='Path para arquivo .env com credenciais do DB')
    args = p.parse_args()

    # --- Argument validation ---
    if args.from_db and args.no_db:
        p.error('--from-db e --no-db sao incompativeis: DB-only mode requer conexao ativa.')
    if not args.from_db and args.candidates_csv is None:
        p.error('candidates_csv obrigatorio sem --from-db.')
    if args.from_db and args.candidates_csv is not None:
        print('[WARN] --from-db ativo: candidates_csv ignorado.', file=sys.stderr)

    # QPA (--phases) needs metadata (Z, M, V) → DB required
    if args.phases is not None and args.no_db:
        print('[WARN] --phases exige DB para quantificacao (Z, M, V); '
              'com --no-db as fracoes em peso serao 0%.', file=sys.stderr)

    db = None
    if not args.no_db:
        try:
            db = DBClient(env_path=args.env)
        except Exception as e:
            print(f"[WARN] DB indisponível: {e}", file=sys.stderr)
            if args.from_db:
                p.error(f'--from-db requer DB ativo. Falha: {e}')
            print('[WARN] Continuando sem metadados (use --no-db para suprimir este aviso)',
                  file=sys.stderr)

    try:
        # --- Candidate loading ---
        if args.from_db:
            top_n = args.prefilter if args.prefilter is not None else 50
            candidates = candidates_from_db(
                args.xye, db, top_n=top_n, wavelength=args.wavelength,
            )
            print(f'[DB] {len(candidates)} candidatos carregados do COD '
                  f'(top-{top_n} por n_matches Hanawalt)',
                  file=sys.stderr)
        else:
            candidates = load_candidates_csv(args.candidates_csv)
            if args.prefilter is not None:
                n_before = len(candidates)
                candidates = prefilter_candidates(
                    args.xye, candidates, top_n=args.prefilter,
                    db_client=db, wavelength=args.wavelength,
                )
                print(f'[Hanawalt] {n_before} candidatos → {len(candidates)} '
                      f'(top-{args.prefilter} por n_matches)',
                      file=sys.stderr)
                if not candidates:
                    print('[WARN] Nenhum candidato sobreviveu ao pré-filtro; '
                          'tente --prefilter maior ou reduzir min_matches.',
                          file=sys.stderr)

        if args.phases is not None:
            mp = run_multiphase(args.xye, candidates, db_client=db, top_k=args.phases)
        else:
            result = run_pipeline(args.xye, candidates, db_client=db)
    finally:
        if db:
            db.close()

    if args.phases is not None:
        print_qpa(mp)
        if args.plot:
            plot_qpa(mp, args.xye, candidates)
    else:
        print_results(result, top=args.top)
        if args.plot:
            plot_best(result, args.xye, candidates)


def print_qpa(mp):
    print(f'\nArquivo : {mp.xye_file}')
    print(f'Pontos  : {mp.n_points}')
    print(f'Modo    : multi-fase ({len(mp.phases)} fases)')
    print(f'Rwp combinado = {mp.Rwp:.4f}   chi2 = {mp.chi2:.2f}   '
          f'(Rwp melhor fase unica = {mp.Rwp_single_best:.4f})\n')

    header = (f"{'cod_id':>10}  {'scale':>11}  {'Z':>4}  {'M(g/mol)':>10}  "
              f"{'V(A^3)':>10}  {'wt%':>7}  {'has_I':>5}  formula")
    print(header)
    print('-' * len(header))
    for ph in sorted(mp.phases, key=lambda p: -p.weight_pct):
        m = ph.metadata
        formula = (m.formula or '?') if m else '-'
        has_i = ('T' if m.has_intensities else 'F') if (m and m.has_intensities is not None) else '?'
        print(f"{str(ph.cod_id):>10}  {ph.scale:>11.4e}  {ph.Z:>4.0f}  {ph.M:>10.2f}  "
              f"{ph.V:>10.2f}  {ph.weight_pct:>7.2f}  {has_i:>5}  {formula}")

    soma = sum(ph.weight_pct for ph in mp.phases)
    print(f"\nSoma wt% = {soma:.2f}  (base 100% cristalina; amorfo nao quantificado)")
    if mp.Rwp_single_best > 0.15 and mp.Rwp < mp.Rwp_single_best:
        print("Nota: melhor fase unica tem Rwp>0.15 — amostra provavelmente e mistura.")


def plot_qpa(mp, xye_path, candidates):
    import matplotlib.pyplot as plt
    from data_loader import parse_xye
    from pattern_calc import build_icalc_unit_absolute

    cand_map = {c.cod_id: c for c in candidates}
    tth, Iobs, sigma = parse_xye(xye_path)

    present = [ph for ph in mp.phases if ph.weight_pct > 0.01]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))

    if present:
        ax1.pie([ph.weight_pct for ph in present],
                labels=[f"{ph.cod_id}\n{ph.weight_pct:.1f}%" for ph in present],
                autopct='%1.1f%%', startangle=90)
        ax1.set_title('Frações em peso (QPA Hill-Howard)')
    else:
        ax1.text(0.5, 0.5, 'Sem wt% (DB inativo)', ha='center', va='center')
        ax1.axis('off')

    ax2.plot(tth, Iobs, 'k-', lw=0.7, label='Iobs')
    for ph in mp.phases:
        if ph.weight_pct > 0.01:
            u = build_icalc_unit_absolute(tth, cand_map[ph.cod_id].reflections)[0]
            ax2.plot(tth, ph.scale * u, lw=0.6, alpha=0.7,
                     label=f'{ph.cod_id} ({ph.weight_pct:.0f}%)')
    ax2.set_xlabel('2 theta (graus)')
    ax2.set_ylabel('Intensidade')
    ax2.legend(fontsize=7)
    ax2.set_title(f'Decomposicao  Rwp={mp.Rwp:.4f}  chi2={mp.chi2:.2f}')
    plt.tight_layout()
    plt.show()


def print_results(result, top=5):
    viable = result.viable()

    print(f'\nArquivo : {result.xye_file}')
    print(f'Pontos  : {result.n_points}')
    print(f'Top {top} de {len(result.candidates)} candidatos\n')

    header = (f"{'#':>3}  {'cod_id':>10}  {'Rwp':>7}  {'Rp':>7}  "
              f"{'chi2':>7}  {'scale':>10}  {'picos':>6}  {'has_I':>5}  "
              f"{'formula':<24}  mineral")
    print(header)
    print('-' * len(header))

    for i, r in enumerate(result.candidates[:top]):
        m = r.metadata
        formula   = (m.formula  or '?') if m else '-'
        mineral   = (m.mineral  or '') if m else '-'
        has_i_str = ('T' if m.has_intensities else 'F') if (m and m.has_intensities is not None) else '?'
        flag = '*' if r in viable else ' '

        # warn if has_intensities=False: intensity_rel missing |F|² -> Rwp unreliable
        warn = ' [!sem F2]' if (m and m.has_intensities is False) else ''

        print(f"{flag}{i+1:>2}  {r.cod_id:>10}  {r.Rwp:>7.4f}  {r.Rp:>7.4f}  "
              f"{r.chi2:>7.2f}  {r.scale:>10.4f}  {r.n_peaks_used:>6}  {has_i_str:>5}  "
              f"{formula:<24}  {mineral}{warn}")

    print()
    if viable:
        print(f"* {len(viable)} candidato(s) viavel(is)  (Rwp<0.15, chi2<3)")
    else:
        print("Nenhum candidato viavel  (Rwp<0.15, chi2<3)")


def plot_best(result, xye_path, candidates):
    import matplotlib.pyplot as plt
    from data_loader import parse_xye
    from pattern_calc import build_icalc_unit
    from linear_fit import linear_fit

    best = result.best()
    cand = next(c for c in candidates if c.cod_id == best.cod_id)

    tth, Iobs, sigma = parse_xye(xye_path)
    Icalc_unit, n_used = build_icalc_unit(tth, cand.reflections)
    _, Icalc = linear_fit(tth, Iobs, sigma, Icalc_unit)
    diff = Iobs - Icalc

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                             gridspec_kw={'height_ratios': [4, 1]})

    axes[0].plot(tth, Iobs,  'k-', lw=0.8, label='Iobs')
    axes[0].plot(tth, Icalc, 'r-', lw=0.8, label=f'Icalc  cod={best.cod_id}', alpha=0.8)
    axes[0].set_ylabel('Intensidade')
    axes[0].legend()

    m = best.metadata
    title = f"cod={best.cod_id}  Rwp={best.Rwp:.4f}  chi2={best.chi2:.2f}  picos={n_used}"
    if m:
        if m.formula:
            title += f"\n{m.formula}"
        if m.mineral:
            title += f"  ({m.mineral})"
        if m.sg_symbol:
            title += f"  SG {m.sg_symbol}"
        if m.has_intensities is False:
            title += "  [sem |F|2]"
    axes[0].set_title(title)

    axes[1].plot(tth, diff, 'g-', lw=0.5)
    axes[1].axhline(0, color='k', lw=0.5)
    axes[1].set_xlabel('2 theta (graus)')
    axes[1].set_ylabel('diff (Iobs-Icalc)')

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
