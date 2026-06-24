"""SVN operations wrapper for COD data retrieval and incremental updates."""
import os
import re
import subprocess
from pathlib import Path


def _run(*args, cwd: Path = None) -> str:
    """Run an SVN command, return stdout. Raises RuntimeError on failure."""
    env = {**os.environ, 'LC_ALL': 'C.UTF-8', 'LANG': 'C.UTF-8'}
    result = subprocess.run(
        ['svn', *args],
        capture_output=True,
        text=True,
        encoding='utf-8',
        errors='replace',
        cwd=str(cwd) if cwd else None,
        env=env,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"svn {' '.join(str(a) for a in args)} failed (exit {result.returncode}):\n"
            f"{result.stderr.strip()}"
        )
    return result.stdout


def get_revision(wc_path: Path) -> int:
    """Return the current SVN revision of a working copy."""
    out = _run('info', '--show-item', 'revision', str(wc_path))
    return int(out.strip())


def checkout_sparse(url: str, local_path: Path, subdirs: list[str]) -> int:
    """
    Sparse SVN checkout: empty root, then expand specified subdirectories to infinity.
    Skips subdirs that already exist in the working copy.
    Returns the revision number after checkout.
    """
    local_path = Path(local_path)

    if not (local_path / '.svn').exists():
        print(f"Sparse checkout {url} → {local_path}")
        _run('checkout', '--depth', 'empty', url, str(local_path))

    for subdir in subdirs:
        target = local_path / subdir
        print(f"Expanding {subdir}/...")
        _run('update', '--set-depth', 'infinity', str(target))

    return get_revision(local_path)


def update(wc_path: Path) -> tuple[int, list[str]]:
    """
    Run `svn update` on the working copy.
    Returns (new_revision, list_of_changed_local_paths).
    Changed paths are absolute strings as reported by SVN.
    """
    out = _run('update', str(wc_path))
    new_rev = get_revision(wc_path)

    changed = []
    for line in out.splitlines():
        # SVN prints: "A  path", "U  path", "D  path", "C  path", "G  path"
        m = re.match(r'^[AUDCG]\s+(.+)$', line)
        if m:
            changed.append(m.group(1).strip())

    return new_rev, changed


def get_changed_files(wc_path: Path, from_rev: int, to_rev: int) -> list[str]:
    """
    List files changed between two revisions using `svn diff --summarize`.
    Returns list of local file paths (absolute).
    """
    out = _run(
        'diff', '--summarize',
        f'--revision={from_rev}:{to_rev}',
        str(wc_path),
    )
    changed = []
    for line in out.splitlines():
        m = re.match(r'^[AMD]\s+(.+)$', line)
        if m:
            changed.append(m.group(1).strip())
    return changed
