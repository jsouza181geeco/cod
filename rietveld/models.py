from dataclasses import dataclass, field


@dataclass
class CandidateInput:
    cod_id: int
    reflections: list[dict]  # [{h,k,l,two_theta,intensity_rel,d_hkl,multiplicity,F_sq}]
    peak_matches: int = 0


@dataclass
class StructureMetadata:
    cod_id: int
    formula: str | None = None
    mineral: str | None = None
    chemname: str | None = None
    sg_number: int | None = None
    sg_symbol: str | None = None
    a: float | None = None
    b: float | None = None
    c: float | None = None
    alpha: float | None = None
    beta: float | None = None
    gamma: float | None = None
    Z: int | None = None
    wavelength: float | None = None
    rad_symbol: str | None = None
    has_intensities: bool | None = None
    authors: str | None = None
    title: str | None = None
    journal: str | None = None
    year: int | None = None
    doi: str | None = None
    method: str | None = None
    status: str | None = None


@dataclass
class CandidateResult:
    cod_id: int
    Rwp: float
    Rp: float
    Rexp: float
    chi2: float
    scale: float
    n_peaks_used: int = 0
    metadata: StructureMetadata | None = None


@dataclass
class RietveldResult:
    xye_file: str
    n_points: int
    candidates: list[CandidateResult] = field(default_factory=list)

    def best(self) -> CandidateResult:
        return self.candidates[0]

    def viable(self, rwp_max: float = 0.15, chi2_max: float = 3.0) -> list[CandidateResult]:
        return [c for c in self.candidates if c.Rwp < rwp_max and c.chi2 < chi2_max]


@dataclass
class PhaseFraction:
    cod_id: int | None
    scale: float              # fitted scale on the ABSOLUTE intensity basis
    weight_pct: float         # Hill-Howard weight fraction (crystalline basis)
    Z: float = 0.0
    M: float = 0.0            # molar mass [g/mol]
    V: float = 0.0            # cell volume [Å³]
    ZMV: float = 0.0
    metadata: StructureMetadata | None = None


@dataclass
class MultiPhaseResult:
    xye_file: str
    n_points: int
    phases: list[PhaseFraction] = field(default_factory=list)
    Rwp: float = 0.0          # combined-model figures of merit
    Rp: float = 0.0
    Rexp: float = 0.0
    chi2: float = 0.0
    Rwp_single_best: float = 0.0   # best single-phase Rwp, for comparison [crit. 6.6]

    def dominant(self) -> PhaseFraction:
        return max(self.phases, key=lambda p: p.weight_pct)
