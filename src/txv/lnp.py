"""Lipid nanoparticle formulation arithmetic.

The quantity that ties the mRNA to the lipid is the **N/P ratio**: moles of
ionisable nitrogen in the ionisable lipid divided by moles of phosphate in the
RNA backbone. One nucleotide contributes exactly one phosphate, so N/P is
computable from sequence length -- and this module computes the RNA molecular
weight from the actual base composition rather than a rule-of-thumb average,
because for a uridine-depleted, m1-pseudouridine-substituted transcript the
rule of thumb drifts by a few percent.

Everything here is stoichiometry and dilution arithmetic. It does not predict
particle size, polydispersity or encapsulation efficiency -- those are
measured, and :class:`FormulationResult` has fields to record them against the
batch.

Molecular weights are nominal values for the free base. PEG-lipids are
polydisperse polymers, so their "MW" is a number-average and the derived masses
carry that uncertainty.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .seqops import clean

# ---------------------------------------------------------------------------
# RNA mass
# ---------------------------------------------------------------------------

#: Residue masses (Da) for internal ribonucleotide monophosphates, i.e. the
#: nucleoside monophosphate minus the water lost on phosphodiester formation.
RNA_RESIDUE_MASS = {"A": 329.21, "G": 345.21, "C": 305.18, "U": 306.17}
_WATER = 18.02
_HPO3 = 79.98
#: N1-methylpseudouridine is an isomer of uridine plus a methyl group.
M1PSI_DELTA = 14.03


def rna_molecular_weight(
    sequence: str,
    five_prime: str = "ppp",
    modified_u_delta: float = 0.0,
) -> float:
    """Molecular weight (Da) of a single-stranded RNA.

    ``sequence`` may be given in DNA or RNA alphabet. ``five_prime`` is
    ``"ppp"`` (uncapped IVT product), ``"p"`` or ``"oh"``. A 5' cap adds only
    a few hundred Da, which is under 0.1% for a transcript of any realistic
    length, so it is not modelled. ``modified_u_delta`` adds a per-uridine mass
    offset -- pass :data:`M1PSI_DELTA` for full m1-pseudouridine substitution.
    """
    rna = clean(sequence).replace("T", "U")
    mass = sum(RNA_RESIDUE_MASS[b] for b in rna) + _WATER
    mass += rna.count("U") * modified_u_delta
    if five_prime == "ppp":
        mass += 2 * _HPO3
    elif five_prime == "oh":
        mass -= _HPO3
    elif five_prime != "p":
        raise ValueError("five_prime must be 'ppp', 'p' or 'oh'")
    return mass


def phosphate_moles(mass_ug: float, sequence: str | None = None,
                    nucleotides: int | None = None, **mw_kwargs) -> float:
    """Moles of backbone phosphate in ``mass_ug`` micrograms of RNA.

    One phosphate per nucleotide. Give either the ``sequence`` (exact) or a
    nucleotide count with an assumed average residue mass (approximate).
    """
    if sequence is not None:
        mw = rna_molecular_weight(sequence, **mw_kwargs)
        n_nt = len(clean(sequence))
    elif nucleotides is not None:
        n_nt = nucleotides
        mw = nucleotides * 330.0
    else:
        raise ValueError("give sequence or nucleotides")
    moles_rna = (mass_ug * 1e-6) / mw
    return moles_rna * n_nt


# ---------------------------------------------------------------------------
# Lipids
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Lipid:
    name: str
    mw: float
    role: str  # ionizable | helper | sterol | peg
    #: Ionisable nitrogens per molecule; only meaningful for the ionisable lipid.
    ionizable_nitrogens: int = 0
    note: str = ""


LIPIDS: dict[str, Lipid] = {
    "ALC-0315": Lipid("ALC-0315", 766.3, "ionizable", 1,
                      "Ionisable aminolipid; apparent pKa ~6.1-6.2."),
    "SM-102": Lipid("SM-102", 710.2, "ionizable", 1,
                    "Ionisable aminolipid; apparent pKa ~6.6-6.8."),
    "DLin-MC3-DMA": Lipid("DLin-MC3-DMA", 642.1, "ionizable", 1,
                          "Ionisable lipid; apparent pKa ~6.4. Non-biodegradable "
                          "linker -- long tissue residence relative to the esters."),
    "DSPC": Lipid("DSPC", 790.1, "helper", 0,
                  "Saturated phosphatidylcholine; structural helper lipid."),
    "DOPE": Lipid("DOPE", 744.0, "helper", 0,
                  "Fusogenic helper lipid; favours the hexagonal phase and can "
                  "improve endosomal escape at the cost of particle stability."),
    "cholesterol": Lipid("cholesterol", 386.65, "sterol", 0,
                         "Packing and stability."),
    "ALC-0159": Lipid("ALC-0159", 2426.8, "peg", 0,
                      "PEG2000 lipid; polydisperse, MW is a number-average."),
    "DMG-PEG2000": Lipid("DMG-PEG2000", 2509.2, "peg", 0,
                         "PEG2000 lipid; polydisperse, MW is a number-average."),
}


@dataclass(frozen=True)
class LipidComposition:
    """Molar percentages of the four lipid components. Must sum to ~100."""

    name: str
    molar_ratios: dict[str, float]
    note: str = ""

    def __post_init__(self) -> None:
        total = sum(self.molar_ratios.values())
        if abs(total - 100.0) > 0.5:
            raise ValueError(
                f"{self.name}: molar ratios sum to {total}, expected 100"
            )
        unknown = sorted(set(self.molar_ratios) - set(LIPIDS))
        if unknown:
            raise ValueError(f"{self.name}: unknown lipids {unknown}")

    @property
    def ionizable(self) -> Lipid:
        for name in self.molar_ratios:
            if LIPIDS[name].role == "ionizable":
                return LIPIDS[name]
        raise ValueError(f"{self.name}: no ionisable lipid in the composition")


COMPOSITIONS: dict[str, LipidComposition] = {
    "sm102_standard": LipidComposition(
        "sm102_standard",
        {"SM-102": 50.0, "DSPC": 10.0, "cholesterol": 38.5, "DMG-PEG2000": 1.5},
        "The 50:10:38.5:1.5 ionisable:helper:sterol:PEG composition.",
    ),
    "alc0315_standard": LipidComposition(
        "alc0315_standard",
        {"ALC-0315": 46.3, "DSPC": 9.4, "cholesterol": 42.7, "ALC-0159": 1.6},
        "The 46.3:9.4:42.7:1.6 composition.",
    ),
    "mc3_standard": LipidComposition(
        "mc3_standard",
        {"DLin-MC3-DMA": 50.0, "DSPC": 10.0, "cholesterol": 38.5, "DMG-PEG2000": 1.5},
        "MC3 benchmark composition; useful as a formulation control.",
    ),
}


@dataclass
class LipidAmount:
    lipid: str
    role: str
    molar_percent: float
    micromoles: float
    micrograms: float
    stock_mg_per_ml: float | None = None
    stock_volume_ul: float | None = None


@dataclass
class FormulationResult:
    construct: str
    mrna_micrograms: float
    mrna_nanomoles: float
    mrna_mw: float
    np_ratio: float
    phosphate_micromoles: float
    composition: str
    lipids: list[LipidAmount]
    total_lipid_micrograms: float
    lipid_to_mrna_mass_ratio: float
    organic_volume_ul: float | None
    aqueous_volume_ul: float | None
    flow_rate_ratio: float | None
    #: Measured after the fact; recorded here so the batch record is one object.
    measured_size_nm: float | None = None
    measured_pdi: float | None = None
    measured_encapsulation: float | None = None
    notes: list[str] = field(default_factory=list)

    def to_text(self) -> str:
        lines = [
            f"LNP formulation: {self.construct}",
            f"  mRNA          {self.mrna_micrograms:.1f} ug "
            f"({self.mrna_nanomoles:.3f} nmol, MW {self.mrna_mw / 1000:.1f} kDa)",
            f"  N/P           {self.np_ratio:.1f} "
            f"({self.phosphate_micromoles:.3f} umol phosphate)",
            f"  composition   {self.composition}",
            f"  {'lipid':<14}{'mol%':>7}{'umol':>10}{'ug':>10}{'stock uL':>11}",
        ]
        for amount in self.lipids:
            stock = f"{amount.stock_volume_ul:.1f}" if amount.stock_volume_ul else "-"
            lines.append(
                f"  {amount.lipid:<14}{amount.molar_percent:>7.1f}"
                f"{amount.micromoles:>10.4f}{amount.micrograms:>10.1f}{stock:>11}"
            )
        lines.append(
            f"  total lipid   {self.total_lipid_micrograms:.1f} ug "
            f"(lipid:mRNA {self.lipid_to_mrna_mass_ratio:.1f}:1 w/w)"
        )
        if self.organic_volume_ul:
            lines.append(
                f"  mixing        {self.aqueous_volume_ul:.0f} uL aqueous : "
                f"{self.organic_volume_ul:.0f} uL ethanol "
                f"(FRR {self.flow_rate_ratio:.0f}:1)"
            )
        for note in self.notes:
            lines.append(f"  note: {note}")
        return "\n".join(lines)

    def to_dict(self) -> dict[str, object]:
        return {
            "construct": self.construct,
            "mrna_ug": self.mrna_micrograms,
            "mrna_nmol": self.mrna_nanomoles,
            "mrna_mw_da": self.mrna_mw,
            "np_ratio": self.np_ratio,
            "composition": self.composition,
            "lipids": [
                {
                    "lipid": a.lipid, "role": a.role, "molar_percent": a.molar_percent,
                    "umol": a.micromoles, "ug": a.micrograms,
                    "stock_volume_ul": a.stock_volume_ul,
                }
                for a in self.lipids
            ],
            "total_lipid_ug": self.total_lipid_micrograms,
            "lipid_to_mrna_w_w": self.lipid_to_mrna_mass_ratio,
            "aqueous_volume_ul": self.aqueous_volume_ul,
            "organic_volume_ul": self.organic_volume_ul,
            "measured": {
                "size_nm": self.measured_size_nm,
                "pdi": self.measured_pdi,
                "encapsulation": self.measured_encapsulation,
            },
        }


def formulate(
    sequence: str,
    mrna_micrograms: float,
    np_ratio: float = 6.0,
    composition: LipidComposition | str = "sm102_standard",
    construct_name: str = "construct",
    modified_u_delta: float = M1PSI_DELTA,
    five_prime: str = "ppp",
    lipid_stocks_mg_per_ml: dict[str, float] | None = None,
    flow_rate_ratio: float | None = 3.0,
) -> FormulationResult:
    """Compute lipid masses and mixing volumes for one LNP batch.

    ``np_ratio`` of about 6 is the usual starting point for mRNA LNPs; it is a
    formulation variable worth titrating (typically 3-9) alongside the PEG
    molar percentage, which is what actually sets particle size.

    ``lipid_stocks_mg_per_ml`` turns the masses into pipetting volumes from
    your ethanolic stocks. ``flow_rate_ratio`` is aqueous:organic for
    microfluidic mixing; the aqueous phase should be the acidic buffer
    (typically citrate, pH ~4) that keeps the ionisable lipid protonated during
    particle formation.
    """
    comp = COMPOSITIONS[composition] if isinstance(composition, str) else composition
    mw = rna_molecular_weight(sequence, five_prime=five_prime,
                              modified_u_delta=modified_u_delta)
    n_nt = len(clean(sequence))
    mrna_moles = (mrna_micrograms * 1e-6) / mw
    p_moles = mrna_moles * n_nt
    p_umol = p_moles * 1e6

    ionizable = comp.ionizable
    if ionizable.ionizable_nitrogens < 1:
        raise ValueError(f"{ionizable.name} declares no ionisable nitrogen")
    ionizable_umol = np_ratio * p_umol / ionizable.ionizable_nitrogens
    ionizable_pct = comp.molar_ratios[ionizable.name]
    total_lipid_umol = ionizable_umol * 100.0 / ionizable_pct

    stocks = lipid_stocks_mg_per_ml or {}
    amounts: list[LipidAmount] = []
    for name, pct in comp.molar_ratios.items():
        lipid = LIPIDS[name]
        umol = total_lipid_umol * pct / 100.0
        ug = umol * lipid.mw
        stock = stocks.get(name)
        # ug / (mg/mL) = ug / (ug/uL) = uL
        volume = (ug / (stock * 1000.0)) * 1000.0 if stock else None
        amounts.append(LipidAmount(name, lipid.role, pct, umol, ug, stock, volume))

    total_ug = sum(a.micrograms for a in amounts)

    organic_ul = aqueous_ul = None
    if stocks and all(a.stock_volume_ul is not None for a in amounts):
        organic_ul = sum(a.stock_volume_ul or 0.0 for a in amounts)
        if flow_rate_ratio:
            aqueous_ul = organic_ul * flow_rate_ratio

    notes = [
        "N/P is computed from the exact base composition of the supplied sequence.",
        "Aqueous phase should be acidic (e.g. citrate pH ~4) so the ionisable "
        "lipid is protonated during mixing; buffer-exchange to a neutral, "
        "cryoprotected buffer afterwards.",
    ]
    if modified_u_delta:
        notes.append(
            f"RNA MW assumes full uridine substitution (+{modified_u_delta} Da per U)."
        )
    if any(LIPIDS[a.lipid].role == "peg" for a in amounts):
        notes.append(
            "PEG-lipid MW is a number-average of a polydisperse polymer; its "
            "molar percentage is the main lever on particle size."
        )

    return FormulationResult(
        construct=construct_name,
        mrna_micrograms=mrna_micrograms,
        mrna_nanomoles=mrna_moles * 1e9,
        mrna_mw=mw,
        np_ratio=np_ratio,
        phosphate_micromoles=p_umol,
        composition=comp.name,
        lipids=amounts,
        total_lipid_micrograms=total_ug,
        lipid_to_mrna_mass_ratio=total_ug / mrna_micrograms,
        organic_volume_ul=organic_ul,
        aqueous_volume_ul=aqueous_ul,
        flow_rate_ratio=flow_rate_ratio if organic_ul else None,
        notes=notes,
    )


def dose_series(
    mg_per_kg: float,
    body_weights_kg: list[float],
    concentration_ug_per_ul: float,
) -> list[dict[str, float]]:
    """Per-subject dose volumes for a given mg/kg and drug-product strength."""
    rows = []
    for weight in body_weights_kg:
        dose_ug = mg_per_kg * weight * 1000.0
        rows.append({
            "body_weight_kg": weight,
            "dose_ug": dose_ug,
            "volume_ul": dose_ug / concentration_ug_per_ul,
        })
    return rows


__all__ = [
    "Lipid", "LIPIDS", "LipidComposition", "COMPOSITIONS", "LipidAmount",
    "FormulationResult", "formulate", "rna_molecular_weight", "phosphate_moles",
    "dose_series", "M1PSI_DELTA",
]
