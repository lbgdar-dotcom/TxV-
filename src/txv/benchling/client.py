"""Push designed constructs into a Benchling registry.

The Benchling SDK is an *optional* dependency (``pip install txv[benchling]``).
Everything in this module works without it via :class:`DryRunRegistrar`, which
produces exactly the payload that would be sent and writes it to disk. Use the
dry run to review a batch -- and in CI, where you do not want a test run
minting registry IDs.

Authentication comes from the environment so that credentials never enter a
design file or a commit:

``BENCHLING_TENANT``
    e.g. ``https://yourcompany.benchling.com``
``BENCHLING_API_KEY``
    A personal API key, or use ``BENCHLING_CLIENT_ID`` / ``BENCHLING_CLIENT_SECRET``
    for an app's OAuth client-credentials flow.
``BENCHLING_FOLDER_ID``
    Destination folder (``lib_...``).
``BENCHLING_SCHEMA_ID``
    DNA sequence schema (``ts_...``) whose fields you want populated.
``BENCHLING_REGISTRY_ID``
    Registry to register entities into (``src_...``); omit to create unregistered.
``BENCHLING_PROJECT_ID``
    Optional project for notebook entries.
"""

from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from ..constructs import Construct
from ..export import to_genbank
from ..qc import QCReport

#: txv feature kind -> Benchling annotation type string.
ANNOTATION_TYPE = {
    "promoter": "promoter",
    "utr5": "misc_feature",
    "utr3": "misc_feature",
    "kozak": "misc_feature",
    "orf": "CDS",
    "signal_peptide": "sig_peptide",
    "trafficking": "misc_feature",
    "linker": "misc_feature",
    "neoepitope": "misc_feature",
    "tumor_associated": "misc_feature",
    "full_length": "misc_feature",
    "helper": "misc_feature",
    "start": "misc_feature",
    "stop": "terminator",
    "polya": "misc_feature",
    "linearization": "misc_feature",
}


@dataclass
class BenchlingConfig:
    tenant: str | None = None
    api_key: str | None = None
    client_id: str | None = None
    client_secret: str | None = None
    folder_id: str | None = None
    schema_id: str | None = None
    registry_id: str | None = None
    project_id: str | None = None
    naming_strategy: str = "NEW_IDS"
    #: Extra schema fields applied to every construct in the batch.
    default_fields: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_env(cls, **overrides: Any) -> "BenchlingConfig":
        config = cls(
            tenant=os.environ.get("BENCHLING_TENANT"),
            api_key=os.environ.get("BENCHLING_API_KEY"),
            client_id=os.environ.get("BENCHLING_CLIENT_ID"),
            client_secret=os.environ.get("BENCHLING_CLIENT_SECRET"),
            folder_id=os.environ.get("BENCHLING_FOLDER_ID"),
            schema_id=os.environ.get("BENCHLING_SCHEMA_ID"),
            registry_id=os.environ.get("BENCHLING_REGISTRY_ID"),
            project_id=os.environ.get("BENCHLING_PROJECT_ID"),
        )
        for key, value in overrides.items():
            if value is not None:
                setattr(config, key, value)
        return config

    @property
    def has_credentials(self) -> bool:
        return bool(self.tenant and (self.api_key or (self.client_id and self.client_secret)))

    def missing(self) -> list[str]:
        gaps = []
        if not self.tenant:
            gaps.append("BENCHLING_TENANT")
        if not (self.api_key or (self.client_id and self.client_secret)):
            gaps.append("BENCHLING_API_KEY or BENCHLING_CLIENT_ID+BENCHLING_CLIENT_SECRET")
        if not self.folder_id:
            gaps.append("BENCHLING_FOLDER_ID")
        return gaps


def construct_payload(
    construct: Construct,
    config: BenchlingConfig,
    qc: QCReport | None = None,
    extra_fields: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build the DNA-sequence creation payload for a construct.

    Benchling annotations are 0-based, end-exclusive, with an explicit strand,
    which matches this package's internal :class:`~txv.constructs.Feature`
    convention -- so the coordinates pass through unchanged.
    """
    summary = construct.summary()
    fields: dict[str, Any] = {
        "Construct type": {"value": "IVT mRNA template"},
        "Antigen count": {"value": len(construct.cassette.antigens)},
        "Antigens": {"value": ", ".join(a.name for a in construct.cassette.antigens)},
        "Transcript length (nt)": {"value": construct.transcript_length},
        "ORF length (nt)": {"value": len(construct.orf)},
        "GC fraction": {"value": summary["gc"]},
        "Uridine fraction": {"value": summary["uridine_fraction"]},
        "CAI": {"value": summary["cai"]},
        "Cassette linker": {"value": construct.cassette.linker},
    }
    if construct.placeholders_used:
        fields["Placeholder parts"] = {"value": ", ".join(construct.placeholders_used)}
    if qc is not None:
        fields["QC status"] = {"value": "PASS" if qc.passed else "FAIL"}
        fields["QC failures"] = {"value": len(qc.failures)}
        fields["QC warnings"] = {"value": len(qc.warnings)}
    fields.update({k: {"value": v} for k, v in config.default_fields.items()})
    fields.update({k: {"value": v} for k, v in (extra_fields or {}).items()})

    annotations = [
        {
            "name": feature.name,
            "type": ANNOTATION_TYPE.get(feature.kind, "misc_feature"),
            "start": feature.start,
            "end": feature.end,
            "strand": 1,
            "color": _color_for(feature.kind),
        }
        for feature in construct.features
    ]

    payload: dict[str, Any] = {
        "name": construct.name,
        "bases": construct.template,
        "isCircular": False,
        "annotations": annotations,
        "fields": fields,
        "customNotes": "\n".join(construct.design_notes),
    }
    if config.folder_id:
        payload["folderId"] = config.folder_id
    if config.schema_id:
        payload["schemaId"] = config.schema_id
    if config.registry_id:
        payload["registryId"] = config.registry_id
        payload["namingStrategy"] = config.naming_strategy
    return payload


_COLORS = {
    "promoter": "#B3E0A6", "utr5": "#C9D7F8", "utr3": "#C9D7F8",
    "kozak": "#F7D08A", "orf": "#F2A0A0", "signal_peptide": "#F5C6E0",
    "trafficking": "#D8C3F0", "linker": "#E3E3E3", "neoepitope": "#9BD1E5",
    "tumor_associated": "#7FB3D5", "helper": "#A3D9C9", "full_length": "#7FB3D5",
    "stop": "#9E9E9E", "polya": "#FFE9A8", "linearization": "#CFCFCF",
    "start": "#9E9E9E",
}


def _color_for(kind: str) -> str:
    return _COLORS.get(kind, "#E3E3E3")


def guard_qc(construct: Construct, qc: QCReport | None) -> None:
    """Raise if ``qc`` failed.

    The registry is the record of what the programme believes is real, and a
    failing construct in it is worse than no construct at all. Pass a report
    that passes, or ``None`` if you are deliberately registering a draft.
    """
    if qc is not None and not qc.passed:
        raise ValueError(
            f"{construct.name}: QC failed ({len(qc.failures)} failure(s)); "
            "refusing to register. Failures: "
            + "; ".join(c.name for c in qc.failures)
        )


class DryRunRegistrar:
    """Writes the payloads it would have sent, and sends nothing.

    This is the default whenever credentials are absent, so an accidental run
    in CI or on a laptop produces files rather than registry entries.
    """

    def __init__(self, config: BenchlingConfig, out_dir: str | Path = "benchling_dry_run") -> None:
        self.config = config
        self.out_dir = Path(out_dir)
        self.sent: list[dict[str, Any]] = []

    def register(
        self,
        construct: Construct,
        qc: QCReport | None = None,
        extra_fields: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        payload = construct_payload(construct, self.config, qc, extra_fields)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        base = self.out_dir / f"{construct.name}"
        base.with_suffix(".payload.json").write_text(json.dumps(payload, indent=2))
        base.with_suffix(".gb").write_text(to_genbank(construct))
        if qc is not None:
            base.with_suffix(".qc.json").write_text(json.dumps(qc.to_dict(), indent=2))
        record = {
            "mode": "dry-run",
            "construct": construct.name,
            "written": str(base.with_suffix(".payload.json")),
            "genbank": str(base.with_suffix(".gb")),
            "timestamp": stamp,
        }
        self.sent.append(record)
        return record

    def register_many(
        self, items: Iterable[tuple[Construct, QCReport | None]]
    ) -> list[dict[str, Any]]:
        return [self.register(construct, qc) for construct, qc in items]


class BenchlingRegistrar:
    """Live Benchling client. Requires ``benchling-sdk``."""

    def __init__(self, config: BenchlingConfig) -> None:
        gaps = config.missing()
        if gaps:
            raise ValueError(
                "Benchling configuration incomplete; missing: " + ", ".join(gaps)
            )
        self.config = config
        self._benchling = self._connect()
        self.sent: list[dict[str, Any]] = []

    def _connect(self):
        try:
            from benchling_sdk.auth.api_key_auth import ApiKeyAuth
            from benchling_sdk.auth.client_credentials_oauth2 import (
                ClientCredentialsOAuth2,
            )
            from benchling_sdk.benchling import Benchling
        except ImportError as exc:  # pragma: no cover - exercised only without the SDK
            raise ImportError(
                "benchling-sdk is not installed. Install it with "
                "`pip install 'txv[benchling]'`, or use DryRunRegistrar."
            ) from exc

        if self.config.api_key:
            auth = ApiKeyAuth(self.config.api_key)
        else:
            auth = ClientCredentialsOAuth2(
                client_id=self.config.client_id,
                client_secret=self.config.client_secret,
            )
        return Benchling(url=self.config.tenant, auth_method=auth)

    def register(
        self,
        construct: Construct,
        qc: QCReport | None = None,
        extra_fields: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Create the DNA sequence in Benchling and return its id and URL.

        Refuses to register a construct whose QC failed (see :func:`guard_qc`).
        """
        from benchling_sdk.models import (  # local import: optional dependency
            AutoAnnotateDnaSequences,
            DnaSequenceCreate,
            Fields,
            SequenceFeatureCreate,
        )

        guard_qc(construct, qc)

        payload = construct_payload(construct, self.config, qc, extra_fields)
        annotations = [
            SequenceFeatureCreate(
                name=a["name"], type=a["type"], start=a["start"],
                end=a["end"], strand=a["strand"], color=a["color"],
            )
            for a in payload["annotations"]
        ]
        create = DnaSequenceCreate(
            name=payload["name"],
            bases=payload["bases"],
            is_circular=False,
            folder_id=payload.get("folderId"),
            schema_id=payload.get("schemaId"),
            annotations=annotations,
            fields=Fields.from_dict(payload["fields"]) if payload.get("fields") else None,
            custom_notes=payload.get("customNotes"),
        )
        created = self._benchling.dna_sequences.create(create)

        record: dict[str, Any] = {
            "mode": "live",
            "construct": construct.name,
            "id": created.id,
            "url": getattr(created, "web_url", None),
        }
        if self.config.registry_id:
            record["registered"] = self._register_entity(created.id)
        self.sent.append(record)
        return record

    def _register_entity(self, entity_id: str) -> bool:
        from benchling_sdk.models import NamingStrategy

        strategy = getattr(NamingStrategy, self.config.naming_strategy)
        self._benchling.registry.register(
            registry_id=self.config.registry_id,
            entity_ids=[entity_id],
            naming_strategy=strategy,
        )
        return True

    def register_many(
        self, items: Iterable[tuple[Construct, QCReport | None]]
    ) -> list[dict[str, Any]]:
        return [self.register(construct, qc) for construct, qc in items]


def make_registrar(
    config: BenchlingConfig | None = None,
    dry_run: bool | None = None,
    out_dir: str | Path = "benchling_dry_run",
):
    """Return a live or dry-run registrar.

    With ``dry_run=None`` (the default) it picks: live if credentials are
    present, dry run otherwise. Pass ``dry_run=True`` to force a dry run even
    when credentials exist -- which is what you want when reviewing a batch.
    """
    config = config or BenchlingConfig.from_env()
    if dry_run is None:
        dry_run = not config.has_credentials
    if dry_run:
        return DryRunRegistrar(config, out_dir)
    return BenchlingRegistrar(config)


__all__ = [
    "ANNOTATION_TYPE", "BenchlingConfig", "BenchlingRegistrar",
    "DryRunRegistrar", "construct_payload", "guard_qc", "make_registrar",
]
