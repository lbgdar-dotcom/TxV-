# Order packages — which one to use

Three packages are kept in step with each other. **Only one is the live route.**

| package | vector | assembly | status |
|---|---|---|---|
| **`pjet12_ivt_units/`** | **pJET1.2/blunt** | **blunt, T4 ligase** | **LIVE — order from this one** |
| `pvax1_ag_panel/` | pVax1_AG | Gibson into the linearised backbone | superseded |
| `pvax1_ag_panel_goldengate/` | minimal IVT vector | BsaI Golden Gate | superseded |

The two superseded packages are regenerated from the same module encodings as
the live one, so they are *consistent*, not abandoned. They are kept because
they are the fallback if pJET1.2 blunt cloning gives trouble, and because the
pVax1_AG maps are the only ones that show the design in a vector already in use
in the lab. They are not what to order today.

## What "in step" means

Every package is built from `txv.pvax1_ag.canonical_encodings()`, so a module
has one DNA sequence across all three. Rebuild all of them after any design
change:

```
python scripts/build_pjet_package.py
python scripts/build_goldengate_package.py
python -m txv.cli order --out orders/pvax1_ag_panel
python -m txv.cli plasmids backbone/pvax1_ag_egfp.gb \
    --prefix pVax1_AG_ --out orders/pvax1_ag_panel/plasmids
```

Then re-run both auditors and the in-silico validation:

```
python audit/audit_pjet.py && python audit/audit_panel.py
python -m pytest -q && python -m pytest audit/ -q
python scripts/validate_in_silico.py
```
