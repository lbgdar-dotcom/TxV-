# Reverse-orientation clones

Blunt cloning has no orientation control, so roughly half your colonies will
carry the insert the other way round. These are the maps of those clones.
**They are not the design** — the maps one level up are.

They are kept because a flipped clone is a real product you may pick, not a
mistake, and because it is useful to be able to check one against a map
rather than guess.

## Why these read the same way as the designed maps

A plasmid is double-stranded, so which strand a map draws is a choice. Shown
on the strand these records happened to be built on, the whole cassette ran
right to left — which is what "the map is upside down" looks like, and it
makes the file useless for checking a design.

Each record here is therefore shown **from its other strand** and rotated so
the insert begins at position 1. The molecule is unchanged; only the strand
drawn and the base it is counted from differ.

The consequence is that the cassette reads 5'→3' from the top in both sets,
and what actually differs between a designed and a flipped clone becomes the
visible difference instead of an apparent mirror image:

| | designed | flipped |
|---|---|---|
| cassette | 1–871, `+` | 1–871, `+` |
| `AmpR` | `-` strand | **`+` strand** |
| `ori` | `-` strand | **`+` strand** |

## Does it matter which one you pick?

For this workflow, no. The BglII digest releases the insert either way, and the
released fragment carries the construct's own T7 promoter, so run-off
transcription gives the same mRNA from either clone.

If you want to know which you have before sequencing: `IVT_F` with the pJET1.2
**reverse** sequencing primer gives a product from a designed clone and
nothing from a flipped one. Expected sizes are in `CLONING_PROTOCOL.md`.
