# Liquid-biopsy-style partial observation maps for hybrid occult residual-disease modes

**Thesis #23.** Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

Given a hybrid (switching) occult residual-disease model, which liquid-biopsy-style sparse observation maps can distinguish mode occupancy versus remaining blind to switches — as an observation-operator question, not as a diagnostic device claim?

The plant is a scalar count with an exogenous mode: quiescence, angiogenic pause, a six-day proliferative dwell, and immune-held latency. Quiescence and angiogenic pause share the count and do not share a turnover. A constant-coefficient map, the cartoon in which shedding per cell ignores the mode, assigns those two modes one law. A two-channel map reads turnover and net growth. On a 28-day clock the six-day dwell contains no sample. On a 2-day clock, at a declared noise level, the same two-channel map recovers the switch times; a louder noise fills the series with false switches.

The mode names are those of a companion architectural thesis. That thesis is about the modes. This one is about partial observers of them. The series are synthetic. They are not a validated assay.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/observer_toy.py` | Seeded hybrid toy and observation maps (seed 20260921) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Latent path, maps, labels, confusion, noise sweep |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/observer_toy.py
python3 build_pdf.py
```

NumPy and Matplotlib are required for the toy. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`.

## Cite

Ogbonna KE. Liquid-biopsy-style partial observation maps for hybrid occult residual-disease modes [Internet]. Thesis #23 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-23-ctdna-partial-observer-occult-modes

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
