# Nature Code and Software Checklist Status

Status reviewed against `docs/nr-software-policy`.

## Completed in the repository

- Source code and manuscript-facing notebooks are included.
- A real demonstration dataset and cached descriptors are included.
- `demo.py` provides a deterministic, timed, end-to-end demo with expected outputs.
- `environment.yml` has been solved from scratch on Windows and the demo passes in that clean environment; the exact resolved versions and installation time are recorded in the README.
- The source code is distributed under the OSI-approved BSD 3-Clause License in the repository-level `LICENSE` file.
- The archived release has the DOI `10.5281/zenodo.21897320`; the README includes a Zenodo badge and the manuscript Code availability statement cites the DOI.
- The README documents system requirements, tested versions, hardware requirements, installation, install-time expectations, demo inputs and outputs, custom-data columns, and manuscript reproduction.
- `docs/ALGORITHM.md` documents the model, fixed seeds, acquisition process, hierarchical filtering, and pseudocode.
- The manuscript describes the transition from the 222,768 nominal combinations to the 205,632 experimentally feasible combinations and includes computational-method and Code availability text.
- The edited Word manuscript passes archive-integrity, Word-to-PDF export, and page-by-page visual inspection; the added Methods and Code availability sections render without overflow or broken pagination.

## Actions remaining before peer review or public release

- Have a colleague unfamiliar with the project install and run the demo in a fresh environment; record their OS, versions, installation time, runtime, and any corrections.
- Replace the provisional Code availability text with the final public repository URL before publication.
- Add the corresponding public GitHub version/tag when the repository is opened. A repository-level `zenodo.json` metadata file remains intentionally deferred.
- Confirm the exact versions of Gaussian, xTB/CREST, and Julia only if the optional descriptor or PySR workflows are claimed as independently reproducible.
- Linux and macOS testing is optional but would strengthen the system-requirements statement.
