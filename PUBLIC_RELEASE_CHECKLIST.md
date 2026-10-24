# Public repository release checklist

The code and derived-data package is being published in the public GitHub repository `Dr3922/huangsw`. This checklist records the remaining versioned-release tasks.

1. Review all files under `data/derived/` for data-sharing and redistribution permissions.
2. Verify that the extracted files are present in the repository (do not upload the ZIP alone as the only repository file).
3. Select an appropriate software license after author approval. If a license is added, update README.md and CITATION.cff accordingly.
4. Run `python -m pip install -r requirements.txt` and `python scripts/smoke_reproduce.py --handoff data/derived`.
5. Create a fixed GitHub Release, e.g. `v1.0.0`, and verify its URL while signed out.
6. Update `CITATION.cff` with the **verified** `repository-code:` URL and `version:` tag. Do not invent a link or a DOI.
7. Insert the fixed release URL in the manuscript, point-by-point response, and editorial system before submission.

The full raw GEO/Visium objects are not packaged. The smoke test uses the supplied derived results and does not constitute a full independent raw-data regeneration.
