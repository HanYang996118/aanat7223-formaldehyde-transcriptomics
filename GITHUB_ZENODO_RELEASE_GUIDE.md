# GitHub and Zenodo release checklist

## GitHub

1. Open https://github.com/new and create a repository named aanat7223-formaldehyde-transcriptomics.
2. Choose Public only after checking the file inventory and institutional release scope.
3. Do not pre-create a README, license, or .gitignore; these are already included in the archive.
4. Extract the ZIP and upload the 47 archive files from its root to the GitHub repository root. The largest source files are small enough for GitHub's browser uploader; GitHub documents a 25 MiB per-file browser limit and a 100-file batch limit.
5. Confirm that LICENSE, CITATION.cff, CODE_CATALOG.md, PATCH_NOTES.md, SHA256SUMS.txt, scripts/, metadata/, outputs/, supplementary/, environment/, and .github/workflows/ appear in the repository.
6. Review the README and license notice, then commit the files.

Official instructions:
- Create a repository: https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository
- Upload files: https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository

## Zenodo

1. Sign in at https://zenodo.org/ using the GitHub account that owns the repository.
2. Connect the GitHub account and authorize Zenodo.
3. In the Zenodo GitHub settings, click Sync now, find the new repository, and switch it On.
4. Only after the repository is enabled, create a GitHub Release from the reviewed version tag (suggested initial tag: v0.1.0).
5. Zenodo should ingest the new release and mint a version-specific DOI. Verify the Zenodo record before using the DOI.
6. Add the verified version DOI and concept DOI to the repository README and manuscript Data Availability statement.

Official instructions:
- Enable repository: https://help.zenodo.org/docs/github/enable-repository/
- GitHub DOI/citation workflow: https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content

Do not insert a placeholder DOI into the manuscript. A DOI exists only after Zenodo has completed the archive record.
