# Release review

Reviewed on 2026-10-16 before the initial Git commit.

- Scope: source code, documentation, configuration examples, filenames and hidden files in this release folder.
- Removed the internal appendix containing server paths, account names, institution labels and original-to-release filename mappings. The private backup is outside this repository.
- Checked for embedded patient/sample identifiers, clinical records, email addresses, phone/national-ID patterns, private network addresses, access tokens, passwords and private keys. No such values were identified in the files selected for this release. Synthetic examples, data schemas and aggregate study results remain.
- Configuration examples contain placeholders. Real inputs, clinical tables, sequence files, sample lists, analysis outputs, local configuration, credentials and caches are excluded from the initial commit; `.gitignore` also protects common future additions.
- Added the MIT license and matching package metadata. External software and datasets retain their own licences and access conditions.
- Validation: all 116 Python source files passed AST syntax parsing in Python 3.8.5; all 43 existing synthetic-input tests passed with pytest 7.2.1. No analysis algorithms were changed for this release.
- Full pipeline reproduction and installation of fresh conda environments were not repeated for this upload. The repository intentionally contains no study data.

This review combines automated pattern checks with manual inspection. It is not a formal guarantee that every possible sensitive value can be detected. Review newly added configuration and data before future commits.
