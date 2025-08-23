# GFFx Changelog

Release v0.3.1:

### Changed
- To avoid short option conflicts in `extract`, renamed:
  - `-f` (feature id) -> `-e`
  - `-F` (feature file) -> `-E`

---

Release v0.3.0 is a **major release** of GFFx, introducing important changes to default behavior, performance improvements, and bug fixes.

  ## Breaking Changes
  - **Default extraction/intersection mode changed**:  
    Non Full-Model Mode is now the default. To preserve the previous behavior (returning the entire gene model), users must explicitly pass the `-F` / `--full-model` flag.
    Non Full-Model Mode and Full-Model Mode have **comparable runtime performance**.

  ## New Features
  - **Non Full-Model Mode for `extract` and `intersect`**  
    Added support for a new mode where these commands return only the features that directly match the queries, instead of always returning the entire gene model.  

  ## Changes
  - **Default Behavior Update**  
    Non Full-Model Mode is now the default (see Breaking Changes).  
  - **Threading Default Update**  
    The default number of threads has been changed from 4 to 12 CPU cores to improve the performance under default settings.
  
  ## Improvements
  - Refactored code to increase cohesion and function reusability.  

  ## Bug Fixes
  - Fixed several minor issues to improve stability and correctness.

---

