# GFFx v0.3.2 Release Notes

Release Note: Architectural Optimization

In this release, we have unified the parsing of various index files under a set of independent and consistently designed structures/modules, significantly improves maintainability, reliability, and usability of the codebase.
This architectural refactoring brings several key benefits:

Improved code reuse and cohesion: Common functionality is encapsulated, reducing duplication and making the system easier to maintain.

Greater testability and evolvability: A standardized interface across different index formats enables more robust testing and smoother future extensions.

Lower learning curve: By consolidating parsing logic into a clear and consistent abstraction, the framework becomes more approachable for researchers and developers who want to use or extend GFFx.

---

## Bug Fixes
- Fixed several issues that could cause the search functionality to malfunction in non-full-model mode.
- Removed modules, functions, and structures that were only useful during development and testing.
