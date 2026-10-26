# Contributing
1. `make cpp && make test` must pass (C++ unit tests, Python unit tests, C++↔Python cross-check).
2. A new mathematical function needs: a docstring stating the theorem it implements, a proof or reference in `docs/theory.md`, and a test that checks it against an independent computation (e.g. dense linear algebra), not against itself.
3. Never commit real-data results without the pinned source SHA and the command that produced them.
