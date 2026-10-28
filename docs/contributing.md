# Contributing to GORILLA_APPLETS

Report bugs and propose changes through the
[issue tracker](https://github.com/itpplasma/GORILLA_APPLETS/issues) and pull
requests against `main`. Run `make test` before opening a pull request.

## Continuous integration and merging

- Every pull request, including each layer of a stacked pull request, runs the
  fast CI tier. Add the label `full-ci` to run the full tier before queueing.
- Review a stack layer by layer. When it is ready, retarget the top layer to
  `main` so that the stack is merged as one pull request; the lower pull
  requests then show as merged.
- Merge with "Merge when ready" (or `gh pr merge --merge --auto`). The pull
  request enters the merge queue, the full CI tier runs on the exact merge
  result, and on success it lands on `main` as a merge commit. `main` is not
  re-tested on push.
- A branch does not need to be up to date with `main`. Merge `main` into a
  pull request only when it conflicts.
