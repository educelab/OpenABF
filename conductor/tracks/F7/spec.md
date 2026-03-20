# F7 — Implement Policy-Based Math Backends

## GitHub Issue
https://github.com/educelab/OpenABF/issues/55

## Summary
Implement a policy-based design for math backends, allowing users to swap out the linear algebra backend (currently Eigen) for alternatives via template policy parameters.

## Acceptance Criteria
- [ ] Policy interface defined for math backend operations
- [ ] Eigen backend implemented as reference/default policy
- [ ] Existing functionality unchanged when using default policy
- [ ] Documentation on implementing custom backends
