# Releases FreeFEM

Currently, the [releases.yml](./workflows/release.yml) workflow creates and
uploads in the FreeFEM release the following assets:
- DEB packages
- EXE packages

## Triggers

Releases are created whenever a **v\*.\*** tag is added. 

---
**NOTE**

This pattern allows to test release candiates on any branch, and adding any suffix to the tag.
For instance v1.0-test triggers the release workflow in the github actions.
---

