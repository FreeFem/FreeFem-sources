# FreeFEM Release Workflow

This document describes the automated release process for FreeFEM using GitHub
Actions.

## Overview

The [`releases.yml`](./workflows/releases.yml) workflow automatically builds and
publishes release assets when version tags are pushed to the repository.

## Release Assets

The workflow generates and uploads the following platform-specific packages:
- **Debian packages** (`.deb`) for configured Ubuntu versions
- **Windows installer** (`.exe`)

## Versioning and Triggers

The workflow triggers on any tag matching the pattern `v*.*` (e.g., `v4.16`,
`v4.16-rc1`, `v4.18-test`).

When triggered, the workflow creates a GitHub Release named after the tag, but
**the actual package versions are parsed from
[`configure.ac`](../configure.ac)**. This deliberate separation allows you to
create test or candidate releases (e.g., `v4.16-test`, `v4.16-rc2`) on any
branch, while keeping the correct semantic version from `configure.ac` in the
assets name.

## Creating a Release

To publish a new official release:

1. Update the version number in [`configure.ac`](../configure.ac) on the
   `develop` branch
2. Merge `develop` into `master`
3. Tag the merge commit with a version tag (see [Tag Patterns](#tag-patterns)
   below)
4. Push the tag to GitHub to trigger the release workflow

## Configuration

To modify which operating systems are targeted by the builds, edit the
`set_versions` job in [`releases.yml`](./workflows/releases.yml):

- **`ubuntu_versions`**: Specifies which Ubuntu versions receive `.deb` packages
- **`windows_versions`**: Specifies which Windows versions receive the `.exe`
  installer
