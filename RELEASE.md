This is a lightweight guide for testing and publishing a new release.

## Test release on TestPyPI with GitHub Actions

Assuming test branch is up to date and synced with remote, tag current commit with a new, unused test version name

```
$ git tag -a v0.1.1b2 -m ""
```

Push tags to remote to initialize TestPyPI action

```
$ git push --tags
```

Visit https://github.com/cchdo/ctdcal/actions to ensure test release was successful.

## Full release on PyPI with GitHub Actions

Checkout and tag commit to be released with appropriate version number (vX.Y.Z)

```
$ git checkout _hash_
$ git tag -a v0.2.0 -m ""
```

Assuming commit is already pushed, now push tags to remote (this will trigger the TestPyPI action above)

```
$ git push --tags
```

[Create a new release](https://github.com/cchdo/ctdcal/releases/new) on GitHub, which should automatically trigger an Action to publish to PyPI (see [.github/workflows/](https://github.com/cchdo/ctdcal/tree/master/.github/workflows)).