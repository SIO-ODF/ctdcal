This is a lightweight guide for testing and publishing a new release.

## Test release on TestPyPI

Tag current commit with a new, unused test version name

```
$ git tag -a v0.1.0b8 -m ""
```

Build distribution

```
$ python -m build
```

Upload using twine

```
$ twine upload --repository testpypi dist/ctdcal-0.1.0b8
```

When prompted:
```
username = __token__
password = TestPyPI API key
```

Visit https://test.pypi.org/project/ctdcal/0.1.0b8/ to ensure test release was successful.

## Full release on PyPI with GitHub Actions

Checkout and tag commit to be released with appropriate version number (vX.Y.Z)

```
$ git checkout _hash_
$ git tag -a v0.2.0 -m ""
```

Assuming commit is already pushed, now push tags to remote

```
$ git push --tags
```

[Create a new release](https://github.com/cchdo/ctdcal/releases/new) on GitHub, which should automatically trigger an Action to publish to PyPI (see [.github/workflows/](https://github.com/cchdo/ctdcal/tree/master/.github/workflows)).