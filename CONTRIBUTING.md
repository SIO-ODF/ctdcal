# Contributing to ctdcal

## Installing the package in developer mode
Clone the package using git:

```
$ git clone https://github.com/cchdo/ctdcal.git
```

(Optional) Create a new virtual environment using the method of your choosing. Using conda:

```
$ conda create -n ctdcal_devel python=3.9
```

Install ctdcal in editable/developer mode (-e) with the [complete] tag to install development, documentation, and testing required packages:

```
$ pip install -e .[complete]
```

Note, for zshell, quotes are required:

```
$ pip install -e ".[complete]"
```

## Developing on new branch
Create feature/bug fix/etc. branch:

```
$ git checkout -b branch_name
```

Make changes and commit with a useful message:

```
$ git commit -m "useful commit message"
```

Push changes to GitHub:

```
$ git push origin branch_name
```

### Rules for pull request approval
* Test coverage must stay the same or increase. If you write new lines of code, you must write tests for it.
* Code must pass tests on entire testing matrix (python 3.8/3.9 and linux/mac/windows)
* Code must be formatted properly (should be handled by `black`/`flake8`/`isort` during pre-commit checks)
* Code must be documented following [numpydoc style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html) to be compatible with Sphinx auto-documentation.
