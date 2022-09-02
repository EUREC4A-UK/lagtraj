# Developing lagtraj

To work on the `lagtraj` codebase yourself you will want to know 1) how to get
and install a local copy of `lagtraj`, 2) how to run the lagtraj tests locally
and 3) how to contribute your changes back to the community. We've covered all
these steps below - welcome in the lagtraj community! :)

Here's a quick TLDR version if you're already familiar with github:

```bash
# clone locally your fork of lagtraj on github
git clone https://github.com/<your-github-username>/lagtraj/

# install this local copy with pip with the necessary development tools. Set up
# pre-commit to automatically run all linting checks on every commit
pip install -e ".[dev]"
pre-commit install

# create a branch for your fix/feature
git checkout -b my-new-feature

# make any modifications you want, make sure to add a test if you're fixing a
# bug and the run the tests
python -m pytest .

# commit these changes locally with git (pre-commit will automatically lint
# your code)
git add .
git commit

# push these changes to your fork on github (providing the name of the branch
# you want to push)
git push origin my-new-feature

# create a pull-request on github to pull in your branch by clicking the link
# the output from the last command, or go to
# https://github.com/EUREC4A-UK/lagtraj/compare
```

## Installing from a local copy

When you're making modifications to the lagtraj codebase you will want to
install this local copy with `pip` instead of installing from github or pypi.
First you will want to [create your own fork]() of lagtraj on github.

```bash
git clone https://github.com/<your-github-username>/lagtraj/
```

You can then install with `pip` using the following command (the `-e`-flag
ensures that pip will pick up any changes you make to the source directory).
Everything needed to get a working developing environment set up is stored
inside `setup.cfg` can be installed with pip by calling:

```bash
python -m pip install -e ".[dev]"
```

Linting is done with [pre-commit](https://pre-commit.com/), run the following
command to have linting run automatically for each git commit:

```bash
pre-commit install
```

## Adding features and fixing bugs

Modifying the lagtraj codebase to add features and fix bugs follows the
common process of branch, fix, test and pull-request. Once you've cloned your
fork of lagtraj locally and installed it with pip (see above) you can start
making your changes!

Start by creating a branch for you feature/fix:

```bash
git checkout -b my-new-feature
```

Then make any changes you'd like and commit them locally:

```bash
git add .
git commit
```

To test that your new feature doesn't break anything you can run the tests
locally (these will also be run automatically on github, see more details on
testing below)

```bash
python -m pytest .
```

Finally, when you're happy with your changes push your branch to github and
make a pull-request. The pull-request will tell the rest of the lagtraj
community that you have a change you'd like to get added to the public version
of lagtraj

```bash
git push origin my-new-feature
```

The last command will give you a link to where you can create a pull-request on
github, otherwise you can also go to
https://github.com/EUREC4A-UK/lagtraj/compare


## Testing lagtraj

`lagtraj` is automatically checked on github with tests that reside in `tests/`.
These are run automatically on all pull-requests against the git
repository at https://github.com/EUREC4A-UK/lagtraj and can be run locally
with `pytest` from the root of the repository:

```bash
pip install pytest
python -m pytest
```

It is useful to have a `ipdb`-debugger open up inline on failing
tests. This can be achieved by first installing `ipdb` and setting the
`PYTEST_ADDOPPS` environment variable:

```bash
export PYTEST_ADDOPTS='--pdb --pdbcls=IPython.terminal.debugger:Pdb'
```

### Speeding up testing

If you are running the tests repeatedly it's a good idea to download the
testdata locally and point lagtraj to where it resides (otherwise lagtraj will
attempt to download the testdata every time the tests are run):

```bash
wget http://gws-access.ceda.ac.uk/public/eurec4auk/testdata/lagtraj.testdata.tar.gz
mkdir /tmp/lagtraj
tar zxvf lagtraj.testdata.tar.gz -C /tmp/lagtraj
export LAGTRAJ_TESTDATA_DIR=/tmp/lagtraj
```

If the computer you are running on has multiple CPUs it can be advantageous to
run the tests in parallel to speed up the testing process. To do this you will
first need to install `pytest-xdist` and run pytest with `-n` to indicate the
number of parallel workers:

```bash
pip install pytest-xdist
pytest -n <n_cpus>
```

You can also speed up testing by reducing the number of test being run (if
you're for example working on fixing just a single breaking test) by using the
`-k` flag which where you provide a regex pattern for the name of tests you
want to run, e.g.

```bash
pytest -k lagrangian
```

### Updating the test data

The testdata used in lagtraj is stored in a .tar.gz-file (currently on JASMIN)
which is generated with the script `tests/make_test_data.py`, which can be run
as normal with python:

```bash
python tests/make_test_data.py
```
