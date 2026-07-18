# Contributing to SemiBin

Thanks for your interest in improving SemiBin! Contributions of all kinds are
welcome: bug reports, documentation fixes, new features, and questions.

## Ways to contribute

- **Report a bug** or request a feature via
  [GitHub issues](https://github.com/BigDataBiology/SemiBin/issues). Reporting
  bugs is actively encouraged, and we will do our best to respond quickly.
- **Ask a question** or start an open-ended discussion on the
  [SemiBin users mailing-list](https://groups.google.com/g/semibin-users).
- **Submit a pull request** with a bugfix, new feature, or documentation
  improvement.

## Reporting bugs

Good bug reports make fixes much faster. Please include:

- The SemiBin version (`SemiBin2 --version`) and how you installed it
  (conda/bioconda, pixi, pip, source).
- Your operating system and Python version.
- The exact command you ran and the full output/error message.
- If possible, a minimal example that reproduces the problem.

## Development setup

SemiBin uses [pixi](https://pixi.sh/) to manage conda dependencies, including
the external command-line tools (samtools, bedtools, hmmer, mmseqs2, prodigal,
fraggenescan).

To run tests inside a managed environment:

```bash
pixi run -e test-py312 pytest
```

### Running the tests

```bash
python -m pytest                                 # all unit tests
python -m pytest test/test_bin.py                # a single test file
python -m pytest test/test_bin.py -k test_name   # a single test
```


Currently, we support Python 3.10 and above. You can run tests in a specific Python version using pixi:

```bash
pixi run -e test-py310 pytest
pixi run -e test-py311 pytest
pixi run -e test-py312 pytest
```

Integration tests require an installed package plus the external tools, but you
can also run them in a specific Python version using pixi:

```bash
pixi run -e test-py311 python integration-tests/generate_data_single_command.py
```

## Pull requests

1. Fork the repository and create a branch for your change.
2. Add or update tests covering your change where it makes sense.
3. Make sure the test suite passes locally.
4. Update the documentation if you change user-facing behaviour.
5. For new features and bugfixes, update **both**:
   - `ChangeLog` — add a line under `Unreleased` summarising the change
     (include the issue number if applicable).
   - `docs/whatsnew.md` — add a bullet under the current unreleased version
     section.
6. Open a pull request against the `main` branch with a clear description of
   what the change does and why.

### Commit message conventions

Prefix commit messages with one of the following tags (they can be combined,
e.g. `BUG+DOC`):

| Tag    | Meaning                          |
|--------|----------------------------------|
| `BUG`  | Bug fix                          |
| `ENH`  | Enhancement / new feature        |
| `MIN`  | Minor change                     |
| `RFCT` | Refactor                         |
| `TST`  | Tests                            |
| `DOC`  | Documentation                    |
| `RLS`  | Release                          |

## Code of conduct

We want SemiBin to be a welcoming and respectful community for everyone. When
participating in the project (through pull requests, github issues, or the
mailing-list)

- Be kind, patient, and considerate towards others.
- Assume good faith
- Respect differing viewpoints and experiences.
- Avoid harassment, discriminatory language, personal attacks, and any other
  disruptive or unwelcome behaviour.

If you experience or witness unacceptable behaviour, please contact the
maintainers at luis@luispedro.org. Reports will be handled confidentially, and
maintainers may take any action they deem appropriate, up to and including
removing a contributor from project spaces.

## License

By contributing to SemiBin, you agree that your contributions will be licensed
under the [MIT License](COPYING.MIT) that covers the project.

