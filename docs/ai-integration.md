# AI integration

This page describes ways to use SemiBin together with AI assistants: an
interactive chat bot that can answer questions about SemiBin, and an agent
*skill* that teaches coding agents how to drive `SemiBin2` correctly.

## Ask the SemiBin chat bot

If you have a question about SemiBin — how to run it, what a particular option
does, how to interpret the output, or how to troubleshoot an error — you can ask
our NotebookLM-powered chat bot, which has been given the SemiBin documentation
as its source material:

[**➜ Chat with the SemiBin assistant**](https://notebooklm.google.com/notebook/e7c2bfae-9756-41d1-99a0-77633c626c58)

The bot answers from the SemiBin documentation and cites the sections it draws
from, so you can follow up in the docs themselves. It is a good first stop for
"how do I…?" questions.

!!! note
    The chat bot is a convenience for exploring the documentation. For bug
    reports, feature requests, or questions that need a human, please use the
    [GitHub issue tracker](https://github.com/BigDataBiology/SemiBin/issues) or
    the [SemiBin users mailing list](https://groups.google.com/g/semibin-users).

## The SemiBin agent skill (for coding agents)

SemiBin ships with an [agent skill](https://docs.claude.com/en/docs/claude-code/skills)
that teaches coding agents (such as Claude Code) how to use `SemiBin2`
correctly: which subcommand to pick for each binning mode, how to validate
inputs, where the output bins are written, and how to resolve common errors.
When the skill is installed, an agent working in your project can pick the right
command and flags on the first try instead of guessing.

### Installing the skill

Use the `install-skills` subcommand. By default it installs into the current
project (`./.claude/skills`), so an agent working in that directory finds it
automatically:

```bash
SemiBin2 install-skills
```

To install it once for your user account (available in every project), install
into `$HOME/.claude/skills` instead:

```bash
SemiBin2 install-skills --user
```

You can also point it at an explicit directory with `--skills-dir`. See the
[`install-skills` subcommand reference](../subcommands/#install-skills) for details.

### What the skill covers

- Choosing between single-sample, co-assembly, and multi-sample binning.
- Building correct command lines (including long-read and pre-trained-model workflows).
- Validating inputs (sorted/indexed BAM files, uniquely named contigs, the multi-sample separator).
- Where the final bins and summary tables are written.
- Fixes for the most common error messages.

The skill is bundled inside the installed SemiBin package, so the copy you
install always matches the version of SemiBin you are running.

## We want your feedback

These AI integrations are new and we are **very interested in hearing how you
use them**. In particular, if the chat bot or an AI agent using the skill ever
gives you **bad, misleading, or outdated advice**, please let us know: that kind
of report is extremely valuable and helps us improve the documentation and the
skill.

Please share your experience (good or bad) through the
[GitHub issue tracker](https://github.com/BigDataBiology/SemiBin/issues) or the
[SemiBin users mailing list](https://groups.google.com/g/semibin-users).

!!! warning "Reporting bad advice"
    When reporting bad advice from the chat bot or an AI agent, it helps to
    include the exact question you asked and the answer you received, so we can
    reproduce and fix the underlying documentation or skill.
