# Contributing

All contributions are welcome. 🔑

Writing documentation and tests is especially important. Below you can find a
brief overview of how to submit your patches to the project. If, instead, you
would like to simply report an issue, feel free to open one
[here](https://github.com/Quantum-for-Life/phase2/issues).

You can start off with forking this repository on GitHub and getting a local
clone:

```bash
git clone https://github.com:[username]/phase2.git
```

where `[username]` is your GitHub account holding the fork.🍴

Check out the main branch and create a topic branch for your work:

```bash
git checkout -b my-topic main
```

It is important that you start `my-topic` from the tip of the current `main`
branch and not from, e.g., another topic branch or `next` (see below).

When your code is ready for review, open a pull request (PR) to merge your
branch into this repository's `main`. A PR on GitHub is the place where your
changes can be discussed in detail, and your work improved and expanded.

In the meantime, you branch will be merged into `next` to see if it plays nicely
with other changes currently to be implemented. You will be able to see if your
code passes our
[test suite](https://github.com/Quantum-for-Life/phase2/actions)
running continuously on GitHub Actions.

Once the changes are accepted, your branch will be merged into `main`. Here's a
few tips how to make this process as smooth as possible:

1. Do not merge `next` into your topic branch. The `next` staging tree may have
   got all the newest hot-rod stuff, but utimately it's a throw-away: it can be
   reset anytime if it gets too messy. (For more info, see e.g.
   [git documentation](https://git-scm.com/docs/gitworkflows); we don't need
   `seen` or `pu` trees here, as our project isn't big enough).
2. If you need other topic branches for you work, merge them in with `--no-ff`
   git option.
3. Do not rebase branches that others have access to. Except for your topic
   branch: you can squash, rebase, twist and contort all you want, even during
   the PR.
4. We will ask you to write unit tests for all the changes you'd like to
   introduce.
5. Trying to keep git history linear and tidy à la
   [Gitflow](https://nvie.com/posts/a-successful-git-branching-model/) isn't
   necessary, though: complex nonlinear commit logs can be easily parsed with
   `git log --first-parent`.
