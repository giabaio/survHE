# survHE: CRAN release checklist

## Introduction

Prepare the current version from `main` on a temporary `cran` branch. r-universe builds `main`, so `main` stays untouched (and keeps the "development version" markers) until CRAN accepts. The code on `cran` is identical to `main`; only the release markers differ.

The submitted commit is kept by a tag (e.g. `v2.0.52`) and a GitHub release. Tags belong to the repository, not to a branch; thus, the tag and release stay visible from `main` (branch/tag menu, Releases page) after `cran` is deleted. The tag also keeps the commit itself in the repository.

- Check that `DESCRIPTION` says `Version: X.X.XX` (i.e. the current development version to move to CRAN) and that this is higher than the version on CRAN.
- If CRAN asks for code changes, make them on `main` first, then merge `main` into `cran`. That way the fix is never only on `cran`.

## 1. Prepare the release on a cran branch

```bash
git checkout main
git pull
git checkout -b cran
```

- [ ] `.github/workflows/check-standard.yaml`: add `cran` to the `push: branches` list, so that pushes to `cran` run the 5 checks. Change this on `cran` only, so nothing is pushed to `main`.

* [ ] `R/zzz.R`: remove the "development version" tag from the startup message
* [ ] `NEWS.md`: change the header `# 2.0.52 (development version)` to `# survHE 2.0.52`
* [ ] `README.Rmd`: update and re-knit `README.md` (both are in `.Rbuildignore`; thus, this does not affect CRAN)
* [ ] `cran-comments.md`: update the test environments and the expected NOTEs (see Step 2) -- this can be created with `usethis::use_cran_comments()` if it does not exist (and added to `.Rbuildignore`).

```bash
git add -A
git commit -m "Prepare for CRAN submission"
git push -u origin cran
```

- [ ] Wait for the 5 GitHub checks on `cran` to pass

## 2. Checks before submitting

The GitHub checks skip the PDF manual and the CRAN incoming checks; thus, run these locally on the release commit.

```r
devtools::check(remote=TRUE,manual=TRUE)  # as CRAN incoming, with the PDF manual
urlchecker::url_check()                   # every URL in the Rd files and DESCRIPTION
devtools::check_win_devel()               # results by email in about 30 minutes
revdepcheck::revdep_check(num_workers=4)  # only if survHE has reverse dependencies
```

- [ ] `check()`: 0 errors and 0 warnings. Read the "CRAN incoming feasibility" section of the output.
- [ ] `url_check()`: no broken links or redirects
- [ ] win-builder: same result as the local check
- [ ] revdeps: no new problems. The `make_data_multi_state()` fix (PR #72) changes results, so watch any reverse dependency that compares against stored numbers.
- [ ] Put the results in `cran-comments.md`

## 3. Submit to CRAN

```r
devtools::submit_cran()
```

This builds the package, asks some final questions and uploads it. It then writes `CRAN-SUBMISSION`, which records the version and the SHA of the commit submitted. That file is already in `.Rbuildignore`.

- [ ] Confirm the submission from the email CRAN sends to the maintainer address
- [ ] Commit and push `CRAN-SUBMISSION`

```bash
git add CRAN-SUBMISSION
git commit -m "Submitted 2.0.52 to CRAN"
git push
```

If CRAN asks for changes: make code fixes on `main` and merge `main` into `cran`; make release-only fixes (e.g. `cran-comments.md`) on `cran`. Then re-run Step 2 and submit again from `cran`. `CRAN-SUBMISSION` is updated with the new SHA each time.

## 4. After acceptance: tag and release

```r
usethis::use_github_release()
```

This reads the SHA from `CRAN-SUBMISSION` and tags exactly that commit as `v2.0.52`. It then creates a GitHub release with the 2.0.52 section of `NEWS.md` as its notes and deletes `CRAN-SUBMISSION`.

Run it from the `cran` branch, where `CRAN-SUBMISSION` lives. I'm not certain whether it complains when you are not on the default branch; if it does, use the manual alternative below.

- [ ] Check that the release and tag are on GitHub

Manual alternative, if you prefer not to use usethis:

```bash
git tag -a v2.0.52 <SHA> -m "survHE 2.0.52 (CRAN)"
git push origin v2.0.52
```

Here `<SHA>` is the commit recorded in `CRAN-SUBMISSION`. Then create the GitHub release from the tag (Releases → Draft a new release) and paste the 2.0.52 NEWS section.

When the tag and release exist, delete the branch. The tag keeps the commit.

```bash
git checkout main
git branch -D cran              # -D: cran is never merged into main
git push origin --delete cran
```

## 5. Back to development (2.0.53)

`main` never received the release edits, so `README.Rmd` and `R/zzz.R` need nothing. Only the version moves on. This push to `main` is also what updates r-universe to 2.0.53.

- [ ] `DESCRIPTION`: `Version: 2.0.53`
- [ ] `NEWS.md`: change `# 2.0.52 (development version)` to `# survHE 2.0.52`, as on `cran`, and add `# 2.0.53 (development version)` above it
- [ ] Copy any edits to `cran-comments.md` from `cran` to `main` if you want to keep them: `git checkout v2.0.52 -- cran-comments.md`

```bash
git commit -am "Start 2.0.53 development"
git push
```

`usethis::use_dev_version()` would do the first two steps, but it uses `2.0.52.9000`-style numbers.

## Expected NOTEs and pitfalls

| Message | Where | What to do |
| --- | --- | --- |
| Suggests not in mainstream repositories (`survHEinla`, `survHEhmc`, `INLA`) | CRAN incoming | Expected. Explain in `cran-comments.md` that they are optional and available from the repositories in `Additional_repositories`. |
| Skipping HTML validation: no command 'tidy'; package 'V8' unavailable | Local check only | Harmless. To remove it: `sudo apt install tidy` and `install.packages("V8")` |
| Code fails when a Suggests package is missing | CRAN machines | Code and tests that use `survHEinla`, `survHEhmc` or `INLA` must check `requireNamespace()` first. The GitHub checks install everything, so they cannot catch this. |
| URL redirects (301) | CRAN incoming | Replace the URL with the final address reported by `urlchecker::url_check()` |
| Snapshot tests | CRAN | `expect_snapshot()` is skipped on CRAN; nothing to do |
