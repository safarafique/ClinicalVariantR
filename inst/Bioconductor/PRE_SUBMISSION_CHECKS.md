# ClinicalVariantR Bioconductor Pre-Submission Checks

Run these steps before opening or updating the Bioconductor submission issue.

## 1. Start From Package Root

Open Windows Command Prompt or PowerShell and move to the package folder:

```cmd
cd /d E:\ACGM\ClinicalVariantR
```

## 2. Clean Old Build Outputs

Remove old tarballs and check folders so results come from the current source:

```cmd
del ClinicalVariantR_0.99.3.tar.gz
rmdir /s /q ClinicalVariantR.Rcheck
rmdir /s /q ClinicalVariantR.BiocCheck
```

It is OK if Windows says a file or folder was not found.

## 3. Run Unit Tests

Use Rscript from the package root:

```cmd
"C:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" -e "devtools::test()"
```

Expected result:

```text
FAIL 0 | WARN 0
```

One skipped Shiny layout test is acceptable if it says the Shiny entry files are not in the current working/install layout.

## 4. Build Source Package

```cmd
"C:\Program Files\R\R-4.6.0\bin\x64\R.exe" CMD build .
```

Expected output:

```text
* building 'ClinicalVariantR_0.99.3.tar.gz'
```

## 5. Run R CMD Check

```cmd
"C:\Program Files\R\R-4.6.0\bin\x64\R.exe" CMD check ClinicalVariantR_0.99.3.tar.gz
```

Review the final summary. Submission target:

```text
Status: OK
```

If there are ERRORs or WARNINGs, fix them before submission. NOTES may be acceptable if they are justified.

## 6. Run BiocCheck

```cmd
"C:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" -e "BiocCheck::BiocCheck('ClinicalVariantR_0.99.3.tar.gz', newPackage=TRUE)"
```

Submission target:

```text
0 ERRORS
```

Warnings and notes should be reviewed. Current acceptable residual notes may include broad style/refactor suggestions, such as long functions, line length, or use of `<<-`.

## 7. Manual Bioconductor Account Checks

Before submitting, confirm:

- You are subscribed to the Bioc-devel mailing list:
  https://stat.ethz.ch/mailman/listinfo/bioc-devel
- You are registered on the Bioconductor Support Site.
- Your Support Site profile has watched tag:

```text
ClinicalVariantR
```

Edit profile:

```text
https://support.bioconductor.org/accounts/edit/profile
```

If you also maintain `GExPipe`, keep both tags in watched tags:

```text
GExPipe
ClinicalVariantR
```

## 8. Check Package Metadata

Confirm `DESCRIPTION` contains:

```text
Package: ClinicalVariantR
Version: 0.99.3
URL: https://github.com/safarafique/ClinicalVariantR
BugReports: https://github.com/safarafique/ClinicalVariantR/issues
```

The GitHub repository name and `Package:` field must match exactly:

```text
ClinicalVariantR
```

## 9. Submit To Bioconductor

Open a new issue:

```text
https://github.com/Bioconductor/BiocContributions/issues/new
```

Use issue title:

```text
ClinicalVariantR
```

Repository URL:

```text
https://github.com/safarafique/ClinicalVariantR
```

After the automated precheck passes, comment:

```text
/accept-policies
```

## 10. After Submission

When reviewers comment:

- Reply point-by-point.
- Fix requested changes in GitHub.
- Bump the package version for each new review build, for example:

```text
0.99.0 -> 0.99.1 -> 0.99.2 -> 0.99.3
```

- Rebuild and rerun `R CMD check` and `BiocCheck` after every change.
