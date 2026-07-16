# Draft: BiocContributions new package issue

**Where to submit:**  
https://github.com/Bioconductor/BiocContributions/issues/new  

**Issue title:** `ClinicalVariantR`  

Use the official **new submission template** without modifying its structure.
Paste your repository URL where indicated:

```text
https://github.com/safarafique/ClinicalVariantR
```

Confirm that you have read:

- Package guidelines: https://contributions.bioconductor.org/  
- Maintainer responsibilities  
- The review process  

After precheck passes, comment **exactly**:

```text
/accept-policies
```

## Suggested short description (if the template asks)

ClinicalVariantR is a Bioconductor / Shiny software package for ACMG/AMP germline variant
classification from VEP-, SnpEff-, or ANNOVAR-annotated VCF files, with streaming
support for large call sets and structured criterion-level evidence export.

## Checklist before clicking Submit

- [ ] Default branch is the package source (`DESCRIPTION` at root).  
- [ ] Version is `0.99.z`.  
- [ ] `R CMD check` and `BiocCheck` cleaned or notes justified in the issue.  
- [ ] Maintainer in `DESCRIPTION` matches the GitHub account that will accept policies.  
- [ ] SSH keys ready for later Bioconductor git access: https://bioconductor.org/developers.html  

## After submit (order of events)

1. Automated precheck  
2. `/accept-policies`  
3. Staging clone + R-universe build report  
4. Fix ERRORs/WARNINGs; bump `0.99.z`; push to **Bioconductor staging** as instructed  
5. Reviewer assigned → point-by-point replies  
6. Acceptance → canonical Bioconductor git + devel manifest  

Do **not** open the issue until `PACKAGE_CONVERSION.md` is complete and checks are clean.
