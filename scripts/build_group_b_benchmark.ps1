# Build Group B benchmark VCF from source testing samples (no R required).
$ErrorActionPreference = "Stop"
$testig = "e:\ACGM\testig\testig"
$outDir = "e:\ACGM\testig\acmgamp_benchmark"
$sampleCol = "GROUP_B_BENCHMARK"

$variants = @(
  @{ file="260324100041.haplotypecaller_VEP.ann.split.vcf"; pattern="55631301" },
  @{ file="260324100041.haplotypecaller_VEP.ann.split.vcf"; pattern="154030906" },
  @{ file="260324100043.clean.vcf"; pattern="62813296" },
  @{ file="260324100039.haplotypecaller_VEP.ann.split.vcf"; pattern="126544289" },
  @{ file="260324100039.haplotypecaller_VEP.ann.split.vcf"; pattern="51687181" },
  @{ file="260324100039.haplotypecaller_VEP.ann.split.vcf"; pattern="132892556" },
  @{ file="260324100041.haplotypecaller_VEP.ann.split.vcf"; pattern="100990864" },
  @{ file="260324100041.haplotypecaller_VEP.ann.split.vcf"; pattern="108292760" },
  @{ file="260324100041.haplotypecaller_VEP.ann.split.vcf"; pattern="23539813" },
  @{ file="260324100041.haplotypecaller_VEP.ann.split.vcf"; pattern="226881976" },
  @{ file="260324100042.haplotypecaller_VEP.ann.split.vcf"; pattern="193614765" },
  @{ file="260324100041.haplotypecaller_VEP.ann.split.vcf"; pattern="55224052" },
  @{ file="260324100039.haplotypecaller_VEP.ann.split.vcf"; pattern="62694680" },
  @{ file="260324100039.haplotypecaller_VEP.ann.split.vcf"; pattern="73223597" },
  @{ file="260324100041.haplotypecaller_VEP.ann.split.vcf"; pattern="110155271" },
  @{ file="260324100041.haplotypecaller_VEP.ann.split.vcf"; pattern="66508138" },
  @{ file="260324100041.haplotypecaller_VEP.ann.split.vcf"; pattern="42914099" },
  @{ file="260324100041.haplotypecaller_VEP.ann.split.vcf"; pattern="55631298" },
  @{ file="260324100041.haplotypecaller_VEP.ann.split.vcf"; pattern="32156037" },
  @{ file="260324100041.haplotypecaller_VEP.ann.split.vcf"; pattern="50623471" }
)

function Get-VcfLine($path, $pattern) {
  $pos = $pattern.Trim("`t")
  $match = Select-String -Path $path -Pattern "`t$([regex]::Escape($pos))`t" | Select-Object -First 1
  if (-not $match) { throw "No match for $pos in $path" }
  return $match.Line
}

function Normalize-SampleColumn($line) {
  $parts = $line -split "`t"
  if ($parts.Count -lt 10) { throw "Invalid VCF line: $($line.Substring(0, [Math]::Min(80, $line.Length)))" }
  $parts[8] = "GT:AD:DP:GQ:PL"
  $parts[9] = "0/1:10,2:12:54:54,0,403"
  if ($parts.Count -gt 10) { $parts = $parts[0..9] }
  return ($parts -join "`t")
}

$headerVcf = Join-Path $testig "260324100041.haplotypecaller_VEP.ann.split.vcf"
$csq = Select-String -Path $headerVcf -Pattern "^##INFO=<ID=CSQ" | Select-Object -First 1
if (-not $csq) { throw "CSQ header not found" }

$header = @(
  "##fileformat=VCFv4.2",
  "##FILTER=<ID=PASS,Description=`"All filters passed`">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=`"Genotype`">",
  "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=`"Allelic depths`">",
  "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=`"Read depth`">",
  "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=`"Genotype quality`">",
  "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=`"Phred-scaled genotype likelihoods`">",
  $csq.Line,
  "#CHROM`tPOS`tID`tREF`tALT`tQUAL`tFILTER`tINFO`tFORMAT`t$sampleCol"
)

$dataLines = foreach ($v in $variants) {
  $path = Join-Path $testig $v.file
  Normalize-SampleColumn (Get-VcfLine $path $v.pattern)
}

$vcfOut = Join-Path $outDir "acmgamp_group_b_benchmark.vcf"
$header + $dataLines | Set-Content -Path $vcfOut -Encoding UTF8
Write-Host "Wrote VCF: $vcfOut ($($dataLines.Count) variants)"
