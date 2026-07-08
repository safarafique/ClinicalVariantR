$p = Join-Path $PSScriptRoot "setup_intervar.sh"
$bytes = [System.IO.File]::ReadAllBytes($p)
$crBefore = ($bytes | Where-Object { $_ -eq 13 }).Count
Write-Output "CR bytes before: $crBefore"
$text = [System.Text.Encoding]::UTF8.GetString($bytes)
$text = $text -replace "`r`n", "`n" -replace "`r", "`n"
$utf8 = New-Object System.Text.UTF8Encoding $false
[System.IO.File]::WriteAllText($p, $text, $utf8)
$bytes2 = [System.IO.File]::ReadAllBytes($p)
$crAfter = ($bytes2 | Where-Object { $_ -eq 13 }).Count
Write-Output "CR bytes after: $crAfter"
