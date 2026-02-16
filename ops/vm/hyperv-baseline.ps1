param(
  [Parameter(Mandatory = $true)]
  [ValidateSet('preflight', 'bringup', 'snapshot', 'reset', 'status')]
  [string]$Action,

  [string]$VmName = 'chem-baseline',
  [string]$SwitchName = 'Default Switch',
  [string]$VhdPath = 'C:\\HyperV\\chem-baseline\\chem-baseline.vhdx',
  [string]$IsoPath = '',
  [string]$SnapshotName = 'baseline',
  [int]$MemoryGb = 8,
  [int]$CpuCount = 4,
  [int]$DiskGb = 120,
  [switch]$StartVm,
  [switch]$ReplaceSnapshot
)

$ErrorActionPreference = 'Stop'

function Assert-Cmdlet {
  param([string]$Name)

  if (-not (Get-Command $Name -ErrorAction SilentlyContinue)) {
    throw "Missing required Hyper-V cmdlet: $Name"
  }
}

function Invoke-Preflight {
  Assert-Cmdlet -Name 'Get-VM'
  Assert-Cmdlet -Name 'Get-VMSwitch'
  Assert-Cmdlet -Name 'Get-VHD'

  if (-not (Get-VMSwitch -Name $SwitchName -ErrorAction SilentlyContinue)) {
    throw "VMSwitch '$SwitchName' was not found. Create it in Hyper-V Manager or pass -SwitchName."
  }

  Write-Host "Preflight OK for VM '$VmName' with switch '$SwitchName'."
}

function Invoke-Bringup {
  Invoke-Preflight

  if ($MemoryGb -le 0 -or $CpuCount -le 0 -or $DiskGb -le 0) {
    throw 'MemoryGb, CpuCount, and DiskGb must be positive integers.'
  }

  $existingVm = Get-VM -Name $VmName -ErrorAction SilentlyContinue
  if ($existingVm) {
    Write-Host "VM '$VmName' already exists. Skipping create."
  }
  else {
    $vhdDir = Split-Path -Path $VhdPath -Parent
    if (-not (Test-Path -LiteralPath $vhdDir)) {
      New-Item -ItemType Directory -Path $vhdDir -Force | Out-Null
    }

    if (-not (Test-Path -LiteralPath $VhdPath)) {
      New-VHD -Path $VhdPath -Dynamic -SizeBytes ($DiskGb * 1GB) | Out-Null
      Write-Host "Created VHD: $VhdPath"
    }
    else {
      Write-Host "Using existing VHD: $VhdPath"
    }

    New-VM `
      -Name $VmName `
      -Generation 2 `
      -MemoryStartupBytes ($MemoryGb * 1GB) `
      -VHDPath $VhdPath `
      -SwitchName $SwitchName | Out-Null

    Set-VMProcessor -VMName $VmName -Count $CpuCount
    Set-VM -Name $VmName -AutomaticStartAction Nothing -AutomaticStopAction ShutDown -CheckpointType Standard | Out-Null

    Write-Host "Created VM '$VmName' (memory=${MemoryGb}GB, cpu=$CpuCount, disk=${DiskGb}GB)."
  }

  if ($IsoPath -ne '') {
    if (-not (Test-Path -LiteralPath $IsoPath)) {
      throw "ISO not found: $IsoPath"
    }

    $dvdDrive = Get-VMDvdDrive -VMName $VmName -ErrorAction SilentlyContinue
    if ($dvdDrive) {
      Set-VMDvdDrive -VMName $VmName -Path $IsoPath | Out-Null
    }
    else {
      Add-VMDvdDrive -VMName $VmName -Path $IsoPath | Out-Null
    }

    Write-Host "Attached ISO: $IsoPath"
  }

  if ($StartVm) {
    Start-VM -Name $VmName | Out-Null
    Write-Host "Started VM '$VmName'."
  }

  Invoke-Status
}

function Invoke-Snapshot {
  Invoke-Preflight

  $vm = Get-VM -Name $VmName -ErrorAction SilentlyContinue
  if (-not $vm) {
    throw "VM '$VmName' not found. Run bringup first."
  }

  $existing = Get-VMSnapshot -VMName $VmName -Name $SnapshotName -ErrorAction SilentlyContinue
  if ($existing) {
    if ($ReplaceSnapshot) {
      Remove-VMSnapshot -VMName $VmName -Name $SnapshotName -Confirm:$false
      Write-Host "Removed existing snapshot '$SnapshotName'."
    }
    else {
      Write-Host "Snapshot '$SnapshotName' already exists. Use -ReplaceSnapshot to refresh it."
      return
    }
  }

  Checkpoint-VM -Name $VmName -SnapshotName $SnapshotName | Out-Null
  Write-Host "Created snapshot '$SnapshotName' for '$VmName'."
}

function Invoke-Reset {
  Invoke-Preflight

  $vm = Get-VM -Name $VmName -ErrorAction SilentlyContinue
  if (-not $vm) {
    throw "VM '$VmName' not found."
  }

  $snapshot = Get-VMSnapshot -VMName $VmName -Name $SnapshotName -ErrorAction SilentlyContinue
  if (-not $snapshot) {
    throw "Snapshot '$SnapshotName' not found for VM '$VmName'."
  }

  if ($vm.State -ne 'Off') {
    Stop-VM -Name $VmName -TurnOff -Force
  }

  Restore-VMSnapshot -VMName $VmName -Name $SnapshotName -Confirm:$false
  Write-Host "Restored '$VmName' to snapshot '$SnapshotName'."

  if ($StartVm) {
    Start-VM -Name $VmName | Out-Null
    Write-Host "Started VM '$VmName'."
  }

  Invoke-Status
}

function Invoke-Status {
  $vm = Get-VM -Name $VmName -ErrorAction SilentlyContinue
  if (-not $vm) {
    Write-Host "VM '$VmName' not found."
    return
  }

  $snapshots = Get-VMSnapshot -VMName $VmName -ErrorAction SilentlyContinue | Select-Object -ExpandProperty Name
  if (-not $snapshots) {
    $snapshotDisplay = '(none)'
  }
  else {
    $snapshotDisplay = ($snapshots -join ', ')
  }

  Write-Host "VM: $($vm.Name)"
  Write-Host "State: $($vm.State)"
  Write-Host "CPU: $($vm.ProcessorCount)"
  Write-Host "MemoryStartupBytes: $($vm.MemoryStartup)"
  Write-Host "Snapshots: $snapshotDisplay"
}

switch ($Action) {
  'preflight' { Invoke-Preflight }
  'bringup' { Invoke-Bringup }
  'snapshot' { Invoke-Snapshot }
  'reset' { Invoke-Reset }
  'status' { Invoke-Status }
  default { throw "Unsupported action: $Action" }
}
