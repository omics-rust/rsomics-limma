#!/bin/bash
set -euo pipefail
export CARGO_HOME=/Volumes/KIOXIA/Developments/cargo-home
export CARGO_TARGET_DIR=/Volumes/KIOXIA/Developments/cargo-target
export TMPDIR=/Volumes/KIOXIA/Developments/tmp
toolchain=/Volumes/KIOXIA/Developments/rustup-home/toolchains/1.91.0-aarch64-apple-darwin/bin
export CARGO="$toolchain/cargo"
export RUSTC="$toolchain/rustc"
export RUSTDOC="$toolchain/rustdoc"
export RUSTFMT="$toolchain/rustfmt"
export CARGO_BUILD_JOBS=2
export PATH="$toolchain:$PATH"
for path in "$CARGO_HOME" "$CARGO_TARGET_DIR" "$TMPDIR"; do
  resolved=$(cd "$path" && pwd -P)
  case "$resolved" in /Volumes/KIOXIA/Developments/*) ;; *) echo "Forbidden build path: $resolved" >&2; exit 1;; esac
  echo "verified build path: $resolved"
done
diskutil info -plist / | /usr/bin/python3 -c 'import plistlib,sys; p=plistlib.loads(sys.stdin.buffer.read()); size=p["APFSContainerSize"]; free=p["APFSContainerFree"]; used=1-free/size; print("Boot APFS occupancy:",used,"free bytes:",free); sys.exit(1 if used>=.8 else 0)'
df -h / /Volumes/KIOXIA
target=$("$CARGO" metadata --no-deps --format-version 1 | /usr/bin/python3 -c 'import json,sys; print(json.load(sys.stdin)["target_directory"])')
test "$target" = "$CARGO_TARGET_DIR"
echo "Cargo resolved target: $target"
if [[ "${1:-}" == "--exec" ]]; then
  shift
  exec "$@"
fi
if [[ "${1:-}" == "clippy" ]]; then
  shift
  exec "$toolchain/cargo-clippy" clippy "$@"
fi
if [[ "${1:-}" == "fmt" ]]; then
  shift
  exec "$toolchain/cargo-fmt" fmt "$@"
fi
exec "$CARGO" "$@"

