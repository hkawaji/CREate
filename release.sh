#!/bin/bash
set -ue
# Usage: ./release.sh <version>   e.g. ./release.sh 0.2.2
# Bump VERSION in CREate.sh, commit that change, and tag v<version>.

ver="${1:-}"
if [ -z "$ver" ]; then
  echo "Usage: $0 <version>   (e.g. $0 0.2.2)" >&2
  exit 1
fi

cd "$(dirname "$0")"

# require a clean working tree so only the version bump is committed
if [ -n "$(git status --porcelain)" ]; then
  echo "Error: working tree not clean. Commit or stash changes first." >&2
  exit 1
fi

cur="$(grep -oE '^VERSION="[^"]*"' CREate.sh | sed -E 's/VERSION="([^"]*)"/\1/')"
if [ "$cur" = "$ver" ]; then
  echo "Error: VERSION is already ${ver}." >&2
  exit 1
fi

sed -i -E "s/^VERSION=\"[^\"]*\"/VERSION=\"${ver}\"/" CREate.sh
grep -q "^VERSION=\"${ver}\"" CREate.sh || { echo "Error: failed to set VERSION." >&2; exit 1; }

git add CREate.sh
git commit -m "bump version to ${ver}"
git tag "v${ver}"

echo "Committed and tagged v${ver}."
echo "Review, then push with:  git push && git push --tags"
