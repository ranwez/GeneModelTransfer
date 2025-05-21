set -euo pipefail

SNAKEFILE_DIR=$1
CONFIG_FILE=$2
OUTPUT_FILE=$3
SENTINEL_FILE=$4

echo -e "\n\n\t\t-----------------------------------------------------------------------------" >>"$OUTPUT_FILE"
echo ">>> DATE AND TIME:" >>"$OUTPUT_FILE"
date >>"$OUTPUT_FILE"
echo "" >>"$OUTPUT_FILE"

cd "$SNAKEFILE_DIR"
if git rev-parse --git-dir >/dev/null 2>&1; then
  echo ">>> COMMIT ID:" >>"$OUTPUT_FILE"
  git rev-parse HEAD >>"$OUTPUT_FILE"
  echo "" >>"$OUTPUT_FILE"
fi
cd - >/dev/null

echo ">>> CONFIG FILE (${CONFIG_FILE}):" >>"$OUTPUT_FILE"
sed 's/#.*//' "$CONFIG_FILE" | grep -vP "^\s*$" >>"$OUTPUT_FILE"

touch "$SENTINEL_FILE"
