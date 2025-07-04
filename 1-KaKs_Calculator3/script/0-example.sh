#!/bin/bash
# Template for running KnKs and KaKs Calculator

# todo 软件所在目录
BASE_DIR="/mnt/e/Scientifc_software/KaKs_Calculator3.0"

# todo 修改输入和输出文件路径
NONCODING_FILE="${BASE_DIR}/examples/noncoding.axt"
CODING_FILE="${BASE_DIR}/examples/coding.axt"
OUTPUT_DIR="${BASE_DIR}/output"
KNKS_OUTPUT="${OUTPUT_DIR}/test.axt.knks"
KAKS_OUTPUT="${OUTPUT_DIR}/test.axt.kaks"
# Make sure output directory exists
mkdir -p "${OUTPUT_DIR}"

# Run KnKs analysis (comparing noncoding and coding sequences)
"${BASE_DIR}/src/KnKs" \
    -i "${NONCODING_FILE}" \
    -j "${CODING_FILE}" \
    -o "${KNKS_OUTPUT}"

# Run KaKs analysis (for coding sequences)
"${BASE_DIR}/src/KaKs" \
    -i "${CODING_FILE}" \
    -o "${KAKS_OUTPUT}"

echo "Analysis complete. Results saved to:"
echo "  ${KNKS_OUTPUT}"
echo "  ${KAKS_OUTPUT}"