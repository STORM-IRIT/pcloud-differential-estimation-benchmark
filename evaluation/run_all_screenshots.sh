#!/bin/bash

CONFIG_DIR="./configs_screenshots"

PYTHON_SCRIPT="run_screenshots.py"

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

if [ ! -d "$CONFIG_DIR" ]; then
    echo -e "${RED}Error: The directory $CONFIG_DIR doesn't exist.${NC}"
    exit 1
fi

if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo -e "${RED}Error: python script $PYTHON_SCRIPT doesn't exist.${NC}"
    exit 1
fi

shopt -s nullglob
files=("$CONFIG_DIR"/*.yaml)

if [ ${#files[@]} -eq 0 ]; then
    echo -e "${YELLOW}No .yaml file founded in $CONFIG_DIR${NC}"
    exit 0
fi

echo -e "${GREEN}=== Begin: ${#files[@]} configuration file(s) founded ===${NC}\n"


count=0
errors=0

for config_file in "${files[@]}"; do
    ((count++))
    filename=$(basename "$config_file")
    
    echo -e "${YELLOW}------------------------------------------------------------${NC}"
    echo -e "${YELLOW}[$count/${#files[@]}] Running: $filename${NC}"
    echo -e "${YELLOW}------------------------------------------------------------${NC}"

    python3 "$PYTHON_SCRIPT" --config "$config_file"

    if [ $? -eq 0 ]; then
        echo -e "${GREEN} Success for $filename${NC}\n"
    else
        echo -e "${RED} Error during $filename${NC}\n"
        ((errors++))
    fi
done

echo -e "${YELLOW}------------------------------------------------------------${NC}"
echo -e "End."
if [ $errors -eq 0 ]; then
    echo -e "${GREEN}Success.${NC}"
else
    echo -e "${RED}Warning : $errors script(s) have encounter errors.${NC}"
fi
