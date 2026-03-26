#!/bin/bash

find ./screenshots -type f -name "*.png" | while read -r file; 
    filename="${file%.png}"
    echo "Convert: $file"
    convert "$file" -background white -flatten -quality 50 "${filename}.jpg"
done

echo "End!"
