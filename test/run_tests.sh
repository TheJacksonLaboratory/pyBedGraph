#!/bin/bash

echo "Running simple tests ..."
python3 test.py

if [ -f "test_files/ENCFF376VCU.bigWig" ] && [ -f "test_files/ENCFF376VCU.bedGraph" ]; then
  echo "Running extensive tests ..."
  python3 extensive_test.py
fi
