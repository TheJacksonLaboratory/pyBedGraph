#!/bin/bash

pip3 install pyBigWig

echo "Running simple tests ..."
python3 test.py

if [ $? -ne 0 ]; then
  exit 1
fi

if [ -f "test_files/ENCFF376VCU.bigWig" ] && [ -f "test_files/ENCFF376VCU.bedGraph" ]; then
  echo "Running extensive tests ..."
  python3 extensive_test.py
  if [ $? -ne 0 ]; then
    exit 1
  fi

fi
