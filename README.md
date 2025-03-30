
# Gene Body Prediction using Hidden Markov Models (HMM)

This repository contains a full Python implementation of a 2-state Hidden Markov Model (HMM) to identify gene body intervals from unlabeled ChIP-seq histone modification data. The project uses the Baum-Welch algorithm (a form of Expectation-Maximization) to learn model parameters in an unsupervised fashion and outputs the 50,000 most probable intervals that overlap protein-coding gene bodies.

---

## Project Background

This project is based on **processed ChIP-seq data** from two histone modifications. The data is anonymized and provided as a sequence of symbolic characters for **200bp genomic intervals**:

| Symbol | Description |
|--------|-------------|
| `x`    | Only histone mark X was above threshold |
| `y`    | Only histone mark Y was above threshold |
| `z`    | Both marks X and Y were above threshold |
| `n`    | Neither mark X nor Y was above threshold |

The goal is to use these inputs to predict intervals that are most likely part of a **protein-coding gene body**, without using any labeled training data.

---

## Input

You will be provided with a file named `input.fasta`, formatted as:

