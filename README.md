
# Gene Body Prediction using Hidden Markov Models (HMM)

This repository contains a full Python implementation of a 2-state Hidden Markov Model (HMM) to identify gene body intervals from unlabeled ChIP-seq histone modification data. The project uses the Baum-Welch algorithm (a form of Expectation-Maximization) to learn model parameters in an unsupervised fashion and outputs the 50,000 most probable intervals that overlap protein-coding gene bodies.

---

## Project Background

This project is based on **processed ChIP-seq data** from two histone modifications. The data is anonymized and provided as a sequence of symbolic characters for genomic intervals (150 - 300 bps):

| Symbol | Description |
|--------|-------------|
| `x`    | Only histone mark X was above threshold |
| `y`    | Only histone mark Y was above threshold |
| `z`    | Both marks X and Y were above threshold |
| `n`    | Neither mark X nor Y was above threshold |

The goal is to use these inputs to predict intervals that are most likely part of a **protein-coding gene body**, without using any labeled training data.

---

## How to Use

You will be provide a file named `input.fasta`, formatted as follows:

<pre>
>seq
n
x
y
z
...
</pre>

Each character corresponds to a genomic interval.

## Function Definitions

In this model, we define two hidden states:

- **State 0**: Background — intervals that are **not** part of a gene body
- **State 1**: Gene Body — intervals that are **likely within** a protein-coding gene

These states are **not observed directly**. Instead, we observe a sequence of symbols (`x`, `y`, `z`, `n`) representing histone modification activity in each 200bp genomic interval. The goal of the model is to infer the most probable hidden state for each interval using these observations.

### `forward_algo(seq, states, trans_prob, emiss_prob, label)`
Computes the **forward probabilities** for each state at every position in the sequence.

- **Parameters:**
  - `seq` *(str)* – The input sequence of symbols (`x`, `y`, `z`, `n`)
  - `states` *(list of int)* – List of hidden states (e.g., `[0, 1]`)
  - `trans_prob` *(2D numpy array)* – Transition probability matrix between hidden states
  - `emiss_prob` *(2D numpy array)* – Emission probability matrix (state × symbol)
  - `label` *(dict)* – Maps symbols (`x`, `y`, `z`, `n`) to column indices

- **Returns:**  
  A normalized matrix of shape `(num_states, sequence_length)` representing the probability of the observed sequence up to each position for each state.

---

### `backward_algo(seq, states, trans_prob, emiss_prob, label)`
Computes the **backward probabilities** for each state at every position in the sequence.

- **Parameters:**
  - `seq`, `states`, `trans_prob`, `emiss_prob`, `label` – Same as `forward_algo`

- **Returns:**  
  A normalized matrix of shape `(num_states, sequence_length)` representing the probability of the remaining sequence from each position for each state.

---

### `baum_welch_algo(seq, states, trans_prob, emiss_prob, label, iterations)`
Trains the HMM using the **Baum-Welch algorithm**, an unsupervised Expectation-Maximization method.

- **Parameters:**
  - `seq` *(str)* – Input sequence of symbols
  - `states` *(list of int)* – Hidden states
  - `trans_prob` *(2D array)* – Initial transition probabilities
  - `emiss_prob` *(2D array)* – Initial emission probabilities
  - `label` *(dict)* – Maps symbols to indices
  - `iterations` *(int)* – Number of training iterations to perform

- **Returns:**
  - `trans_prob` – Updated transition probability matrix
  - `emiss_prob` – Updated emission probability matrix
  - `var` – Posterior probability matrix (`var[1, i]` = probability that interval *i* is part of a gene body)

---

## Output

After training the HMM and calculating posterior probabilities, the model generates a file called `predictions.csv`, which contains your top 50,000 predicted genomic intervals that are most likely to be part of a **protein-coding gene body**.





